/***
Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.

Copyright 2016. UChicago Argonne, LLC. This software was produced
under U.S. Government contract DE-AC02-06CH11357 for Argonne National
Laboratory (ANL), which is operated by UChicago Argonne, LLC for the
U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
modified to produce derivative works, such modified software should
be clearly marked, so as not to confuse it with the version available
from ANL.

Additionally, redistribution and use in source and binary forms, with
or without modification, are permitted provided that the following
conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in
      the documentation and/or other materials provided with the
      distribution.

    * Neither the name of UChicago Argonne, LLC, Argonne National
      Laboratory, ANL, the U.S. Government, nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago
Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
***/

/// Initial Author <2026>: Arthur Glowacki


#include "fitting/routines/gpu/metal_nnls_backend.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <cmath>
#include <cstring>
#include <vector>

namespace fitting
{
namespace routines
{
namespace gpu
{

// ----------------------------------------------------------------------------

/**
 * Threads per threadgroup for the solver kernel. Must be a power of two: the
 * threadgroup reductions halve the active thread count each round.
 */
static const int TG_THREADS = 64;

// ----------------------------------------------------------------------------

/**
 * The compute kernels, compiled at load time with newLibraryWithSource. Building
 * from source rather than shipping a .metallib keeps the shaders in the repository
 * next to the host code and avoids having to locate a resource bundle at runtime.
 * The cost is a one time compile of a few milliseconds during initialize().
 *
 * Everything here is float: the Metal Shading Language has no double type.
 */
static NSString* const METAL_NNLS_SOURCE = @R"METAL(
#include <metal_stdlib>
using namespace metal;

#define TG_THREADS 64
#define MAX_N 64

// Sum a value across the threadgroup. Every thread receives the total.
// Requires a power of two thread count and uniform control flow.
inline float tg_reduce_sum(threadgroup float* scratch, uint tid, uint tcount, float val)
{
    scratch[tid] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint off = tcount >> 1; off > 0; off >>= 1)
    {
        if (tid < off)
        {
            scratch[tid] += scratch[tid + off];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    float result = scratch[0];
    // Hold everyone here so a later reduction cannot overwrite scratch[0]
    // before every thread has read it.
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return result;
}

inline float tg_reduce_max(threadgroup float* scratch, uint tid, uint tcount, float val)
{
    scratch[tid] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint off = tcount >> 1; off > 0; off >>= 1)
    {
        if (tid < off)
        {
            scratch[tid] = max(scratch[tid], scratch[tid + off]);
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    float result = scratch[0];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    return result;
}

// c = Ascaled^T * b, one thread per (component, spectrum) pair.
// This is the only kernel that walks the full channel dimension, and it is where
// the bulk of the memory traffic goes.
kernel void nnls_project_rhs(device const float*  A         [[buffer(0)]],
                             device const float*  B         [[buffer(1)]],
                             device float*        C         [[buffer(2)]],
                             constant int&        m         [[buffer(3)]],
                             constant int&        n         [[buffer(4)]],
                             constant int&        n_spectra [[buffer(5)]],
                             uint2 gid [[thread_position_in_grid]])
{
    int j = int(gid.x);
    int s = int(gid.y);
    if (j >= n || s >= n_spectra)
    {
        return;
    }
    device const float* col = A + (size_t)j * (size_t)m;
    device const float* b   = B + (size_t)s * (size_t)m;
    float acc = 0.0f;
    for (int i = 0; i < m; i++)
    {
        acc = fma(col[i], b[i], acc);
    }
    C[(size_t)s * (size_t)n + j] = acc;
}

// bnorm2[s] = ||b_s||^2, needed for the objective value the descent test compares.
kernel void nnls_bnorm2(device const float* B         [[buffer(0)]],
                        device float*       BN        [[buffer(1)]],
                        constant int&       m         [[buffer(2)]],
                        constant int&       n_spectra [[buffer(3)]],
                        uint gid [[thread_position_in_grid]])
{
    int s = int(gid);
    if (s >= n_spectra)
    {
        return;
    }
    device const float* b = B + (size_t)s * (size_t)m;
    float acc = 0.0f;
    for (int i = 0; i < m; i++)
    {
        acc = fma(b[i], b[i], acc);
    }
    BN[s] = acc;
}

// One threadgroup per spectrum. Projected gradient with Barzilai-Borwein steps and
// a non monotone descent test, mirroring nsNNLS::nnls, but driven by the normal
// equations so each iteration is O(n^2) in threadgroup memory rather than O(m n)
// against the full spectrum.
//
//   gradient = G y - c          (G = Ascaled^T Ascaled, c = Ascaled^T b)
//   2 * obj  = y^T G y - 2 y^T c + ||b||^2
//
// G stays in the constant address space: it is identical for every spectrum, so
// every thread reads the same addresses and the hardware caches it well. Keeping it
// out of threadgroup memory leaves room for high occupancy.
kernel void nnls_solve_batch(constant float*      G         [[buffer(0)]],
                             device const float*  C         [[buffer(1)]],
                             device const float*  BNORM2    [[buffer(2)]],
                             constant float*      COLNORM   [[buffer(3)]],
                             device float*        X         [[buffer(4)]],
                             device int*          ITERS     [[buffer(5)]],
                             device float*        NPG       [[buffer(6)]],
                             constant int&        n         [[buffer(7)]],
                             constant int&        maxit     [[buffer(8)]],
                             constant float&      pgtol     [[buffer(9)]],
                             constant int&        mcheck    [[buffer(10)]],
                             constant float&      sigma     [[buffer(11)]],
                             constant float&      decay     [[buffer(12)]],
                             uint tgid   [[threadgroup_position_in_grid]],
                             uint tid    [[thread_position_in_threadgroup]],
                             uint tcount [[threads_per_threadgroup]])
{
    threadgroup float sc[MAX_N];        // right hand side c
    threadgroup float sy[MAX_N];        // current iterate, in column scaled units
    threadgroup float sgrad[MAX_N];     // gradient
    threadgroup float soldy[MAX_N];
    threadgroup float soldg[MAX_N];
    threadgroup float sydelta[MAX_N];
    threadgroup float sgdelta[MAX_N];
    threadgroup float srefy[MAX_N];     // reference iterate for the descent test
    threadgroup float st[MAX_N];        // scratch for G*y
    threadgroup float sred[TG_THREADS]; // reduction scratch

    const uint s = tgid;
    const float bnorm2 = BNORM2[s];
    device const float* c_in = C + (size_t)s * (size_t)n;

    // ---- initialize -------------------------------------------------------
    // The CPU solver starts from x = 0.5. In column scaled units that same point
    // is y = 0.5 * colnorm, which keeps the two solvers starting together.
    for (uint i = tid; i < (uint)n; i += tcount)
    {
        sc[i]    = c_in[i];
        sy[i]    = 0.5f * COLNORM[i];
        soldy[i] = 0.0f;
        soldg[i] = -c_in[i];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint i = tid; i < (uint)n; i += tcount)
    {
        float acc = 0.0f;
        for (int j = 0; j < n; j++)
        {
            acc = fma(G[i * (uint)n + j], sy[j], acc);
        }
        sgrad[i] = acc - sc[i];
        srefy[i] = sy[i];
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    float beta = 1.0f;
    float obj = 0.0f;
    float obj_prev_check = 0.0f;
    bool  have_prev_check = false;
    int   iter = -1;
    int   term = 0;
    float npg = 0.0f;

    // The loop body runs once more after the termination flag is raised, matching
    // the CPU solver, which sets term at the top of the body and only tests it at
    // the bottom. Keeping that quirk keeps iteration counts comparable.
    while (term == 0)
    {
        iter++;

        if (iter >= maxit)
        {
            term = 1;
        }

        // Inf norm of the projected gradient, reported in unscaled units so that
        // pgtol means the same thing it does on the CPU. Scaling the gradient back
        // is a multiply by colnorm, since grad_x = diag(colnorm) * grad_y.
        float local_max = 0.0f;
        for (uint i = tid; i < (uint)n; i += tcount)
        {
            float g  = sgrad[i];
            float pg = (sy[i] > 0.0f) ? g : min(0.0f, g);
            local_max = max(local_max, fabs(pg * COLNORM[i]));
        }
        npg = tg_reduce_max(sred, tid, tcount, local_max);
        if (npg < pgtol)
        {
            term = 2;
        }

        // Fixed set, and the deltas the BB step is built from. A variable sitting at
        // the lower bound with a gradient pushing it further down contributes nothing.
        for (uint i = tid; i < (uint)n; i += tcount)
        {
            bool fixed  = (sy[i] == 0.0f && sgrad[i] > 0.0f);
            sydelta[i]  = fixed ? 0.0f : (sy[i] - soldy[i]);
            sgdelta[i]  = fixed ? 0.0f : (sgrad[i] - soldg[i]);
            soldy[i]    = sy[i];
            soldg[i]    = sgrad[i];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);

        // Barzilai-Borwein step, alternating the two variants by iteration parity.
        float nr_local = 0.0f;
        float dr_local = 0.0f;
        if (iter % 2)
        {
            for (uint i = tid; i < (uint)n; i += tcount)
            {
                nr_local = fma(sydelta[i], sydelta[i], nr_local);
                dr_local = fma(sydelta[i], sgdelta[i], dr_local);
            }
        }
        else
        {
            for (uint i = tid; i < (uint)n; i += tcount)
            {
                nr_local = fma(sydelta[i], sgdelta[i], nr_local);
                dr_local = fma(sgdelta[i], sgdelta[i], dr_local);
            }
        }
        float nr = tg_reduce_sum(sred, tid, tcount, nr_local);
        float dr = tg_reduce_sum(sred, tid, tcount, dr_local);
        // Guard the denominator as well as the numerator. The CPU version only
        // checks nr, which lets a zero dr produce an infinite step; in float that
        // turns the whole iterate into NaN on the next update.
        float step = (nr == 0.0f || dr == 0.0f) ? 0.0f : beta * (nr / dr);

        for (uint i = tid; i < (uint)n; i += tcount)
        {
            float v = fma(-step, sgrad[i], sy[i]);
            sy[i] = (v < 0.0f) ? 0.0f : v;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);

        // Objective and gradient from the normal equations.
        for (uint i = tid; i < (uint)n; i += tcount)
        {
            float acc = 0.0f;
            for (int j = 0; j < n; j++)
            {
                acc = fma(G[i * (uint)n + j], sy[j], acc);
            }
            st[i] = acc;
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);

        float ytgy_local = 0.0f;
        float ytc_local  = 0.0f;
        for (uint i = tid; i < (uint)n; i += tcount)
        {
            ytgy_local = fma(sy[i], st[i], ytgy_local);
            ytc_local  = fma(sy[i], sc[i], ytc_local);
            sgrad[i]   = st[i] - sc[i];
        }
        float ytgy = tg_reduce_sum(sred, tid, tcount, ytgy_local);
        float ytc  = tg_reduce_sum(sred, tid, tcount, ytc_local);
        obj = 0.5f * (ytgy - 2.0f * ytc + bnorm2);
        threadgroup_barrier(mem_flags::mem_threadgroup);

        // Non monotone descent test, every mcheck iterations. On failure the step
        // is damped; on success the reference iterate advances.
        if (iter % mcheck == 0)
        {
            float d_local = 0.0f;
            for (uint i = tid; i < (uint)n; i += tcount)
            {
                d_local = fma(sgrad[i], srefy[i] - sy[i], d_local);
            }
            float d = sigma * tg_reduce_sum(sred, tid, tcount, d_local);

            d = have_prev_check ? (obj_prev_check - obj - d) : (obj - d);

            if (d < 0.0f)
            {
                beta *= decay;
            }
            else
            {
                for (uint i = tid; i < (uint)n; i += tcount)
                {
                    srefy[i] = sy[i];
                }
                threadgroup_barrier(mem_flags::mem_threadgroup);
            }
            obj_prev_check  = obj;
            have_prev_check = true;
        }
    }

    // Undo the column scaling: x = y / colnorm.
    for (uint i = tid; i < (uint)n; i += tcount)
    {
        X[(size_t)s * (size_t)n + i] = sy[i] / COLNORM[i];
    }
    if (tid == 0)
    {
        ITERS[s] = iter;
        NPG[s]   = npg;
    }
}
)METAL";

/**
 * The SNIP kernel, kept in its own source string because it is compiled with fast math
 * DISABLED while METAL_NNLS_SOURCE is compiled with it enabled.
 *
 * Two reasons for the split. The expand step relies on isfinite() to reject the overflow
 * the CPU routine also guards against, and fast math lets the compiler assume no Inf or
 * NaN ever occurs, which permits folding that test away. And unlike the solver, which is
 * iterative and convergence tested, this kernel is a straight-line port whose whole job
 * is to reproduce data_struct::snip_background channel for channel, so relaxed transcendentals
 * would cost accuracy against the CPU for no speed that matters: the kernel is dominated by
 * the serial clipping sweep, not by the log/exp passes.
 */
static NSString* const METAL_SNIP_SOURCE = @R"METAL(
#include <metal_stdlib>
using namespace metal;

// One threadgroup per spectrum. Computes the SNIP background over the full spectrum, then
// writes the background-subtracted, clamped fit-range slice to OUTB (the solver's b) and
// the fit-range background to OUT_BKG. This is a faithful float port of
// data_struct::snip_background (spectra.h): a boxcar(5) convolution, a log-log compress,
// n_first + n_decay clipping passes with the width schedule, then an exp-exp expand.
//
// The clipping recurrence updates the background in place while sweeping channels
// ascending (background[k] depends on an already-updated lo_index < k), so it is
// inherently sequential within a pass and is run by a single thread over one threadgroup
// working buffer. The device parallelism is across the thousands of spectra in the batch
// (one threadgroup each), not across channels within a spectrum.
//
// The working buffer holds full_len floats and is allocated by the host with
// setThreadgroupMemoryLength; the host guarantees full_len <= METAL_SNIP_MAX_CHANNELS and
// that the allocation fits the device before encoding.
kernel void snip_background_kernel(device const float*  FULL            [[buffer(0)]],
                                   device const float*  BASE_WIDTH      [[buffer(1)]],
                                   constant int&        full_len        [[buffer(2)]],
                                   constant int&        n_first_passes  [[buffer(3)]],
                                   constant int&        n_decay_passes  [[buffer(4)]],
                                   constant int&        lo_limit        [[buffer(5)]],
                                   constant int&        up_limit        [[buffer(6)]],
                                   constant int&        fit_min         [[buffer(7)]],
                                   constant int&        n_fit           [[buffer(8)]],
                                   device float*        OUTB            [[buffer(9)]],
                                   device float*        OUT_BKG         [[buffer(10)]],
                                   constant int&        has_bkg_out     [[buffer(11)]],
                                   threadgroup float*   cur             [[threadgroup(0)]],
                                   uint tgid   [[threadgroup_position_in_grid]],
                                   uint tid    [[thread_position_in_threadgroup]],
                                   uint tcount [[threads_per_threadgroup]])
{
    const uint s = tgid;
    device const float* spec = FULL + (size_t)s * (size_t)full_len;

    // ---- boxcar(5) convolution, matching convolve1d (spectra.h:310) -----------------
    // convolve1d centers a length-5 boxcar: output channel k is the mean of the 5 input
    // channels centered on k, and the first/last 2 channels are left zero (the valid
    // region is [2, full_len-3]). Reproduce that exactly so the edges agree with the CPU.
    for (uint k = tid; k < (uint)full_len; k += tcount)
    {
        float v = 0.0f;
        if (k >= 2u && k + 2u < (uint)full_len)
        {
            v = (spec[k - 2] + spec[k - 1] + spec[k] + spec[k + 1] + spec[k + 2]) * 0.2f;
        }
        // log(log(x+1)+1)
        cur[k] = log(log(v + 1.0f) + 1.0f);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // ---- clipping passes ------------------------------------------------------------
    // Pass p uses width[k] * (1/sqrt2)^(decay index). The first n_first_passes run at the
    // base width (decay index 0); each subsequent pass divides the width by sqrt2, matching
    // the CPU's two fixed passes followed by the while loop that halves current_width.
    //
    // The CPU sweep (spectra.h:401-433) updates background[k] IN PLACE while iterating k
    // ascending, so background[k] reads a lo_index < k that was already lowered earlier in
    // the same pass. That within-pass dependency is sequential by construction, so the
    // sweep is done by a single thread per threadgroup (tid 0) over the shared buffer `cur`.
    // The parallel parts (convolution above, expand below) still use the whole threadgroup;
    // only this inherently serial recurrence is single threaded. One threadgroup still
    // handles one spectrum, so the device parallelism is across the batch.
    const float inv_sqrt2 = 0.70710678118654752440f;
    const int total_passes = n_first_passes + n_decay_passes;

    if (tid == 0)
    {
        for (int p = 0; p < total_passes; p++)
        {
            // decay index: 0 for the first passes, then 0,1,2,... over the decay passes.
            const int decay_idx = (p < n_first_passes) ? 0 : (p - n_first_passes);
            float scale = 1.0f;
            for (int d = 0; d < decay_idx; d++)
            {
                scale *= inv_sqrt2;
            }

            for (int k = 0; k < full_len; k++)
            {
                float w = BASE_WIDTH[k] * scale;
                // Reproduce the CPU index math (spectra.h:405-408), including its
                // asymmetry: lo_index is clamped to [lo_limit, up_limit]; hi_index is only
                // capped at up_limit (the lo_limit floor on line 407 is overwritten by 408).
                int lo_index = (int)((float)k - w);
                if (lo_index < lo_limit) lo_index = lo_limit;
                if (lo_index > up_limit) lo_index = up_limit;
                int hi_index = (int)((float)k + w);
                if (hi_index > up_limit) hi_index = up_limit;
                if (hi_index < 0) hi_index = 0;

                float temp = (cur[lo_index] + cur[hi_index]) * 0.5f;
                if (cur[k] > temp)
                {
                    cur[k] = temp;   // in place, exactly like the CPU
                }
            }
        }
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    // ---- expand back: exp(exp(x)-1)-1, finite-guarded (spectra.h:435-436) ------------
    // Then subtract from the raw spectrum over the fit range, clamp >= 0, and emit.
    for (uint t = tid; t < (uint)n_fit; t += tcount)
    {
        const int k = fit_min + (int)t;
        float bg = exp(exp(cur[k]) - 1.0f) - 1.0f;
        if (!isfinite(bg)) bg = 0.0f;

        float sub = spec[k] - bg;
        if (sub < 0.0f) sub = 0.0f;

        OUTB[(size_t)s * (size_t)n_fit + t] = sub;
        if (has_bkg_out != 0)
        {
            OUT_BKG[(size_t)s * (size_t)n_fit + t] = bg;
        }
    }
}
)METAL";

// ----------------------------------------------------------------------------

struct Metal_NNLS_Backend::Impl
{
    id<MTLDevice>               device          = nil;
    id<MTLCommandQueue>         queue           = nil;
    id<MTLComputePipelineState> pipe_rhs        = nil;
    id<MTLComputePipelineState> pipe_bnorm2     = nil;
    id<MTLComputePipelineState> pipe_solve      = nil;
    id<MTLComputePipelineState> pipe_snip       = nil;

    // Per fit matrix, rebuilt by set_fit_matrix.
    id<MTLBuffer> buf_A       = nil;    // m x n, column major, columns unit norm
    id<MTLBuffer> buf_G       = nil;    // n x n, row major
    id<MTLBuffer> buf_colnorm = nil;    // n

    // Per batch, grown on demand and reused across calls.
    id<MTLBuffer> buf_B       = nil;
    id<MTLBuffer> buf_C       = nil;
    id<MTLBuffer> buf_bnorm2  = nil;
    id<MTLBuffer> buf_X       = nil;
    id<MTLBuffer> buf_iters   = nil;
    id<MTLBuffer> buf_npg     = nil;

    // For the on-device SNIP path, grown on demand alongside the per batch buffers.
    id<MTLBuffer> buf_full    = nil;    // n_spectra x full_len, raw spectra
    id<MTLBuffer> buf_bkg     = nil;    // n_spectra x n_fit, returned background
    id<MTLBuffer> buf_width   = nil;    // full_len, base width in channels

    int batch_capacity = 0;
    int full_capacity  = 0;             // n_spectra the buf_full / buf_bkg buffers hold
    int full_len_cap   = 0;             // full_len the buf_full / buf_width buffers were sized for
    int fit_len_cap    = 0;             // n_fit the buf_bkg buffer was sized for

    int tg_threads = TG_THREADS;

    int snip_tg_threads = TG_THREADS;

    // Bytes of threadgroup memory the device allows one kernel, from
    // maxThreadgroupMemoryLength. Caps the full spectrum length the SNIP kernel accepts.
    int max_tg_memory = 0;

    /**
     * Grow the per batch buffers to hold n_spectra. Shared storage mode keeps them
     * directly addressable from the CPU, which on Apple Silicon's unified memory
     * means the host writes land where the GPU reads with no staging copy.
     */
    bool ensure_batch_capacity(int n_spectra, int m, int n)
    {
        if (n_spectra <= batch_capacity && buf_B != nil)
        {
            return true;
        }
        const MTLResourceOptions opts = MTLResourceStorageModeShared;
        buf_B      = [device newBufferWithLength:sizeof(float) * (size_t)n_spectra * (size_t)m options:opts];
        buf_C      = [device newBufferWithLength:sizeof(float) * (size_t)n_spectra * (size_t)n options:opts];
        buf_bnorm2 = [device newBufferWithLength:sizeof(float) * (size_t)n_spectra options:opts];
        buf_X      = [device newBufferWithLength:sizeof(float) * (size_t)n_spectra * (size_t)n options:opts];
        buf_iters  = [device newBufferWithLength:sizeof(int)   * (size_t)n_spectra options:opts];
        buf_npg    = [device newBufferWithLength:sizeof(float) * (size_t)n_spectra options:opts];

        if (buf_B == nil || buf_C == nil || buf_bnorm2 == nil || buf_X == nil || buf_iters == nil || buf_npg == nil)
        {
            batch_capacity = 0;
            return false;
        }
        batch_capacity = n_spectra;
        return true;
    }

    /**
     * Grow the SNIP staging buffers (raw full spectra in, fit-range background out) to
     * hold n_spectra at the given dimensions, reusing them across calls.
     */
    bool ensure_snip_capacity(int n_spectra, int full_len, int n_fit)
    {
        if (n_spectra <= full_capacity && full_len <= full_len_cap && n_fit <= fit_len_cap
            && buf_full != nil && buf_bkg != nil && buf_width != nil)
        {
            return true;
        }
        const MTLResourceOptions opts = MTLResourceStorageModeShared;
        buf_full  = [device newBufferWithLength:sizeof(float) * (size_t)n_spectra * (size_t)full_len options:opts];
        buf_bkg   = [device newBufferWithLength:sizeof(float) * (size_t)n_spectra * (size_t)n_fit options:opts];
        buf_width = [device newBufferWithLength:sizeof(float) * (size_t)full_len options:opts];

        if (buf_full == nil || buf_bkg == nil || buf_width == nil)
        {
            full_capacity = 0;
            full_len_cap  = 0;
            fit_len_cap   = 0;
            return false;
        }
        full_capacity = n_spectra;
        full_len_cap  = full_len;
        fit_len_cap   = n_fit;
        return true;
    }
};

// ----------------------------------------------------------------------------

bool Metal_NNLS_Backend::is_available()
{
    @autoreleasepool
    {
        id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
        return dev != nil;
    }
}

// ----------------------------------------------------------------------------

Metal_NNLS_Backend::Metal_NNLS_Backend() : _impl(nullptr), _m(0), _n(0)
{

}

// ----------------------------------------------------------------------------

Metal_NNLS_Backend::~Metal_NNLS_Backend()
{
    delete _impl;
    _impl = nullptr;
}

// ----------------------------------------------------------------------------

bool Metal_NNLS_Backend::initialize()
{
    if (_impl != nullptr)
    {
        return true;
    }

    @autoreleasepool
    {
        Impl* impl = new Impl();

        impl->device = MTLCreateSystemDefaultDevice();
        if (impl->device == nil)
        {
            _last_error = "No Metal capable device found";
            delete impl;
            return false;
        }
        _device_name = std::string([[impl->device name] UTF8String]);

        impl->queue = [impl->device newCommandQueue];
        if (impl->queue == nil)
        {
            _last_error = "Could not create Metal command queue";
            delete impl;
            return false;
        }

        NSError* err = nil;
        MTLCompileOptions* opts = [[MTLCompileOptions alloc] init];
        // The solver is iterative and convergence tested, not bit reproducible, so
        // relaxed math is a fair trade for the speed.
        [opts setFastMathEnabled:YES];

        id<MTLLibrary> lib = [impl->device newLibraryWithSource:METAL_NNLS_SOURCE options:opts error:&err];
        if (lib == nil)
        {
            _last_error = "Metal shader compile failed: ";
            if (err != nil)
            {
                _last_error += std::string([[err localizedDescription] UTF8String]);
            }
            delete impl;
            return false;
        }

        // The SNIP kernel is built from its own library with fast math off; see
        // METAL_SNIP_SOURCE for why. A failure here is not fatal: the routine simply
        // keeps doing SNIP on the host and calling solve_batch.
        MTLCompileOptions* snip_opts = [[MTLCompileOptions alloc] init];
        [snip_opts setFastMathEnabled:NO];
        id<MTLLibrary> snip_lib = [impl->device newLibraryWithSource:METAL_SNIP_SOURCE
                                                             options:snip_opts
                                                               error:nil];

        struct { const char* name; __strong id<MTLComputePipelineState>* dest; id<MTLLibrary> lib; bool required; } kernels[] = {
            { "nnls_project_rhs",      &impl->pipe_rhs,    lib,      true },
            { "nnls_bnorm2",           &impl->pipe_bnorm2, lib,      true },
            { "nnls_solve_batch",      &impl->pipe_solve,  lib,      true },
            { "snip_background_kernel", &impl->pipe_snip,  snip_lib, false },
        };

        for (auto& k : kernels)
        {
            id<MTLFunction> fn = (k.lib == nil) ? nil
                                 : [k.lib newFunctionWithName:[NSString stringWithUTF8String:k.name]];
            if (fn != nil)
            {
                err = nil;
                *k.dest = [impl->device newComputePipelineStateWithFunction:fn error:&err];
            }
            if (*k.dest == nil && k.required)
            {
                _last_error = std::string("Could not build pipeline for ") + k.name;
                if (err != nil)
                {
                    _last_error += ": " + std::string([[err localizedDescription] UTF8String]);
                }
                delete impl;
                return false;
            }
        }

        // The reductions need a power of two thread count that the solver pipeline
        // will actually accept. TG_THREADS is the target; halve it if the pipeline
        // reports a lower ceiling.
        int max_threads = (int)[impl->pipe_solve maxTotalThreadsPerThreadgroup];
        int threads = TG_THREADS;
        while (threads > max_threads && threads > 1)
        {
            threads >>= 1;
        }
        impl->tg_threads = threads;

        // The SNIP kernel does no reductions, so it only needs a count its own pipeline
        // accepts; it still gets the same target so the two dispatches look alike.
        if (impl->pipe_snip != nil)
        {
            int snip_max = (int)[impl->pipe_snip maxTotalThreadsPerThreadgroup];
            int snip_threads = TG_THREADS;
            while (snip_threads > snip_max && snip_threads > 1)
            {
                snip_threads >>= 1;
            }
            impl->snip_tg_threads = snip_threads;
            impl->max_tg_memory = (int)[impl->device maxThreadgroupMemoryLength];
        }

        _impl = impl;
    }

    return true;
}

// ----------------------------------------------------------------------------

bool Metal_NNLS_Backend::set_fit_matrix(const double* const fitmatrix_col_major, int m, int n)
{
    if (_impl == nullptr)
    {
        _last_error = "set_fit_matrix called before initialize";
        return false;
    }
    if (m <= 0 || n <= 0)
    {
        _last_error = "Fit matrix has a zero dimension";
        return false;
    }
    if (n > METAL_NNLS_MAX_COMPONENTS)
    {
        _last_error = "Too many components for the Metal solver";
        return false;
    }

    _m = m;
    _n = n;

    // Column norms, in double. A zero column (an element model that is identically
    // zero over the range) gets a norm of one so the scaling stays well defined; the
    // column contributes nothing either way.
    std::vector<double> colnorm((size_t)n, 1.0);
    for (int j = 0; j < n; j++)
    {
        const double* col = fitmatrix_col_major + (size_t)j * (size_t)m;
        double acc = 0.0;
        for (int i = 0; i < m; i++)
        {
            acc += col[i] * col[i];
        }
        double nrm = std::sqrt(acc);
        colnorm[(size_t)j] = (nrm > 0.0) ? nrm : 1.0;
    }

    const MTLResourceOptions opts = MTLResourceStorageModeShared;

    @autoreleasepool
    {
        _impl->buf_A       = [_impl->device newBufferWithLength:sizeof(float) * (size_t)m * (size_t)n options:opts];
        _impl->buf_G       = [_impl->device newBufferWithLength:sizeof(float) * (size_t)n * (size_t)n options:opts];
        _impl->buf_colnorm = [_impl->device newBufferWithLength:sizeof(float) * (size_t)n options:opts];

        if (_impl->buf_A == nil || _impl->buf_G == nil || _impl->buf_colnorm == nil)
        {
            _last_error = "Could not allocate Metal buffers for the fit matrix";
            return false;
        }

        float* dst_A  = (float*)[_impl->buf_A contents];
        float* dst_G  = (float*)[_impl->buf_G contents];
        float* dst_cn = (float*)[_impl->buf_colnorm contents];

        for (int j = 0; j < n; j++)
        {
            dst_cn[j] = (float)colnorm[(size_t)j];
            const double* col = fitmatrix_col_major + (size_t)j * (size_t)m;
            float* dcol = dst_A + (size_t)j * (size_t)m;
            const double inv = 1.0 / colnorm[(size_t)j];
            for (int i = 0; i < m; i++)
            {
                dcol[i] = (float)(col[i] * inv);
            }
        }

        // G = Ascaled^T Ascaled, accumulated in double then narrowed. Symmetric, so
        // only the upper triangle is computed and then mirrored.
        for (int a = 0; a < n; a++)
        {
            const double* ca = fitmatrix_col_major + (size_t)a * (size_t)m;
            const double inv_a = 1.0 / colnorm[(size_t)a];
            for (int b = a; b < n; b++)
            {
                const double* cb = fitmatrix_col_major + (size_t)b * (size_t)m;
                const double inv_b = 1.0 / colnorm[(size_t)b];
                double acc = 0.0;
                for (int i = 0; i < m; i++)
                {
                    acc += ca[i] * cb[i];
                }
                float v = (float)(acc * inv_a * inv_b);
                dst_G[(size_t)a * (size_t)n + b] = v;
                dst_G[(size_t)b * (size_t)n + a] = v;
            }
        }
    }

    return true;
}

// ----------------------------------------------------------------------------

int Metal_NNLS_Backend::recommended_batch_size() const
{
    // One threadgroup per spectrum, so the batch wants to be at least a few times
    // the core count to keep every core fed, while staying small enough that the
    // staging buffers remain a modest fraction of memory.
    return 4096;
}

// ----------------------------------------------------------------------------

bool Metal_NNLS_Backend::solve_batch(const float* const spectra_batch,
                                     int n_spectra,
                                     float* out_x,
                                     int* out_num_iter,
                                     float* out_npg,
                                     int max_iter,
                                     float pgtol)
{
    if (_impl == nullptr)
    {
        _last_error = "solve_batch called before initialize";
        return false;
    }
    if (_impl->buf_A == nil)
    {
        _last_error = "solve_batch called before set_fit_matrix";
        return false;
    }
    if (n_spectra <= 0)
    {
        return true;
    }

    const int m = _m;
    const int n = _n;

    // Mirror the CPU solver's convergence controls so the two behave alike.
    const int   mcheck = 100;
    const float sigma  = 0.01f;
    const float decay  = 0.9f;

    @autoreleasepool
    {
        if (false == _impl->ensure_batch_capacity(n_spectra, m, n))
        {
            _last_error = "Could not allocate Metal batch buffers";
            return false;
        }

        std::memcpy([_impl->buf_B contents], spectra_batch, sizeof(float) * (size_t)n_spectra * (size_t)m);

        id<MTLCommandBuffer> cmd = [_impl->queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

        // c = Ascaled^T b
        [enc setComputePipelineState:_impl->pipe_rhs];
        [enc setBuffer:_impl->buf_A offset:0 atIndex:0];
        [enc setBuffer:_impl->buf_B offset:0 atIndex:1];
        [enc setBuffer:_impl->buf_C offset:0 atIndex:2];
        [enc setBytes:&m         length:sizeof(int) atIndex:3];
        [enc setBytes:&n         length:sizeof(int) atIndex:4];
        [enc setBytes:&n_spectra length:sizeof(int) atIndex:5];
        [enc dispatchThreads:MTLSizeMake((NSUInteger)n, (NSUInteger)n_spectra, 1)
              threadsPerThreadgroup:MTLSizeMake((NSUInteger)((n < 8) ? n : 8), 8, 1)];

        // ||b||^2
        [enc setComputePipelineState:_impl->pipe_bnorm2];
        [enc setBuffer:_impl->buf_B      offset:0 atIndex:0];
        [enc setBuffer:_impl->buf_bnorm2 offset:0 atIndex:1];
        [enc setBytes:&m         length:sizeof(int) atIndex:2];
        [enc setBytes:&n_spectra length:sizeof(int) atIndex:3];
        [enc dispatchThreads:MTLSizeMake((NSUInteger)n_spectra, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];

        // The solve, one threadgroup per spectrum.
        [enc setComputePipelineState:_impl->pipe_solve];
        [enc setBuffer:_impl->buf_G       offset:0 atIndex:0];
        [enc setBuffer:_impl->buf_C       offset:0 atIndex:1];
        [enc setBuffer:_impl->buf_bnorm2  offset:0 atIndex:2];
        [enc setBuffer:_impl->buf_colnorm offset:0 atIndex:3];
        [enc setBuffer:_impl->buf_X       offset:0 atIndex:4];
        [enc setBuffer:_impl->buf_iters   offset:0 atIndex:5];
        [enc setBuffer:_impl->buf_npg     offset:0 atIndex:6];
        [enc setBytes:&n        length:sizeof(int)   atIndex:7];
        [enc setBytes:&max_iter length:sizeof(int)   atIndex:8];
        [enc setBytes:&pgtol    length:sizeof(float) atIndex:9];
        [enc setBytes:&mcheck   length:sizeof(int)   atIndex:10];
        [enc setBytes:&sigma    length:sizeof(float) atIndex:11];
        [enc setBytes:&decay    length:sizeof(float) atIndex:12];
        [enc dispatchThreadgroups:MTLSizeMake((NSUInteger)n_spectra, 1, 1)
              threadsPerThreadgroup:MTLSizeMake((NSUInteger)_impl->tg_threads, 1, 1)];

        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];

        if ([cmd error] != nil)
        {
            _last_error = std::string("Metal dispatch failed: ") +
                          std::string([[[cmd error] localizedDescription] UTF8String]);
            return false;
        }

        std::memcpy(out_x, [_impl->buf_X contents], sizeof(float) * (size_t)n_spectra * (size_t)n);
        if (out_num_iter != nullptr)
        {
            std::memcpy(out_num_iter, [_impl->buf_iters contents], sizeof(int) * (size_t)n_spectra);
        }
        if (out_npg != nullptr)
        {
            std::memcpy(out_npg, [_impl->buf_npg contents], sizeof(float) * (size_t)n_spectra);
        }
    }

    return true;
}

// ----------------------------------------------------------------------------

bool Metal_NNLS_Backend::snip_supported(int full_len) const
{
    if (_impl == nullptr || _impl->pipe_snip == nil)
    {
        return false;
    }
    if (full_len <= 0 || full_len > METAL_SNIP_MAX_CHANNELS)
    {
        return false;
    }
    // The kernel keeps one full-length float working buffer in threadgroup memory.
    return (int)(sizeof(float) * (size_t)full_len) <= _impl->max_tg_memory;
}

// ----------------------------------------------------------------------------

bool Metal_NNLS_Backend::solve_batch_snip(const float* const full_spectra,
                                          int n_spectra,
                                          int full_len,
                                          const float* const base_width_chan,
                                          int n_first_passes,
                                          int n_decay_passes,
                                          int lo_limit,
                                          int up_limit,
                                          int fit_min,
                                          int n_fit,
                                          float* out_x,
                                          float* out_background,
                                          int* out_num_iter,
                                          float* out_npg,
                                          int max_iter,
                                          float pgtol)
{
    if (_impl == nullptr)
    {
        _last_error = "solve_batch_snip called before initialize";
        return false;
    }
    if (_impl->buf_A == nil)
    {
        _last_error = "solve_batch_snip called before set_fit_matrix";
        return false;
    }
    if (n_spectra <= 0)
    {
        return true;
    }
    if (false == snip_supported(full_len))
    {
        _last_error = "Full spectrum length exceeds the Metal SNIP threadgroup-memory limit";
        return false;
    }
    if (n_fit != _m)
    {
        _last_error = "solve_batch_snip n_fit does not match the fit matrix rows";
        return false;
    }
    if (fit_min < 0 || fit_min + n_fit > full_len)
    {
        _last_error = "solve_batch_snip fit range lies outside the full spectrum";
        return false;
    }

    const int m = _m;   // == n_fit
    const int n = _n;

    // Mirror the CPU solver's convergence controls so the two behave alike.
    const int   mcheck = 100;
    const float sigma  = 0.01f;
    const float decay  = 0.9f;

    const int has_bkg_out = (out_background != nullptr) ? 1 : 0;

    @autoreleasepool
    {
        if (false == _impl->ensure_batch_capacity(n_spectra, m, n))
        {
            _last_error = "Could not allocate Metal batch buffers";
            return false;
        }
        if (false == _impl->ensure_snip_capacity(n_spectra, full_len, n_fit))
        {
            _last_error = "Could not allocate Metal SNIP buffers";
            return false;
        }

        // Unified memory, so these land where the GPU reads them with no staging copy.
        // The width is tiny and constant across the run; re-writing it per batch keeps
        // the API stateless.
        std::memcpy([_impl->buf_full contents], full_spectra,
                    sizeof(float) * (size_t)n_spectra * (size_t)full_len);
        std::memcpy([_impl->buf_width contents], base_width_chan,
                    sizeof(float) * (size_t)full_len);

        id<MTLCommandBuffer> cmd = [_impl->queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

        // SNIP background + subtract + clamp, one threadgroup per spectrum. The kernel
        // writes the fit-range solver input into buf_B and the fit-range background into
        // buf_bkg. Threadgroup allocations are rounded up to 16 bytes.
        {
            const NSUInteger tg_mem = ((sizeof(float) * (size_t)full_len) + 15u) & ~((size_t)15u);
            [enc setComputePipelineState:_impl->pipe_snip];
            [enc setBuffer:_impl->buf_full  offset:0 atIndex:0];
            [enc setBuffer:_impl->buf_width offset:0 atIndex:1];
            [enc setBytes:&full_len       length:sizeof(int) atIndex:2];
            [enc setBytes:&n_first_passes length:sizeof(int) atIndex:3];
            [enc setBytes:&n_decay_passes length:sizeof(int) atIndex:4];
            [enc setBytes:&lo_limit       length:sizeof(int) atIndex:5];
            [enc setBytes:&up_limit       length:sizeof(int) atIndex:6];
            [enc setBytes:&fit_min        length:sizeof(int) atIndex:7];
            [enc setBytes:&n_fit          length:sizeof(int) atIndex:8];
            [enc setBuffer:_impl->buf_B   offset:0 atIndex:9];
            [enc setBuffer:_impl->buf_bkg offset:0 atIndex:10];
            [enc setBytes:&has_bkg_out    length:sizeof(int) atIndex:11];
            [enc setThreadgroupMemoryLength:tg_mem atIndex:0];
            [enc dispatchThreadgroups:MTLSizeMake((NSUInteger)n_spectra, 1, 1)
                  threadsPerThreadgroup:MTLSizeMake((NSUInteger)_impl->snip_tg_threads, 1, 1)];
        }

        // c = Ascaled^T b
        [enc setComputePipelineState:_impl->pipe_rhs];
        [enc setBuffer:_impl->buf_A offset:0 atIndex:0];
        [enc setBuffer:_impl->buf_B offset:0 atIndex:1];
        [enc setBuffer:_impl->buf_C offset:0 atIndex:2];
        [enc setBytes:&m         length:sizeof(int) atIndex:3];
        [enc setBytes:&n         length:sizeof(int) atIndex:4];
        [enc setBytes:&n_spectra length:sizeof(int) atIndex:5];
        [enc dispatchThreads:MTLSizeMake((NSUInteger)n, (NSUInteger)n_spectra, 1)
              threadsPerThreadgroup:MTLSizeMake((NSUInteger)((n < 8) ? n : 8), 8, 1)];

        // ||b||^2
        [enc setComputePipelineState:_impl->pipe_bnorm2];
        [enc setBuffer:_impl->buf_B      offset:0 atIndex:0];
        [enc setBuffer:_impl->buf_bnorm2 offset:0 atIndex:1];
        [enc setBytes:&m         length:sizeof(int) atIndex:2];
        [enc setBytes:&n_spectra length:sizeof(int) atIndex:3];
        [enc dispatchThreads:MTLSizeMake((NSUInteger)n_spectra, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(64, 1, 1)];

        // The solve, one threadgroup per spectrum.
        [enc setComputePipelineState:_impl->pipe_solve];
        [enc setBuffer:_impl->buf_G       offset:0 atIndex:0];
        [enc setBuffer:_impl->buf_C       offset:0 atIndex:1];
        [enc setBuffer:_impl->buf_bnorm2  offset:0 atIndex:2];
        [enc setBuffer:_impl->buf_colnorm offset:0 atIndex:3];
        [enc setBuffer:_impl->buf_X       offset:0 atIndex:4];
        [enc setBuffer:_impl->buf_iters   offset:0 atIndex:5];
        [enc setBuffer:_impl->buf_npg     offset:0 atIndex:6];
        [enc setBytes:&n        length:sizeof(int)   atIndex:7];
        [enc setBytes:&max_iter length:sizeof(int)   atIndex:8];
        [enc setBytes:&pgtol    length:sizeof(float) atIndex:9];
        [enc setBytes:&mcheck   length:sizeof(int)   atIndex:10];
        [enc setBytes:&sigma    length:sizeof(float) atIndex:11];
        [enc setBytes:&decay    length:sizeof(float) atIndex:12];
        [enc dispatchThreadgroups:MTLSizeMake((NSUInteger)n_spectra, 1, 1)
              threadsPerThreadgroup:MTLSizeMake((NSUInteger)_impl->tg_threads, 1, 1)];

        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];

        if ([cmd error] != nil)
        {
            _last_error = std::string("Metal dispatch failed: ") +
                          std::string([[[cmd error] localizedDescription] UTF8String]);
            return false;
        }

        std::memcpy(out_x, [_impl->buf_X contents], sizeof(float) * (size_t)n_spectra * (size_t)n);
        if (out_background != nullptr)
        {
            std::memcpy(out_background, [_impl->buf_bkg contents],
                        sizeof(float) * (size_t)n_spectra * (size_t)n_fit);
        }
        if (out_num_iter != nullptr)
        {
            std::memcpy(out_num_iter, [_impl->buf_iters contents], sizeof(int) * (size_t)n_spectra);
        }
        if (out_npg != nullptr)
        {
            std::memcpy(out_npg, [_impl->buf_npg contents], sizeof(float) * (size_t)n_spectra);
        }
    }

    return true;
}

// ----------------------------------------------------------------------------

} //namespace gpu
} //namespace routines
} //namespace fitting
