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


#include "fitting/routines/gpu/cuda_nnls_backend.h"

#include <cuda_runtime.h>

#include <cmath>
#include <cstring>
#include <string>
#include <vector>

namespace fitting
{
namespace routines
{
namespace gpu
{

// ----------------------------------------------------------------------------

/**
 * Threads per block for the solver kernel. Must be a power of two: the block
 * reductions halve the active thread count each round. Mirrors the Metal backend's
 * TG_THREADS so the two solvers reduce in the same order.
 */
static const int TG_THREADS = 64;

/**
 * Must match CUDA_NNLS_MAX_COMPONENTS; used to size the shared memory working set.
 */
#define MAX_N 64

// ----------------------------------------------------------------------------
// Device kernels, ported 1:1 from the Metal Shading Language source in
// metal_nnls_backend.mm. Everything is float to match that backend; results agree
// with the CPU solver only to single precision.
// ----------------------------------------------------------------------------

// Sum a value across the block. Every thread receives the total.
// Requires a power of two thread count and uniform control flow.
__device__ inline float tg_reduce_sum(float* scratch, unsigned tid, unsigned tcount, float val)
{
    scratch[tid] = val;
    __syncthreads();
    for (unsigned off = tcount >> 1; off > 0; off >>= 1)
    {
        if (tid < off)
        {
            scratch[tid] += scratch[tid + off];
        }
        __syncthreads();
    }
    float result = scratch[0];
    // Hold everyone here so a later reduction cannot overwrite scratch[0]
    // before every thread has read it.
    __syncthreads();
    return result;
}

__device__ inline float tg_reduce_max(float* scratch, unsigned tid, unsigned tcount, float val)
{
    scratch[tid] = val;
    __syncthreads();
    for (unsigned off = tcount >> 1; off > 0; off >>= 1)
    {
        if (tid < off)
        {
            scratch[tid] = fmaxf(scratch[tid], scratch[tid + off]);
        }
        __syncthreads();
    }
    float result = scratch[0];
    __syncthreads();
    return result;
}

// c = Ascaled^T * b, one thread per (component, spectrum) pair.
// This is the only kernel that walks the full channel dimension, and it is where
// the bulk of the memory traffic goes.
__global__ void nnls_project_rhs(const float* __restrict__ A,
                                 const float* __restrict__ B,
                                 float* __restrict__ C,
                                 int m,
                                 int n,
                                 int n_spectra)
{
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int s = blockIdx.y * blockDim.y + threadIdx.y;
    if (j >= n || s >= n_spectra)
    {
        return;
    }
    const float* col = A + (size_t)j * (size_t)m;
    const float* b   = B + (size_t)s * (size_t)m;
    float acc = 0.0f;
    for (int i = 0; i < m; i++)
    {
        acc = fmaf(col[i], b[i], acc);
    }
    C[(size_t)s * (size_t)n + j] = acc;
}

// bnorm2[s] = ||b_s||^2, needed for the objective value the descent test compares.
__global__ void nnls_bnorm2(const float* __restrict__ B,
                            float* __restrict__ BN,
                            int m,
                            int n_spectra)
{
    int s = blockIdx.x * blockDim.x + threadIdx.x;
    if (s >= n_spectra)
    {
        return;
    }
    const float* b = B + (size_t)s * (size_t)m;
    float acc = 0.0f;
    for (int i = 0; i < m; i++)
    {
        acc = fmaf(b[i], b[i], acc);
    }
    BN[s] = acc;
}

// One block per spectrum. Projected gradient with Barzilai-Borwein steps and a non
// monotone descent test, mirroring nsNNLS::nnls, but driven by the normal equations
// so each iteration is O(n^2) in shared memory rather than O(m n) against the full
// spectrum.
//
//   gradient = G y - c          (G = Ascaled^T Ascaled, c = Ascaled^T b)
//   2 * obj  = y^T G y - 2 y^T c + ||b||^2
//
// G is read through the normal global path (__restrict__): it is identical for every
// spectrum, so every thread reads the same addresses and the cache serves it well.
// Metal keeps G in the constant address space; CUDA constant memory would cap at
// 64 KB and re-introduce a hard size limit, so a cached global read is used instead.
// Numerically identical.
__global__ void nnls_solve_batch(const float* __restrict__ G,
                                 const float* __restrict__ C,
                                 const float* __restrict__ BNORM2,
                                 const float* __restrict__ COLNORM,
                                 float* __restrict__ X,
                                 int* __restrict__ ITERS,
                                 float* __restrict__ NPG,
                                 int n,
                                 int maxit,
                                 float pgtol,
                                 int mcheck,
                                 float sigma,
                                 float decay)
{
    __shared__ float sc[MAX_N];         // right hand side c
    __shared__ float sy[MAX_N];         // current iterate, in column scaled units
    __shared__ float sgrad[MAX_N];      // gradient
    __shared__ float soldy[MAX_N];
    __shared__ float soldg[MAX_N];
    __shared__ float sydelta[MAX_N];
    __shared__ float sgdelta[MAX_N];
    __shared__ float srefy[MAX_N];      // reference iterate for the descent test
    __shared__ float st[MAX_N];         // scratch for G*y
    __shared__ float sred[TG_THREADS];  // reduction scratch

    const unsigned s      = blockIdx.x;
    const unsigned tid    = threadIdx.x;
    const unsigned tcount = blockDim.x;

    const float bnorm2 = BNORM2[s];
    const float* c_in  = C + (size_t)s * (size_t)n;

    // ---- initialize -------------------------------------------------------
    // The CPU solver starts from x = 0.5. In column scaled units that same point
    // is y = 0.5 * colnorm, which keeps the two solvers starting together.
    for (unsigned i = tid; i < (unsigned)n; i += tcount)
    {
        sc[i]    = c_in[i];
        sy[i]    = 0.5f * COLNORM[i];
        soldy[i] = 0.0f;
        soldg[i] = -c_in[i];
    }
    __syncthreads();

    for (unsigned i = tid; i < (unsigned)n; i += tcount)
    {
        float acc = 0.0f;
        for (int j = 0; j < n; j++)
        {
            acc = fmaf(G[i * (unsigned)n + j], sy[j], acc);
        }
        sgrad[i] = acc - sc[i];
        srefy[i] = sy[i];
    }
    __syncthreads();

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
        for (unsigned i = tid; i < (unsigned)n; i += tcount)
        {
            float g  = sgrad[i];
            float pg = (sy[i] > 0.0f) ? g : fminf(0.0f, g);
            local_max = fmaxf(local_max, fabsf(pg * COLNORM[i]));
        }
        npg = tg_reduce_max(sred, tid, tcount, local_max);
        if (npg < pgtol)
        {
            term = 2;
        }

        // Fixed set, and the deltas the BB step is built from. A variable sitting at
        // the lower bound with a gradient pushing it further down contributes nothing.
        for (unsigned i = tid; i < (unsigned)n; i += tcount)
        {
            bool fixed  = (sy[i] == 0.0f && sgrad[i] > 0.0f);
            sydelta[i]  = fixed ? 0.0f : (sy[i] - soldy[i]);
            sgdelta[i]  = fixed ? 0.0f : (sgrad[i] - soldg[i]);
            soldy[i]    = sy[i];
            soldg[i]    = sgrad[i];
        }
        __syncthreads();

        // Barzilai-Borwein step, alternating the two variants by iteration parity.
        float nr_local = 0.0f;
        float dr_local = 0.0f;
        if (iter % 2)
        {
            for (unsigned i = tid; i < (unsigned)n; i += tcount)
            {
                nr_local = fmaf(sydelta[i], sydelta[i], nr_local);
                dr_local = fmaf(sydelta[i], sgdelta[i], dr_local);
            }
        }
        else
        {
            for (unsigned i = tid; i < (unsigned)n; i += tcount)
            {
                nr_local = fmaf(sydelta[i], sgdelta[i], nr_local);
                dr_local = fmaf(sgdelta[i], sgdelta[i], dr_local);
            }
        }
        float nr = tg_reduce_sum(sred, tid, tcount, nr_local);
        float dr = tg_reduce_sum(sred, tid, tcount, dr_local);
        // Guard the denominator as well as the numerator. The CPU version only
        // checks nr, which lets a zero dr produce an infinite step; in float that
        // turns the whole iterate into NaN on the next update.
        float step = (nr == 0.0f || dr == 0.0f) ? 0.0f : beta * (nr / dr);

        for (unsigned i = tid; i < (unsigned)n; i += tcount)
        {
            float v = fmaf(-step, sgrad[i], sy[i]);
            sy[i] = (v < 0.0f) ? 0.0f : v;
        }
        __syncthreads();

        // Objective and gradient from the normal equations.
        for (unsigned i = tid; i < (unsigned)n; i += tcount)
        {
            float acc = 0.0f;
            for (int j = 0; j < n; j++)
            {
                acc = fmaf(G[i * (unsigned)n + j], sy[j], acc);
            }
            st[i] = acc;
        }
        __syncthreads();

        float ytgy_local = 0.0f;
        float ytc_local  = 0.0f;
        for (unsigned i = tid; i < (unsigned)n; i += tcount)
        {
            ytgy_local = fmaf(sy[i], st[i], ytgy_local);
            ytc_local  = fmaf(sy[i], sc[i], ytc_local);
            sgrad[i]   = st[i] - sc[i];
        }
        float ytgy = tg_reduce_sum(sred, tid, tcount, ytgy_local);
        float ytc  = tg_reduce_sum(sred, tid, tcount, ytc_local);
        obj = 0.5f * (ytgy - 2.0f * ytc + bnorm2);
        __syncthreads();

        // Non monotone descent test, every mcheck iterations. On failure the step
        // is damped; on success the reference iterate advances.
        if (iter % mcheck == 0)
        {
            float d_local = 0.0f;
            for (unsigned i = tid; i < (unsigned)n; i += tcount)
            {
                d_local = fmaf(sgrad[i], srefy[i] - sy[i], d_local);
            }
            float d = sigma * tg_reduce_sum(sred, tid, tcount, d_local);

            d = have_prev_check ? (obj_prev_check - obj - d) : (obj - d);

            if (d < 0.0f)
            {
                beta *= decay;
            }
            else
            {
                for (unsigned i = tid; i < (unsigned)n; i += tcount)
                {
                    srefy[i] = sy[i];
                }
                __syncthreads();
            }
            obj_prev_check  = obj;
            have_prev_check = true;
        }
    }

    // Undo the column scaling: x = y / colnorm.
    for (unsigned i = tid; i < (unsigned)n; i += tcount)
    {
        X[(size_t)s * (size_t)n + i] = sy[i] / COLNORM[i];
    }
    if (tid == 0)
    {
        ITERS[s] = iter;
        NPG[s]   = npg;
    }
}

// One block per spectrum. Computes the SNIP background over the full spectrum, then
// writes the background-subtracted, clamped fit-range slice to OUTB (the solver's b) and
// the fit-range background to OUT_BKG. This is a faithful float port of
// data_struct::snip_background (spectra.h): a boxcar(5) convolution, a log-log compress,
// n_first + n_decay clipping passes with the width schedule, then an exp-exp expand.
//
// The clipping recurrence updates the background in place while sweeping channels
// ascending (background[k] depends on an already-updated lo_index < k), so it is
// inherently sequential within a pass and is run by a single thread over one shared
// working buffer. The device parallelism is across the thousands of spectra in the batch
// (one block each), not across channels within a spectrum.
//
// Shared memory holds full_len floats. The host guarantees full_len <=
// CUDA_SNIP_MAX_CHANNELS before launching, and passes the size in as a dynamic shared
// allocation (full_len * sizeof(float)).
extern __shared__ float snip_smem[];

__global__ void snip_background_kernel(const float* __restrict__ FULL,
                                       const float* __restrict__ BASE_WIDTH,
                                       int full_len,
                                       int n_first_passes,
                                       int n_decay_passes,
                                       int lo_limit,
                                       int up_limit,
                                       int fit_min,
                                       int n_fit,
                                       float* __restrict__ OUTB,      // n_spectra x n_fit (solver b)
                                       float* __restrict__ OUT_BKG)   // n_spectra x n_fit, may be null
{
    float* cur = snip_smem;                 // working background, full_len

    const unsigned s      = blockIdx.x;
    const unsigned tid    = threadIdx.x;
    const unsigned tcount = blockDim.x;

    const float* spec = FULL + (size_t)s * (size_t)full_len;

    // ---- boxcar(5) convolution, matching convolve1d (spectra.h:310) -----------------
    // convolve1d centers a length-5 boxcar: output channel k is the mean of the 5 input
    // channels centered on k, and the first/last 2 channels are left zero (the valid
    // region is [2, full_len-3]). Reproduce that exactly so the edges agree with the CPU.
    for (unsigned k = tid; k < (unsigned)full_len; k += tcount)
    {
        float v = 0.0f;
        if (k >= 2u && k + 2u < (unsigned)full_len)
        {
            v = (spec[k - 2] + spec[k - 1] + spec[k] + spec[k + 1] + spec[k + 2]) * 0.2f;
        }
        // log(log(x+1)+1)
        cur[k] = logf(logf(v + 1.0f) + 1.0f);
    }
    __syncthreads();

    // ---- clipping passes ------------------------------------------------------------
    // Pass p uses width[k] * (1/sqrt2)^(decay index). The first n_first_passes run at the
    // base width (decay index 0); each subsequent pass divides the width by sqrt2, matching
    // the CPU's two fixed passes followed by the while loop that halves current_width.
    //
    // The CPU sweep (spectra.h:401-433) updates background[k] IN PLACE while iterating k
    // ascending, so background[k] reads a lo_index < k that was already lowered earlier in
    // the same pass. That within-pass dependency is sequential by construction, so the
    // sweep is done by a single thread per block (tid 0) over the shared buffer `cur`. The
    // parallel parts (convolution above, expand below) still use the whole block; only this
    // inherently serial recurrence is single threaded. One block still handles one
    // spectrum, so the device parallelism is across the thousands of spectra in the batch.
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
    __syncthreads();

    // ---- expand back: exp(exp(x)-1)-1, finite-guarded (spectra.h:435-436) ------------
    // Then subtract from the raw spectrum over the fit range, clamp >= 0, and emit.
    for (unsigned t = tid; t < (unsigned)n_fit; t += tcount)
    {
        const int k = fit_min + (int)t;
        float bg = expf(expf(cur[k]) - 1.0f) - 1.0f;
        if (!isfinite(bg)) bg = 0.0f;

        float sub = spec[k] - bg;
        if (sub < 0.0f) sub = 0.0f;

        OUTB[(size_t)s * (size_t)n_fit + t] = sub;
        if (OUT_BKG != nullptr)
        {
            OUT_BKG[(size_t)s * (size_t)n_fit + t] = bg;
        }
    }
}

// ----------------------------------------------------------------------------

struct Cuda_NNLS_Backend::Impl
{
    int          device = 0;
    cudaStream_t stream = nullptr;

    // Per fit matrix, rebuilt by set_fit_matrix.
    float* d_A       = nullptr;     // m x n, column major, columns unit norm
    float* d_G       = nullptr;     // n x n, row major
    float* d_colnorm = nullptr;     // n

    // Per batch, grown on demand and reused across calls.
    float* d_B       = nullptr;
    float* d_C       = nullptr;
    float* d_bnorm2  = nullptr;
    float* d_X       = nullptr;
    int*   d_iters   = nullptr;
    float* d_npg     = nullptr;

    // For the on-device SNIP path, grown on demand alongside the per batch buffers.
    float* d_full    = nullptr;     // n_spectra x full_len, raw spectra
    float* d_bkg     = nullptr;     // n_spectra x n_fit, returned background
    float* d_width   = nullptr;     // full_len, base width in channels (per fit matrix)

    int batch_capacity = 0;
    int full_capacity  = 0;         // n_spectra the d_full / d_bkg buffers hold
    int full_len_cap   = 0;         // full_len the d_full buffer was sized for
    int fit_len_cap    = 0;         // n_fit the d_bkg buffer was sized for

    int tg_threads = TG_THREADS;

    ~Impl()
    {
        free_fit_matrix();
        free_batch();
        free_snip();
        if (stream != nullptr)
        {
            cudaStreamDestroy(stream);
            stream = nullptr;
        }
    }

    void free_fit_matrix()
    {
        if (d_A != nullptr)       { cudaFree(d_A);       d_A = nullptr; }
        if (d_G != nullptr)       { cudaFree(d_G);       d_G = nullptr; }
        if (d_colnorm != nullptr) { cudaFree(d_colnorm); d_colnorm = nullptr; }
        if (d_width != nullptr)   { cudaFree(d_width);   d_width = nullptr; }
    }

    void free_snip()
    {
        if (d_full != nullptr) { cudaFree(d_full); d_full = nullptr; }
        if (d_bkg != nullptr)  { cudaFree(d_bkg);  d_bkg = nullptr; }
        full_capacity = 0;
        full_len_cap  = 0;
        fit_len_cap   = 0;
    }

    void free_batch()
    {
        if (d_B != nullptr)      { cudaFree(d_B);      d_B = nullptr; }
        if (d_C != nullptr)      { cudaFree(d_C);      d_C = nullptr; }
        if (d_bnorm2 != nullptr) { cudaFree(d_bnorm2); d_bnorm2 = nullptr; }
        if (d_X != nullptr)      { cudaFree(d_X);      d_X = nullptr; }
        if (d_iters != nullptr)  { cudaFree(d_iters);  d_iters = nullptr; }
        if (d_npg != nullptr)    { cudaFree(d_npg);    d_npg = nullptr; }
        batch_capacity = 0;
    }

    // Grow the per batch device buffers to hold n_spectra, reusing them across calls.
    bool ensure_batch_capacity(int n_spectra, int m, int n)
    {
        if (n_spectra <= batch_capacity && d_B != nullptr)
        {
            return true;
        }
        free_batch();
        bool ok = true;
        ok = ok && cudaMalloc((void**)&d_B,      sizeof(float) * (size_t)n_spectra * (size_t)m) == cudaSuccess;
        ok = ok && cudaMalloc((void**)&d_C,      sizeof(float) * (size_t)n_spectra * (size_t)n) == cudaSuccess;
        ok = ok && cudaMalloc((void**)&d_bnorm2, sizeof(float) * (size_t)n_spectra)             == cudaSuccess;
        ok = ok && cudaMalloc((void**)&d_X,      sizeof(float) * (size_t)n_spectra * (size_t)n) == cudaSuccess;
        ok = ok && cudaMalloc((void**)&d_iters,  sizeof(int)   * (size_t)n_spectra)             == cudaSuccess;
        ok = ok && cudaMalloc((void**)&d_npg,    sizeof(float) * (size_t)n_spectra)             == cudaSuccess;
        if (false == ok)
        {
            free_batch();
            return false;
        }
        batch_capacity = n_spectra;
        return true;
    }

    // Grow the SNIP staging buffers (raw full spectra in, fit-range background out) to
    // hold n_spectra at the given dimensions, reusing them across calls.
    bool ensure_snip_capacity(int n_spectra, int full_len, int n_fit)
    {
        if (n_spectra <= full_capacity && full_len <= full_len_cap && n_fit <= fit_len_cap
            && d_full != nullptr && d_bkg != nullptr)
        {
            return true;
        }
        free_snip();
        bool ok = true;
        ok = ok && cudaMalloc((void**)&d_full, sizeof(float) * (size_t)n_spectra * (size_t)full_len) == cudaSuccess;
        ok = ok && cudaMalloc((void**)&d_bkg,  sizeof(float) * (size_t)n_spectra * (size_t)n_fit)    == cudaSuccess;
        if (false == ok)
        {
            free_snip();
            return false;
        }
        full_capacity = n_spectra;
        full_len_cap  = full_len;
        fit_len_cap   = n_fit;
        return true;
    }
};

// ----------------------------------------------------------------------------

bool Cuda_NNLS_Backend::is_available()
{
    int count = 0;
    cudaError_t err = cudaGetDeviceCount(&count);
    // Clear any sticky error the probe may have set so it does not leak into later
    // calls. A machine with no driver or no GPU returns an error here, not a throw.
    cudaGetLastError();
    return (err == cudaSuccess) && (count > 0);
}

// ----------------------------------------------------------------------------

Cuda_NNLS_Backend::Cuda_NNLS_Backend() : _impl(nullptr), _m(0), _n(0)
{

}

// ----------------------------------------------------------------------------

Cuda_NNLS_Backend::~Cuda_NNLS_Backend()
{
    delete _impl;
    _impl = nullptr;
}

// ----------------------------------------------------------------------------

bool Cuda_NNLS_Backend::initialize()
{
    if (_impl != nullptr)
    {
        return true;
    }

    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count <= 0)
    {
        cudaGetLastError();
        _last_error = "No CUDA capable device found";
        return false;
    }

    Impl* impl = new Impl();
    impl->device = 0;

    if (cudaSetDevice(impl->device) != cudaSuccess)
    {
        _last_error = "Could not select CUDA device 0";
        delete impl;
        return false;
    }

    cudaDeviceProp prop;
    if (cudaGetDeviceProperties(&prop, impl->device) != cudaSuccess)
    {
        _last_error = "Could not query CUDA device properties";
        delete impl;
        return false;
    }
    _device_name = std::string(prop.name);

    if (cudaStreamCreate(&impl->stream) != cudaSuccess)
    {
        _last_error = "Could not create CUDA stream";
        delete impl;
        return false;
    }

    // The reductions need a power of two thread count the device will accept.
    // TG_THREADS is the target; halve it if the device reports a lower ceiling.
    int max_threads = prop.maxThreadsPerBlock;
    int threads = TG_THREADS;
    while (threads > max_threads && threads > 1)
    {
        threads >>= 1;
    }
    impl->tg_threads = threads;

    _impl = impl;
    return true;
}

// ----------------------------------------------------------------------------

bool Cuda_NNLS_Backend::set_fit_matrix(const double* const fitmatrix_col_major, int m, int n)
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
    if (n > CUDA_NNLS_MAX_COMPONENTS)
    {
        _last_error = "Too many components for the CUDA solver";
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

    // Build the host staging copies: column scaled A (float), the scaled Gram matrix
    // (accumulated in double, narrowed), and the column norms.
    std::vector<float> host_A((size_t)m * (size_t)n);
    std::vector<float> host_G((size_t)n * (size_t)n);
    std::vector<float> host_cn((size_t)n);

    for (int j = 0; j < n; j++)
    {
        host_cn[(size_t)j] = (float)colnorm[(size_t)j];
        const double* col = fitmatrix_col_major + (size_t)j * (size_t)m;
        float* dcol = host_A.data() + (size_t)j * (size_t)m;
        const double inv = 1.0 / colnorm[(size_t)j];
        for (int i = 0; i < m; i++)
        {
            dcol[i] = (float)(col[i] * inv);
        }
    }

    // G = Ascaled^T Ascaled, accumulated in double then narrowed. Symmetric, so only
    // the upper triangle is computed and then mirrored.
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
            host_G[(size_t)a * (size_t)n + b] = v;
            host_G[(size_t)b * (size_t)n + a] = v;
        }
    }

    // (Re)allocate the device buffers for this fit matrix and upload.
    _impl->free_fit_matrix();

    bool ok = true;
    ok = ok && cudaMalloc((void**)&_impl->d_A,       sizeof(float) * (size_t)m * (size_t)n) == cudaSuccess;
    ok = ok && cudaMalloc((void**)&_impl->d_G,       sizeof(float) * (size_t)n * (size_t)n) == cudaSuccess;
    ok = ok && cudaMalloc((void**)&_impl->d_colnorm, sizeof(float) * (size_t)n)             == cudaSuccess;
    if (false == ok)
    {
        _impl->free_fit_matrix();
        _last_error = "Could not allocate CUDA buffers for the fit matrix";
        return false;
    }

    ok = ok && cudaMemcpy(_impl->d_A, host_A.data(), sizeof(float) * (size_t)m * (size_t)n,
                          cudaMemcpyHostToDevice) == cudaSuccess;
    ok = ok && cudaMemcpy(_impl->d_G, host_G.data(), sizeof(float) * (size_t)n * (size_t)n,
                          cudaMemcpyHostToDevice) == cudaSuccess;
    ok = ok && cudaMemcpy(_impl->d_colnorm, host_cn.data(), sizeof(float) * (size_t)n,
                          cudaMemcpyHostToDevice) == cudaSuccess;
    if (false == ok)
    {
        _impl->free_fit_matrix();
        _last_error = "Could not upload the fit matrix to the CUDA device";
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

int Cuda_NNLS_Backend::recommended_batch_size() const
{
    // One block per spectrum, so the batch wants to be at least a few times the SM
    // count to keep every SM fed, while staying small enough that the staging buffers
    // remain a modest fraction of memory. Matches the Metal backend's default.
    return 4096;
}

// ----------------------------------------------------------------------------

bool Cuda_NNLS_Backend::solve_batch(const float* const spectra_batch,
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
    if (_impl->d_A == nullptr)
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

    if (false == _impl->ensure_batch_capacity(n_spectra, m, n))
    {
        _last_error = "Could not allocate CUDA batch buffers";
        return false;
    }

    cudaStream_t stream = _impl->stream;

    if (cudaMemcpyAsync(_impl->d_B, spectra_batch,
                        sizeof(float) * (size_t)n_spectra * (size_t)m,
                        cudaMemcpyHostToDevice, stream) != cudaSuccess)
    {
        _last_error = "Could not upload spectra batch to the CUDA device";
        return false;
    }

    // c = Ascaled^T b, 2D grid over (component, spectrum).
    {
        dim3 block((unsigned)((n < 8) ? n : 8), 8, 1);
        dim3 grid((unsigned)((n + (int)block.x - 1) / (int)block.x),
                  (unsigned)((n_spectra + (int)block.y - 1) / (int)block.y), 1);
        nnls_project_rhs<<<grid, block, 0, stream>>>(_impl->d_A, _impl->d_B, _impl->d_C,
                                                     m, n, n_spectra);
    }

    // ||b||^2, one thread per spectrum.
    {
        dim3 block(64, 1, 1);
        dim3 grid((unsigned)((n_spectra + 63) / 64), 1, 1);
        nnls_bnorm2<<<grid, block, 0, stream>>>(_impl->d_B, _impl->d_bnorm2, m, n_spectra);
    }

    // The solve, one block per spectrum.
    {
        dim3 block((unsigned)_impl->tg_threads, 1, 1);
        dim3 grid((unsigned)n_spectra, 1, 1);
        nnls_solve_batch<<<grid, block, 0, stream>>>(_impl->d_G, _impl->d_C, _impl->d_bnorm2,
                                                     _impl->d_colnorm, _impl->d_X, _impl->d_iters,
                                                     _impl->d_npg, n, max_iter, pgtol,
                                                     mcheck, sigma, decay);
    }

    cudaError_t launch_err = cudaGetLastError();
    if (launch_err != cudaSuccess)
    {
        _last_error = std::string("CUDA kernel launch failed: ") +
                      std::string(cudaGetErrorString(launch_err));
        return false;
    }

    // Copy the results back once the solve completes.
    cudaMemcpyAsync(out_x, _impl->d_X, sizeof(float) * (size_t)n_spectra * (size_t)n,
                    cudaMemcpyDeviceToHost, stream);
    if (out_num_iter != nullptr)
    {
        cudaMemcpyAsync(out_num_iter, _impl->d_iters, sizeof(int) * (size_t)n_spectra,
                        cudaMemcpyDeviceToHost, stream);
    }
    if (out_npg != nullptr)
    {
        cudaMemcpyAsync(out_npg, _impl->d_npg, sizeof(float) * (size_t)n_spectra,
                        cudaMemcpyDeviceToHost, stream);
    }

    cudaError_t err = cudaStreamSynchronize(stream);
    if (err != cudaSuccess)
    {
        _last_error = std::string("CUDA dispatch failed: ") +
                      std::string(cudaGetErrorString(err));
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

bool Cuda_NNLS_Backend::snip_supported(int full_len) const
{
    // The SNIP kernel keeps two full-length float buffers in dynamic shared memory.
    return (full_len > 0) && (full_len <= CUDA_SNIP_MAX_CHANNELS);
}

// ----------------------------------------------------------------------------

bool Cuda_NNLS_Backend::solve_batch_snip(const float* const full_spectra,
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
    if (_impl->d_A == nullptr)
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
        _last_error = "Full spectrum length exceeds the CUDA SNIP shared-memory limit";
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

    if (false == _impl->ensure_batch_capacity(n_spectra, m, n))
    {
        _last_error = "Could not allocate CUDA batch buffers";
        return false;
    }
    if (false == _impl->ensure_snip_capacity(n_spectra, full_len, n_fit))
    {
        _last_error = "Could not allocate CUDA SNIP buffers";
        return false;
    }
    if (_impl->d_width == nullptr)
    {
        if (cudaMalloc((void**)&_impl->d_width, sizeof(float) * (size_t)full_len) != cudaSuccess)
        {
            _impl->d_width = nullptr;
            _last_error = "Could not allocate CUDA SNIP width buffer";
            return false;
        }
    }

    cudaStream_t stream = _impl->stream;

    // Upload the raw full spectra and the per-channel base width (width is tiny and
    // constant across the run; re-uploading it per batch keeps the API stateless).
    if (cudaMemcpyAsync(_impl->d_full, full_spectra,
                        sizeof(float) * (size_t)n_spectra * (size_t)full_len,
                        cudaMemcpyHostToDevice, stream) != cudaSuccess)
    {
        _last_error = "Could not upload full spectra to the CUDA device";
        return false;
    }
    if (cudaMemcpyAsync(_impl->d_width, base_width_chan,
                        sizeof(float) * (size_t)full_len,
                        cudaMemcpyHostToDevice, stream) != cudaSuccess)
    {
        _last_error = "Could not upload SNIP width to the CUDA device";
        return false;
    }

    // SNIP background + subtract + clamp, one block per spectrum. The kernel writes the
    // fit-range solver input into d_B and (optionally) the fit-range background into d_bkg.
    {
        const int threads = _impl->tg_threads;
        const size_t shmem = sizeof(float) * (size_t)full_len;
        snip_background_kernel<<<(unsigned)n_spectra, (unsigned)threads, shmem, stream>>>(
            _impl->d_full, _impl->d_width, full_len, n_first_passes, n_decay_passes,
            lo_limit, up_limit, fit_min, n_fit,
            _impl->d_B, (out_background != nullptr) ? _impl->d_bkg : nullptr);
    }

    // c = Ascaled^T b, 2D grid over (component, spectrum).
    {
        dim3 block((unsigned)((n < 8) ? n : 8), 8, 1);
        dim3 grid((unsigned)((n + (int)block.x - 1) / (int)block.x),
                  (unsigned)((n_spectra + (int)block.y - 1) / (int)block.y), 1);
        nnls_project_rhs<<<grid, block, 0, stream>>>(_impl->d_A, _impl->d_B, _impl->d_C,
                                                     m, n, n_spectra);
    }

    // ||b||^2, one thread per spectrum.
    {
        dim3 block(64, 1, 1);
        dim3 grid((unsigned)((n_spectra + 63) / 64), 1, 1);
        nnls_bnorm2<<<grid, block, 0, stream>>>(_impl->d_B, _impl->d_bnorm2, m, n_spectra);
    }

    // The solve, one block per spectrum.
    {
        dim3 block((unsigned)_impl->tg_threads, 1, 1);
        dim3 grid((unsigned)n_spectra, 1, 1);
        nnls_solve_batch<<<grid, block, 0, stream>>>(_impl->d_G, _impl->d_C, _impl->d_bnorm2,
                                                     _impl->d_colnorm, _impl->d_X, _impl->d_iters,
                                                     _impl->d_npg, n, max_iter, pgtol,
                                                     mcheck, sigma, decay);
    }

    cudaError_t launch_err = cudaGetLastError();
    if (launch_err != cudaSuccess)
    {
        _last_error = std::string("CUDA kernel launch failed: ") +
                      std::string(cudaGetErrorString(launch_err));
        return false;
    }

    // Copy the results back once the solve completes.
    cudaMemcpyAsync(out_x, _impl->d_X, sizeof(float) * (size_t)n_spectra * (size_t)n,
                    cudaMemcpyDeviceToHost, stream);
    if (out_background != nullptr)
    {
        cudaMemcpyAsync(out_background, _impl->d_bkg,
                        sizeof(float) * (size_t)n_spectra * (size_t)n_fit,
                        cudaMemcpyDeviceToHost, stream);
    }
    if (out_num_iter != nullptr)
    {
        cudaMemcpyAsync(out_num_iter, _impl->d_iters, sizeof(int) * (size_t)n_spectra,
                        cudaMemcpyDeviceToHost, stream);
    }
    if (out_npg != nullptr)
    {
        cudaMemcpyAsync(out_npg, _impl->d_npg, sizeof(float) * (size_t)n_spectra,
                        cudaMemcpyDeviceToHost, stream);
    }

    cudaError_t err = cudaStreamSynchronize(stream);
    if (err != cudaSuccess)
    {
        _last_error = std::string("CUDA dispatch failed: ") +
                      std::string(cudaGetErrorString(err));
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

} //namespace gpu
} //namespace routines
} //namespace fitting
