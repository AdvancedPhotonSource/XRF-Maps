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


#ifndef Cuda_NNLS_Backend_H
#define Cuda_NNLS_Backend_H

#include <string>

namespace fitting
{
namespace routines
{
namespace gpu
{

/**
 * @brief Largest number of fitted components the batched CUDA solver accepts.
 *
 * The solver keeps the whole per-spectrum working set (roughly eight n-vectors plus
 * a reduction scratch) in shared memory so the iteration never touches global memory.
 * At n = 64 that is a few KB, well inside the 48 KB of shared memory every CUDA
 * architecture guarantees per block. Problems wider than this fall back to the CPU
 * solver rather than silently spilling. Kept equal to METAL_NNLS_MAX_COMPONENTS so
 * the two GPU backends accept exactly the same problems.
 */
const int CUDA_NNLS_MAX_COMPONENTS = 64;

/**
 * @brief Largest full-spectrum length the on-device SNIP kernel accepts.
 *
 * The SNIP kernel keeps one full-length float working buffer in shared memory. At this
 * cap that is 4096 * 4 B = 16 KB, well inside the 48 KB per-block shared memory every
 * supported architecture guarantees. Spectra longer than this fall back to the host
 * SNIP path.
 */
const int CUDA_SNIP_MAX_CHANNELS = 4096;

/**
 * @brief Batched non negative least squares solver running on CUDA.
 *
 * This is the NVIDIA counterpart to Metal_NNLS_Backend, with an identical interface
 * and identical numerics. All device arithmetic is single precision to match the
 * Metal backend; callers templated on double must narrow on the way in and widen on
 * the way out, and should expect results to agree with the CPU solver only to single
 * precision. See GPU_NNLS_Fit_Routine for how that is handled.
 *
 * The solver works on the normal equations. Given the fit matrix A (m x n) it
 * precomputes the Gram matrix G = A^T A once in set_fit_matrix, so each iteration
 * costs O(n^2) in shared memory instead of the O(m n) the CPU solver pays walking the
 * full spectrum. Because n (element count, tens) is far smaller than m (channel count,
 * thousands), that alone is a large win; the GPU then solves a whole batch of spectra
 * at once, one thread block per spectrum.
 *
 * Forming the normal equations squares the condition number, which matters in single
 * precision. To hold that in check the columns of A are scaled to unit norm before G
 * is formed, and the solution is unscaled on the way out. The scaling is diagonal and
 * positive, so it preserves the non negativity constraint exactly.
 *
 * This class is not thread safe. Each instance owns one CUDA stream; use one instance
 * per thread, or serialise calls.
 */
class Cuda_NNLS_Backend
{
public:

    /**
     * @brief is_available : True when a usable CUDA device is present. Safe to call on
     *                       a machine with no GPU or no driver, and on a build without
     *                       CUDA.
     */
    static bool is_available();

    Cuda_NNLS_Backend();

    ~Cuda_NNLS_Backend();

    Cuda_NNLS_Backend(const Cuda_NNLS_Backend&) = delete;
    Cuda_NNLS_Backend& operator=(const Cuda_NNLS_Backend&) = delete;

    /**
     * @brief initialize : Select the device and create the stream.
     * @return false if CUDA is unavailable, in which case last_error() explains why and
     *         the caller should use the CPU path.
     */
    bool initialize();

    /**
     * @brief set_fit_matrix : Upload the fit matrix and precompute the Gram matrix.
     *                         Call once per energy range / element list, not per spectrum.
     * Takes double regardless of the caller's precision: the column scaling and the
     * Gram matrix are accumulated in double on the CPU and only then narrowed to float
     * for the device. This is the one place where the extra precision is cheap (it
     * happens once per energy range) and where it buys the most, since forming A^T A is
     * what squares the condition number.
     *
     * @param fitmatrix_col_major : A, m x n, column major (Eigen's default storage).
     * @param m : rows, the number of energy channels being fit.
     * @param n : columns, the number of components being fit. Must be <= CUDA_NNLS_MAX_COMPONENTS.
     */
    bool set_fit_matrix(const double* const fitmatrix_col_major, int m, int n);

    /**
     * @brief solve_batch : Solve min ||A x - b|| subject to x >= 0 for a batch of spectra.
     * @param spectra_batch : n_spectra x m, row major. Each row is one background
     *                        subtracted spectrum, already clamped to >= 0 by the caller.
     * @param n_spectra : number of spectra in this batch.
     * @param out_x : n_spectra x n, row major, filled with the solutions. Caller allocates.
     * @param out_num_iter : per spectrum iteration count. Caller allocates, may be null.
     * @param out_npg : per spectrum inf-norm of the projected gradient at exit, the same
     *                  residual measure the CPU solver reports. Caller allocates, may be null.
     * @param max_iter : iteration cap, mirroring nsNNLS::nnls maxit.
     * @param pgtol : projected gradient tolerance, mirroring nsNNLS::nnls pgtol.
     * @return false on dispatch failure; last_error() explains.
     */
    bool solve_batch(const float* const spectra_batch,
                     int n_spectra,
                     float* out_x,
                     int* out_num_iter,
                     float* out_npg,
                     int max_iter,
                     float pgtol);

    /**
     * @brief snip_supported : True when solve_batch_snip can run for the given full
     *                         spectrum length, i.e. the per-block SNIP working set fits
     *                         in shared memory. When false the caller must subtract the
     *                         background on the host and use solve_batch instead.
     */
    bool snip_supported(int full_len) const;

    /**
     * @brief solve_batch_snip : The full GPU pipeline for the SNIP-background case.
     *                           Computes the SNIP background, subtracts it over the fit
     *                           range, clamps to >= 0, and solves the NNLS problem, all on
     *                           the device. This exists to keep the expensive SNIP pass off
     *                           the single host thread the batched routine runs on.
     *
     * The spectrum-independent SNIP inputs (the per-channel base width in channels, the
     * limits, and the pass counts) are computed once on the host from the energy
     * calibration and uploaded, exactly matching data_struct::snip_background so the
     * device result tracks the CPU routine to single precision.
     *
     * @param full_spectra   : n_spectra x full_len, row major, the raw spectra (not yet
     *                          background subtracted). SNIP runs on the whole spectrum.
     * @param n_spectra      : number of spectra in this batch.
     * @param full_len       : channels per spectrum (spectra->size()).
     * @param base_width_chan: length full_len, current_width in channels after the
     *                          width * fwhm / energy_linear conversion (spectra.h). The
     *                          kernel scales this by (1/sqrt2)^pass for the decay passes.
     * @param n_first_passes : fixed "first snipping" passes at full width (== 2).
     * @param n_decay_passes : while-loop passes, host computed from the width schedule.
     * @param lo_limit       : clamp lower bound for the snip indices (== max(xmin,0)).
     * @param up_limit       : clamp upper bound (== min(xmax, full_len-1)).
     * @param fit_min        : index into the full spectrum where the fit range starts
     *                          (energy_range.min); the device slices [fit_min, fit_min+n_fit).
     * @param n_fit          : fit-range channel count (== m used by set_fit_matrix).
     * @param out_x          : n_spectra x n, row major, the solutions. Caller allocates.
     * @param out_background : n_spectra x n_fit, row major, the SNIP background over the
     *                          fit range, so the host can accumulate the integrated
     *                          background. Caller allocates, may be null.
     * @param out_num_iter   : per spectrum iteration count. Caller allocates, may be null.
     * @param out_npg        : per spectrum projected gradient inf-norm. May be null.
     * @param max_iter       : iteration cap, mirroring nsNNLS::nnls maxit.
     * @param pgtol          : projected gradient tolerance.
     * @return false on dispatch failure; last_error() explains.
     */
    bool solve_batch_snip(const float* const full_spectra,
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
                          float pgtol);

    /**
     * @brief recommended_batch_size : Spectra per dispatch that keeps the device busy
     *                                 without oversizing the staging buffers.
     */
    int recommended_batch_size() const;

    const std::string& last_error() const { return _last_error; }

    const std::string& device_name() const { return _device_name; }

private:

    struct Impl;

    Impl* _impl;

    std::string _last_error;

    std::string _device_name;

    int _m;

    int _n;

};

} //namespace gpu

} //namespace routines

} //namespace fitting

#endif // Cuda_NNLS_Backend_H
