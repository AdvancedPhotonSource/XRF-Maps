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


#include "gpu_nnls_fit_routine.h"

#include <cmath>

#ifdef _BUILD_WITH_METAL
#include "fitting/routines/gpu/metal_nnls_backend.h"
#endif

#ifdef _BUILD_WITH_CUDA
#include "fitting/routines/gpu/cuda_nnls_backend.h"
#endif

// The maximum component count is a property of whichever backend is active. Both
// backends define the same value (64), kept in sync on purpose.
#if defined(_BUILD_WITH_METAL)
#define GPU_NNLS_MAX_COMPONENTS gpu::METAL_NNLS_MAX_COMPONENTS
#elif defined(_BUILD_WITH_CUDA)
#define GPU_NNLS_MAX_COMPONENTS gpu::CUDA_NNLS_MAX_COMPONENTS
#endif

namespace fitting
{
namespace routines
{

// ----------------------------------------------------------------------------

template<typename T_real>
GPU_NNLS_Fit_Routine<T_real>::GPU_NNLS_Fit_Routine() : NNLS_Fit_Routine<T_real>()
{
    _backend = nullptr;
    _gpu_ready = false;
    _max_batch_size = 4096;
}

// ----------------------------------------------------------------------------

template<typename T_real>
GPU_NNLS_Fit_Routine<T_real>::GPU_NNLS_Fit_Routine(size_t max_iter) : NNLS_Fit_Routine<T_real>(max_iter)
{
    _backend = nullptr;
    _gpu_ready = false;
    _max_batch_size = 4096;
}

// ----------------------------------------------------------------------------

template<typename T_real>
GPU_NNLS_Fit_Routine<T_real>::~GPU_NNLS_Fit_Routine()
{
    // Under _BUILD_WITH_GPU, _backend is a unique_ptr and frees itself.
    _backend = nullptr;
    _gpu_ready = false;
}

// ----------------------------------------------------------------------------

template<typename T_real>
std::string GPU_NNLS_Fit_Routine<T_real>::gpu_device_name() const
{
#ifdef _BUILD_WITH_GPU
    if (_backend != nullptr)
    {
        return _backend->device_name();
    }
#endif
    return "";
}

// ----------------------------------------------------------------------------

template<typename T_real>
int GPU_NNLS_Fit_Routine<T_real>::preferred_batch_size() const
{
#ifdef _BUILD_WITH_GPU
    if (_gpu_ready && _backend != nullptr)
    {
        int rec = _backend->recommended_batch_size();
        return (rec < _max_batch_size) ? rec : _max_batch_size;
    }
#endif
    return 1;
}

// ----------------------------------------------------------------------------

template<typename T_real>
bool GPU_NNLS_Fit_Routine<T_real>::_setup_gpu()
{
    _gpu_ready = false;

#ifdef _BUILD_WITH_GPU

    const int m = (int)this->_fitmatrix.rows();
    const int n = (int)this->_fitmatrix.cols();

    if (m <= 0 || n <= 0)
    {
        logW << "GPU NNLS: fit matrix is empty, using CPU solver\n";
        return false;
    }

    if (n > GPU_NNLS_MAX_COMPONENTS)
    {
        logW << "GPU NNLS: " << n << " components exceeds the device limit of "
             << GPU_NNLS_MAX_COMPONENTS << ", using CPU solver\n";
        return false;
    }

    if (false == GPU_NNLS_BACKEND_TYPE::is_available())
    {
        logW << "GPU NNLS: no GPU device available, using CPU solver\n";
        return false;
    }

    if (_backend == nullptr)
    {
        _backend = std::make_unique<GPU_NNLS_BACKEND_TYPE>();
    }

    if (false == _backend->initialize())
    {
        logW << "GPU NNLS: " << _backend->last_error() << ", using CPU solver\n";
        return false;
    }

    // Hand the backend the fit matrix in double, column major, whatever T_real is.
    // The scaling and the Gram matrix are accumulated in double on the way in; see
    // the active backend's set_fit_matrix.
    std::vector<double> host_a((size_t)m * (size_t)n);
    for (int j = 0; j < n; j++)
    {
        for (int i = 0; i < m; i++)
        {
            host_a[(size_t)j * (size_t)m + (size_t)i] = (double)this->_fitmatrix(i, j);
        }
    }

    if (false == _backend->set_fit_matrix(host_a.data(), m, n))
    {
        logW << "GPU NNLS: " << _backend->last_error() << ", using CPU solver\n";
        return false;
    }

    logI << "GPU NNLS: using " << _backend->device_name() << " for " << m
         << " channels x " << n << " components\n";

    _gpu_ready = true;

#else

    logI << "GPU NNLS: built without GPU support, using CPU solver\n";

#endif

    return _gpu_ready;
}

// ----------------------------------------------------------------------------

template<typename T_real>
void GPU_NNLS_Fit_Routine<T_real>::initialize(models::Base_Model<T_real>& model,
                                              const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                              const struct Range energy_range,
                                              ArrayTr<T_real>* custom_background)
{
    NNLS_Fit_Routine<T_real>::initialize(model, elements_to_fit, energy_range, custom_background);
    _setup_gpu();
}

// ----------------------------------------------------------------------------

template<typename T_real>
void GPU_NNLS_Fit_Routine<T_real>::initialize_mp(models::Base_Model<T_real>& model,
                                                 const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                 const struct Range energy_range)
{
    NNLS_Fit_Routine<T_real>::initialize_mp(model, elements_to_fit, energy_range);
    _setup_gpu();
}

// ----------------------------------------------------------------------------

template<typename T_real>
void GPU_NNLS_Fit_Routine<T_real>::fit_spectra_batch(const models::Base_Model<T_real>& model,
                                                     const std::vector<const Spectra<T_real>*>& spectra_batch,
                                                     const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                     std::vector<std::unordered_map<std::string, T_real>>& out_counts)
{
    const size_t n_in = spectra_batch.size();

    out_counts.clear();
    out_counts.resize(n_in);

    if (n_in == 0)
    {
        return;
    }

    if (false == _gpu_ready)
    {
        // Nothing to batch onto, so fall back to the inherited CPU solve. Identical
        // results to selecting the plain NNLS routine.
        for (size_t s = 0; s < n_in; s++)
        {
            if (spectra_batch[s] == nullptr || spectra_batch[s]->size() == 0)
            {
                continue;
            }
            NNLS_Fit_Routine<T_real>::fit_spectra(model, spectra_batch[s], elements_to_fit, out_counts[s]);
        }
        return;
    }

#ifdef _BUILD_WITH_GPU

    const int m = this->_energy_range.count();
    const int n = (int)this->_fitmatrix.cols();

    Fit_Parameters<T_real> fit_params = model.fit_parameters();

    // A custom background does not depend on the spectrum, so resolve it once.
    const bool has_custom_bkg = (this->_custom_background != nullptr);
    const bool has_snip = (false == has_custom_bkg) && fit_params.contains(STR_SNIP_WIDTH);

    ArrayTr<T_real> custom_bkg;
    if (has_custom_bkg)
    {
        custom_bkg = this->_custom_background->segment(this->_energy_range.min, this->_energy_range.count());
    }

    // The full spectrum length (SNIP runs on the whole spectrum, not just the fit range).
    // All spectra in a scan share a length; take it from the first usable one.
    size_t full_len = 0;
    for (size_t s = 0; s < n_in; s++)
    {
        if (spectra_batch[s] != nullptr && spectra_batch[s]->size() > 0)
        {
            full_len = (size_t)spectra_batch[s]->size();
            break;
        }
    }

    // The on-device SNIP path. Both backends implement it; it applies to the SNIP
    // background case when the full spectrum fits the kernel's threadgroup/shared-memory
    // working set. It moves the expensive per-spectrum SNIP off the single host thread
    // this routine runs on (the batch path has no thread pool), which is what makes the
    // GPU routine faster than the CPU routine rather than slower. Every other case keeps
    // the host loop below.
    bool use_gpu_snip = false;
#ifdef _BUILD_WITH_GPU
    use_gpu_snip = has_snip && full_len > 0 && _backend->snip_supported((int)full_len);
#endif

    // Rows are packed (null/empty spectra skipped), so keep a map back to caller indices.
    std::vector<size_t> row_to_src;
    row_to_src.reserve(n_in);

    ArrayTr<T_real> sum_background = ArrayTr<T_real>::Zero(m);

    std::vector<float> host_x;
    std::vector<int>   host_iters;
    std::vector<float> host_npg;

#ifdef _BUILD_WITH_GPU
    if (use_gpu_snip)
    {
        // Compute the spectrum-independent SNIP inputs once, matching data_struct::
        // snip_background (spectra.h) exactly so the device result tracks the CPU routine.
        const T_real e_off  = fit_params.value(STR_ENERGY_OFFSET);
        const T_real e_lin  = fit_params.value(STR_ENERGY_SLOPE);
        const T_real e_quad = fit_params.value(STR_ENERGY_QUADRATIC);
        const T_real width  = fit_params.value(STR_SNIP_WIDTH);

        ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced((Eigen::Index)full_len, 0, (T_real)(full_len - 1));
        energy = e_off + (energy * e_lin) + (Eigen::pow(energy, (T_real)2.0) * e_quad);
        ArrayTr<T_real> tmp = std::pow((e_off / (T_real)2.3548), (T_real)2.0) + energy * (T_real)3.58 * e_lin;
        tmp = tmp.unaryExpr([](T_real r) { return r < 0.0 ? (T_real)0.0 : r; });
        ArrayTr<T_real> current_width = (T_real)2.35 * Eigen::sqrt(tmp);
        current_width = width * current_width / e_lin;   // in channels

        const int lo_limit = std::max<int>((int)this->_energy_range.min, 0);
        const int up_limit = std::min<int>((int)this->_energy_range.max, (int)full_len - 1);

        // Fixed two "first snipping" passes, then the decaying while loop; the trip count
        // is data independent, so compute it here once (mirrors the CPU while condition).
        const int n_first_passes = 2;
        int n_decay_passes = 0;
        {
            T_real mx = current_width.maxCoeff();
            while (mx >= (T_real)0.5)
            {
                n_decay_passes++;
                mx = mx / (T_real)M_SQRT2;
            }
        }

        std::vector<float> base_width((size_t)full_len);
        for (size_t k = 0; k < full_len; k++)
        {
            base_width[k] = (float)current_width[(Eigen::Index)k];
        }

        // Pack the raw full spectra for every usable pixel.
        std::vector<float> host_full;
        host_full.reserve(n_in * full_len);
        for (size_t s = 0; s < n_in; s++)
        {
            const Spectra<T_real>* const spectra = spectra_batch[s];
            if (spectra == nullptr || spectra->size() == 0)
            {
                continue;
            }
            for (size_t k = 0; k < full_len; k++)
            {
                host_full.push_back((float)(*spectra)[(Eigen::Index)k]);
            }
            row_to_src.push_back(s);
        }

        const int n_rows = (int)row_to_src.size();
        if (n_rows == 0)
        {
            return;
        }

        host_x.assign((size_t)n_rows * (size_t)n, 0.f);
        host_iters.assign((size_t)n_rows, 0);
        host_npg.assign((size_t)n_rows, 0.f);
        std::vector<float> host_bkg((size_t)n_rows * (size_t)m, 0.f);

        if (false == _backend->solve_batch_snip(host_full.data(),
                                                n_rows,
                                                (int)full_len,
                                                base_width.data(),
                                                n_first_passes,
                                                n_decay_passes,
                                                lo_limit,
                                                up_limit,
                                                (int)this->_energy_range.min,
                                                m,
                                                host_x.data(),
                                                host_bkg.data(),
                                                host_iters.data(),
                                                host_npg.data(),
                                                (int)this->_max_iter,
                                                (float)1e-3))
        {
            logE << "GPU NNLS: " << _backend->last_error() << ", falling back to CPU for this batch\n";
            for (size_t s = 0; s < n_in; s++)
            {
                if (spectra_batch[s] == nullptr || spectra_batch[s]->size() == 0)
                {
                    continue;
                }
                NNLS_Fit_Routine<T_real>::fit_spectra(model, spectra_batch[s], elements_to_fit, out_counts[s]);
            }
            return;
        }

        // Accumulate the integrated background from the per-row backgrounds the device
        // returned over the fit range.
        for (int r = 0; r < n_rows; r++)
        {
            const float* const bkr = host_bkg.data() + (size_t)r * (size_t)m;
            for (int i = 0; i < m; i++)
            {
                sum_background[i] += (T_real)bkr[i];
            }
        }
    }
    else
#endif
    {
        // Host preprocessing path: subtract the background per spectrum on the CPU, then
        // run the fit-range-only device solve. Used for the custom-background and no-SNIP
        // cases, and for spectra too long for the SNIP kernel.
        std::vector<float> host_b;
        host_b.reserve(n_in * (size_t)m);

        for (size_t s = 0; s < n_in; s++)
        {
            const Spectra<T_real>* const spectra = spectra_batch[s];
            if (spectra == nullptr || spectra->size() == 0)
            {
                continue;
            }

            ArrayTr<T_real> background;
            if (has_custom_bkg)
            {
                background = custom_bkg;
            }
            else if (has_snip)
            {
                ArrayTr<T_real> bkg = snip_background<T_real>(spectra,
                                                              fit_params.value(STR_ENERGY_OFFSET),
                                                              fit_params.value(STR_ENERGY_SLOPE),
                                                              fit_params.value(STR_ENERGY_QUADRATIC),
                                                              fit_params.value(STR_SNIP_WIDTH),
                                                              this->_energy_range.min,
                                                              this->_energy_range.max);
                background = bkg.segment(this->_energy_range.min, this->_energy_range.count());
            }
            else
            {
                background.setZero(m);
            }

            ArrayTr<T_real> sub = spectra->segment(this->_energy_range.min, this->_energy_range.count());
            sub -= background;
            sub = sub.unaryExpr([](T_real v) { return v > 0.0 ? v : (T_real)0.0; });

            sum_background += background;

            for (int i = 0; i < m; i++)
            {
                host_b.push_back((float)sub[i]);
            }
            row_to_src.push_back(s);
        }

        const int n_rows = (int)row_to_src.size();
        if (n_rows == 0)
        {
            return;
        }

        host_x.assign((size_t)n_rows * (size_t)n, 0.f);
        host_iters.assign((size_t)n_rows, 0);
        host_npg.assign((size_t)n_rows, 0.f);

        if (false == _backend->solve_batch(host_b.data(),
                                           n_rows,
                                           host_x.data(),
                                           host_iters.data(),
                                           host_npg.data(),
                                           (int)this->_max_iter,
                                           (float)1e-3))
        {
            logE << "GPU NNLS: " << _backend->last_error() << ", falling back to CPU for this batch\n";
            for (size_t s = 0; s < n_in; s++)
            {
                if (spectra_batch[s] == nullptr || spectra_batch[s]->size() == 0)
                {
                    continue;
                }
                NNLS_Fit_Routine<T_real>::fit_spectra(model, spectra_batch[s], elements_to_fit, out_counts[s]);
            }
            return;
        }
    }

    const int n_rows = (int)row_to_src.size();

    // The integrated fitted spectrum is a sum over the batch of (background + A x),
    // and escape_peak is linear in the spectrum it is given, so the whole batch needs
    // one matrix-vector product and one escape peak rather than a reconstruction per
    // spectrum. Summing the solution vectors first is what makes that possible.
    Eigen::Matrix<T_real, Eigen::Dynamic, 1> sum_x = Eigen::Matrix<T_real, Eigen::Dynamic, 1>::Zero(n);

    for (int r = 0; r < n_rows; r++)
    {
        const size_t src = row_to_src[(size_t)r];
        std::unordered_map<std::string, T_real>& counts = out_counts[src];
        const float* const xr = host_x.data() + (size_t)r * (size_t)n;

        for (const auto& itr : this->_element_models)
        {
            const int idx = this->_element_row_index[itr.first];
            const T_real v = (T_real)xr[idx];
            if (std::isfinite(v))
            {
                counts[itr.first] = v;
                sum_x(idx) += v;
            }
            else
            {
                // Same contract as the CPU routine: a non finite amplitude is reported
                // as zero and contributes nothing to the fitted model.
                counts[itr.first] = (T_real)0.0;
            }
        }

        counts[STR_NUM_ITR] = (T_real)host_iters[(size_t)r];
        counts[STR_RESIDUAL] = (T_real)host_npg[(size_t)r];
    }

    Spectra<T_real> batch_model(m);
    batch_model.setZero(m);
    {
        Eigen::Matrix<T_real, Eigen::Dynamic, 1> fitted = this->_fitmatrix * sum_x;
        for (int i = 0; i < m; i++)
        {
            batch_model[i] = sum_background[i] + fitted(i);
        }
    }

    if (fit_params.contains(STR_SI_ESCAPE)
        && fit_params.at(STR_SI_ESCAPE).value > 0.0)
    {
        ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(this->_energy_range.count(),
                                                            this->_energy_range.min,
                                                            this->_energy_range.max);
        ArrayTr<T_real> ev = fit_params.value(STR_ENERGY_OFFSET)
                             + (energy * fit_params.value(STR_ENERGY_SLOPE))
                             + (pow(energy, (T_real)2.0) * fit_params.value(STR_ENERGY_QUADRATIC));

        batch_model += model.escape_peak(batch_model,
                                          ev,
                                          fit_params.value(STR_SI_ESCAPE),
                                          fit_params.value_or(STR_ESCAPE_LINEAR, (T_real)0.0));
    }

    {
        std::lock_guard<std::mutex> lock(this->_int_spec_mutex);

        if (false == this->_integrated_fitted_spectra.add(batch_model))
        {
            logW << "Could not add model spectra to _integrated_fitted_spectra\n";
        }
        if (false == this->_integrated_background.add(sum_background))
        {
            logW << "Could not add background spectra to _integrated_background\n";
        }
    }

#endif

}

// ----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT GPU_NNLS_Fit_Routine<float>;
TEMPLATE_CLASS_DLL_EXPORT GPU_NNLS_Fit_Routine<double>;

} //namespace routines
} //namespace fitting
