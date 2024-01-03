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

/// Initial Author <2016>: Arthur Glowacki



#ifndef Gaussian_Model_H
#define Gaussian_Model_H

#include "fitting/models/base_model.h"
#include "data_struct/fit_parameters.h"

namespace fitting
{
namespace models
{

using namespace data_struct;

template<typename T_real>
class DLL_EXPORT Gaussian_Model: public Base_Model<T_real>
{
public:
    Gaussian_Model();

    virtual ~Gaussian_Model();

    virtual const Fit_Parameters<T_real>& fit_parameters() const { return _fit_parameters; }

    // single threaded
    virtual const Spectra<T_real> model_spectrum(const Fit_Parameters<T_real>* const fit_params,
                                                    const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                    std::unordered_map<std::string, ArrayTr<T_real>>* labeled_spectras,
                                                    const struct Range energy_range);

    const std::tuple<std::vector<std::string>, std::vector<ArrayTr<T_real>>> model_spectrum_info(const Fit_Parameters<T_real>* const fit_params,
                                                    const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                    std::unordered_map<std::string, ArrayTr<T_real>>* labeled_spectras,
                                                    const struct Range energy_range);

    // multi threaded
    virtual const Spectra<T_real> model_spectrum_mp(const Fit_Parameters<T_real>* const fit_params,
                                                        const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                        const struct Range energy_range);

    virtual const Spectra<T_real> model_spectrum_element(const Fit_Parameters<T_real>* const fitp,
                                                            const Fit_Element_Map<T_real>* const element_to_fit,
                                                            const ArrayTr<T_real> &ev,
                                                            std::unordered_map<std::string, ArrayTr<T_real>>* labeled_spectras);

    void set_fit_params_preset(Fit_Params_Preset lock_macro);

    /**
     * @brief gauss_peak :  models a gaussian fluorescence peak, see also van espen, spectrum evaluation,
                            in van grieken, handbook of x-ray spectrometry, 2nd ed, page 182 ff
                            gain / (sigma * sqrt( 2.0 * M_PI) ) * exp( -0.5 * ( (delta_energy / sigma) ** 2 )
     * @param gain
     * @param sigma
     * @param delta_energy
     * @return
     */
    virtual const ArrayTr<T_real> peak(T_real gain, T_real sigma, const ArrayTr<T_real>& delta_energy) const;

    /**
     * @brief gauss_step : gain / 2.0 /  peak_E * erfc(delta_energy/(M_SQRT2 * sigma));
     * @param gain
     * @param sigma
     * @param delta_energy
     * @param peak_E
     * @return
     */
    virtual const ArrayTr<T_real> step(T_real gain, T_real sigma, const ArrayTr<T_real>& delta_energy, T_real peak_E) const;

    virtual const ArrayTr<T_real> tail(T_real gain, T_real sigma, ArrayTr<T_real> delta_energy, T_real gamma) const;

    virtual const ArrayTr<T_real> elastic_peak(const Fit_Parameters<T_real>* const fitp, const ArrayTr<T_real>& ev, T_real gain) const;

    virtual const ArrayTr<T_real> compton_peak(const Fit_Parameters<T_real>* const fitp, const ArrayTr<T_real>& ev, T_real gain) const;

    virtual const ArrayTr<T_real> escape_peak(const ArrayTr<T_real>& spectra, const ArrayTr<T_real>& ev, T_real escape_factor) const;

    virtual void reset_to_default_fit_params() { _fit_parameters = _generate_default_fit_parameters(); }

    virtual void update_fit_params_values(const Fit_Parameters<T_real> *fit_params) { _fit_parameters.update_values(fit_params); }

    void update_and_add_fit_params_values_gt_zero(Fit_Parameters<T_real>* fit_params) { _fit_parameters.update_and_add_values_gt_zero(fit_params); }

protected:

    Fit_Parameters<T_real> _generate_default_fit_parameters();

    Fit_Parameters<T_real> _fit_parameters;

};


template<typename T_real>
DLL_EXPORT ArrayTr<T_real> generate_ev_array(Range energy_range, Fit_Parameters<T_real>& fit_params);

template<typename T_real>
DLL_EXPORT ArrayTr<T_real> generate_ev_array(Range energy_range, T_real energy_offset, T_real energy_slope, T_real energy_quad);

} //namespace models

} //namespace fitting

#endif // Gaussian_Model_H
