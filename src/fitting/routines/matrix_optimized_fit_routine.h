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



#ifndef Matrix_Optimized_Fit_Routine_H
#define Matrix_Optimized_Fit_Routine_H

#include <mutex>

#include "fitting/routines/param_optimized_fit_routine.h"
#include "data_struct/fit_parameters.h"

namespace fitting
{
namespace routines
{

using namespace data_struct;

/**
 * @brief The Matrix_Optimized_Fit_Routine class : Matrix fit model
 */
template<typename T_real>
class DLL_EXPORT Matrix_Optimized_Fit_Routine: public Param_Optimized_Fit_Routine<T_real>
{
public:
    Matrix_Optimized_Fit_Routine();

    virtual ~Matrix_Optimized_Fit_Routine();

    virtual optimizers::OPTIMIZER_OUTCOME fit_spectra(const models::Base_Model<T_real>* const model,
                                          const Spectra<T_real>* const spectra,
                                          const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                          std::unordered_map<std::string, T_real>& out_counts);

    virtual std::string get_name() { return STR_FIT_GAUSS_MATRIX; }

    virtual void initialize(models::Base_Model<T_real>* const model,
                            const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                            const struct Range energy_range);

    void initialize_mp(models::Base_Model<T_real>* const model,
                        const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                        const struct Range energy_range);


    virtual void model_spectrum(const Base_Model<T_real>* const model,
                                const Fit_Parameters<T_real>* const fit_params,
                                const struct Range * const energy_range,
                                Spectra<T_real>* spectra_model);

    const Spectra<T_real>& fitted_integrated_spectra() {return _integrated_fitted_spectra;}

    const Spectra<T_real>& fitted_integrated_background() { return _integrated_background; }

	const Spectra<T_real>& max_integrated_spectra() { return _max_channels_spectra; }

	const Spectra<T_real>& max_10_integrated_spectra() { return _max_10_channels_spectra; }

    void set_use_weights(bool val) {_use_weights = val;}

protected:

    std::unordered_map<std::string, Spectra<T_real>> _generate_element_models(models::Base_Model<T_real>* const model,
                                                            const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                            struct Range energy_range);

    std::unordered_map<std::string, Spectra<T_real>> _generate_element_models_mp(models::Base_Model<T_real>* const model,
        const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
        struct Range energy_range);


	data_struct::Spectra<T_real> _integrated_fitted_spectra;
    data_struct::Spectra<T_real> _integrated_background;
	data_struct::Spectra<T_real> _max_channels_spectra;
	data_struct::Spectra<T_real> _max_10_channels_spectra;

    bool _use_weights;

    std::unordered_map<std::string, Spectra<T_real>> _element_models;

    static std::mutex _int_spec_mutex;

};

} //namespace routines

} //namespace fitting

#endif // Matrix_Optimized_Fit_Routine_H
