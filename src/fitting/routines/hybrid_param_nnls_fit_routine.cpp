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

/// Initial Author <2022>: Arthur Glowacki



#include "hybrid_param_nnls_fit_routine.h"

#include <iostream>
#include <algorithm>
#include <math.h>

#include <string.h>

using namespace data_struct;

namespace fitting
{
namespace routines
{

// ----------------------------------------------------------------------------

template<typename T_real>
Hybrid_Param_NNLS_Fit_Routine<T_real>::Hybrid_Param_NNLS_Fit_Routine() : NNLS_Fit_Routine<T_real>()
{
    _model = nullptr;
    _elements_to_fit = nullptr;
    _spectra = nullptr;
    this->_max_iter = 4000;
}

// ----------------------------------------------------------------------------

template<typename T_real>
Hybrid_Param_NNLS_Fit_Routine<T_real>::~Hybrid_Param_NNLS_Fit_Routine()
{

}

// ----------------------------------------------------------------------------

template<typename T_real>
void Hybrid_Param_NNLS_Fit_Routine<T_real>::model_spectrum(const Fit_Parameters<T_real>* const fit_params,
                                                    const struct Range* const energy_range,
                                                    Spectra<T_real>* spectra_model)
{
    if (_model != nullptr && _elements_to_fit != nullptr)
    {
        _model->update_fit_params_values(fit_params);
        this->initialize(_model, _elements_to_fit, *energy_range);
        this->fit_spectrum_model(_spectra, &_background, _elements_to_fit, spectra_model);
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME Hybrid_Param_NNLS_Fit_Routine<T_real>::fit_spectra_parameters(const models::Base_Model<T_real>* const model,
                                                                    const Spectra<T_real>* const spectra,
                                                                    const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                                    Fit_Parameters<T_real>& out_fit_params,
                                                                    Callback_Func_Status_Def* status_callback)
{
    OPTIMIZER_OUTCOME ret_val = OPTIMIZER_OUTCOME::FAILED;

    Fit_Parameters<T_real> fit_params = model->fit_parameters();
    
    fit_params.add_parameter(Fit_Param<T_real>(STR_NUM_ITR, 0.0));
    fit_params.add_parameter(Fit_Param<T_real>(STR_RESIDUAL, 0.0));
    fit_params.add_parameter(Fit_Param<T_real>(STR_CHISQUARE, 0.0));
    fit_params.add_parameter(Fit_Param<T_real>(STR_CHISQRED, 0.0));
    fit_params.add_parameter(Fit_Param<T_real>(STR_FREE_PARS, 0.0));

    if (fit_params.contains(STR_COMPTON_AMPLITUDE))
    {
        fit_params[STR_COMPTON_AMPLITUDE].bound_type = E_Bound_Type::FIXED;
    }
    if (fit_params.contains(STR_COHERENT_SCT_AMPLITUDE))
    {
        fit_params[STR_COHERENT_SCT_AMPLITUDE].bound_type = E_Bound_Type::FIXED;
    }

    out_fit_params.append_and_update(fit_params);

    //If the sum of the spectra we are trying to fit to is zero then set out counts to -10.0 == log(0.0000000001)
    if(spectra->sum() == 0)
    {
        fit_params.set_all_value(-10.0, E_Bound_Type::FIT);
		ret_val = OPTIMIZER_OUTCOME::FOUND_ZERO;
    }
    else
    {
        if(this->_optimizer != nullptr)
        {
            //ret_val = _optimizer->minimize(&fit_params, spectra, elements_to_fit, model, _energy_range, status_callback);
            std::function<void(const Fit_Parameters<T_real>* const, const  Range* const, Spectra<T_real>*)> gen_func = std::bind(&Hybrid_Param_NNLS_Fit_Routine<T_real>::model_spectrum, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

            if (fit_params.contains(STR_SNIP_WIDTH))
            {
                ArrayTr<T_real> bkg = snip_background<T_real>(spectra,
                    fit_params.value(STR_ENERGY_OFFSET),
                    fit_params.value(STR_ENERGY_SLOPE),
                    fit_params.value(STR_ENERGY_QUADRATIC),
                    fit_params.value(STR_SNIP_WIDTH),
                    this->_energy_range.min,
                    this->_energy_range.max);

                _background = bkg.segment(this->_energy_range.min, this->_energy_range.count());
            }
            else
            {
                _background.setZero(this->_energy_range.count());
            }
            _model = (models::Base_Model<T_real>*)model;
            _elements_to_fit = elements_to_fit;
            _spectra = spectra;

            ret_val = this->_optimizer->minimize_func(&fit_params, spectra, this->_energy_range, &_background, gen_func);

            _model->update_fit_params_values(&fit_params);
            this->initialize(_model, elements_to_fit, this->_energy_range);
            std::unordered_map<std::string, T_real> out_counts;
            this->fit_spectra(_model, spectra, elements_to_fit, out_counts);

            out_fit_params.append_and_update(fit_params);
            for (const auto& itr : *elements_to_fit)
            {
                out_fit_params.add_parameter(Fit_Param<T_real>(itr.first, log10(out_counts[itr.first])));
            }
        }
    }
    return ret_val;
}

// ----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT Hybrid_Param_NNLS_Fit_Routine<float>;
TEMPLATE_CLASS_DLL_EXPORT Hybrid_Param_NNLS_Fit_Routine<double>;

} //namespace routines
} //namespace fitting
