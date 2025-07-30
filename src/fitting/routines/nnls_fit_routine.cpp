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


#include "nnls_fit_routine.h"

//debug
#include <iostream>

namespace fitting
{
namespace routines
{

template<typename T_real>
std::mutex Matrix_Optimized_Fit_Routine<T_real>::_int_spec_mutex;

template<typename T_real>
NNLS_Fit_Routine<T_real>::NNLS_Fit_Routine() : Matrix_Optimized_Fit_Routine<T_real>()
{

    _max_iter = 200;

}

// ----------------------------------------------------------------------------

template<typename T_real>
NNLS_Fit_Routine<T_real>::NNLS_Fit_Routine(size_t max_iter) : Matrix_Optimized_Fit_Routine<T_real>()
{

    _max_iter = max_iter;

}

// ----------------------------------------------------------------------------

template<typename T_real>
NNLS_Fit_Routine<T_real>::~NNLS_Fit_Routine()
{
    _element_row_index.clear();
}

// ----------------------------------------------------------------------------

template<typename T_real>
void NNLS_Fit_Routine<T_real>::_generate_fitmatrix()
{

    _element_row_index.clear();
    _fitmatrix.resize(this->_energy_range.count(), this->_element_models.size());

    int i = 0;
    for(const auto& itr : this->_element_models)
    {
        //Spectra element_model = itr.second;
        for (int j=0; j<itr.second.size(); j++)
        {
            _fitmatrix(j,i) = itr.second[j];
        }
        //save element index for later
        _element_row_index[itr.first] = i;
        i++;
    }

}

// ----------------------------------------------------------------------------

template<typename T_real>
void NNLS_Fit_Routine<T_real>::fit_spectrum_model(const Spectra<T_real>* const spectra,
                                          const ArrayTr<T_real>* const background,
                                          const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                          Spectra<T_real>* spectra_model)
{
    //spectra_model->setZero(this->_energy_range.count());
    spectra_model->setZero();

    data_struct::ArrayTr<T_real>* result;
    int num_iter;
    T_real npg;

    ArrayTr<T_real> spectra_sub_background = spectra->segment(this->_energy_range.min, this->_energy_range.count());
    spectra_sub_background -= *background;
    spectra_sub_background = spectra_sub_background.unaryExpr([](T_real v) { return v > 0.0 ? v : (T_real)0.0; });
    nsNNLS::nnls<T_real> solver(&_fitmatrix, &spectra_sub_background, _max_iter);

    solver.optimize(num_iter, npg);
    //logI << "NNLS num iter: " << num_iter << " : npg : " << npg << "\n";
    if (num_iter < 0)
    {
        logE << "Num iter < 0" << "\n";
    }

    result = solver.getSolution();

    for (const auto& itr : *elements_to_fit)
    {
        if (std::isfinite((*result)[_element_row_index[itr.first]]))
        {
            for (int j = 0; j < this->_energy_range.count(); j++)
            {
                T_real val = _fitmatrix(j, _element_row_index[itr.first]) * (*result)[_element_row_index[itr.first]];
                if (std::isfinite(val))
                {
                    (*spectra_model)[j] += val;
                }
            }
        }
        else
        {
            for (int j = 0; j < this->_energy_range.count(); j++)
            {
                (*spectra_model)[j] = std::numeric_limits<T_real>::quiet_NaN();
            }
            break;
        }
    }

}

// ----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME NNLS_Fit_Routine<T_real>::fit_spectra(const models::Base_Model<T_real>* const model,
                                                const Spectra<T_real>* const spectra,
                                                const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                std::unordered_map<std::string, T_real>& out_counts)
{
	data_struct::ArrayTr<T_real>* result;
    int num_iter;
    T_real npg;
    Fit_Parameters<T_real> fit_params = model->fit_parameters();
    fit_params.add_parameter(Fit_Param<T_real>(STR_RESIDUAL, 0.0));
    ArrayTr<T_real> background;

    if(this->_custom_background == nullptr)
    {
        if (fit_params.contains(STR_SNIP_WIDTH))
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
            background.setZero(this->_energy_range.count());
        }
    }
    else
    {
        background = this->_custom_background->segment(this->_energy_range.min, this->_energy_range.count());
    }

    ArrayTr<T_real> spectra_sub_background = spectra->segment(this->_energy_range.min, this->_energy_range.count());
    spectra_sub_background -= background;
    spectra_sub_background = spectra_sub_background.unaryExpr([](T_real v) { return v>0.0 ? v : (T_real)0.0; });
    nsNNLS::nnls<T_real> solver(&_fitmatrix, &spectra_sub_background, _max_iter);

    Spectra<T_real> spectra_model = background;

    solver.optimize(num_iter, npg);
    if (num_iter < 0)
    {
        logE<<"num_iter < 0"<<"\n";
    }

    result = solver.getSolution();

    for(const auto& itr : this->_element_models)
    {
        if (std::isfinite((*result)[_element_row_index[itr.first]]))
        {
            out_counts[itr.first] = (*result)[_element_row_index[itr.first]];

            for (int j = 0; j < this->_energy_range.count(); j++)
            {
                T_real val = _fitmatrix(j, _element_row_index[itr.first]) * (*result)[_element_row_index[itr.first]];
                if (std::isfinite(val))
                {
                    spectra_model[j] += val;
                }
            }
        }
        else
        {
            out_counts[itr.first] = 0.;
        }
    }
    // add escape peak
    if (fit_params.at(STR_SI_ESCAPE).value > 0.0 && model != nullptr)
    {
        ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(this->_energy_range.count(), this->_energy_range.min, this->_energy_range.max);
        ArrayTr<T_real> ev = fit_params.value(STR_ENERGY_OFFSET) + (energy * fit_params.value(STR_ENERGY_SLOPE)) + (pow(energy, (T_real)2.0) * fit_params.value(STR_ENERGY_QUADRATIC));

        spectra_model += model->escape_peak(spectra_model, ev, fit_params.value(STR_SI_ESCAPE));
    }


    out_counts[STR_NUM_ITR] = static_cast<T_real>(num_iter);
    out_counts[STR_RESIDUAL] = npg;

	//lock and integrate results
	{
		std::lock_guard<std::mutex> lock(this->_int_spec_mutex);
		this->_integrated_fitted_spectra.add(spectra_model);
        this->_integrated_background.add(background);
	}

    if (num_iter == solver.getMaxit())
    {
        return OPTIMIZER_OUTCOME::EXHAUSTED;
    }
    return OPTIMIZER_OUTCOME::CONVERGED;

}

// ----------------------------------------------------------------------------

template<typename T_real>
void NNLS_Fit_Routine<T_real>::initialize(models::Base_Model<T_real>* const model,
                                  const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                  const struct Range energy_range,
                                  ArrayTr<T_real>* custom_background)
{
    Matrix_Optimized_Fit_Routine<T_real>::initialize(model, elements_to_fit, energy_range, custom_background);
    _generate_fitmatrix();
}

// ----------------------------------------------------------------------------

template<typename T_real>
void NNLS_Fit_Routine<T_real>::initialize_mp(models::Base_Model<T_real>* const model,
    const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
    const struct Range energy_range)
{
    Matrix_Optimized_Fit_Routine<T_real>::initialize_mp(model, elements_to_fit, energy_range);
    _generate_fitmatrix();
}
// ----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT NNLS_Fit_Routine<float>;
TEMPLATE_CLASS_DLL_EXPORT NNLS_Fit_Routine<double>;

} //namespace routines
} //namespace fitting
