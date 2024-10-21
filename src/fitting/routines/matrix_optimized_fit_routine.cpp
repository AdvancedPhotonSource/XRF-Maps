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


#include "matrix_optimized_fit_routine.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace fitting
{
namespace routines
{

template<typename T_real>
std::mutex Matrix_Optimized_Fit_Routine<T_real>::_int_spec_mutex;

// ----------------------------------------------------------------------------

template<typename T_real>
Matrix_Optimized_Fit_Routine<T_real>::Matrix_Optimized_Fit_Routine() : Param_Optimized_Fit_Routine<T_real>()
{
    _use_weights = true;
}

// ----------------------------------------------------------------------------

template<typename T_real>
Matrix_Optimized_Fit_Routine<T_real>::~Matrix_Optimized_Fit_Routine()
{

    //logD<<"******** destroy element models *******"<<"\n";
    _element_models.clear();

}

// --------------------------------------------------------------------------------------------------------------------

template<typename T_real>
void Matrix_Optimized_Fit_Routine<T_real>::model_spectrum(const Base_Model<T_real>* const model,
                                                        const Fit_Parameters<T_real>* const fit_params,
                                                        const struct Range * const energy_range,
                                                        Spectra<T_real>* spectra_model)
{
	spectra_model->setZero();

    for(const auto& itr : _element_models)
    {
        if(fit_params->contains(itr.first))
        {
            Fit_Param<T_real> param = fit_params->at(itr.first);
            (*spectra_model) += (pow((T_real)10.0, param.value) * itr.second);
        }
    }
    if (model != nullptr && fit_params->value(STR_SI_ESCAPE) > 0.0)
    {
        ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(energy_range->count(), energy_range->min, energy_range->max);
        ArrayTr<T_real> ev = fit_params->value(STR_ENERGY_OFFSET) + (energy * fit_params->value(STR_ENERGY_SLOPE)) + (pow(energy, (T_real)2.0) * fit_params->value(STR_ENERGY_QUADRATIC));

        (*spectra_model) += model->escape_peak((*spectra_model), ev, fit_params->value(STR_SI_ESCAPE));
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
std::unordered_map<std::string, Spectra<T_real>> Matrix_Optimized_Fit_Routine<T_real>::_generate_element_models(models::Base_Model<T_real>* const model,
    const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
    struct Range energy_range)
{
    // fitmatrix(energy_range.count(), elements_to_fit->size()+2); //+2 for compton and elastic //n_pileup)
    std::unordered_map<std::string, Spectra<T_real>> element_spectra;

    Fit_Parameters<T_real> fit_parameters = model->fit_parameters();
    //set all fit parameters to be fixed. We only want to fit element counts
    fit_parameters.set_all(E_Bound_Type::FIXED);

    T_real energy_offset = fit_parameters.value(STR_ENERGY_OFFSET);
    T_real energy_slope = fit_parameters.value(STR_ENERGY_SLOPE);
    T_real energy_quad = fit_parameters.value(STR_ENERGY_QUADRATIC);

    ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    ArrayTr<T_real> ev = energy_offset + (energy * energy_slope) + (pow(energy, (T_real)2.0) * energy_quad);


    for (const auto& itr : (*elements_to_fit))
    {
        Fit_Element_Map<T_real>* element = itr.second;
        // Set value to 0.0 . This is the pre_faktor in gauss_tails_model. we do 10.0 ^ pre_faktor = 1.0
        if (false == fit_parameters.contains(itr.first))
        {
            Fit_Param<T_real> fp(itr.first, MIN_COUNTS_LIMIT_LOG, MAX_COUNTS_LIMIT_LOG, (T_real)0.0, STEP_COUNTS_LIMIT_LOG, E_Bound_Type::LIMITED_LO_HI);
            fit_parameters[itr.first] = fp;
        }
        else
        {
            fit_parameters[itr.first].value = (T_real)0.0;
        }
        element_spectra[itr.first] = model->model_spectrum_element(&fit_parameters, element, ev, nullptr);
    }
    //i = elements_to_fit->size();
    // scattering:
    // elastic peak

    Spectra<T_real> elastic_model(energy_range.count());
    // Set value to 0 because log10(0) = 1.0
    fit_parameters[STR_COHERENT_SCT_AMPLITUDE].value = 0.0;
    elastic_model += model->elastic_peak(&fit_parameters, ev, fit_parameters.at(STR_ENERGY_SLOPE).value);
    element_spectra[STR_COHERENT_SCT_AMPLITUDE] = elastic_model;
    //Set it so we fit coherent amp in fit params
    ///(*fit_params)[STR_COHERENT_SCT_AMPLITUDE].bound_type = data_struct::E_Bound_Type::FIT;


    // compton peak
    Spectra<T_real> compton_model(energy_range.count());
    // Set value to 0 because log10(0) = 1.0
    fit_parameters[STR_COMPTON_AMPLITUDE].value = 0.0;
    compton_model += model->compton_peak(&fit_parameters, ev, fit_parameters.at(STR_ENERGY_SLOPE).value);
    element_spectra[STR_COMPTON_AMPLITUDE] = compton_model;
    //Set it so we fit STR_COMPTON_AMPLITUDE  in fit params
    ///(*fit_params)[STR_COMPTON_AMPLITUDE].bound_type = data_struct::FIT;

    return element_spectra;

}


// ----------------------------------------------------------------------------

template<typename T_real>
std::unordered_map<std::string, Spectra<T_real>> Matrix_Optimized_Fit_Routine<T_real>::_generate_element_models_mp(models::Base_Model<T_real>* const model,
                                                                                      const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                                                      struct Range energy_range)
{
    // fitmatrix(energy_range.count(), elements_to_fit->size()+2); //+2 for compton and elastic //n_pileup)
    std::unordered_map<std::string, Spectra<T_real>> element_spectra;

    Fit_Parameters<T_real> fit_parameters = model->fit_parameters();
    //set all fit parameters to be fixed. We only want to fit element counts
    fit_parameters.set_all(E_Bound_Type::FIXED);

    T_real energy_offset = fit_parameters.value(STR_ENERGY_OFFSET);
    T_real energy_slope = fit_parameters.value(STR_ENERGY_SLOPE);
    T_real energy_quad = fit_parameters.value(STR_ENERGY_QUADRATIC);

    ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    ArrayTr<T_real> ev = energy_offset + (energy * energy_slope) + (pow(energy, (T_real)2.0) * energy_quad);


#ifdef _OPENMP
    std::vector<std::string> keys;
    for (const auto& itr : (*elements_to_fit))
    {
        if (false == fit_parameters.contains(itr.first))
        {
            Fit_Param<T_real> fp(itr.first, MIN_COUNTS_LIMIT_LOG, MAX_COUNTS_LIMIT_LOG, (T_real)0.0, STEP_COUNTS_LIMIT_LOG, E_Bound_Type::LIMITED_LO_HI);
            fit_parameters[itr.first] = fp;
        }
        else
        {
            fit_parameters[itr.first].value = (T_real)0.0;
        }

        if (itr.first == STR_COHERENT_SCT_AMPLITUDE || itr.first == STR_COMPTON_AMPLITUDE)
        {
            continue;
        }
        else
        {
            keys.push_back(itr.first);
        }
    }
#pragma omp parallel for
    for (int i = 0; i < (int)keys.size(); i++)
    {
        Fit_Element_Map<T_real>* element = elements_to_fit->at(keys[i]);
        // Set value to 0.0 . This is the pre_faktor in gauss_tails_model. we do 10.0 ^ pre_faktor = 1.0
        data_struct::Spectra<T_real> spec = model->model_spectrum_element(&fit_parameters, element, ev, nullptr);
#pragma omp critical
        element_spectra[keys[i]] = spec;
    }
#else
    for(const auto& itr : (*elements_to_fit))
    {
        Fit_Element_Map<T_real>* element = itr.second;
        // Set value to 0.0 . This is the pre_faktor in gauss_tails_model. we do 10.0 ^ pre_faktor = 1.0
        if( false == fit_parameters.contains(itr.first) )
        {
            Fit_Param<T_real> fp(itr.first, MIN_COUNTS_LIMIT_LOG, MAX_COUNTS_LIMIT_LOG, (T_real)0.0, STEP_COUNTS_LIMIT_LOG, E_Bound_Type::LIMITED_LO_HI);
            fit_parameters[itr.first] = fp;
        }
        else
        {
            fit_parameters[itr.first].value = (T_real)0.0;
        }
        element_spectra[itr.first] = model->model_spectrum_element(&fit_parameters, element, ev, nullptr);
    }
#endif
    //i = elements_to_fit->size();
    // scattering:
    // elastic peak

    Spectra<T_real> elastic_model(energy_range.count());
    // Set value to 0 because log10(0) = 1.0
    fit_parameters[STR_COHERENT_SCT_AMPLITUDE].value = 0.0;
    elastic_model += model->elastic_peak(&fit_parameters, ev, fit_parameters.at(STR_ENERGY_SLOPE).value);
    element_spectra[STR_COHERENT_SCT_AMPLITUDE] = elastic_model;
    //Set it so we fit coherent amp in fit params
    ///(*fit_params)[STR_COHERENT_SCT_AMPLITUDE].bound_type = data_struct::E_Bound_Type::FIT;


    // compton peak
    Spectra<T_real> compton_model(energy_range.count());
    // Set value to 0 because log10(0) = 1.0
    fit_parameters[STR_COMPTON_AMPLITUDE].value = 0.0;
    compton_model += model->compton_peak(&fit_parameters, ev, fit_parameters.at(STR_ENERGY_SLOPE).value);
    element_spectra[STR_COMPTON_AMPLITUDE] = compton_model;
    //Set it so we fit STR_COMPTON_AMPLITUDE  in fit params
    ///(*fit_params)[STR_COMPTON_AMPLITUDE].bound_type = data_struct::FIT;

    return element_spectra;

}

// ----------------------------------------------------------------------------

template<typename T_real>
void Matrix_Optimized_Fit_Routine<T_real>::initialize(models::Base_Model<T_real>* const model,
                                              const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                              const struct Range energy_range)
{

    this->_energy_range = energy_range;
    _element_models.clear();
    //logI<<"-------- Generating element models ---------"<<"\n";
    _element_models = _generate_element_models(model, elements_to_fit, energy_range);

    {
        std::lock_guard<std::mutex> lock(_int_spec_mutex);
        _integrated_fitted_spectra.setZero(energy_range.count());
        _integrated_background.setZero(energy_range.count());
    }

}

// ----------------------------------------------------------------------------

template<typename T_real>
void Matrix_Optimized_Fit_Routine<T_real>::initialize_mp(models::Base_Model<T_real>* const model,
    const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
    const struct Range energy_range)
{

    this->_energy_range = energy_range;
    _element_models.clear();
    //logI<<"-------- Generating element models ---------"<<"\n";
    _element_models = _generate_element_models_mp(model, elements_to_fit, energy_range);

    {
        std::lock_guard<std::mutex> lock(_int_spec_mutex);
        _integrated_fitted_spectra.setZero(energy_range.count());
        _integrated_background.setZero(energy_range.count());
    }

}

// ----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME Matrix_Optimized_Fit_Routine<T_real>:: fit_spectra(const models::Base_Model<T_real>* const model,
                                                            const Spectra<T_real>* const spectra,
                                                            const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                                            std::unordered_map<std::string, T_real>& out_counts)
{

    Fit_Parameters<T_real> fit_params = model->fit_parameters();
    //Add fit param for number of iterations
    fit_params.add_parameter(Fit_Param<T_real>(STR_NUM_ITR, 0.0));
    fit_params.add_parameter(Fit_Param<T_real>(STR_RESIDUAL, 0.0));
    this->_add_elements_to_fit_parameters(&fit_params, spectra, elements_to_fit);
    this->_calc_and_update_coherent_amplitude(&fit_params, spectra);
    OPTIMIZER_OUTCOME ret_val = OPTIMIZER_OUTCOME::FAILED;

    if(this->_optimizer != nullptr)
    {
        //todo : snip background here and pass to optimizer, then add to integrated background to save in h5
        
        ArrayTr<T_real> background;
        
        
        if(fit_params.contains(STR_SNIP_WIDTH))
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

        std::function<void(const models::Base_Model<T_real>* const model, const Fit_Parameters<T_real>* const, const  Range* const, Spectra<T_real>*)> gen_func = std::bind(&Matrix_Optimized_Fit_Routine<T_real>::model_spectrum, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);

        //set num iter to 300;
        std::unordered_map<std::string, T_real> opt_options{ {STR_OPT_MAXITER, 2000.} };
        std::unordered_map<std::string, T_real> saved_options = this->_optimizer->get_options();
        this->_optimizer->set_options(opt_options);


        ret_val = this->_optimizer->minimize_func(model, &fit_params, spectra, this->_energy_range, &background, gen_func, _use_weights);
        //Save the counts from fit parameters into fit count dict for each element
        for (auto el_itr : *elements_to_fit)
        {
            T_real value =  fit_params.at(el_itr.first).value;
            //convert from log10
            value = std::pow((T_real)10.0, value);
            out_counts[el_itr.first] = value;
        }

        out_counts[STR_NUM_ITR] = fit_params.at(STR_NUM_ITR).value;
        out_counts[STR_RESIDUAL] = fit_params.at(STR_RESIDUAL).value;

		//get max and top 10 max channels
        std::vector<std::pair<int, T_real> > max_map;
		data_struct::Spectra<T_real> max_vals = *spectra;
		typename data_struct::Spectra<T_real>::Index idx;
		for (int i = 0; i < 10; i++)
		{
			T_real max_val = max_vals.maxCoeff(&idx);
			max_map.push_back({ idx, max_val });
			max_vals[idx] = 0;
		}

		//model fit spectra
        Spectra<T_real> model_spectra(this->_energy_range.count());
        this->model_spectrum(model, &fit_params, &this->_energy_range, &model_spectra);
        
        model_spectra += background;
        model_spectra = (ArrayTr<T_real>)model_spectra.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });

		//lock and integrate results
		{
            std::lock_guard<std::mutex> lock(_int_spec_mutex);
            _integrated_fitted_spectra.add(model_spectra);
            _integrated_background.add(background);

			//we don't know the spectra size during initlaize() will have to resize here
			if (_max_channels_spectra.size() < spectra->size())
			{
				_max_channels_spectra.setZero(spectra->size());
			}
			if (_max_10_channels_spectra.size() < spectra->size())
			{
				_max_10_channels_spectra.setZero(spectra->size());
			}

			_max_channels_spectra[max_map[0].first] += max_map[0].second;

			for (auto &itr : max_map)
			{
				_max_10_channels_spectra[itr.first] += itr.second;
			}
        }

        this->_optimizer->set_options(saved_options);
    }

    return ret_val;

}

// ----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT Matrix_Optimized_Fit_Routine<float>;
TEMPLATE_CLASS_DLL_EXPORT Matrix_Optimized_Fit_Routine<double>;

} //namespace routines
} //namespace fitting
