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

namespace fitting
{
namespace routines
{

std::mutex Matrix_Optimized_Fit_Routine::_int_spec_mutex;

// ----------------------------------------------------------------------------

Matrix_Optimized_Fit_Routine::Matrix_Optimized_Fit_Routine() : Param_Optimized_Fit_Routine()
{

}

// ----------------------------------------------------------------------------

Matrix_Optimized_Fit_Routine::~Matrix_Optimized_Fit_Routine()
{

    //logD<<"******** destroy element models *******"<<"\n";
    _element_models.clear();

}

// --------------------------------------------------------------------------------------------------------------------

void Matrix_Optimized_Fit_Routine::model_spectrum(const Fit_Parameters * const fit_params,
                                                  const struct Range * const energy_range,
												  Spectra* spectra_model)
{
	spectra_model->setZero();

    for(const auto& itr : _element_models)
    {
        if(fit_params->contains(itr.first))
        {
            Fit_Param param = fit_params->at(itr.first);
            (*spectra_model) += (pow((real_t)10.0, param.value) * itr.second);
        }
    }

/*
    if (np.sum(this->add_matrixfit_pars[3:6]) >= 0.)
    {
        ev = this->add_matrixfit_pars[keywords.energy_pos[0]] + energy * this->add_matrixfit_pars[keywords.energy_pos[1]] + (energy)**2 * this->add_matrixfit_pars[keywords.energy_pos[2]];
        counts_escape = counts.copy();
        counts_escape[:] = 0.0;
        if (this->add_matrixfit_pars[3] > 0.0)
        {
            real_t escape_E = 1.73998;
            wo = np.where(ev > escape_E+ev[0]);

            escape_factor = np.abs(p[len(p)-3] + p[len(p)-1] * ev);
            if (len(wo[0]) > 0)
            {
                for (size_t ii=0; ii<(len(wo[0]); ii++)
                {
                    counts_escape[ii] = counts[wo[0][ii]]*np.amax(np.append(escape_factor[wo[0][ii]],0.0));
                }
            }
            counts = counts + counts_escape;
        }
    }
*/

}

// ----------------------------------------------------------------------------

unordered_map<string, Spectra> Matrix_Optimized_Fit_Routine::_generate_element_models(models::Base_Model * const model,
                                                                                      const Fit_Element_Map_Dict * const elements_to_fit,
                                                                                      struct Range energy_range)
{
    // fitmatrix(energy_range.count(), elements_to_fit->size()+2); //+2 for compton and elastic //n_pileup)
    unordered_map<string, Spectra> element_spectra;

    //n_pileup = 9
    //valarray<real_t> value(0.0, energy_range.count());
    //real_t start_val = (real_t)0.0;
    //Spectra counts(energy_range.count());

    Fit_Parameters fit_parameters = model->fit_parameters();
    //set all fit parameters to be fixed. We only want to fit element counts
    fit_parameters.set_all(E_Bound_Type::FIXED);

    real_t energy_offset = fit_parameters.value(STR_ENERGY_OFFSET);
    real_t energy_slope = fit_parameters.value(STR_ENERGY_SLOPE);
    real_t energy_quad = fit_parameters.value(STR_ENERGY_QUADRATIC);

	ArrayXr energy = ArrayXr::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    ArrayXr ev = energy_offset + (energy * energy_slope) + (pow(energy, (real_t)2.0) * energy_quad);

    for(const auto& itr : (*elements_to_fit))
    {
        Fit_Element_Map* element = itr.second;
        // Set value to 0.0 . This is the pre_faktor in gauss_tails_model. we do 10.0 ^ pre_faktor = 1.0
        if( false == fit_parameters.contains(itr.first) )
        {
            Fit_Param fp(itr.first, (real_t)-100.0, std::numeric_limits<real_t>::max(), 0.0, (real_t)0.00001, E_Bound_Type::FIT);
            fit_parameters[itr.first] = fp;
        }
        else
        {
            fit_parameters[itr.first].value = 0.0;
        }
        element_spectra[itr.first] = model->model_spectrum_element(&fit_parameters, element, ev, nullptr);
    }

    //i = elements_to_fit->size();
    // scattering:
    // elastic peak

    Spectra elastic_model(energy_range.count());
    // Set value to 0 because log10(0) = 1.0
    fit_parameters[STR_COHERENT_SCT_AMPLITUDE].value = 0.0;
    elastic_model += model->elastic_peak(&fit_parameters, ev, fit_parameters.at(STR_ENERGY_SLOPE).value);
    element_spectra[STR_COHERENT_SCT_AMPLITUDE] = elastic_model;
    //Set it so we fit coherent amp in fit params
    ///(*fit_params)[STR_COHERENT_SCT_AMPLITUDE].bound_type = data_struct::E_Bound_Type::FIT;


    // compton peak
    Spectra compton_model(energy_range.count());
    // Set value to 0 because log10(0) = 1.0
    fit_parameters[STR_COMPTON_AMPLITUDE].value = 0.0;
    compton_model += model->compton_peak(&fit_parameters, ev, fit_parameters.at(STR_ENERGY_SLOPE).value);
    element_spectra[STR_COMPTON_AMPLITUDE] = compton_model;
    //Set it so we fit STR_COMPTON_AMPLITUDE  in fit params
    ///(*fit_params)[STR_COMPTON_AMPLITUDE].bound_type = data_struct::FIT;

    /*
    //int this_i = i + 2;
        i = np.amax(keywords.mele_pos)-np.amin(keywords.kele_pos)+1+ii;
        if (add_pars[i, j].energy <= 0.0)
        {
            continue;
        }
        delta_energy = ev.copy() - (add_pars[i, j].energy);
        faktor = add_pars[i, j].ratio;
        counts = faktor * this->model_gauss_peak(fit_parameters.at(STR_ENERGY_SLOPE).value, sigma[i, j], delta_energy);

        //fitmatrix[:, this_i+ii] = fitmatrix[:, this_i+ii]+counts[:];
        fitmatrix.row(this_i + ii) = fitmatrix.row(this_i + ii) + counts;
        counts = 0.0;
    }
    */
    //return fitmatrix;
    return element_spectra;

}

// ----------------------------------------------------------------------------

void Matrix_Optimized_Fit_Routine::initialize(models::Base_Model * const model,
                                              const Fit_Element_Map_Dict * const elements_to_fit,
                                              const struct Range energy_range)
{

    _energy_range = energy_range;
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

OPTIMIZER_OUTCOME Matrix_Optimized_Fit_Routine:: fit_spectra(const models::Base_Model * const model,
                                                            const Spectra * const spectra,
                                                            const Fit_Element_Map_Dict * const elements_to_fit,
                                                            std::unordered_map<std::string, real_t>& out_counts)
{

    Fit_Parameters fit_params = model->fit_parameters();
    //Add fit param for number of iterations
    fit_params.add_parameter(Fit_Param(STR_NUM_ITR, 0.0));
    fit_params.add_parameter(Fit_Param(STR_RESIDUAL, 0.0));
    _add_elements_to_fit_parameters(&fit_params, spectra, elements_to_fit);
    _calc_and_update_coherent_amplitude(&fit_params, spectra);
    OPTIMIZER_OUTCOME ret_val = OPTIMIZER_OUTCOME::FAILED;

    if(_optimizer != nullptr)
    {
        //todo : snip background here and pass to optimizer, then add to integrated background to save in h5
        
        ArrayXr background;
        
        
        if(fit_params.contains(STR_SNIP_WIDTH))
        {
            real_t spectral_binning = 0.0;
            ArrayXr bkg = snip_background(spectra,
                                         fit_params.value(STR_ENERGY_OFFSET),
                                         fit_params.value(STR_ENERGY_SLOPE),
                                         fit_params.value(STR_ENERGY_QUADRATIC),
                                         spectral_binning,
                                         fit_params.value(STR_SNIP_WIDTH),
                                         _energy_range.min,
                                         _energy_range.max);
            background = bkg.segment(_energy_range.min, _energy_range.count());
        }
        else
        {
            background.setZero(_energy_range.count());
        }

        std::function<void(const Fit_Parameters* const, const  Range* const, Spectra*)> gen_func = std::bind(&Matrix_Optimized_Fit_Routine::model_spectrum, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3);

        //set num iter to 300;
        unordered_map<string, real_t> opt_options{ {STR_OPT_MAXITER, 300.}, {STR_OPT_FTOL, 1.0e-11 }, {STR_OPT_GTOL, 1.0e-11 } };
        unordered_map<string, real_t> saved_options = _optimizer->get_options();        
        _optimizer->set_options(opt_options);


        ret_val = _optimizer->minimize_func(&fit_params, spectra, _energy_range, &background, gen_func);
        //Save the counts from fit parameters into fit count dict for each element
        for (auto el_itr : *elements_to_fit)
        {
            real_t value =  fit_params.at(el_itr.first).value;
            //convert from log10
            value = std::pow((real_t)10.0, value);
            out_counts[el_itr.first] = value;
        }

        out_counts[STR_NUM_ITR] = fit_params.at(STR_NUM_ITR).value;
        out_counts[STR_RESIDUAL] = fit_params.at(STR_RESIDUAL).value;

		//get max and top 10 max channels
		vector<pair<int, real_t> > max_map;
		data_struct::Spectra max_vals = *spectra;
		data_struct::Spectra::Index idx;
		for (int i = 0; i < 10; i++)
		{
			real_t max_val = max_vals.maxCoeff(&idx);
			max_map.push_back({ idx, max_val });
			max_vals[idx] = 0;
		}

		//model fit spectra
        Spectra model_spectra(_energy_range.count());
        this->model_spectrum(&fit_params, &_energy_range, &model_spectra);
        
        model_spectra += background;
        model_spectra = (ArrayXr)model_spectra.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });

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

        _optimizer->set_options(saved_options);
    }

    return ret_val;

}

// ----------------------------------------------------------------------------

} //namespace routines
} //namespace fitting
