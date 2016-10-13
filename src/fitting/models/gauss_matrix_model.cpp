/***

Copyright (c) 2016 Arthur Glowacki

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

***/

#include "gauss_matrix_model.h"

//debug
#include <iostream>

namespace fitting
{
namespace models
{


Gauss_Matrix_Model::Gauss_Matrix_Model() : Gauss_Tails_Model()
{
    //TODO: maybe not use pointer?
    _element_models = nullptr;
    _counts_log_10 = true;
    _update_element_guess_value = true;
}

// ----------------------------------------------------------------------------

Gauss_Matrix_Model::~Gauss_Matrix_Model()
{

    std::cout<<"******** destroy element models *******"<<std::endl;
    if(_element_models != nullptr)
        delete _element_models;
    _element_models = nullptr;

}

// -----------------------------------------------------------------------------

Spectra Gauss_Matrix_Model::model_spectrum(const Fit_Parameters * const fit_params,
                                           const Spectra * const spectra,
                                           const Calibration_Standard * const calibration,
                                           const Fit_Element_Map_Dict * const elements_to_fit,
                                           const struct Range energy_range)
{
    Spectra spectra_model(energy_range.count());

//    valarray<real_t> energy((real_t)0.0, energy_range.count());
//    real_t e_val = energy_range.min;
//    for(int i=0; i < (energy_range.max - energy_range.min )+1; i++)
//    {
//        energy[i] = e_val;
//        e_val += 1.0;
//    }

//    real_t gain = calibration->slope();
//    valarray<real_t> ev = calibration->offset() + energy * calibration->slope() + pow(energy, (real_t)2.0) * calibration->quad();
/*
    if( _snip_background )
    {
        real_t spectral_binning = 0.0;
        spectra->snip_background(_background_counts, calibration->offset(), calibration->slope(), calibration->quad(), spectral_binning, fit_params->at(STR_SNIP_WIDTH).value, energy_range.min, energy_range.max);
    }
*/
/*
    if (keywords.spectral_binning > 0)
    {
        ind = energy/keywords.spectral_binning;
        counts_background = keywords.background[ind.astype(int)];
    }
    else
    {
        counts_background = keywords.background[energy];
    }
*/
    if(_element_models != nullptr)
    {
        for(const auto& itr : *_element_models)
        {
            if(fit_params->contains(itr.first))
            {
                Fit_Param param = fit_params->at(itr.first);
                real_t va = pow(10.0, param.value);
                spectra_model += pow((real_t)10.0, param.value) * itr.second;
            }
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
//   *counts += _background_counts;

    return spectra_model;
}

// ----------------------------------------------------------------------------

void Gauss_Matrix_Model::initialize(Fit_Parameters *fit_params,
                                    const Calibration_Standard * const calibration,
                                    const Fit_Element_Map_Dict * const elements_to_fit,
                                    const struct Range energy_range)
{

    Base_Model::initialize(fit_params, calibration, elements_to_fit, energy_range);

    if(_element_models == nullptr)// || _element_models->size() != elements_to_fit->size())
    {
        _element_models = new unordered_map<string, Spectra>();
    }
    else
    {
        _element_models->clear();
    }

    std::cout<<"-------- Generating element models ---------"<<std::endl;
    *_element_models = _generate_element_models(fit_params, calibration, elements_to_fit, energy_range);

}

// ----------------------------------------------------------------------------

unordered_map<string, Spectra> Gauss_Matrix_Model::_generate_element_models(Fit_Parameters *fit_params,
                                                                            const Calibration_Standard * const calibration,
                                                                            const Fit_Element_Map_Dict * const elements_to_fit,
                                                                            struct Range energy_range)
{
    //Eigen::MatrixXd fitmatrix(energy_range.count(), elements_to_fit->size()+2); //+2 for compton and elastic //n_pileup)
    unordered_map<string, Spectra> element_spectra;

    //n_pileup = 9
    //valarray<real_t> value(0.0, energy_range.count());
    valarray<real_t> counts(0.0, energy_range.count());

    Fit_Parameters fit_parameters = *fit_params;

    valarray<real_t> energy((real_t)0.0, energy_range.count());
    real_t e_val = energy_range.min;
    for(int i=0; i < (energy_range.max - energy_range.min )+1; i++)
    {
        energy[i] = e_val;
        e_val += 1.0;
    }

    real_t gain = calibration->slope();
    valarray<real_t> ev = calibration->offset() + energy * calibration->slope() + pow(energy, (real_t)2.0) * calibration->quad();

    int i = 0;
    for(const auto& itr : (*elements_to_fit))
    {
        Fit_Element_Map* element = itr.second;
        // Set value to 0.0 . This is the pre_faktor in gauss_tails_model. we do 10.0 ^ pre_faktor = 1.0
        if( false == fit_parameters.contains(itr.first) )
        {
            data_struct::xrf::Fit_Param fp(itr.first, (real_t)-100.0, std::numeric_limits<real_t>::max(), 0.0, (real_t)0.00001, data_struct::xrf::E_Bound_Type::FIT);
            fit_parameters[itr.first] = fp;
        }
        else
        {
            fit_parameters[itr.first].value = 0.0;
        }
        element_spectra[itr.first] = model_spectrum_element(&fit_parameters, element, calibration, energy);
    }

    //i = elements_to_fit->size();
    // scattering:
    // elastic peak

    Spectra elastic_model(energy_range.count());
    // Set value to 0 because log10(0) = 1.0
    fit_parameters[STR_COHERENT_SCT_AMPLITUDE].value = 0.0;
    elastic_model += elastic_peak(&fit_parameters, ev, gain);
    element_spectra[STR_COHERENT_SCT_AMPLITUDE] = elastic_model;
    //Set it so we fit coherent amp in fit params
    (*fit_params)[STR_COHERENT_SCT_AMPLITUDE].bound_type = data_struct::xrf::E_Bound_Type::FIT;


    // compton peak
    Spectra compton_model(energy_range.count());
    // Set value to 0 because log10(0) = 1.0
    fit_parameters[STR_COMPTON_AMPLITUDE].value = 0.0;
    compton_model += compton_peak(&fit_parameters, ev, gain);
    element_spectra[STR_COMPTON_AMPLITUDE] = compton_model;
    //Set it so we fit STR_COMPTON_AMPLITUDE  in fit params
    (*fit_params)[STR_COMPTON_AMPLITUDE].bound_type = data_struct::xrf::FIT;

    /*
    //int this_i = i + 2;
        i = np.amax(keywords.mele_pos)-np.amin(keywords.kele_pos)+1+ii;
        if (add_pars[i, j].energy <= 0.0)
        {
            continue;
        }
        delta_energy = ev.copy() - (add_pars[i, j].energy);
        faktor = add_pars[i, j].ratio;
        counts = faktor * this->model_gauss_peak(gain, sigma[i, j], delta_energy);

        //fitmatrix[:, this_i+ii] = fitmatrix[:, this_i+ii]+counts[:];
        fitmatrix.row(this_i + ii) = fitmatrix.row(this_i + ii) + counts;
        counts = 0.0;
    }
    */
    //return fitmatrix;
    return element_spectra;

}


} //namespace models
} //namespace fitting
