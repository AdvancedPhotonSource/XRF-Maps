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

#include "base_model.h"
#include <functional>
#include <numeric>
#include <iostream>

namespace fitting
{
namespace models
{


Range get_energy_range(real_t min_energy, real_t max_energy, size_t spectra_size, const Calibration_Standard * const calibration)
{
    //real_t MIN_ENERGY_TO_FIT = 1.0;
    //real_t MAX_ENERGY_TO_FIT = 11.0;

    real_t MIN_ENERGY_TO_FIT = min_energy;
    real_t MAX_ENERGY_TO_FIT = max_energy;


    struct Range energy_range;
    energy_range.min = (int)ceil( (MIN_ENERGY_TO_FIT - calibration->offset()) / calibration->slope() );
    energy_range.max = (int)ceil( (MAX_ENERGY_TO_FIT - calibration->offset()) / calibration->slope() );
    //if (xmax > used_chan - 1) or (xmax <= np.amin([xmin, used_chan / 20.])):
    if ( (energy_range.max > spectra_size - 1) || (energy_range.max <= energy_range.min) )
    {
        energy_range.max = spectra_size - 1;
    }
    if (energy_range.min < 0 || energy_range.min > energy_range.max)
    {
        energy_range.min = 0;
    }
    return energy_range;

}

// ====================================================================================================================

Base_Model::Base_Model()
{
    _update_element_guess_value = true;
    _counts_log_10 = false;
    _save_counts_per_sec = true;
}

// --------------------------------------------------------------------------------------------------------------------

Base_Model::~Base_Model()
{

}

// --------------------------------------------------------------------------------------------------------------------

Fit_Parameters Base_Model::fit_spectra(const Fit_Parameters fit_params,
                                       const Spectra * const spectra,
                                       const Calibration_Standard * const calibration,
                                       const Fit_Element_Map_Dict * const elements_to_fit,
                                       Fit_Count_Dict *out_counts_dic,
                                       size_t row_idx,
                                       size_t col_idx)
{

    Fit_Parameters local_fit_params = fit_params;

    _pre_process(&local_fit_params, spectra, calibration, elements_to_fit);

    _fit_spectra(&local_fit_params, spectra, calibration, elements_to_fit);

    _post_process(&local_fit_params, spectra, calibration, elements_to_fit, out_counts_dic, row_idx, col_idx);

    return local_fit_params;

}

// --------------------------------------------------------------------------------------------------------------------

void Base_Model::fit_spectra_volume(const Fit_Parameters fit_params,
                                    const Spectra_Volume * const spectra_volume,
                                    const Calibration_Standard * const calibration,
                                    const Fit_Element_Map_Dict * const elements_to_fit,
                                    Fit_Count_Dict *out_counts_dic,
                                    size_t row_idx,
                                    size_t col_idx)
{
/*
    Fit_Parameters local_fit_params = fit_params;

    _pre_process(&local_fit_params, spectra_volume, calibration, elements_to_fit, row_idx, col_idx);

    _fit_spectra(&local_fit_params, spectra_volume, calibration, elements_to_fit, row_idx, col_idx);

    _post_process(&local_fit_params, spectra_volume, calibration, elements_to_fit, row_idx, col_idx);
*/
}

// --------------------------------------------------------------------------------------------------------------------

void Base_Model::initialize(Fit_Parameters *fit_params,
                            const Calibration_Standard * const calibration,
                            const Fit_Element_Map_Dict * const elements_to_fit,
                            const struct Range energy_range)
{

    _add_elements_to_fit_parameters(fit_params, nullptr, calibration, elements_to_fit);

}

// --------------------------------------------------------------------------------------------------------------------

void Base_Model::_pre_process(Fit_Parameters *fit_params,
                              const Spectra * const spectra,
                              const Calibration_Standard * const calibration,
                              const Fit_Element_Map_Dict * const elements_to_fit)
{

    if (_update_element_guess_value)
    {
        _update_elements_guess(fit_params, spectra, calibration, elements_to_fit);
    }

}

// --------------------------------------------------------------------------------------------------------------------

void Base_Model::_post_process(Fit_Parameters *fit_params,
                               const Spectra * const spectra,
                               const Calibration_Standard * const calibration,
                               const Fit_Element_Map_Dict * const elements_to_fit,
                               Fit_Count_Dict * out_counts_dic,
                               size_t row_idx,
                               size_t col_idx)
{

    for (auto el_itr : *elements_to_fit)
    {
        //data_struct::xrf::Fit_Element_Map *element = el_itr.second;
        real_t value =  fit_params->at(el_itr.first).value;
        // save values counts / second
        if (_counts_log_10)
        {
            value = std::pow(10.0, value);
        }
        if(_save_counts_per_sec)
        {
            real_t counts_per_sec = value / spectra->elapsed_lifetime();
            // if val is not nan then save it , otherwise leave it as 0.
            if(false == std::isnan(counts_per_sec))
            {
                out_counts_dic->at(el_itr.first).counts[row_idx][col_idx] = counts_per_sec;
            }
        }
        else
        {
            out_counts_dic->at(el_itr.first).counts[row_idx][col_idx] = value;
        }
        //(*element)[row_idx][col_idx] = value / spectra->elapsed_lifetime();
    }
    //std::cout<<std::endl;

}

// --------------------------------------------------------------------------------------------------------------------

void Base_Model::_add_elements_to_fit_parameters(Fit_Parameters *fit_params,
                                                 const Spectra * const spectra,
                                                 const Calibration_Standard * const calibration,
                                                 const Fit_Element_Map_Dict * const elements_to_fit)
{

    real_t this_factor = (real_t)10.0;

    for (auto el_itr : *elements_to_fit)
    {
        real_t e_guess = (real_t)1.0e-10;
        if(false == fit_params->contains(el_itr.first))
        {
            data_struct::xrf::Fit_Element_Map *element = el_itr.second;
            //Set Si detector
            element->init_energy_ratio_for_detector_element( Element_Info_Map::inst()->get_element("Si") );
            if(element->energy_ratios().size() > 0)
            {
                std::vector<Element_Energy_Ratio> energies = element->energy_ratios();
                if( false == (*fit_params).contains(el_itr.first) )
                {

                    if (spectra != nullptr && energies.size() > 0)
                    {
                        real_t e_energy = element->energy_ratios()[0].energy;
                        real_t min_e =  e_energy - (real_t)0.1;
                        real_t max_e =  e_energy + (real_t)0.1;

                        struct Range energy_range = fitting::models::get_energy_range(min_e, max_e, spectra->size(), calibration);
                        real_t sum = (*spectra)[std::slice(energy_range.min, energy_range.count(), 1)].sum();
                        sum /= energy_range.count();
                        e_guess = std::max( sum * this_factor + (real_t)0.01, (real_t)1.0);
                        //e_guess = std::max( (spectra->mean(energy_range.min, energy_range.max + 1) * this_factor + (real_t)0.01), 1.0);
                    }
                    e_guess = std::log10(e_guess);
                    //                             name                   min               max                               val       step size       fit
                    //data_struct::xrf::Fit_Param fp(element->full_name(), (real_t)1.0e-10, std::numeric_limits<real_t>::max(), e_guess, (real_t)0.00001, data_struct::xrf::E_Bound_Type::LIMITED_LO);
                    data_struct::xrf::Fit_Param fp(element->full_name(), (real_t)-11.0, 300, e_guess, (real_t)0.00001, data_struct::xrf::E_Bound_Type::FIT);
                    (*fit_params)[el_itr.first] = fp;
                    //std::cout<<"add "<<el_itr.first<<" val "<<e_guess<<std::endl;
                }
            }
        }
    }
    std::cout<<std::endl;
}

// --------------------------------------------------------------------------------------------------------------------

void Base_Model::_update_elements_guess(Fit_Parameters *fit_params,
                                        const Spectra * const spectra,
                                        const Calibration_Standard * const calibration,
                                        const Fit_Element_Map_Dict * const elements_to_fit)
{

    real_t this_factor = (real_t)14.0;

    for (auto el_itr : *elements_to_fit)
    {
        real_t e_guess = (real_t)1.0e-10;
        data_struct::xrf::Fit_Element_Map *element = el_itr.second;

        std::vector<Element_Energy_Ratio> energies = element->energy_ratios();
        if( fit_params->contains(el_itr.first) )
        {
            if(spectra != nullptr  && energies.size() > 0)
            {

                real_t e_energy = element->energy_ratios()[0].energy;
                real_t min_e =  e_energy - (real_t)0.1;
                real_t max_e =  e_energy + (real_t)0.1;

                struct Range energy_range = fitting::models::get_energy_range(min_e, max_e, spectra->size(), calibration);

                real_t sum = (*spectra)[std::slice(energy_range.min, energy_range.count(), 1)].sum();
                sum /= energy_range.count();
                e_guess = std::max( sum * this_factor + (real_t)0.01, (real_t)1.0);
                //e_guess = std::max( (spectra->mean(energy_range.min, energy_range.max + 1) * this_factor + (real_t)0.01), 1.0);
                e_guess = std::log10(e_guess);

                (*fit_params)[el_itr.first].value = e_guess;
                //std::cout<<"if "<<el_itr.first<<" val "<<e_guess<<std::endl;
            }
            else
            {
                e_guess = std::log10(e_guess);
                (*fit_params)[el_itr.first].value = e_guess;
                //std::cout<<"el "<<el_itr.first<<" val "<<e_guess<<std::endl;
            }
        }
    }

}


} //namespace models
} //namespace fitting
