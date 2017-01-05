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



#include "param_optimized_fit_routine.h"

#include <iostream>
#include <algorithm>
#include <math.h>

#include "Faddeeva.hh"
#include <string.h>

#define SQRT_2xPI (real_t)2.506628275 // sqrt ( 2.0 * M_PI )

using namespace data_struct::xrf;


namespace fitting
{
namespace routines
{

// ----------------------------------------------------------------------------

Param_Optimized_Fit_Routine::Param_Optimized_Fit_Routine() : Base_Fit_Routine()
{

    _optimizer = nullptr;

}

// ----------------------------------------------------------------------------

Param_Optimized_Fit_Routine::~Param_Optimized_Fit_Routine()
{

}

// ----------------------------------------------------------------------------

void Param_Optimized_Fit_Routine::_add_elements_to_fit_parameters(Fit_Parameters *fit_params,
                                                                  const Spectra * const spectra,
                                                                  const Fit_Element_Map_Dict * const elements_to_fit)
{


    real_t this_factor = (real_t)14.0; //used to be 10.0

    for (auto el_itr : *elements_to_fit)
    {
        if( false == fit_params->contains(el_itr.first) )
        {
            real_t e_guess = (real_t)1.0e-10;
            data_struct::xrf::Fit_Element_Map *element = el_itr.second;

            //element->init_energy_ratio_for_detector_element( detector->get_element() );
            std::vector<Element_Energy_Ratio> energies = element->energy_ratios();

            //if element counts is not in fit params structure, add it
            if( false == fit_params->contains(el_itr.first) )
            {
                data_struct::xrf::Fit_Param fp(element->full_name(), (real_t)-11.0, 300, e_guess, (real_t)0.00001, data_struct::xrf::E_Bound_Type::FIT);
                (*fit_params)[el_itr.first] = fp;
            }

            if(spectra != nullptr  && energies.size() > 0)
            {
                real_t e_energy = element->energy_ratios()[0].energy;
                real_t min_e =  e_energy - (real_t)0.1;
                real_t max_e =  e_energy + (real_t)0.1;

                struct Range energy_range = data_struct::xrf::get_energy_range(min_e, max_e, spectra->size(), fit_params->at(STR_ENERGY_OFFSET).value, fit_params->at(STR_ENERGY_SLOPE).value);

                real_t sum = (*spectra)[std::slice(energy_range.min, energy_range.count(), 1)].sum();
                sum /= energy_range.count();
                e_guess = std::max( sum * this_factor + (real_t)0.01, (real_t)1.0);
                //e_guess = std::max( (spectra->mean(energy_range.min, energy_range.max + 1) * this_factor + (real_t)0.01), 1.0);
                e_guess = std::log10(e_guess);

                (*fit_params)[el_itr.first].value = e_guess;
            }
            else
            {
                e_guess = std::log10(e_guess);
                (*fit_params)[el_itr.first].value = e_guess;
            }
        }
    }

    if( false == fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        //add number of iteration it took
        data_struct::xrf::Fit_Param fp(data_struct::xrf::STR_NUM_ITR, (real_t)-1.0, 999999, 0.0, (real_t)0.00001, data_struct::xrf::E_Bound_Type::FIXED);
        (*fit_params)[data_struct::xrf::STR_NUM_ITR] = fp;
    }
    (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = 0.0;
}

// ----------------------------------------------------------------------------

void Param_Optimized_Fit_Routine::_calc_and_update_coherent_amplitude(Fit_Parameters* fitp,
                                                                      const Spectra * const spectra)
{
    //STR_COHERENT_SCT_ENERGY
    //STR_COHERENT_SCT_AMPLITUDE
    real_t min_e = fitp->at(STR_COHERENT_SCT_ENERGY).value - (real_t)0.4;
    real_t max_e = fitp->at(STR_COHERENT_SCT_ENERGY).value + (real_t)0.4;
    real_t this_factor = (real_t)35.0; //was 8.0 in MAPS, this gets closer though
    fitting::models::Range energy_range = fitting::models::get_energy_range(min_e, max_e, spectra->size(), fitp->at(STR_ENERGY_OFFSET).value, fitp->at(STR_ENERGY_SLOPE).value);
    size_t e_size = (energy_range.max + 1) - energy_range.min;
    real_t sum = (*spectra)[std::slice(energy_range.min, e_size, 1)].sum();
    sum /= energy_range.count();
    real_t e_guess = std::max(sum * this_factor + (real_t)0.01, (real_t)1.0);
    real_t logval = std::log10(e_guess);
    (*fitp)[STR_COMPTON_AMPLITUDE].value = logval;
    (*fitp)[STR_COHERENT_SCT_AMPLITUDE].value = logval;

}

// ----------------------------------------------------------------------------

void Param_Optimized_Fit_Routine::_save_counts(Fit_Parameters *fit_params,
                                               const Spectra * const spectra,
                                               const Fit_Element_Map_Dict * const elements_to_fit,
                                               Fit_Count_Dict * out_counts_dic,
                                               size_t row_idx,
                                               size_t col_idx)
{
    //Save the counts from fit parameters into fit count dict for each element
    for (auto el_itr : *elements_to_fit)
    {
        //data_struct::xrf::Fit_Element_Map *element = el_itr.second;
        real_t value =  fit_params->at(el_itr.first).value;
        // save values counts / second

        //convert from log10
        value = std::pow(10.0, value);
/*
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
        */
            out_counts_dic->at(el_itr.first).counts[row_idx][col_idx] = value;
       // }

    }

    //check if we are saving the number of iterations and save if so
    if(fit_params->contains(data_struct::xrf::STR_NUM_ITR) && out_counts_dic->count(data_struct::xrf::STR_NUM_ITR) > 0)
    {
        out_counts_dic->at(data_struct::xrf::STR_NUM_ITR).counts[row_idx][col_idx] = fit_params->at(data_struct::xrf::STR_NUM_ITR).value;
    }

}

// ----------------------------------------------------------------------------

std::unordered_map<std::string, real_t> Param_Optimized_Fit_Routine::fit_spectra(const models::Base_Model * const model,
                                                                                 const Spectra * const spectra,
                                                                                 const Fit_Element_Map_Dict * const elements_to_fit)
{
    //int xmin = np.argmin(abs(x - (fitp.g.xmin - fitp.s.val[keywords.energy_pos[0]]) / fitp.s.val[keywords.energy_pos[1]]));
    //int xmax = np.argmin(abs(x - (fitp.g.xmax - fitp.s.val[keywords.energy_pos[0]]) / fitp.s.val[keywords.energy_pos[1]]));
    // fitp.g.xmin = MIN_ENERGY_TO_FIT
    // fitp.g.xmax = MAX_ENERGY_TO_FIT
    /*
    Range energy_range = get_energy_range(1, 11,spectra_volume->samples_size(), detector);
    */

    std::unordered_map<std::string, real_t> counts_dict;
    Fit_Parameters fit_params = model->fit_parameters();
    //Add fit param for number of iterations
    fit_params.add_parameter(data_struct::xrf::STR_NUM_ITR, Fit_Param(data_struct::xrf::STR_NUM_ITR));
    _add_elements_to_fit_parameters(&fit_params, spectra, elements_to_fit);
    _calc_and_update_coherent_amplitude(&fit_params, spectra);

    //If the sum of the spectra we are trying to fit to is zero then set out counts to -10.0 == log(0.0000000001)
    if(spectra->sum() == 0)
    {

        fit_params.set_all_value(-10.0, E_Bound_Type::FIT);

        for (auto el_itr : *elements_to_fit)
        {
            counts_dict[el_itr.first] = -10.0;
        }
        return counts_dict;
    }

    if(_optimizer != nullptr)
    {
        _optimizer->minimize(&fit_params, spectra, elements_to_fit, model);

        //Save the counts from fit parameters into fit count dict for each element
        for (auto el_itr : *elements_to_fit)
        {
            real_t value =  fit_params.at(el_itr.first).value;
            //convert from log10
            value = std::pow(10.0, value);
            counts_dict[el_itr.first] = value;
        }

        //model->update_fit_params_values(fit_params);

        //check if we are saving the number of iterations and save if so
        if(fit_params.contains(data_struct::xrf::STR_NUM_ITR))
        {
            counts_dict[data_struct::xrf::STR_NUM_ITR] = fit_params.at(data_struct::xrf::STR_NUM_ITR).value;
        }
    }

    return counts_dict;
}

// ----------------------------------------------------------------------------
//TODO: FIX SO this isn't a repeat function
Fit_Parameters Param_Optimized_Fit_Routine::fit_spectra_parameters(const models::Base_Model * const model,
                                                        const Spectra * const spectra,
                                                        const Fit_Element_Map_Dict * const elements_to_fit)
{
    //int xmin = np.argmin(abs(x - (fitp.g.xmin - fitp.s.val[keywords.energy_pos[0]]) / fitp.s.val[keywords.energy_pos[1]]));
    //int xmax = np.argmin(abs(x - (fitp.g.xmax - fitp.s.val[keywords.energy_pos[0]]) / fitp.s.val[keywords.energy_pos[1]]));
    // fitp.g.xmin = MIN_ENERGY_TO_FIT
    // fitp.g.xmax = MAX_ENERGY_TO_FIT
    /*
    Range energy_range = get_energy_range(1, 11,spectra_volume->samples_size(), detector);
    */

    std::unordered_map<std::string, real_t> counts_dict;
    Fit_Parameters fit_params = model->fit_parameters();
    //Add fit param for number of iterations
    fit_params.add_parameter(data_struct::xrf::STR_NUM_ITR, Fit_Param(data_struct::xrf::STR_NUM_ITR));
    _add_elements_to_fit_parameters(&fit_params, spectra, elements_to_fit);
    _calc_and_update_coherent_amplitude(&fit_params, spectra);

    //If the sum of the spectra we are trying to fit to is zero then set out counts to -10.0 == log(0.0000000001)
    if(spectra->sum() == 0)
    {

        fit_params.set_all_value(-10.0, E_Bound_Type::FIT);

        for (auto el_itr : *elements_to_fit)
        {
            counts_dict[el_itr.first] = -10.0;
        }
    }
    else
    {
        if(_optimizer != nullptr)
        {
            _optimizer->minimize(&fit_params, spectra, elements_to_fit, model);

            //Save the counts from fit parameters into fit count dict for each element
            for (auto el_itr : *elements_to_fit)
            {
                real_t value =  fit_params.at(el_itr.first).value;
                //convert from log10
                value = std::pow(10.0, value);
                counts_dict[el_itr.first] = value;
            }

            //model->update_fit_params_values(fit_params);

            //check if we are saving the number of iterations and save if so
            if(fit_params.contains(data_struct::xrf::STR_NUM_ITR))
            {
                counts_dict[data_struct::xrf::STR_NUM_ITR] = fit_params.at(data_struct::xrf::STR_NUM_ITR).value;
            }
        }
    }

    return fit_params;
}


// ----------------------------------------------------------------------------

void Param_Optimized_Fit_Routine::initialize(models::Base_Model * const model,
                                             const Fit_Element_Map_Dict * const elements_to_fit,
                                             const struct Range energy_range)
{

}

// ----------------------------------------------------------------------------

void Param_Optimized_Fit_Routine::set_optimizer(Optimizer *optimizer)
{

    _optimizer = optimizer;

}

// ----------------------------------------------------------------------------

} //namespace routines
} //namespace fitting
