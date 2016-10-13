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

#include "fit_element_map.h"
#include <iostream>

namespace data_struct
{
namespace xrf
{

Fit_Element_Map::Fit_Element_Map(std::string name, Element_Info* element_info)
{

    _full_name = name;

    _element_info = element_info;

    _width_multi = 1.0;

    int idx = _full_name.find_last_of("_") + 1;
    std::string shell_type = _full_name.substr(idx);

    size_t num_ratios = 1;
    if (idx == 0) // K line
    {
        num_ratios = 4;
    }
    else if (shell_type == "L")
    {
        num_ratios = 12;
    }
    for(size_t i=0; i<num_ratios; i++)
    {
        _energy_ratio_custom_multipliers.push_back(1.0);
    }

}

Fit_Element_Map::~Fit_Element_Map()
{

    _energy_ratios.clear();
    _energy_ratio_custom_multipliers.clear();

}

void Fit_Element_Map::init_energy_ratio_for_detector_element(Element_Info* detector_element)
{

    if(_element_info == nullptr)
    {
        //TODO: LOG Something
        return;
    }
    int idx = _full_name.find_last_of("_") + 1;
    std::string shell_type = _full_name.substr(idx);

    //real_t weight =  Element_Weight.at(element_info->number);

    _energy_ratios.clear();

    if (idx == 0) // K line
    {
        _center = _element_info->xrf["ka1"];

        generate_energy_ratio(_element_info->xrf["ka1"], 1.0, Element_Param_Type::Ka_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["ka2"], (_element_info->xrf_abs_yield["ka2"] / _element_info->xrf_abs_yield["ka1"]) * _energy_ratio_custom_multipliers[1], Element_Param_Type::Ka_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["kb1"], (_element_info->xrf_abs_yield["kb1"] / _element_info->xrf_abs_yield["ka1"]) * _energy_ratio_custom_multipliers[2], Element_Param_Type::Kb_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["kb2"], (_element_info->xrf_abs_yield["kb2"] / _element_info->xrf_abs_yield["ka1"]) * _energy_ratio_custom_multipliers[3], Element_Param_Type::Kb_Line, detector_element);
    }
    else if (shell_type == "L")
    {
        _center = _element_info->xrf["la1"];

        generate_energy_ratio(_element_info->xrf["la1"], 1.0, Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["la2"], (_element_info->xrf_abs_yield["la2"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[1], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lb1"], (_element_info->xrf_abs_yield["lb1"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[2], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lb2"], (_element_info->xrf_abs_yield["lb2"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[3], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lb3"], (_element_info->xrf_abs_yield["lb3"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[4], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lb4"], (_element_info->xrf_abs_yield["lb4"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[5], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lg1"], (_element_info->xrf_abs_yield["lg1"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[6], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lg2"], (_element_info->xrf_abs_yield["lg2"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[7], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lg3"], (_element_info->xrf_abs_yield["lg3"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[8], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lg4"], (_element_info->xrf_abs_yield["lg4"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[9], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["ll"], (_element_info->xrf_abs_yield["ll"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[10], Element_Param_Type::L_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["ln"], (_element_info->xrf_abs_yield["ln"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[11], Element_Param_Type::L_Line, detector_element);
    }
    else if(shell_type == "M")
    {
        //TODO: finish adding M info
        generate_energy_ratio(_element_info->xrf["ma1"], 1.0, Element_Param_Type::M_Line, detector_element);
    }

    _width = int( std::sqrt( std::pow(ENERGY_RES_OFFSET, 2) + std::pow( (_center * ENERGY_RES_SQRT), 2)  ) );

}

void Fit_Element_Map::generate_energy_ratio(real_t energy, real_t ratio, Element_Param_Type et, Element_Info* detector_element)
{
    if(_element_info == nullptr)
    {
        //TODO: LOG Something
        return;
    }

    real_t e_tmp = energy * (real_t)1000.0;
    real_t molecules_per_cc;
    real_t mu_fraction;
    real_t density = Henke_Compound_Density_Map.at(detector_element->name);

    if (energy == (real_t)0.0)
    {
        return;
    }

    if (_element_info->mass != 0.0)
    {
        molecules_per_cc = density * AVOGADRO / _element_info->mass;
    }
    else
    {
        molecules_per_cc = (real_t)0.0;
    }

    real_t wavelength_angstroms = HC_ANGSTROMS / e_tmp;
    // This constant has wavelength in angstroms and then
    // they are converted to centimeters

    real_t constant = RE * ((real_t)1.0e-16 * wavelength_angstroms * wavelength_angstroms) * molecules_per_cc / ((real_t)2.0 * M_PI);

    size_t hi_e_ind;
    size_t lo_e_ind;
    real_t lower_energy;
    real_t high_energy;

    _element_info->get_energies_between(e_tmp, &lower_energy, &high_energy, &lo_e_ind, &hi_e_ind);

    real_t ln_lower_energy = std::log(lower_energy);
    real_t ln_higher_energy = std::log(high_energy);
    real_t fraction = (std::log(e_tmp)-ln_lower_energy)/(ln_higher_energy-ln_lower_energy);
    real_t ln_f2_lower = std::log(std::abs(detector_element->f2_atomic_scattering_imaginary[lo_e_ind]));
    real_t ln_f2_higher = std::log(std::abs(detector_element->f2_atomic_scattering_imaginary[hi_e_ind]));
    real_t f2 = std::exp(ln_f2_lower + fraction * (ln_f2_higher - ln_f2_lower));
    real_t beta = constant * f2;

    mu_fraction = (e_tmp * (real_t)4.0 * (real_t)M_PI * beta) / (density * (real_t)1.239852 ) * (real_t)10000.0;

    _energy_ratios.push_back(Element_Energy_Ratio(energy, ratio, mu_fraction, et));
}

void Fit_Element_Map::set_custom_multiply_ratio(unsigned int idx, real_t multi)
{
    if (idx > 0 && idx < _energy_ratio_custom_multipliers.size())
    {
        _energy_ratio_custom_multipliers[idx] = multi;
    }
}

} //namespace xrf
} //namespace data_struct
