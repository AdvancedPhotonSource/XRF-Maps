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


#include "fit_element_map.h"
#include <iostream>

namespace data_struct
{

//-----------------------------------------------------------------------------

Fit_Element_Map::Fit_Element_Map(std::string name, Element_Info* element_info)
{

    _full_name = name;

    _element_info = element_info;

    _pileup_element_info = nullptr;

    _width_multi = 1.0;

    size_t num_ratios = 1;

    int idx = _full_name.find_last_of("_") + 1;
    if(idx == 0)
    {
        _shell_type = "K";
    }
    else
    {
        _shell_type = _full_name.substr(idx);
    }


    if (_shell_type == "K") // K line
    {
        num_ratios = 4;
    }
    else if (_shell_type == "L")
    {
        num_ratios = 12;
    }
    for(size_t i=0; i<num_ratios; i++)
    {
        _energy_ratio_custom_multipliers.push_back(1.0);
    }

}

//-----------------------------------------------------------------------------

Fit_Element_Map::~Fit_Element_Map()
{

    _energy_ratios.clear();
    _energy_ratio_custom_multipliers.clear();

}

//-----------------------------------------------------------------------------

void Fit_Element_Map::init_energy_ratio_for_detector_element(const Element_Info * const detector_element)
{

    if(_element_info == nullptr)
    {
        //Don't log because we can have non elements such as Compton that doesn't have element properties
        //logE<<"Element info was not properly loaded. Variable is null! Can't initialize energy ratios!\n";
        return;
    }
    //real_t weight =  Element_Weight.at(element_info->number);

    _energy_ratios.clear();

    if (_shell_type == "K") // K line
    {
        if(_pileup_element_info != nullptr)
        {
            _center = _element_info->xrf["ka1"] + _pileup_element_info->xrf["ka1"];

            generate_energy_ratio(_element_info->xrf["ka1"] + _pileup_element_info->xrf["ka1"], 1.0, Element_Param_Type::Ka1_Line, detector_element);
            generate_energy_ratio(_element_info->xrf["kb1"] + _pileup_element_info->xrf["kb1"], (_element_info->xrf_abs_yield["kb1"] / _element_info->xrf_abs_yield["ka1"]) * (_pileup_element_info->xrf_abs_yield["kb1"] / _pileup_element_info->xrf_abs_yield["ka1"]), Element_Param_Type::Kb1_Line, detector_element);
            generate_energy_ratio(_element_info->xrf["ka1"] + _pileup_element_info->xrf["kb1"], (_pileup_element_info->xrf_abs_yield["kb1"] / _pileup_element_info->xrf_abs_yield["ka1"]), Element_Param_Type::Kb1_Line, detector_element);

            if(_element_info != _pileup_element_info)
            {
                generate_energy_ratio(_element_info->xrf["kb1"] + _pileup_element_info->xrf["ka1"], (_element_info->xrf_abs_yield["kb1"] / _element_info->xrf_abs_yield["ka1"]), Element_Param_Type::Kb1_Line, detector_element);
            }

        }
        else
        {
            _center = _element_info->xrf["ka1"];

            generate_energy_ratio(_element_info->xrf["ka1"], 1.0, Element_Param_Type::Ka1_Line, detector_element);
            generate_energy_ratio(_element_info->xrf["ka2"], (_element_info->xrf_abs_yield["ka2"] / _element_info->xrf_abs_yield["ka1"]) * _energy_ratio_custom_multipliers[1], Element_Param_Type::Ka2_Line, detector_element);
            generate_energy_ratio(_element_info->xrf["kb1"], (_element_info->xrf_abs_yield["kb1"] / _element_info->xrf_abs_yield["ka1"]) * _energy_ratio_custom_multipliers[2], Element_Param_Type::Kb1_Line, detector_element);
            generate_energy_ratio(_element_info->xrf["kb2"], (_element_info->xrf_abs_yield["kb2"] / _element_info->xrf_abs_yield["ka1"]) * _energy_ratio_custom_multipliers[3], Element_Param_Type::Kb2_Line, detector_element);
        }


    }
    else if (_shell_type == "L")
    {
        _center = _element_info->xrf["la1"];

        generate_energy_ratio(_element_info->xrf["la1"], 1.0, Element_Param_Type::La1_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["la2"], (_element_info->xrf_abs_yield["la2"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[1], Element_Param_Type::La2_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lb1"], (_element_info->xrf_abs_yield["lb1"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[2], Element_Param_Type::Lb1_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lb2"], (_element_info->xrf_abs_yield["lb2"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[3], Element_Param_Type::Lb2_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lb3"], (_element_info->xrf_abs_yield["lb3"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[4], Element_Param_Type::Lb3_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lb4"], (_element_info->xrf_abs_yield["lb4"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[5], Element_Param_Type::Lb4_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lg1"], (_element_info->xrf_abs_yield["lg1"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[6], Element_Param_Type::Lg1_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lg2"], (_element_info->xrf_abs_yield["lg2"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[7], Element_Param_Type::Lg2_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lg3"], (_element_info->xrf_abs_yield["lg3"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[8], Element_Param_Type::Lg3_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["lg4"], (_element_info->xrf_abs_yield["lg4"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[9], Element_Param_Type::Lg4_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["ll"], (_element_info->xrf_abs_yield["ll"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[10], Element_Param_Type::Ll_Line, detector_element);
        generate_energy_ratio(_element_info->xrf["ln"], (_element_info->xrf_abs_yield["ln"] / _element_info->xrf_abs_yield["la1"]) * _energy_ratio_custom_multipliers[11], Element_Param_Type::Ln_Line, detector_element);
    }
    else if(_shell_type == "M")
    {
        //TODO: finish adding M info
        generate_energy_ratio(_element_info->xrf["ma1"], 1.0, Element_Param_Type::Ma1_Line, detector_element);
        //generate_energy_ratio(_element_info->xrf["ma2"], 1.0, Element_Param_Type::Ma2_Line, detector_element);
        //generate_energy_ratio(_element_info->xrf["mb"], 1.0, Element_Param_Type::Mb_Line, detector_element);
        //generate_energy_ratio(_element_info->xrf["mg"], 1.0, Element_Param_Type::Mg_Line, detector_element);
    }

    _width = int( std::sqrt( std::pow(ENERGY_RES_OFFSET, 2) + std::pow( (_center * ENERGY_RES_SQRT), 2)  ) );
    // _center position for peaks/rois is in keV, widths of ROIs is in eV
    _width /= 2000.0; // divide by 2 and divide by 1000

}

//-----------------------------------------------------------------------------

void Fit_Element_Map::generate_energy_ratio(real_t energy, real_t ratio, Element_Param_Type et, const Element_Info * const detector_element)
{
    if(_element_info == nullptr)
    {
        //Don't log because we can have non elements such as Compton that doesn't have element properties
        //logE<<"Element info was not properly loaded. Variable is null! Can't initialize energy ratios!\n";
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

    real_t constant = RE * ((real_t)1.0e-16 * wavelength_angstroms * wavelength_angstroms) * molecules_per_cc / ((real_t)2.0 * (real_t)M_PI);

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

//-----------------------------------------------------------------------------

void Fit_Element_Map::set_custom_multiply_ratio(unsigned int idx, real_t multi)
{
    if (idx > 0 && idx < _energy_ratio_custom_multipliers.size())
    {
        _energy_ratio_custom_multipliers[idx] = multi;
    }
}

//-----------------------------------------------------------------------------

void Fit_Element_Map::set_as_pileup(std::string name, Element_Info* element_info)
{
    if(_pileup_element_info == nullptr)
    {
        int idx = name.find_last_of("_") + 1;
        if(idx == 0)
        {
            _pileup_shell_type = "K";
        }
        else
        {
            _pileup_shell_type = name.substr(idx);
        }
        _full_name += "_"+name;
        _pileup_element_info = element_info;
    }
    else
    {
        logW<<"Could not add second pileup to "<<_full_name<<"\n";
    }
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

data_struct::Fit_Element_Map* gen_element_map(std::string element_symb)
{
	data_struct::Element_Info_Map * element_info_map = data_struct::Element_Info_Map::inst();
	data_struct::Fit_Element_Map* fit_map = nullptr;

	element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());

	// check if element_symb contains '_'
	std::string base_element_symb = element_symb.substr(0, element_symb.find_last_of("_"));

	//logI<<element_symb<<" : "<<base_element_symb<<"\n";

	data_struct::Element_Info* e_info = element_info_map->get_element(base_element_symb);
	if (e_info == nullptr)
	{
		logW << "Can not find element " << base_element_symb << "\n";
	}
	else
	{
		fit_map = new data_struct::Fit_Element_Map(element_symb, e_info);
		fit_map->init_energy_ratio_for_detector_element(data_struct::Element_Info_Map::inst()->get_element("Si"));
	}

	return fit_map;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} //namespace data_struct
