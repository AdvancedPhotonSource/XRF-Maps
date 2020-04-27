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



#include "quantification_model.h"

#include <string>
#include <cmath>
#include <algorithm>

//debug
#include <iostream>

namespace quantification
{
namespace models
{

//-----------------------------------------------------------------------------

Quantification_Model::Quantification_Model()
{

}

//-----------------------------------------------------------------------------

Quantification_Model::~Quantification_Model()
{

}

//-----------------------------------------------------------------------------
/*
std::unordered_map<std::string, Element_Quant> Quantification_Model::generate_quant_map(real_t incident_energy,
                                                                                           Element_Info* detector_element,
                                                                                           Electron_Shell shell,
                                                                                           bool airpath,
                                                                                           real_t detector_chip_thickness,
                                                                                           real_t beryllium_window_thickness,
                                                                                           real_t germanium_dead_layer,
                                                                                           size_t start_z,
                                                                                           size_t end_z)
{
    std::unordered_map<std::string, Element_Quant> element_quant_map;
    start_z = std::max(start_z, (size_t)1);
    end_z = std::min(end_z, (size_t)92);

    for (int i=start_z; i <= end_z; i ++)
    {
        Element_Quant element_quant = generate_element_quant(incident_energy,
                                                             detector_element,
                                                             shell,
                                                             airpath,
                                                             detector_chip_thickness,
                                                             beryllium_window_thickness,
                                                             germanium_dead_layer,
                                                             i);
        Element_Info * element_info = Element_Info_Map::inst()->get_element(i);
        element_quant_map.emplace( std::pair <std::string, Element_Quant>(element_info->name, element_quant) );
    }

    return element_quant_map;
}
*/

std::vector<Element_Quant> Quantification_Model::generate_quant_vec(real_t incident_energy,
                                                                   Element_Info* detector_element,
                                                                   Electron_Shell shell,
                                                                   real_t airpath,
                                                                   real_t detector_chip_thickness,
                                                                   real_t beryllium_window_thickness,
                                                                   real_t germanium_dead_layer,
                                                                   size_t start_z,
                                                                   size_t end_z)
{
    start_z = std::max(start_z, (size_t)1);
    end_z = std::min(end_z, (size_t)92);
    std::vector<Element_Quant> element_quant_vec( end_z - start_z + 1 );

    size_t idx = 0;
    for (size_t i=start_z; i <= end_z; i ++)
    {
        element_quant_vec[idx] = generate_element_quant(incident_energy,
                                                         detector_element,
                                                         shell,
                                                         airpath,
                                                         detector_chip_thickness,
                                                         beryllium_window_thickness,
                                                         germanium_dead_layer,
                                                         i);
        idx++;
    }

    return element_quant_vec;
}

//-----------------------------------------------------------------------------

Element_Quant Quantification_Model::generate_element_quant(real_t incident_energy,
                                                              Element_Info* detector_element,
                                                              Electron_Shell shell,
                                                              real_t airpath,
                                                              real_t detector_chip_thickness,
                                                              real_t beryllium_window_thickness,
                                                              real_t germanium_dead_layer,
                                                              size_t z_number)
{

    //incident_E(incident_energy) == COHERENT_SCT_ENERGY: Maps_fit_params
    //fit_t_be == BE_WINDOW_THICKNESS * 1000
    //fit_t_ge == GE_DEAD_LAYER * 1000 if detector is not Si
    //add_float['a'] == DET_CHIP_THICKNESS   * 1000

    real_t beta;
    real_t shell_factor = 0.0;
    real_t ev;
    real_t jump_factor = 0.0;
    real_t total_jump_factor = 0.0;
    Element_Quant element_quant;
    Element_Info* element_info = Element_Info_Map::inst()->get_element(z_number);
    if(element_info == nullptr)
    {
        return element_quant;
    }
    switch(shell)
    {
    case K_SHELL:
        ev = element_info->xrf.at("ka1") * (real_t)1000.0;
        element_quant.yield = element_info->yieldD.at("k");
        if( incident_energy > element_info->bindingE.at("K") )
        {
            jump_factor = element_info->jump.at("K");
        }
        break;
    case L_SHELL:
        ev = element_info->xrf.at("la1") * (real_t)1000.0;
        element_quant.yield = element_info->xrf_abs_yield.at("la1");
        jump_factor = element_info->jump.at("L3");
        if( incident_energy > element_info->bindingE.at("L2") )
        {
            total_jump_factor = element_info->jump.at("L2");
        }
        if( incident_energy > element_info->bindingE.at("L1") )
        {
            total_jump_factor = total_jump_factor * element_info->jump.at("L1");
        }
        break;
    case M_SHELL:
        ev = element_info->xrf.at("ma1") * (real_t)1000.0;
        element_quant.yield = element_info->xrf_abs_yield.at("ma1");
        jump_factor = element_info->jump.at("M5");
        total_jump_factor = element_info->jump.at("M1") * element_info->jump.at("M2") * element_info->jump.at("M3") * element_info->jump.at("M4");
        break;
    default:
        throw "Unsupported shell. Currently only support for K, L, and M ";
        break;
    }

    if( jump_factor != 0.0 )
    {
        if( total_jump_factor == 0.0 )
        {
            shell_factor = (jump_factor - (real_t)1.0) / jump_factor;
        }
        else
        {
            shell_factor = (jump_factor - (real_t)1.0) / jump_factor / total_jump_factor;
        }
    }
    // replace straight henke routines, with those
    // that take the absorption edges into account
    // make sure we are a bit above the absorption edge to make sure that for calibration purposes we do not eoncouner any weird things.
    beta = element_info->calc_beta(element_info->density, (incident_energy + (real_t)0.1) * (real_t)1000.0);
    // stds in microgram/cm2
    // density rho = g/cm3 = 1 microgram/cm2 /1000/1000/cm = 1 microgram/cm2 /1000/1000/*10*1000/um = 1 microgram/cm2 /100/um
    // thickness for 1 ugr/cm2
    // =1/(density[g/cm3]/10)

    real_t thickness = (real_t)1.0 / (element_info->density * (real_t)10.0) * (real_t)1000.0;
    ////aux_arr[mm, 0] = self.absorption(thickness, beta, 1239.852/((self.maps_conf.incident_E+0.1)*1000.), shell_factor=shell_factor)
    element_quant.absorption = absorption(thickness, beta, (real_t)1239.852 / ((incident_energy + (real_t)0.1) * (real_t)1000.0), shell_factor);

    if(ev > 0)
    {
        beta  = Element_Info_Map::inst()->calc_beta("Be", (real_t)1.848, ev);
        ////aux_arr[mm, 1] = self.transmission(self.maps_conf.fit_t_be, beta, 1239.852/ev)
        element_quant.transmission_Be = transmission(beryllium_window_thickness, beta, (real_t)1239.852 / ev);

        beta  = Element_Info_Map::inst()->calc_beta("Ge", (real_t)5.323, ev);
        ////aux_arr[mm, 2] = self.transmission(self.maps_conf.fit_t_ge, beta, 1239.852/ev)
        element_quant.transmission_Ge = transmission(germanium_dead_layer, beta, (real_t)1239.852 / ev);
    }
    ////aux_arr[mm, 3] = yieldd
    //element_quant.yield = element_info->yieldD["K"]; //yieldd === newrel_yield * info_elements[element_temp].yieldD['k']

    if (detector_element->name == "Si" && detector_chip_thickness > 0,0 && ev > 0) //  (self.maps_conf.add_long['a'] == 1)
    {
        beta  = Element_Info_Map::inst()->calc_beta("Si", (real_t)2.3, ev);
        element_quant.transmission_through_Si_detector = transmission(detector_chip_thickness, beta, (real_t)1239.852 / ev);
    }
    ////aux_arr[mm, 4] = self.transmission(self.maps_conf.add_float['a'], beta, 1239.852/ev)
    else // ( (self.maps_conf.add_float['a'] == 0.) || (self.maps_conf.add_long['a'] != 1) )
    {
        ////aux_arr[mm, 4] = 0.
        element_quant.transmission_through_Si_detector = 0.0;
    }


    if( airpath > 0 && ev > 0)
    {
        //density = 1.0
        real_t density = (real_t)0.00117;
        //air_ele = 'N78.08O20.95Ar0.93'
        //density = 1.2047e-3
        //f1, f2, delta, beta, graze_mrad, reflect, inverse_mu, atwt = Chenke.get_henke_single('air', density, ev)
        beta = Element_Info_Map::inst()->calc_beta("N:78.08,O:20.95,Ar:0.93", density, ev); // air
        ////aux_arr[mm, 5] = self.transmission(airpath*1000., beta, 1239.852/ev)  // airpath is read in microns, transmission function expects nm
        // airpath is read in microns, transmission function expects nm
        element_quant.transmission_through_air = transmission( airpath * (real_t)1000.0, beta, (real_t)1239.852 / ev);
    }
    else
    {
        ////aux_arr[mm, 5] = 1.
        element_quant.transmission_through_air = 1.0;
    }

    return element_quant;
}

//-----------------------------------------------------------------------------

real_t Quantification_Model::transmission(real_t thickness, real_t beta, real_t llambda) const
{
    return std::abs( std::exp( ( (real_t)-4.0 * (real_t)M_PI * thickness * beta / llambda ) ) );
}

//-----------------------------------------------------------------------------

real_t Quantification_Model::absorption(real_t thickness, real_t beta, real_t llambda, real_t shell_factor) const
{
    // shell factor <1 is to determine how much is
    // absorbed by a subshell, and is essentially the
    // ratio of jump factor -1 / jump factor

    return ( 1 - std::abs( std::exp( ((real_t)-4.0 * (real_t)M_PI * thickness * shell_factor * beta / llambda) ) ) );
}

//-----------------------------------------------------------------------------

std::unordered_map<std::string, real_t> Quantification_Model::model_calibrationcurve(std::unordered_map<std::string, Element_Quant> quant_map, real_t p)
{
    // aux_arr[mm, 0] = absorption
    // aux_arr[mm, 1] = transmission, Be
    // aux_arr[mm, 2] = transmission, Ge or Si dead layer
    // aux_arr[mm, 3] = yield
    // aux_arr[mm, 4] = transmission through Si detector
    // aux_arr[mm, 5] = transmission through  air (N2)
//aux array[92, 6] (92 elements)
    //aux_arr = self.aux_arr;
//returns array size 3
    //z_prime is array size 3 of element index of calibraion elements
    //std::vector<real_t> result(aux_arr.size);
    std::unordered_map<std::string, real_t> result_map;
    for(auto& itr : quant_map)
    {
        real_t val = p * itr.second.absorption * itr.second.transmission_Be * itr.second.transmission_Ge * itr.second.yield * ((real_t)1. - itr.second.transmission_through_Si_detector) * itr.second.transmission_through_air;
        result_map.emplace(std::pair<std::string, real_t>(itr.first, val));
    }

    return result_map;
    //return p[0] * aux_arr[z_prime, 0] * aux_arr[z_prime, 1] * aux_arr[z_prime, 2] * aux_arr[z_prime, 3] * ( 1. - aux_arr[z_prime, 4]) * aux_arr[z_prime, 5];
}

//-----------------------------------------------------------------------------

std::vector<real_t> Quantification_Model::model_calibrationcurve(std::vector<Element_Quant> quant_vec, real_t p)
{
    std::vector<real_t> result_vec(quant_vec.size());
    for(size_t i=0; i < quant_vec.size(); i++)
    {
        result_vec[i] = p * quant_vec[i].absorption * quant_vec[i].transmission_Be * quant_vec[i].transmission_Ge * quant_vec[i].yield * ( 1. - quant_vec[i].transmission_through_Si_detector) * quant_vec[i].transmission_through_air;
    }

    return result_vec;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

Electron_Shell get_shell_by_name(std::string element_name)
{
    int idx = element_name.find_last_of("_") + 1;
    std::string shell_type =  element_name.substr(idx);
    if(idx == 0)
    {
        return quantification::models::K_SHELL;
    }
    else
    {
        if(shell_type == "L")
        {
            return quantification::models::L_SHELL;
        }
        if(shell_type == "M")
        {
            return quantification::models::M_SHELL;
        }
    }


    return quantification::models::K_SHELL;
}

} //namespace models
} //namespace quantification
