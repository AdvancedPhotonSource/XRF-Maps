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

template<typename T_real>
Quantification_Model<T_real>::Quantification_Model()
{

}

//-----------------------------------------------------------------------------

template<typename T_real>
Quantification_Model<T_real>::~Quantification_Model()
{

}

//-----------------------------------------------------------------------------

template<typename T_real>
void Quantification_Model<T_real>::init_element_quant(Element_Quant<T_real>& element_quant,
                                            T_real incident_energy,
                                            Element_Info<T_real>* detector_element,
                                            data_struct::Electron_Shell shell,
                                            T_real airpath,
                                            T_real detector_chip_thickness,
                                            T_real beryllium_window_thickness,
                                            T_real germanium_dead_layer,
                                            size_t z_number)
{
    //incident_E(incident_energy) == COHERENT_SCT_ENERGY: Maps_fit_params
    //fit_t_be == BE_WINDOW_THICKNESS * 1000
    //fit_t_ge == GE_DEAD_LAYER * 1000 if detector is not Si
    //add_float['a'] == DET_CHIP_THICKNESS   * 1000

    T_real beta;
    T_real shell_factor = 0.0;
    T_real ev;
    T_real jump_factor = 0.0;
    T_real total_jump_factor = 0.0;

    

    Element_Info<T_real>* element_info = Element_Info_Map<T_real>::inst()->get_element(z_number);
    if(element_info == nullptr)
    {
        return;
    }

    element_quant.name = element_info->name;
    element_quant.Z = z_number;

    switch(shell)
    {
    case data_struct::Electron_Shell::K_SHELL:
        ev = element_info->xrf.at("ka1") * (T_real)1000.0;
        element_quant.yield = element_info->yieldD.at("k");
        if( incident_energy > element_info->bindingE.at("K") )
        {
            jump_factor = element_info->jump.at("K");
        }
        break;
    case data_struct::Electron_Shell::L_SHELL:
        ev = element_info->xrf.at("la1") * (T_real)1000.0;
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
    case data_struct::Electron_Shell::M_SHELL:
        ev = element_info->xrf.at("ma1") * (T_real)1000.0;
        element_quant.yield = element_info->xrf_abs_yield.at("ma1");
        jump_factor = element_info->jump.at("M5");
        total_jump_factor = element_info->jump.at("M1") * element_info->jump.at("M2") * element_info->jump.at("M3") * element_info->jump.at("M4");
        break;
    default:
        throw "Unsupported shell. Currently only support for K, L, and M ";
        break;
    }

    ////aux_arr[mm, 3] = yieldd
 //element_quant.yield = element_info->yieldD["K"]; //yieldd === newrel_yield * info_elements[element_temp].yieldD['k']

    if( jump_factor != 0.0 )
    {
        if( total_jump_factor == 0.0 )
        {
            shell_factor = (jump_factor - (T_real)1.0) / jump_factor;
        }
        else
        {
            shell_factor = (jump_factor - (T_real)1.0) / jump_factor / total_jump_factor;
        }
    }
    // replace straight henke routines, with those
    // that take the absorption edges into account
    // make sure we are a bit above the absorption edge to make sure that for calibration purposes we do not eoncouner any weird things.
    beta = element_info->calc_beta(element_info->density, (incident_energy + (T_real)0.1) * (T_real)1000.0);
    // stds in microgram/cm2
    // density rho = g/cm3 = 1 microgram/cm2 /1000/1000/cm = 1 microgram/cm2 /1000/1000/*10*1000/um = 1 microgram/cm2 /100/um
    // thickness for 1 ugr/cm2
    // =1/(density[g/cm3]/10)

    T_real thickness = (T_real)1.0 / (element_info->density * (T_real)10.0) * (T_real)1000.0;
    ////aux_arr[mm, 0] = self.absorption(thickness, beta, 1239.852/((self.maps_conf.incident_E+0.1)*1000.), shell_factor=shell_factor)
    element_quant.absorption = absorption(thickness, beta, (T_real)1239.852 / ((incident_energy + (T_real)0.1) * (T_real)1000.0), shell_factor);

    if( false == std::isfinite(element_quant.absorption) )
    {
        logW<<element_quant.name<<" absorption =  "<<element_quant.absorption<<"; Setting to 0.0\n";
        element_quant.absorption = 0.0;
    }
    if(ev > 0)
    {
        beta  = Element_Info_Map<T_real>::inst()->calc_beta("Be", (T_real)1.848, ev);
        ////aux_arr[mm, 1] = self.transmission(self.maps_conf.fit_t_be, beta, 1239.852/ev)
        element_quant.transmission_Be = transmission(beryllium_window_thickness, beta, (T_real)1239.852 / ev);
        // get density (second arg) from henke lookup
        beta  = Element_Info_Map<T_real>::inst()->calc_beta("Ge", (T_real)5.323, ev);
        ////aux_arr[mm, 2] = self.transmission(self.maps_conf.fit_t_ge, beta, 1239.852/ev)
        element_quant.transmission_Ge = transmission(germanium_dead_layer, beta, (T_real)1239.852 / ev);
    }
 

    if (detector_element->name.length() > 0 && data_struct::Henke_Compound_Density_Map<T_real>.count(detector_element->name) > 0 && detector_chip_thickness > 0.0 && ev > 0) //  (self.maps_conf.add_long['a'] == 1)
    {
        beta  = Element_Info_Map<T_real>::inst()->calc_beta(detector_element->name, data_struct::Henke_Compound_Density_Map<T_real>.at(detector_element->name), ev);
        element_quant.transmission_through_Si_detector = transmission(detector_chip_thickness, beta, (T_real)1239.852 / ev);
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
        //T_real density = (T_real)0.00117;
        T_real density = (T_real)1.2047e-3;
        //air_ele = 'N78.08O20.95Ar0.93'
        //density = 1.2047e-3
        //f1, f2, delta, beta, graze_mrad, reflect, inverse_mu, atwt = Chenke.get_henke_single('air', density, ev)
        beta = Element_Info_Map<T_real>::inst()->calc_beta("N:78.08,O:20.95,Ar:0.93", density, ev); // air
        ////aux_arr[mm, 5] = self.transmission(airpath*1000., beta, 1239.852/ev)  // airpath is read in microns, transmission function expects nm
        // airpath is read in microns, transmission function expects nm
        element_quant.transmission_through_air = transmission( airpath , beta, (T_real)1239.852 / ev);
    }
    else
    {
        ////aux_arr[mm, 5] = 1.
        element_quant.transmission_through_air = 1.0;
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
T_real Quantification_Model<T_real>::transmission(T_real thickness, T_real beta, T_real llambda) const
{
    return std::abs( std::exp( ( (T_real)-4.0 * (T_real)M_PI * thickness * beta / llambda ) ) );
}

//-----------------------------------------------------------------------------

template<typename T_real>
T_real Quantification_Model<T_real>::absorption(T_real thickness, T_real beta, T_real llambda, T_real shell_factor) const
{
    // shell factor <1 is to determine how much is
    // absorbed by a subshell, and is essentially the
    // ratio of jump factor -1 / jump factor

    return ( 1 - std::abs( std::exp( ((T_real)-4.0 * (T_real)M_PI * thickness * shell_factor * beta / llambda) ) ) );
}

//-----------------------------------------------------------------------------

template<typename T_real>
std::unordered_map<std::string, T_real> Quantification_Model<T_real>::model_calibrationcurve(std::unordered_map<std::string, Element_Quant<T_real>> quant_map, T_real p)
{
    // aux_arr[mm, 0] = absorption
    // aux_arr[mm, 1] = transmission, Be
    // aux_arr[mm, 2] = transmission, detector dead layer
    // aux_arr[mm, 3] = yield
    // aux_arr[mm, 4] = transmission through detector element
    // aux_arr[mm, 5] = transmission through  air (N2)
//aux array[92, 6] (92 elements)
    //aux_arr = self.aux_arr;
//returns array size 3
    //z_prime is array size 3 of element index of calibraion elements
    //std::vector<T_real> result(aux_arr.size);
    std::unordered_map<std::string, T_real> result_map;
    for(auto& itr : quant_map)
    {
        T_real val = p * itr.second.absorption * itr.second.transmission_Be * itr.second.transmission_Ge * itr.second.yield * ((T_real)1. - itr.second.transmission_through_Si_detector) * itr.second.transmission_through_air;
        
        if (false == std::isfinite(val))
        {
            logW << itr.second.name<<" val = " << val << 
            "\np = "<<p<<
            "\nabsorption = "<< itr.second.absorption <<
            "\ntransmission_Be = "<<itr.second.transmission_Be <<
            "\ntransmission_Ge = "<<itr.second.transmission_Ge <<
            "\nyield = "<<itr.second.yield <<
            "\ntransmission_through_Si_detector = "<<itr.second.transmission_through_Si_detector <<
            "\ntransmission_through_air = "<<itr.second.transmission_through_air << "\n";
            val = std::numeric_limits<T_real>::max() * .9;  // 90% of max real
        }
        
        result_map.emplace(std::pair<std::string, T_real>(itr.first, val));
    }

    return result_map;
    //return p[0] * aux_arr[z_prime, 0] * aux_arr[z_prime, 1] * aux_arr[z_prime, 2] * aux_arr[z_prime, 3] * ( 1. - aux_arr[z_prime, 4]) * aux_arr[z_prime, 5];
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Quantification_Model<T_real>::model_calibrationcurve(std::vector<Element_Quant<T_real>> *quant_vec, T_real p)
{
    for(auto &itr : *quant_vec)
    {
        T_real val = p * itr.absorption * itr.transmission_Be * itr.transmission_Ge * itr.yield * ((T_real)1. - itr.transmission_through_Si_detector) * itr.transmission_through_air;
        
        if (false == std::isfinite(val))
        {
            logW << itr.name<<" val = " << val << 
            "\np = "<<p<<
            "\nabsorption = "<< itr.absorption <<
            "\ntransmission_Be = "<<itr.transmission_Be <<
            "\ntransmission_Ge = "<<itr.transmission_Ge <<
            "\nyield = "<<itr.yield <<
            "\ntransmission_through_Si_detector = "<<itr.transmission_through_Si_detector <<
            "\ntransmission_through_air = "<<itr.transmission_through_air << "\n";
            val = std::numeric_limits<T_real>::max() * .9;  // 90% of max real
        }

        itr.calib_curve_val = val;
    }
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


TEMPLATE_CLASS_DLL_EXPORT Quantification_Model<float>;
TEMPLATE_CLASS_DLL_EXPORT Quantification_Model<double>;

} //namespace models
} //namespace quantification
