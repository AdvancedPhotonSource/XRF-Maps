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



#include "quantification_standard.h"

#include <string>
#include <cmath>


namespace data_struct
{
namespace xrf
{

//-----------------------------------------------------------------------------

Quantification_Standard::Quantification_Standard()
{

    _sr_current = 1.0;
    _IC_US = 1.0;
    _IC_DS = 1.0;

}

//-----------------------------------------------------------------------------

Quantification_Standard::~Quantification_Standard()
{

}

//-----------------------------------------------------------------------------

void Quantification_Standard::append_element(std::string name, real_t weight)
{
    _element_quants.emplace(std::pair<std::string, Element_Quant>(name, Element_Quant(weight)));
}

//-----------------------------------------------------------------------------

std::unordered_map<std::string, Element_Quant> Quantification_Standard::generate_quant_map(real_t incident_energy, Element_Info* detector_element)
{

    std::unordered_map<std::string, Element_Quant> element_quant_map;
    /*
    //incident_E(incident_energy) == COHERENT_SCT_ENERGY: Maps_fit_params
    //fit_t_be == BE_WINDOW_THICKNESS * 1000
    //fit_t_ge == GE_DEAD_LAYER * 1000 if detector is not Si
    //add_float['a'] == DET_CHIP_THICKNESS   * 1000
    for (int i=0; i < 92; i ++)
    {
        real_t beta;
        real_t shell_factor;
        Element_Quant element_quant;
        Element_Info* element_info = Element_Info_Map::inst()->get_element(i);
        if(element_info == nullptr)
        {
            continue;
        }
        // replace straight henke routines, with those
        // that take the absorption edges into account
        // make sure we are a bit above the absorption edge to make sure that for calibration purposes we do not eoncouner any weird things.
        beta = element_info->calc_beta(element_info->density, (incident_energy + 0.1)*1000.0);
        // stds in microgram/cm2
        // density rho = g/cm3 = 1 microgram/cm2 /1000/1000/cm = 1 microgram/cm2 /1000/1000/*10*1000/um = 1 microgram/cm2 /100/um
        // thickness for 1 ugr/cm2
        // =1/(density[g/cm3]/10)

        real_t thickness = 1.0 / (element_info->density * 10.0) * 1000.0;
        ////aux_arr[mm, 0] = self.absorption(thickness, beta, 1239.852/((self.maps_conf.incident_E+0.1)*1000.), shell_factor=shell_factor)
        element_quant.absorption = absorption(thickness, beta, 1239.852 / ((incident_energy + 0.1) * 1000.0), shell_factor);

        beta  = Element_Info_Map::inst()->calc_beta('Be', 1.848, ev);
        ////aux_arr[mm, 1] = self.transmission(self.maps_conf.fit_t_be, beta, 1239.852/ev)
        element_quant.transmission_Be = transmission(self.maps_conf.fit_t_be, beta, 1239.852 / ev);

        beta  = Element_Info_Map::inst()->calc_beta( 'Ge', 5.323, ev);
        ////aux_arr[mm, 2] = self.transmission(self.maps_conf.fit_t_ge, beta, 1239.852/ev)
        element_quant.transmission_Ge = transmission(self.maps_conf.fit_t_be, beta, 1239.852 / ev);

        ////aux_arr[mm, 3] = yieldd
        element_quant.yield = yieldd; //yieldd === newrel_yield * info_elements[element_temp].yieldD['k']

        if (self.maps_conf.add_long['a'] == 1)
        {
            beta  = Element_Info_Map::inst()->calc_beta( 'Si', 2.3, ev);
        }
        ////aux_arr[mm, 4] = self.transmission(self.maps_conf.add_float['a'], beta, 1239.852/ev)
        element_quant.transmission_through_Si_detector = self.transmission(self.maps_conf.add_float['a'], beta, 1239.852 / ev);
        if ( (self.maps_conf.add_float['a'] == 0.) || (self.maps_conf.add_long['a'] != 1) )
        {
            ////aux_arr[mm, 4] = 0.
            element_quant.transmission_through_Si_detector = 0.0;
        }
        if( airpath > 0)
        {
            //density = 1.0
            density = 0.00117;
            //air_ele = 'N78.08O20.95Ar0.93'
            //density = 1.2047e-3
            //f1, f2, delta, beta, graze_mrad, reflect, inverse_mu, atwt = Chenke.get_henke_single('air', density, ev)
            beta = Element_Info_Map::inst()->calc_beta('air', density, ev);
            ////aux_arr[mm, 5] = self.transmission(airpath*1000., beta, 1239.852/ev)  // airpath is read in microns, transmission function expects nm
            // airpath is read in microns, transmission function expects nm
            element_quant.transmission_through_Si_detector = transmission( airpath * 1000.0, beta, 1239.852 / ev);
        }
        else
        {
            ////aux_arr[mm, 5] = 1.
            element_quant.transmission_through_air = 1.0;
        }
        element_quant_map.emplace( std::pair <std::string, Element_Quant>(element_info->name, element_quant) );
    }
*/
    return element_quant_map;
}

//-----------------------------------------------------------------------------

real_t Quantification_Standard::transmission(real_t thickness, real_t beta, real_t llambda) const
{
    return std::abs( std::exp( ( -4.0 * M_PI * thickness * beta / llambda ) ) );
}

//-----------------------------------------------------------------------------

real_t Quantification_Standard::absorption(real_t thickness, real_t beta, real_t llambda, real_t shell_factor) const
{
    // make sure shell_factor is defined, and if not, set it to 1
    // shell factor <1 is to determine how much is
    // absorbed by a subshell, and is essentially the
    // ratio of jump factor -1 / jump factor
    if (shell_factor == 0)
        shell_factor = 1;

    return ( 1 - std::abs( std::exp( (-4.0 * M_PI * thickness * shell_factor * beta / llambda) ) ) );
}

//-----------------------------------------------------------------------------

std::unordered_map<std::string, real_t> Quantification_Standard::fit_calibrationcurve(std::unordered_map<std::string, Element_Quant> quant_map, real_t p)
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
        real_t val = p * itr.second.absorption * itr.second.transmission_Be * itr.second.transmission_Ge * itr.second.yield * ( 1. - itr.second.transmission_through_Si_detector) * itr.second.transmission_through_air;
        result_map.emplace(std::pair<std::string, real_t>(itr.first, val));
    }

    return result_map;
    //return p[0] * aux_arr[z_prime, 0] * aux_arr[z_prime, 1] * aux_arr[z_prime, 2] * aux_arr[z_prime, 3] * ( 1. - aux_arr[z_prime, 4]) * aux_arr[z_prime, 5];
}
/*
//-----------------------------------------------------------------------------

void Quantification_Standard::residuals(std::valarray<real_t> p, std::valarray<real_t> y, std::valarray<real_t> x) const
{
    //y is array of elements standards
    //x is indexes of elements in standard
    //p is array size 2 but seems only first index is used
    return (y - this->fit_calibrationcurve(x, p));
}
*/
//-----------------------------------------------------------------------------

} //namespace xrf
} //namespace data_struct
