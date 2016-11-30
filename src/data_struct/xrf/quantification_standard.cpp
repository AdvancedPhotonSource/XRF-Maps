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

#include "quantification_model.h"

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

bool Quantification_Standard::quantifiy(fitting::optimizers::Optimizer * optimizer,
                                        data_struct::xrf::Fit_Count_Dict  *element_counts,
                                        real_t incident_energy,
                                        Element_Info* detector_element,
                                        bool airpath,
                                        real_t detector_chip_thickness,
                                        real_t beryllium_window_thickness,
                                        real_t germanium_dead_layer)
{

    quantification::models::Quantification_Model quantification_model;
    std::vector<quantification::models::Electron_Shell> shells_to_quant = {quantification::models::K_SHELL, quantification::models::L_SHELL, quantification::models::M_SHELL};
    std::unordered_map<std::string, real_t> quant_list = {std::pair<std::string, real_t>("current", 10.0), std::pair<std::string, real_t>("us_ic", 1.0), std::pair<std::string, real_t>("ds_ic", 1.0) };

    for(auto shell : shells_to_quant)
    {
        for(auto &itr : _element_quants)
        {

            Element_Info* element = Element_Info_Map::inst()->get_element(itr.first);
            if (element == nullptr)
            {
                continue;
            }
            Element_Quant element_quant = quantification_model.generate_element_quant(incident_energy,
                                                                               detector_element,
                                                                               shell,
                                                                               airpath,
                                                                               detector_chip_thickness,
                                                                               beryllium_window_thickness,
                                                                               germanium_dead_layer,
                                                                               element->number);
            itr.second.absorption = element_quant.absorption;
            itr.second.transmission_Be = element_quant.transmission_Be;
            itr.second.transmission_Ge = element_quant.transmission_Ge;
            itr.second.transmission_through_air = element_quant.transmission_through_air;
            itr.second.transmission_through_Si_detector = element_quant.transmission_through_Si_detector;
            itr.second.yield = element_quant.yield;

        }


        for (auto& quant_itr : quant_list)
        {

            for(auto& element_itr : _element_quants)
            {
                //factor = quant_itr.second;
                real_t e_cal_factor = (element_itr.second.weight * quant_itr.second);
                real_t e_cal = e_cal_factor / element_counts->at(element_itr.first)[0][0];
                element_itr.second.e_cal_ratio = 1.0 / e_cal;
                //initial guess: parinfo_value[0] = 100000.0 / factor
            }
            Fit_Parameters fit_params;
            fit_params.add_parameter("quantifier", Fit_Param("quantifier", 0.0, 0.0, 1.0, 0.001, FIT));
            fit_params["quantifier"].value = 100000.0 / quant_itr.second;
            optimizer->minimize_quantification(&fit_params, &_element_quants, &quantification_model);
            quant_itr.second = fit_params["quantifier"].value;
        }

    }

//    if l == 0:
//        factor = std.current
//    if l == 1:
//        factor = std.us_ic
//    if l == 2:
//        factor = std.ds_ic
//   e_cal_factor = (weight_ugr_cm * factor)
//   e_cal = e_cal_factor / counts
//   y = 1.0 / e_cal
//initial guess: parinfo_value[0] = 100000.0 / factor
    //optimizer->minimize_quantification(_element_quants, &quantification_model);

    return false;

}

//-----------------------------------------------------------------------------

} //namespace xrf
} //namespace data_struct
