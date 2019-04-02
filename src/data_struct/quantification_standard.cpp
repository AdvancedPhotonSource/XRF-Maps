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

#include "quantification/models/quantification_model.h"

namespace data_struct
{

//-----------------------------------------------------------------------------

Quantification_Standard::Quantification_Standard()
{

    _sr_current = 1.0;
    _US_IC = 1.0;
    _DS_IC = 1.0;
    _processed = false;

}

//-----------------------------------------------------------------------------

Quantification_Standard::~Quantification_Standard()
{

}

//-----------------------------------------------------------------------------

void Quantification_Standard::append_element(string proc_type, string name, real_t weight)
{
	if (_element_quants.count(proc_type) > 0)
	{
		_element_quants[proc_type][name] = Element_Quant(weight);
	}
	else
	{
		unordered_map<string, Element_Quant> e_map;
		e_map[name] = Element_Quant(weight);
		_element_quants[proc_type] = e_map;
	}
}

//-----------------------------------------------------------------------------

void Quantification_Standard::element_counts(string proc_type_str, unordered_map<string, real_t> map)
{
    _element_counts[proc_type_str].clear();
    for(auto itr: map)
    {
          _element_counts[proc_type_str][itr.first] = itr.second;
    }
}

//-----------------------------------------------------------------------------

std::string get_shell_element_label(int shell, size_t l)
{
    std::string shell_str = "";
    switch(shell)
    {
    case quantification::models::K_SHELL:
        //shell_str = "_K";
        break;
    case quantification::models::L_SHELL:
        shell_str = "_L";
        break;
    case quantification::models::M_SHELL:
        shell_str = "_M";
        break;
    case quantification::models::N_SHELL:
        shell_str = "_N";
        break;
    case quantification::models::O_SHELL:
        shell_str = "_O";
        break;
    case quantification::models::P_SHELL:
        shell_str = "_P";
        break;
    case quantification::models::Q_SHELL:
        shell_str = "_Q";
        break;
    default:
        break;
    }

    return Element_Symbols[l]+shell_str;
}

bool Quantification_Standard::quantifiy(fitting::optimizers::Optimizer * optimizer,
                                        string proc_type_str,
                                        unordered_map<string, real_t>  *element_counts,
                                        real_t incident_energy,
                                        Element_Info* detector_element,
                                        bool airpath,
                                        real_t detector_chip_thickness,
                                        real_t beryllium_window_thickness,
                                        real_t germanium_dead_layer)
{

    quantification::models::Quantification_Model quantification_model;
    vector<quantification::models::Electron_Shell> shells_to_quant = {quantification::models::K_SHELL, quantification::models::L_SHELL, quantification::models::M_SHELL};
    unordered_map<size_t, real_t*> quant_list =
    {
        {Quantifiers::CURRENT, &_sr_current},
        {Quantifiers::US_IC, &_US_IC},
        {Quantifiers::DS_IC, &_DS_IC}
    };
    /*{
        {Quantifiers::CURRENT, 101.94},
        {Quantifiers::US_IC, 268303.0},
        {Quantifiers::DS_IC, 134818.0}
        };
*/
    _element_counts[proc_type_str] = *element_counts;

    calibration_curves.emplace(pair<string, Quantifiers>(proc_type_str, Quantifiers(93)) );
    Quantifiers *quantifiers = &calibration_curves.at(proc_type_str);

    for (auto& quant_itr : quant_list)
    {

        for(auto &itr : _element_quants[proc_type_str])
        {

            //TODO: Parse the name to see if it is K, L, or M
            quantification::models::Electron_Shell shell = quantification::models::K_SHELL;

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

            //factor = quant_itr.second;
            real_t e_cal_factor = (itr.second.weight * (*quant_itr.second));
            real_t e_cal = e_cal_factor / element_counts->at(itr.first);
            itr.second.e_cal_ratio = (real_t)1.0 / e_cal;
            //initial guess: parinfo_value[0] = 100000.0 / factor

        }

        Fit_Parameters fit_params;
        fit_params.add_parameter(Fit_Param("quantifier", 0.0, 0.0, 1.0, 0.001, E_Bound_Type::FIT));
        fit_params["quantifier"].value = (real_t)100000.0 / (*quant_itr.second);
        optimizer->minimize_quantification(&fit_params, &_element_quants[proc_type_str], &quantification_model);
        real_t val = fit_params["quantifier"].value;

        for(auto shell : shells_to_quant)
        {
            //generate calibration curve for all elements
            vector<Element_Quant> element_quant_vec = quantification_model.generate_quant_vec(incident_energy,
                                                                                           detector_element,
                                                                                           shell,
                                                                                           airpath,
                                                                                           detector_chip_thickness,
                                                                                           beryllium_window_thickness,
                                                                                           germanium_dead_layer,
                                                                                           1,
                                                                                           92);

            quantifiers->calib_curves[quant_itr.first].shell_curves[shell] = quantification_model.model_calibrationcurve(element_quant_vec, val);
            //change nan's to zeros
            quantifiers->calib_curves[quant_itr.first].shell_curves_labels[shell].resize(quantifiers->calib_curves[quant_itr.first].shell_curves[shell].size());
            for(size_t l = 0; l<quantifiers->calib_curves[quant_itr.first].shell_curves[shell].size(); l++)
            {
                quantifiers->calib_curves[quant_itr.first].shell_curves_labels[shell][l] = get_shell_element_label(shell, l);
                if( std::isnan(quantifiers->calib_curves[quant_itr.first].shell_curves[shell][l]) )
                {
                    quantifiers->calib_curves[quant_itr.first].shell_curves[shell][l] = 0.0;
                }
            }
        }

    }
    _processed = true;

    return true;

}

//-----------------------------------------------------------------------------

} //namespace data_struct
