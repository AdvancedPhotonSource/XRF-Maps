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

/// Initial Author <2017>: Arthur Glowacki



#include "data_struct/detector.h"

namespace data_struct
{

Detector::Detector(unsigned int number)
{
    model = nullptr;
    detector_element = data_struct::Element_Info_Map::inst()->get_element("Si");
    beryllium_window_thickness = 0.0; //24.000000
    germanium_dead_layer = 0.0; //350.0
    detector_chip_thickness = 0.0;
    incident_energy = 10.0;
    airpath = 0.0;
    _number = number;
}

Detector::~Detector()
{
    if (model != nullptr)
    {
        delete model;
        model = nullptr;
    }
    for (auto& itr : fit_routines)
    {
        fitting::routines::Base_Fit_Routine* fit_routine = itr.second;
        if (fit_routine != nullptr)
        {
            delete fit_routine;
        }
    }
    fit_routines.clear();
    quantification_standards.clear();
    for (auto& itr : all_element_quants)
    {
        for (auto& itr2 : itr.second)
        {
            itr2.second.clear();
        }
        itr.second.clear();
    }
    all_element_quants.clear();
}

//-----------------------------------------------------------------------------

void Detector::append_element(Fitting_Routines routine, string quant_scaler, string name, real_t weight)
{
    if (fitting_quant_map.count(routine) == 0)
    {
        fitting_quant_map[routine] = Fitting_Quantification_Struct();
    }

    Element_Info* element = Element_Info_Map::inst()->get_element(name);
    if (element != nullptr)
    {
        Electron_Shell shell = get_shell_by_name(name);
        fitting_quant_map.at(routine).update_weight(shell, element->number, weight);

        if (fitting_quant_map.at(routine).quant_scaler_map.count(quant_scaler) > 0)
        {
            all_element_quants[routine][quant_scaler][name] = &(fitting_quant_map.at(routine).quant_scaler_map.at(quant_scaler).curve_quant_map[shell][element->number - 1]);
        }
    }
    else
    {
        logW << "Could not add element " << name << ". Not found in Element_Info_Map\n";
    }
}

//-----------------------------------------------------------------------------

void Detector::update_element_quants(Fitting_Routines routine,
                                    string quantifier_scaler,
                                    Quantification_Standard* standard,
                                    Quantification_Model* quantification_model,
                                    real_t ic_quantifier)
{
    if (fitting_quant_map.count(routine) > 0)
    {
        if (fitting_quant_map.at(routine).quant_scaler_map.count(quantifier_scaler) > 0)
        {
            for (const auto& shell_itr : Shells_Quant_List)
            {
                for (auto& eq_itr : fitting_quant_map.at(routine).quant_scaler_map.at(quantifier_scaler).curve_quant_map.at(shell_itr))
                {
                    Element_Info* element = Element_Info_Map::inst()->get_element(eq_itr.Z);
                    if (element == nullptr)
                    {
                        continue;
                    }
                    quantification_model->init_element_quant(eq_itr,
                        incident_energy,
                        detector_element,
                        shell_itr,
                        airpath,
                        detector_chip_thickness,
                        beryllium_window_thickness,
                        germanium_dead_layer,
                        element->number);

                    // if we have weight for this element, update e_cal_ratio
                    if (eq_itr.weight > 0.0)
                    {
                        string name = eq_itr.name;
                        if (shell_itr == Electron_Shell::L_SHELL)
                        {
                            name += "_L";
                        }
                        if (shell_itr == Electron_Shell::M_SHELL)
                        {
                            name += "_M";
                        }
                        // with 2 standards we can have weights from another standard so we have to check if we have counts
                        if (standard->element_counts.at(routine).count(name) > 0)
                        {
                            if (ic_quantifier == 0.0)
                            {
                                ic_quantifier = 1.0;
                            }

                            real_t counts = standard->element_counts.at(routine).at(name);
                            real_t e_cal_factor = (eq_itr.weight * (ic_quantifier));
                            if (counts > 0.)
                            {
                                real_t e_cal = e_cal_factor / counts;
                                if (eq_itr.e_cal_ratio == 0.)
                                {
                                    eq_itr.e_cal_ratio = (real_t)1.0 / e_cal;
                                }
                                else // else avg the two
                                {
                                    real_t second_cal_ratio = (real_t)1.0 / e_cal;
                                    eq_itr.e_cal_ratio += second_cal_ratio;
                                    eq_itr.e_cal_ratio *= 0.5;
                                }
                            }
                            else
                            {
                                if (eq_itr.e_cal_ratio == 0.)
                                {
                                    eq_itr.e_cal_ratio = 1.0e-10;
                                }
                            }
                        }
                    }
                }
            }
        }
        else
        {
            logW << "Could not find quantifier scalers : " << quantifier_scaler << " .\n";
        }
    }
    else
    {
        logW << "Could not find fitting routine " << Fitting_Routine_To_Str.at(routine) << " .\n";
    }
}

//-----------------------------------------------------------------------------

void Detector::generage_avg_quantification_scalers()
{
    real_t avg_sr_current = 0.0;
    real_t avg_US_IC = 0.0;
    real_t avg_DS_IC = 0.0;

    real_t crnt_cnt = 0.0;
    real_t us_cnt = 0.0;
    real_t ds_cnt = 0.0;

    //average quantification scalers
    for (const auto& itr : quantification_standards)
    {
        if (itr.second.sr_current > 0.0)
        {
            avg_sr_current += itr.second.sr_current;
            crnt_cnt += 1.0;
        }
        if (itr.second.US_IC > 0.0)
        {
            avg_US_IC += itr.second.US_IC;
            us_cnt += 1.0;
        }
        if (itr.second.DS_IC > 0.0)
        {
            avg_DS_IC += itr.second.DS_IC;
            ds_cnt += 1.0;
        }
    }
    
    if (avg_sr_current == 0.0 && avg_US_IC == 0.0 && avg_DS_IC == 0.0)
    {
        logE << "Could not find SR_Current, US_IC, and DS_IC. Not going to perform quantification\n";
    }

    if (crnt_cnt != us_cnt && crnt_cnt != ds_cnt)
    {
        logE << "Averaging encounded one of the values to be 0 so it will not be correct!\n";
    }


    if (avg_sr_current == 0.0)
    {
        logW"SR_Current is 0. Probably couldn't find it in the dataset. Setting it to 1. Quantification will be incorrect.\n";
        avg_sr_current = 100.0;
    }
    else
    {
        avg_sr_current /= crnt_cnt;
    }
    if (avg_US_IC == 0.0)
    {
        logW"US_IC is 0. Probably couldn't find it in the dataset. Setting it to 1. Quantification will be incorrect.\n";
        avg_US_IC = 1.0;
    }
    else
    {
        avg_US_IC /= us_cnt;
    }
    if (avg_DS_IC == 0.0)
    {
        logW"DS_IC is 0. Probably couldn't find it in the dataset. Setting it to 1. Quantification will be incorrect.\n";
        avg_DS_IC = 1.0;
    }
    else
    {
        avg_DS_IC /= ds_cnt;
    }

    avg_quantification_scaler_map[STR_SR_CURRENT] = avg_sr_current;
    avg_quantification_scaler_map[STR_US_IC] = avg_US_IC;
    avg_quantification_scaler_map[STR_DS_IC] = avg_DS_IC;
}

//-----------------------------------------------------------------------------

void Detector::update_calibration_curve(Fitting_Routines routine,
                                        string quantifier_scaler,
                                        Quantification_Model* quantification_model,
                                        real_t val)
{
    if (fitting_quant_map.count(routine) > 0)
    {
        if (fitting_quant_map.at(routine).quant_scaler_map.count(quantifier_scaler) > 0)
        {
            for (const auto& shell_itr : Shells_Quant_List)
            {
                vector<Element_Quant>* quant_vec = &(fitting_quant_map.at(routine).quant_scaler_map.at(quantifier_scaler).curve_quant_map.at(shell_itr));
                quantification_model->model_calibrationcurve(quant_vec, val);
            }
        }
    }
}

//-----------------------------------------------------------------------------

void Detector::update_from_fit_paramseters()
{
    //Parameters for calibration curve
    if (fit_params_override_dict.detector_element.length() > 0)
    {
        // Get the element info class                                           // detector element as string "Si" or "Ge" usually
        detector_element = (data_struct::Element_Info_Map::inst()->get_element(fit_params_override_dict.detector_element));
    }
    if (fit_params_override_dict.be_window_thickness.length() > 0)
    {
        beryllium_window_thickness = std::stof(fit_params_override_dict.be_window_thickness) * 1000.0;
    }
    if (fit_params_override_dict.ge_dead_layer.length() > 0)
    {
        germanium_dead_layer = std::stof(fit_params_override_dict.ge_dead_layer) * 1000.0;
    }
    if (fit_params_override_dict.det_chip_thickness.length() > 0)
    {
        detector_chip_thickness = std::stof(fit_params_override_dict.det_chip_thickness) * 1000.0;
    }
    if (fit_params_override_dict.airpath.length() > 0)
    {
        airpath = std::stof(fit_params_override_dict.airpath) * 1000.0;
    }
    if (fit_params_override_dict.fit_params.contains(STR_COHERENT_SCT_ENERGY))
    {
        incident_energy = (fit_params_override_dict.fit_params.at(STR_COHERENT_SCT_ENERGY).value);
    }
}

//-----------------------------------------------------------------------------

} //namespace data_struct
