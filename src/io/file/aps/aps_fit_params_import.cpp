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


#include "aps_fit_params_import.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>
#include <chrono>
#include <ctime>


namespace io
{
namespace file
{
namespace aps
{

/***
 * Translation map from APS file tags to internal tags
 */
const std::unordered_map<std::string, std::string> FILE_TAGS_TRANSLATION = {
	{"CAL_OFFSET_[E_OFFSET]", STR_ENERGY_OFFSET},
    {"CAL_OFFSET_[E_OFFSET]_MAX", STR_ENERGY_OFFSET},
    {"CAL_OFFSET_[E_OFFSET]_MIN", STR_ENERGY_OFFSET},
    {"CAL_SLOPE_[E_LINEAR]", STR_ENERGY_SLOPE},
    {"CAL_SLOPE_[E_LINEAR]_MAX", STR_ENERGY_SLOPE},
    {"CAL_SLOPE_[E_LINEAR]_MIN", STR_ENERGY_SLOPE},
    {"CAL_QUAD_[E_QUADRATIC]", STR_ENERGY_QUADRATIC},
    {"CAL_QUAD_[E_QUADRATIC]_MAX", STR_ENERGY_QUADRATIC},
    {"CAL_QUAD_[E_QUADRATIC]_MIN", STR_ENERGY_QUADRATIC},
    {"FWHM_OFFSET", STR_FWHM_OFFSET},
    {"FWHM_FANOPRIME", STR_FWHM_FANOPRIME},
    {"COHERENT_SCT_ENERGY", STR_COHERENT_SCT_ENERGY},
    {"COHERENT_SCT_ENERGY_MAX", STR_COHERENT_SCT_ENERGY},
    {"COHERENT_SCT_ENERGY_MIN", STR_COHERENT_SCT_ENERGY},
    {"COMPTON_ANGLE", STR_COMPTON_ANGLE},
    {"COMPTON_ANGLE_MAX", STR_COMPTON_ANGLE},
    {"COMPTON_ANGLE_MIN", STR_COMPTON_ANGLE},
    {"COMPTON_FWHM_CORR", STR_COMPTON_FWHM_CORR},
    {"COMPTON_STEP", STR_COMPTON_F_STEP},
    {"COMPTON_F_TAIL", STR_COMPTON_F_TAIL},
    {"COMPTON_GAMMA", STR_COMPTON_GAMMA},
    {"COMPTON_HI_F_TAIL", STR_COMPTON_HI_F_TAIL},
    {"COMPTON_HI_GAMMA", STR_COMPTON_HI_GAMMA},
    {"STEP_OFFSET", STR_F_STEP_OFFSET},
    {"STEP_LINEAR", STR_F_STEP_LINEAR},
    {"STEP_QUADRATIC", STR_F_STEP_QUADRATIC},
    {"F_TAIL_OFFSET", STR_F_TAIL_OFFSET},
    {"F_TAIL_LINEAR", STR_F_TAIL_LINEAR},
    {"F_TAIL_QUADRATIC", STR_F_TAIL_QUADRATIC},
    {"KB_F_TAIL_OFFSET", STR_KB_F_TAIL_OFFSET},
    {"KB_F_TAIL_LINEAR", STR_KB_F_TAIL_LINEAR},
    {"KB_F_TAIL_QUADRATIC", STR_KB_F_TAIL_QUADRATIC},
    {"GAMMA_OFFSET", STR_GAMMA_OFFSET},
    {"GAMMA_LINEAR", STR_GAMMA_LINEAR},
    {"GAMMA_QUADRATIC", STR_GAMMA_QUADRATIC},
    {"SNIP_WIDTH", STR_SNIP_WIDTH},
    {"SI_ESCAPE_FACTOR", STR_SI_ESCAPE},
    {"GE_ESCAPE_FACTOR", STR_GE_ESCAPE},
    {"ESCAPE_LINEAR", STR_ESCAPE_LINEAR},
	{"MAX_ENERGY_TO_FIT", STR_MAX_ENERGY_TO_FIT },
	{"MIN_ENERGY_TO_FIT", STR_MIN_ENERGY_TO_FIT }
};

const std::vector<std::string> Updatable_TAGS = {
                                        "CAL_OFFSET_[E_OFFSET]",
                                        "CAL_SLOPE_[E_LINEAR]",
                                        "CAL_QUAD_[E_QUADRATIC]",
                                        "FWHM_OFFSET",
                                        "FWHM_FANOPRIME",
                                        "COHERENT_SCT_ENERGY",
                                        "COMPTON_ANGLE",
                                        "COMPTON_FWHM_CORR",
                                        "COMPTON_STEP",
                                        "COMPTON_F_TAIL",
                                        "COMPTON_GAMMA",
                                        "COMPTON_HI_F_TAIL",
                                        "COMPTON_HI_GAMMA",
                                        "STEP_OFFSET",
                                        "STEP_LINEAR",
                                        "STEP_QUADRATIC",
                                        "F_TAIL_OFFSET",
                                        "F_TAIL_LINEAR",
                                        "F_TAIL_QUADRATIC",
                                        "KB_F_TAIL_OFFSET",
                                        "KB_F_TAIL_LINEAR",
                                        "KB_F_TAIL_QUADRATIC",
                                        "GAMMA_OFFSET",
                                        "GAMMA_LINEAR",
                                        "GAMMA_QUADRATIC",
                                        "SNIP_WIDTH",
                                        "SI_ESCAPE_FACTOR",
                                        "GE_ESCAPE_FACTOR",
                                        "ESCAPE_LINEAR"
};

const std::vector<std::string> Updatable_detector_dependand_TAGS = {
                                                                "ELT1",
                                                                "ERT1",
                                                                "ICR1",
                                                                "OCR1"
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

APS_Fit_Params_Import::APS_Fit_Params_Import()
{


}

//-----------------------------------------------------------------------------

APS_Fit_Params_Import::~APS_Fit_Params_Import()
{


}

//-----------------------------------------------------------------------------

bool APS_Fit_Params_Import::load(std::string path,
                                 Element_Info_Map * element_info_map,
                                 Params_Override *params_override)
{

    std::ifstream paramFileStream(path);


    if (paramFileStream.is_open() )
    {
        //paramFileStream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        paramFileStream.exceptions(std::ifstream::failbit);
        //std::string line;
        std::string tag;
        try
        {
            for (std::string line; std::getline(paramFileStream, line, '\n'); )
            //while(std::getline(paramFileStream, line))
            {
                std::istringstream strstream(line);
                std::getline(strstream, tag, ':');
                //logit<<"tag : "<<tag<<"\n";
                if (tag == "VERSION" || tag == "DATE")
                {
                    //logit << line << "\n";
                }
                else if (tag == "DETECTOR_ELEMENTS")
                {

                }
                else if (tag == "ELEMENTS_TO_FIT")
                {
                    std::string element_symb;
                    while(std::getline(strstream, element_symb, ','))
                    {
                        element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());

                        // check if element_symb contains '_'
                        std::string base_element_symb = element_symb.substr(0, element_symb.find_last_of("_"));

                        //logit<<"Element : "<<element_symb<<" : "<<base_element_symb<<"\n";

                        Element_Info* e_info = element_info_map->get_element(base_element_symb);
                        if(e_info == nullptr)
                        {
                            logit<<"Can not find element "<<base_element_symb<<"\n";
                        }
                        else
                        {
                            Fit_Element_Map* fit_map;
                            if(params_override->elements_to_fit.count(element_symb) > 0)
                            {
                                fit_map = params_override->elements_to_fit[element_symb];
                            }
                            else
                            {
                                fit_map = new Fit_Element_Map(element_symb, e_info);
                                params_override->elements_to_fit[element_symb] = fit_map;
                            }
                        }
                    }

                }
                else if (tag == "ELEMENTS_WITH_PILEUP")
                {
                    /* compile, need to update logic
                    std::string element_symb;
                    while(std::getline(strstream, element_symb, ','))
                    {
                        element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());
                        logit<<"Element with pileup : "<<element_symb<<"\n";
                        //Element_Param* element_param = new Element_Param();
                        //element_param->name = element_symb;
                        //params_override->fit_params.append_element(element_param);
                    }
                    */
                }
                else if(FILE_TAGS_TRANSLATION.count(tag)> 0)
                {
                    std::string tag_name = FILE_TAGS_TRANSLATION.at(tag);
                    if( false == params_override->fit_params.contains(tag_name))
                    {
                        params_override->fit_params.add_parameter(Fit_Param(tag_name));
                    }

                    std::string str_value;
                    std::getline(strstream, str_value, ':');
                    float fvalue = std::stof(str_value);

                    if (tag.find("_MAX") != std::string::npos)
                    {
                        params_override->fit_params[tag_name].max_val = fvalue;
                    }
                    else if (tag.find("_MIN") != std::string::npos)
                    {
                        params_override->fit_params[tag_name].min_val = fvalue;
                    }
                    else
                    {
                        params_override->fit_params[tag_name].value = fvalue;
                    }
                }
                else if ( tag == "BRANCHING_FAMILY_ADJUSTMENT_L" || tag == "BRANCHING_RATIO_ADJUSTMENT_L" || tag == "BRANCHING_RATIO_ADJUSTMENT_K")
                {
                    unsigned int cnt = 0;

                    if (tag == "BRANCHING_FAMILY_ADJUSTMENT_L")
                    {
                        cnt = 3;
                    }
                    else if (tag == "BRANCHING_RATIO_ADJUSTMENT_K")
                    {
                        cnt = 4;
                    }
                    else if (tag == "BRANCHING_RATIO_ADJUSTMENT_L")
                    {
                        cnt = 12;
                    }

                    std::string element_symb;
                    std::string str_value;

                    std::getline(strstream, element_symb, ',');
                    element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());

                    Fit_Element_Map* fit_map;
                    if(params_override->elements_to_fit.count(element_symb) > 0)
                    {
                        fit_map = params_override->elements_to_fit[element_symb];
                        for (unsigned int i = 0; i<cnt; i++)
                        {
                            float factor = 1.0;
                            std::getline(strstream, str_value, ',');
                            factor = std::stof(str_value);
                            fit_map->set_custom_multiply_ratio(i, factor);
                        }
                    }
                }
                else if (tag == STR_FIT_SNIP_WIDTH)
                {
                    std::string str_value;
                    std::getline(strstream, str_value, ':');
                    str_value.erase(std::remove(str_value.begin(), str_value.end(), '\n'), str_value.end());
                    str_value.erase(std::remove(str_value.begin(), str_value.end(), '\r'), str_value.end());
                    str_value.erase(std::remove(str_value.begin(), str_value.end(), ' '), str_value.end());
                    float fvalue = std::stof(str_value);
                    params_override->fit_snip_width = fvalue;

                    if( false == params_override->fit_params.contains(STR_SNIP_WIDTH))
                    {
                        params_override->fit_params.add_parameter(Fit_Param(STR_SNIP_WIDTH));
                    }

                    if(fvalue > 0.0)
                        params_override->fit_params[STR_SNIP_WIDTH].bound_type = E_Bound_Type::FIT;
                    else
                        params_override->fit_params[STR_SNIP_WIDTH].bound_type = E_Bound_Type::FIXED;
                }
                else if (tag == "TIME_SCALER_PV")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->time_scaler = value;
                }
                else if (tag == "TIME_SCALER_CLOCK")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->time_scaler_clock = value;
                }
                else if (tag == "TIME_NORMALIZED_SCALER")
                {
                    std::string name, value;
                    std::getline(strstream, name, ';');
                    std::getline(strstream, value);
                    name.erase(std::remove(name.begin(), name.end(), '\n'), name.end());
                    name.erase(std::remove(name.begin(), name.end(), '\r'), name.end());
                    name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->time_normalized_scalers.insert(std::pair<std::string, std::string>(name, value));
                }
                else if (tag == "DETECTOR_MATERIAL") // =  0 = Germanium, 1 = Si
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());


                    if(value == "0")
                    {
                        params_override->detector_element = "Ge";
                    }
                    else if(value == "1")
                    {
                        params_override->detector_element = "Si";
                    }
                    else
                    {
                        logit<<"Error: Unknown detector element enumeration : "<<value<<"\n";
                    }
                }
                else if (tag == "ELT1")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->elt_pv = value;
                }
                else if (tag == "ERT1")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->ert_pv = value;
                }
                else if (tag == "ICR1")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->in_cnt_pv = value;
                }
                else if (tag == "OCR1")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->out_cnt_pv = value;
                }
                else if (tag == "US_AMP_SENS_NUM_PV")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->us_amp_sens_num_pv = value;
                }
                else if (tag == "US_AMP_SENS_UNIT_PV")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->us_amp_sens_unit_pv = value;
                }
                else if (tag == "DS_AMP_SENS_NUM_PV")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->ds_amp_sens_num_pv = value;
                }
                else if (tag == "DS_AMP_SENS_UNIT_PV")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->ds_amp_sens_unit_pv = value;
                }
                else if (tag == "US_AMP_SENS_NUM")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->us_amp_sens_num = str_to_real(value);
                }
                else if (tag == "US_AMP_SENS_UNIT")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->us_amp_sens_unit = str_to_real(value);
                }
                else if (tag == "DS_AMP_SENS_NUM")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->ds_amp_sens_num = str_to_real(value);
                }
                else if (tag == "DS_AMP_SENS_UNIT")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->ds_amp_sens_unit = str_to_real(value);
                }
                else if (tag == "BE_WINDOW_THICKNESS")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->be_window_thickness = value;
                }
                else if (tag == "DET_CHIP_THICKNESS")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->det_chip_thickness = value;
                }
                else if (tag == "GE_DEAD_LAYER")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->ge_dead_layer = value;
                }
                else if (tag == "SI_ESCAPE_ENABLE")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());


                    if(value == "1")
                    {
                        params_override->si_escape_enabled = true;
                    }
                    else
                    {
                        params_override->si_escape_enabled = false;
                    }
                }
                else if (tag == "GE_ESCAPE_ENABLE")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());


                    if(value == "1")
                    {
                        params_override->ge_escape_enabled = true;
                    }
                    else
                    {
                        params_override->ge_escape_enabled = false;
                    }
                }
				/*
                else if (tag == "MAX_ENERGY_TO_FIT")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->max_energy = std::stof(value);
                }
                else if (tag == "MIN_ENERGY_TO_FIT")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->min_energy = std::stof(value);
                }
				*/
                else if (tag == "LINEAR_ESCAPE_FACTOR")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    //params_override->ge_dead_layer = std::stoi(value);
                }
                else if (tag == "TAIL_FRACTION_ADJUST_SI")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    //params_override->ge_dead_layer = std::stoi(value);
                }
                else if (tag == "TAIL_WIDTH_ADJUST_SI")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    //params_override->ge_dead_layer = std::stoi(value);
                }
                else if (tag == "IDENTIFYING_NAME_[WHATEVERE_YOU_LIKE]")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    //params_override->ge_dead_layer = std::stoi(value);
                }
                else
                {
                    if (tag.length() > 0 && tag[0] != ' ' && (line.find(":") != std::string::npos))
                    {
                        std::string value;
                        std::getline(strstream, value);
                        value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                        value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                        value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                        params_override->scaler_pvs[tag] = value;
                    }
                }

            }
        }
        catch(std::exception e)
        {
            if (paramFileStream.eof() == 0 && (paramFileStream.bad() || paramFileStream.fail()) )
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << paramFileStream.fail()
                    << "\neofbit: " << paramFileStream.eof()
                    << "\nbadbit: " << paramFileStream.bad() << "\n";
            }
        }

        paramFileStream.close();
        return true;
    }
    return false;

}

//-----------------------------------------------------------------------------

bool APS_Fit_Params_Import::save(std::string path,
                                 Fit_Parameters fit_params,
                                 int detector_num)
{

    std::ifstream in_stream(path);
    std::string save_path = path+std::to_string(detector_num);
    std::ofstream out_stream(save_path);

    logit<<save_path<<"\n";

    std::chrono::system_clock::time_point today = std::chrono::system_clock::now();
    std::time_t tt;
    tt = std::chrono::system_clock::to_time_t(today);

    if ( in_stream.is_open() && out_stream.is_open() )
    {
        in_stream.exceptions(std::ifstream::failbit);
        //std::string line;
        std::string tag;
        try
        {
            for (std::string line; std::getline(in_stream, line, '\n'); )
            {
                std::istringstream strstream(line);
                std::getline(strstream, tag, ':');
                if(tag == "DATE")
                {
                    out_stream<<tag<<": "<<std::ctime(&tt)<<"\n";
                }
                else if( std::find(Updatable_TAGS.begin(), Updatable_TAGS.end(), tag) != Updatable_TAGS.end() )
                {
                    std::string trans_str = FILE_TAGS_TRANSLATION.at(tag);
                    if(fit_params.contains( trans_str ) )
                    {
                        std::string str_val = std::to_string(fit_params.at( trans_str ).value);
                        out_stream<<tag<<": "<<str_val<<"\n";
                    }
                    else
                    {
                        out_stream<<line<<"\n";
                    }
                }
                else if( std::find(Updatable_detector_dependand_TAGS.begin(), Updatable_detector_dependand_TAGS.end(), tag) != Updatable_detector_dependand_TAGS.end() )
                {
                    /*
                    ELT1: dxpXMAP2xfm3:mca1.ELTM
                    ERT1: dxpXMAP2xfm3:mca1.ERTM
                    ICR1: dxpXMAP2xfm3:dxp1:InputCountRate
                    OCR1: dxpXMAP2xfm3:dxp1:OutputCountRate
                    */
                    std::string str_det_num = std::to_string(detector_num + 1);
                    int idx = line.find(":mca");
                    if(idx > -1)
                    {
                        line[idx+4] = str_det_num[0];
                    }

                    idx = line.find(":dxp");
                    if(idx > -1)
                    {

                        line[idx+4] = str_det_num[0];
                    }
                    out_stream<<line<<"\n";
                }
                else
                {
                    out_stream<<line<<"\n";
                }
            }
        }

        catch(std::exception e)
        {
            if (in_stream.eof() == 0 && (in_stream.bad() || in_stream.fail()) )
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << in_stream.fail()
                    << "\neofbit: " << in_stream.eof()
                    << "\nbadbit: " << in_stream.bad() << "\n";
            }
        }

        in_stream.close();
        out_stream.close();
        return true;

    }

    logit<<"Error: APS_Fit_Params_Import::save() : couldn't opening file "<<path<<"\n";
    return false;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} //end namespace aps
} //end namespace file
}// end namespace io
