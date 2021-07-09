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
//-----------------------------------------------------------------------------

bool load_parameters_override(std::string path, Params_Override *params_override)
{

    data_struct::Element_Info_Map *element_info_map = data_struct::Element_Info_Map::inst();
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
                //tag.erase(std::remove_if(tag.begin(), tag.end(), ::isspace), tag.end());
                //logD<<"tag : "<<tag<<"\n";
                if (tag == "VERSION" || tag == "DATE")
                {
                    //logD << line << "\n";
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

                        //logD<<"Element : "<<element_symb<<" : "<<base_element_symb<<"\n";

                        Element_Info* e_info = element_info_map->get_element(base_element_symb);
                        if(e_info == nullptr)
                        {
                            logW<<"Can not find element "<<base_element_symb<<"\n";
                        }
                        else
                        {
                            Fit_Element_Map* fit_map;
                            if(params_override->elements_to_fit.count(element_symb) < 1)
                            {
                                fit_map = new Fit_Element_Map(element_symb, e_info);
                                params_override->elements_to_fit[element_symb] = fit_map;
                            }
                        }
                    }

                }
                else if (tag == "ELEMENTS_WITH_PILEUP")
                {
                    std::string element_symb;
                    while(std::getline(strstream, element_symb, ','))
                    {
                        Element_Info* e_info1 = nullptr;
                        Element_Info* e_info2 = nullptr;
                        std::string efull_name1;
                        std::string efull_name2;

                        element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());
                        //logI<<"Element with pileup : "<<element_symb<<"\n";
                        std::string orig_el_symb = element_symb;

                        std::vector<std::string> string_list;
                        std::size_t found = element_symb.find("_");
                        size_t prev = 0;
                        while(found!=std::string::npos)
                        {
                            string_list.push_back(element_symb.substr(prev, found));
                            element_symb = element_symb.substr(found+1, element_symb.size());
                            prev = found;
                            found = element_symb.find("_");
                        }
                        if(element_symb.size() > 0)
                        {
                            string_list.push_back(element_symb);
                        }


                        if(string_list.size() == 4)
                        {
                            if((string_list[1] == "L" || string_list[1] == "M") && (string_list[3] == "L" || string_list[3] == "M"))
                            {
                                e_info1 = element_info_map->get_element(string_list[0]);
                                e_info2 = element_info_map->get_element(string_list[2]);
                                efull_name1 = string_list[0] + "_" + string_list[1];
                                efull_name2 = string_list[2] + "_" + string_list[3];
                            }
                        }
                        else if(string_list.size() == 3)
                        {
                            if(string_list[1] == "L" || string_list[1] == "M")
                            {
                                e_info1 = element_info_map->get_element(string_list[0]);
                                e_info2 = element_info_map->get_element(string_list[2]);
                                efull_name1 = string_list[0] + "_" + string_list[1];
                                efull_name2 = string_list[2];
                            }
                            else if(string_list[2] == "L" || string_list[2] == "M")
                            {
                                e_info1 = element_info_map->get_element(string_list[0]);
                                e_info2 = element_info_map->get_element(string_list[1]);
                                efull_name1 = string_list[0];
                                efull_name2 = string_list[1] + "_" + string_list[2];
                            }
                        }
                        else if(string_list.size() == 2)
                        {
                            e_info1 = element_info_map->get_element(string_list[0]);
                            e_info2 = element_info_map->get_element(string_list[1]);
                            efull_name1 = string_list[0];
                            efull_name2 = string_list[1];
                        }


                        if(e_info1 != nullptr && e_info2 != nullptr)
                        {
                            Fit_Element_Map* fit_map;
                            if(params_override->elements_to_fit.count(orig_el_symb) < 1)
                            {
                                fit_map = new Fit_Element_Map(efull_name1, e_info1);
                                fit_map->set_as_pileup(efull_name2, e_info2);
                                params_override->elements_to_fit[orig_el_symb] = fit_map;
                            }
                        }
                        else
                        {
                            logW<<"Could not parse pileup string: "<<orig_el_symb<<".\n";
                        }
                    }
                }
                else if(FILE_TAGS_TRANSLATION.count(tag)> 0)
                {
                    //ignore quadratic min because we don't want it to be negative so we default min to 0
                    if (tag == "CAL_QUAD_[E_QUADRATIC]_MIN")
                    {
                        continue;
                    }

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
                        params_override->branching_family_L.push_back(line);
                        cnt = 3;
                    }
                    else if (tag == "BRANCHING_RATIO_ADJUSTMENT_K")
                    {
                        params_override->branching_ratio_K.push_back(line);
                        cnt = 4;
                    }
                    else if (tag == "BRANCHING_RATIO_ADJUSTMENT_L")
                    {
                        params_override->branching_ratio_L.push_back(line);
                        cnt = 12;
                    }

                    std::string element_symb;
                    std::string str_value;

                    std::getline(strstream, element_symb, ',');
                    element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());

                    Fit_Element_Map* fit_map;
					if (params_override->elements_to_fit.count(element_symb) > 0)
					{
						fit_map = params_override->elements_to_fit[element_symb];
						if (cnt == 3) // family
						{
							float factor = 1.0;

							// 1
							std::getline(strstream, str_value, ',');
							factor = std::stof(str_value);
							fit_map->multiply_custom_multiply_ratio(4, factor);
							fit_map->multiply_custom_multiply_ratio(5, factor);
							fit_map->multiply_custom_multiply_ratio(7, factor);
							fit_map->multiply_custom_multiply_ratio(8, factor);
							fit_map->multiply_custom_multiply_ratio(9, factor);

							// 2
							std::getline(strstream, str_value, ',');
							factor = std::stof(str_value);
							fit_map->multiply_custom_multiply_ratio(2, factor);
							fit_map->multiply_custom_multiply_ratio(6, factor);
							fit_map->multiply_custom_multiply_ratio(11, factor);

							//3
							std::getline(strstream, str_value, ',');
							factor = std::stof(str_value);
							fit_map->multiply_custom_multiply_ratio(0, factor);
							fit_map->multiply_custom_multiply_ratio(1, factor);
							fit_map->multiply_custom_multiply_ratio(3, factor);
							fit_map->multiply_custom_multiply_ratio(10, factor);
						}
						else // ratio's
						{
							for (unsigned int i = 0; i < cnt; i++)
							{
								float factor = 1.0;
								std::getline(strstream, str_value, ',');
								factor = std::stof(str_value);
								fit_map->multiply_custom_multiply_ratio(i, factor);
							}
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
                        logE<<"Unknown detector element enumeration : "<<value<<"\n";
                    }
                }
                else if (tag == "ELT" || tag == "ELT1" || tag == "ELT2" || tag == "ELT3" || tag == "ELT4")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->elt_pv = value;
                }
                else if (tag == "ERT" || tag == "ERT1" || tag == "ERT2" || tag == "ERT3" || tag == "ERT4")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->ert_pv = value;
                }
                else if (tag == "ICR" || tag == "ICR1" || tag == "ICR2" || tag == "ICR3" || tag == "ICR4")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->in_cnt_pv = value;
                }
                else if (tag == "OCR" || tag == "OCR1" || tag == "OCR2" || tag == "OCR3" || tag == "OCR4")
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
                else if (tag == "AIRPATH")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->airpath = value;
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
                else if (tag == "THETA_PV")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    params_override->theta_pv = value;
                }
                else if (tag == "SUMMED_SCALER")
                {
                    data_struct::Summed_Scaler s_scaler;
                    s_scaler.normalize_by_time = false;
                    std::string value;
                    std::string scaler_name;
                    std::getline(strstream, value, ':');
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    // split scalers names by ','
                    s_scaler.scaler_name = value;
                    std::string last_scaler;
                    std::getline(strstream, scaler_name, ',');
                    while(last_scaler != scaler_name)
                    {
                        scaler_name.erase(std::remove(scaler_name.begin(), scaler_name.end(), '\n'), scaler_name.end());
                        scaler_name.erase(std::remove(scaler_name.begin(), scaler_name.end(), '\r'), scaler_name.end());
                        scaler_name.erase(std::remove(scaler_name.begin(), scaler_name.end(), ' '), scaler_name.end());
                        last_scaler = scaler_name;
                        // add scaler name and set mda_idx to -1, we will search for the index later and unpdate
						if (scaler_name.length() > 0)
						{
							s_scaler.scalers_to_sum.push_back(scaler_name);
						}
                        std::getline(strstream, scaler_name, ',');
                    }
                    params_override->summed_scalers.push_back(s_scaler);
                }
                /* // we don't need this anymore since the time normalization will happen before summing.
                else if (tag == "TIME_NORMALIZED_SUMMED_SCALER")
                {
                    data_struct::Summed_Scaler s_scaler;
                    s_scaler.normalize_by_time = true;
                    std::string value;
                    std::string scaler_name;
                    std::getline(strstream, value, ':');
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    // split scalers names by ','
                    s_scaler.scaler_name = value;
                    std::string last_scaler;
                    std::getline(strstream, scaler_name, ',');
                    while(last_scaler != scaler_name)
                    {
                        scaler_name.erase(std::remove(scaler_name.begin(), scaler_name.end(), '\n'), scaler_name.end());
                        scaler_name.erase(std::remove(scaler_name.begin(), scaler_name.end(), '\r'), scaler_name.end());
                        scaler_name.erase(std::remove(scaler_name.begin(), scaler_name.end(), ' '), scaler_name.end());
                        last_scaler = scaler_name;
                        // add scaler name and set mda_idx to -1, we will search for the index later and unpdate
                        s_scaler.scalers_to_sum.push_back(scaler_name);
                        std::getline(strstream, scaler_name, ',');
                    }
                    params_override->summed_scalers.push_back(s_scaler);
                }
                */
                else
                {
                    if (tag.length() > 0 && tag[0] != ' ' && tag[0] != '\t' && (line.find(":") != std::string::npos))
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
        catch(std::exception& e)
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

		params_override->parse_and_gen_branching_ratios();

        paramFileStream.close();
        return true;
    }
    return false;

}

//-----------------------------------------------------------------------------

bool save_parameters_override(std::string path, Params_Override *params_override)
{

    /*
    else if (tag == "AIRPATH")

    else if (tag == "TAIL_FRACTION_ADJUST_SI")

    else if (tag == "TAIL_WIDTH_ADJUST_SI")
    */
	std::vector<std::string> branching_ratios_updated;

    if (params_override == nullptr)
    {
        logE << "No parameters to save in " << path << "\n";
        return false;
    }

    std::ofstream out_stream(path);

    logI << path << "\n";

    if (out_stream.is_open())
    {
        out_stream << "   This file will override default fit settings for the maps program for a 3 element detector remove : removeme_ * elementdetector_to make it work. \n";
        out_stream << "   NOTE : the filename MUST be maps_fit_parameters_override.txt\n";
        out_stream << "VERSION: 5.00000\n";
        out_stream << "DATE: \n";
        out_stream << "   put below the number of detectors that were used to acquire spectra. IMPORTANT: \n";
        out_stream << "   this MUST come after VERSION, and before all other options!\n";
        out_stream << "DETECTOR_ELEMENTS:       1\n";
        out_stream << "   give this file an internal name, whatever you like\n";
        out_stream << "IDENTIFYING_NAME_[WHATEVERE_YOU_LIKE]:automatic\n";
        out_stream << "   list the elements that you want to be fit. For K lines, just use the element\n";
        out_stream << "   name, for L lines add _L, e.g., Au_L, for M lines add _M\n";
        out_stream << "ELEMENTS_TO_FIT: ";
        for (const auto& itr : params_override->elements_to_fit)
        {
            // if not pileup
            if (itr.second->pileup_element() == nullptr && (itr.first != "COHERENT_SCT_AMPLITUDE" || itr.first != "COMPTON_AMPLITUDE"))
            {
                out_stream << itr.first << " , ";
            }
        }
        out_stream << "\n";
        out_stream << "   list the element combinations you want to fit for pileup, e.g., Si_Si, Si_Si_Si, Si_Cl, etc\n";
        out_stream << "ELEMENTS_WITH_PILEUP: ";
        for (const auto& itr : params_override->elements_to_fit)
        {
            // if not pileup
            if (itr.second->pileup_element() != nullptr)
            {
                out_stream << itr.first << " , ";
            }
        }
        out_stream << "\n";
        out_stream << "   offset of energy calibration, in kev\n";
        out_stream << "CAL_OFFSET_[E_OFFSET]:   "<< params_override->fit_params.at(STR_ENERGY_OFFSET).value <<"\n";
        out_stream << "CAL_OFFSET_[E_OFFSET]_MIN:   "<< params_override->fit_params.at(STR_ENERGY_OFFSET).min_val <<"\n";
        out_stream << "CAL_OFFSET_[E_OFFSET]_MAX:   "<< params_override->fit_params.at(STR_ENERGY_OFFSET).max_val <<"\n";
        out_stream << "   slope of energy calibration, in leV / channel\n";
        out_stream << "CAL_SLOPE_[E_LINEAR]:   "<< params_override->fit_params.at(STR_ENERGY_SLOPE).value <<"\n";
        out_stream << "CAL_SLOPE_[E_LINEAR]_MIN:   "<< params_override->fit_params.at(STR_ENERGY_SLOPE).min_val <<"\n";
        out_stream << "CAL_SLOPE_[E_LINEAR]_MAX:   "<< params_override->fit_params.at(STR_ENERGY_SLOPE).max_val <<"\n";
        out_stream << "   quadratic correction for energy calibration, unless you know exactly what you are doing, please leave it at 0.\n";
        out_stream << "CAL_QUAD_[E_QUADRATIC]:   "<< params_override->fit_params.at(STR_ENERGY_QUADRATIC).value <<"\n";
        out_stream << "CAL_QUAD_[E_QUADRATIC]_MIN:   "<< params_override->fit_params.at(STR_ENERGY_QUADRATIC).min_val <<"\n";
        out_stream << "CAL_QUAD_[E_QUADRATIC]_MAX:   "<< params_override->fit_params.at(STR_ENERGY_QUADRATIC).max_val <<"\n";
        out_stream << "    energy_resolution at 0keV\n";
        out_stream << "FWHM_OFFSET: "<< params_override->fit_params.at(STR_FWHM_OFFSET).value << "\n";
        out_stream << "    energy dependence of the energy resolution\n";
        out_stream << "FWHM_FANOPRIME: " << params_override->fit_params.at(STR_FWHM_FANOPRIME).value << "\n";
        out_stream << "    incident energy\n";
        out_stream << "COHERENT_SCT_ENERGY: " << params_override->fit_params.at(STR_COHERENT_SCT_ENERGY).value << "\n";
        out_stream << "    upper contstraint for the incident energy\n";
        out_stream << "COHERENT_SCT_ENERGY_MAX: " << params_override->fit_params.at(STR_COHERENT_SCT_ENERGY).max_val << "\n";
        out_stream << "    lower contstraint for the incident energy\n";
        out_stream << "COHERENT_SCT_ENERGY_MIN: " << params_override->fit_params.at(STR_COHERENT_SCT_ENERGY).min_val << "\n";
        out_stream << "    angle for the compton scatter (in degrees)\n";
        out_stream << "COMPTON_ANGLE: " << params_override->fit_params.at(STR_COMPTON_ANGLE).value << "\n";
        out_stream << "COMPTON_ANGLE_MAX: " << params_override->fit_params.at(STR_COMPTON_ANGLE).max_val << "\n";
        out_stream << "COMPTON_ANGLE_MIN: " << params_override->fit_params.at(STR_COMPTON_ANGLE).min_val << "\n";
        out_stream << "    additional width of the compton\n";
        out_stream << "COMPTON_FWHM_CORR: " << params_override->fit_params.at(STR_COMPTON_FWHM_CORR).value << "\n";
        out_stream << "COMPTON_STEP: " << params_override->fit_params.at(STR_COMPTON_F_STEP).value << "\n";
        out_stream << "COMPTON_F_TAIL: " << params_override->fit_params.at(STR_COMPTON_F_TAIL).value << "\n";
        out_stream << "COMPTON_GAMMA: " << params_override->fit_params.at(STR_COMPTON_GAMMA).value << "\n";
        out_stream << "COMPTON_HI_F_TAIL: " << params_override->fit_params.at(STR_COMPTON_HI_F_TAIL).value << "\n";
        out_stream << "COMPTON_HI_GAMMA: " << params_override->fit_params.at(STR_COMPTON_HI_GAMMA).value << "\n";
        out_stream << "    tailing parameters, see also Grieken, Markowicz, Handbook of X-ray spectrometry\n";
        out_stream << "    2nd ed, van Espen spectrum evaluation page 287.  _A corresponds to f_S, _B to\n";
        out_stream << "    f_T and _C to gamma\n";
        out_stream << "STEP_OFFSET: " << params_override->fit_params.at(STR_F_STEP_OFFSET).value << "\n";
        out_stream << "STEP_LINEAR: " << params_override->fit_params.at(STR_F_STEP_LINEAR).value << "\n";
        out_stream << "STEP_QUADRATIC: " << params_override->fit_params.at(STR_F_STEP_QUADRATIC).value << "\n";
        out_stream << "F_TAIL_OFFSET: " << params_override->fit_params.at(STR_F_TAIL_OFFSET).value << "\n";
        out_stream << "F_TAIL_LINEAR: " << params_override->fit_params.at(STR_F_TAIL_LINEAR).value << "\n";
        out_stream << "F_TAIL_QUADRATIC: " << params_override->fit_params.at(STR_F_TAIL_QUADRATIC).value << "\n";
        out_stream << "KB_F_TAIL_OFFSET: " << params_override->fit_params.at(STR_KB_F_TAIL_OFFSET).value << "\n";
        out_stream << "KB_F_TAIL_LINEAR: " << params_override->fit_params.at(STR_KB_F_TAIL_LINEAR).value << "\n";
        out_stream << "KB_F_TAIL_QUADRATIC: " << params_override->fit_params.at(STR_KB_F_TAIL_QUADRATIC).value << "\n";
        out_stream << "GAMMA_OFFSET: " << params_override->fit_params.at(STR_GAMMA_OFFSET).value << "\n";
        out_stream << "GAMMA_LINEAR: " << params_override->fit_params.at(STR_GAMMA_LINEAR).value << "\n";
        out_stream << "GAMMA_QUADRATIC: " << params_override->fit_params.at(STR_GAMMA_QUADRATIC).value << "\n";
        out_stream << "    snip width is the width used for estimating background. 0.5 is typically a good start\n";
        out_stream << "SNIP_WIDTH: " << params_override->fit_params.at(STR_SNIP_WIDTH).value << "\n";
        out_stream << "    set FIT_SNIP_WIDTH to 1 to fit the width of the snipping for background estimate, set to 0 not to. Only use if you know what it is doing!\n";
        out_stream << "FIT_SNIP_WIDTH: " << params_override->fit_snip_width << "\n";
        out_stream << "    detector material: 0= Germanium, 1 = Si\n";
        out_stream << "DETECTOR_MATERIAL: ";
        if(params_override->detector_element == "Si")
            out_stream << "1\n";
        else
            out_stream << "0\n";
        out_stream << "    beryllium window thickness, in micrometers, typically 8 or 24\n";
        out_stream << "BE_WINDOW_THICKNESS: " << params_override->be_window_thickness << "\n";
        out_stream << "    thickness of the detector chip, e.g., 350 microns for an SDD\n";
        out_stream << "DET_CHIP_THICKNESS: " << params_override->det_chip_thickness << "\n";
        out_stream << "    thickness of the Germanium detector dead layer, in microns, for the purposes of the NBS calibration\n";
        out_stream << "GE_DEAD_LAYER: " << params_override->ge_dead_layer << "\n";
        out_stream << "    maximum energy value to fit up to [keV]\n";
        out_stream << "MAX_ENERGY_TO_FIT: " << params_override->fit_params.at(STR_MAX_ENERGY_TO_FIT).value << "\n";
        out_stream << "    minimum energy value [keV]\n";
        out_stream << "MIN_ENERGY_TO_FIT: " << params_override->fit_params.at(STR_MIN_ENERGY_TO_FIT).value << "\n";


		for (const auto& itr : params_override->elements_to_fit)
		{
			for (float ratio : itr.second->energy_ratio_multipliers())
			{
				if (ratio != 1.0)
				{
					branching_ratios_updated.push_back(itr.first+",");
					if (itr.second->shell_type_as_string() == "K")
					{
						out_stream << "BRANCHING_RATIO_ADJUSTMENT_K: " << itr.first;
					}
					else if (itr.second->shell_type_as_string() == "L")
					{
						out_stream << "BRANCHING_RATIO_ADJUSTMENT_L: " << itr.first;
					}
					for (float r : itr.second->energy_ratio_multipliers())
					{
						out_stream << "," << r;
					}
					out_stream << "\n";
					break;
				}
			}
		}
		// we still need these to be saved if not included above
		out_stream << "    this allows manual adjustment of the branhcing ratios between the different lines of L1, L2, and L3.\n";
		out_stream << "    note, the numbers that are put in should be RELATIVE modifications, i.e., a 1 will correspond to exactly the literature value,\n";
		out_stream << "    0.8 will correspond to to 80% of that, etc.\n";
        for (const auto& itr : params_override->branching_family_L)
        {
			bool found = false;
			for (std::string element_name : branching_ratios_updated)
			{
				if (itr.find(element_name) != std::string::npos)
				{
					found = true;
					break;
				}
			}
			if (false == found)
			{
				out_stream << itr << "\n";
			}
        }
		out_stream << "    this allows manual adjustment of the branhcing ratios between the different L lines, such as La 1, la2, etc.\n";
		out_stream << "    Please note, these are all RELATIVE RELATIVE modifications, i.e., a 1 will correspond to exactly the literature value, etc.\n";
		out_stream << "    all will be normalized to the La1 line, and the values need to be in the following order:\n";
		out_stream << "    La1, La2, Lb1, Lb2, Lb3, Lb4, Lg1, Lg2, Lg3, Lg4, Ll, Ln\n";
		out_stream << "    please note, the first value (la1) MUST BE A 1. !!!\n";
        for (const auto& itr : params_override->branching_ratio_L)
        {
			bool found = false;
			for (std::string element_name : branching_ratios_updated)
			{
				if (itr.find(element_name) != std::string::npos)
				{
					found = true;
					break;
				}
			}
			if (false == found)
			{
				out_stream << itr << "\n";
			}
        }
		out_stream << "    this allows manual adjustment of the branhcing ratios between the different K lines, such as Ka1, Ka2, Kb1, Kb2\n";
		out_stream << "    Please note, these are all RELATIVE RELATIVE modifications, i.e., a 1 will correspond to exactly the literature value, etc.\n";
		out_stream << "    all will be normalized to the Ka1 line, and the values need to be in the following order:\n";
		out_stream << "    Ka1, Ka2, Kb1(+3), Kb2\n";
		out_stream << "    please note, the first value (Ka1) MUST BE A 1. !!!\n";
        for (const auto& itr : params_override->branching_ratio_K)
        {
			bool found = false;
			for (std::string element_name : branching_ratios_updated)
			{
				if (itr.find(element_name) != std::string::npos)
				{
					found = true;
					break;
				}
			}
			if (false == found)
			{
				out_stream << itr << "\n";
			}
        }
		

        out_stream << "    the parameter adds the escape peaks (offset) to the fit if larger than 0. You should not enable Si and Ge at the same time, ie, one of these two values should be zero\n";
        out_stream << "SI_ESCAPE_FACTOR: " << params_override->si_escape_factor << "\n";
        out_stream << "GE_ESCAPE_FACTOR: " << params_override->ge_escape_factor << "\n";
        out_stream << "    this parameter adds a component to the escape peak that depends linear on energy\n";
        out_stream << "LINEAR_ESCAPE_FACTOR: 0.0\n";
        out_stream << "    the parameter enables fitting of the escape peak strengths. set 1 to enable, set to 0 to disable. (in matrix fitting always disabled)\n";
        out_stream << "SI_ESCAPE_ENABLE: " << params_override->si_escape_enabled << "\n";
        out_stream << "GE_ESCAPE_ENABLE: " << params_override->ge_escape_enabled << "\n";
        out_stream << "    the lines below(if any) give backup description of IC amplifier sensitivity, in case it cannot be found in the mda file\n";
        out_stream << "    for the amps, the _NUM value should be between 0 and 8 where 0 = 1, 1 = 2, 2 = 5, 3 = 10, 4 = 20, 5 = 50, 6 = 100, 7 = 200, 8 = 500\n";
        out_stream << "    for the amps, the _UNIT value should be between 0 and 3 where 0 = pa / v, 1 = na / v, 2 = ua / v 3 = ma / v\n";
        out_stream << "US_AMP_SENS_NUM: " << params_override->us_amp_sens_num << "\n";
        out_stream << "US_AMP_SENS_UNIT: " << params_override->us_amp_sens_unit << "\n";
        out_stream << "DS_AMP_SENS_NUM: " << params_override->ds_amp_sens_num << "\n";
        out_stream << "DS_AMP_SENS_UNIT: "<< params_override->ds_amp_sens_unit <<"\n";
        out_stream << "THETA_PV: " << params_override->theta_pv << "\n";
        out_stream << "    the lines (if any) below will override the detector names built in to maps. please modify only if you are sure you understand the effect\n";
        for (const auto& itr : params_override->scaler_pvs)
        {
            out_stream << itr.first << ":" << itr.second << "\n";
        }
        out_stream << "TIME_SCALER_PV: " << params_override->time_scaler << "\n";
        out_stream << "TIME_SCALER_CLOCK: " << params_override->time_scaler_clock << "\n";
        for (const auto& itr : params_override->time_normalized_scalers)
        {
            out_stream <<"TIME_NORMALIZED_SCALER: "<< itr.first << ";" << itr.second << "\n";
        }

        for (const auto& itr : params_override->summed_scalers)
        {
            if (itr.normalize_by_time)
            {
                out_stream << "TIME_NORMALIZED_SUMMED_SCALER: " << itr.scaler_name << ":";
                for (const auto& inner_itr : itr.scalers_to_sum)
                {
                    out_stream << inner_itr << ",";
                }
                out_stream << "\n";
            }
            else
            {
                out_stream << "SUMMED_SCALER: "<< itr.scaler_name<<":";
                for (const auto& inner_itr : itr.scalers_to_sum)
                {
                    out_stream << inner_itr << ",";
                }
                out_stream << "\n";
            }
        }

        out_stream.close();
        return true;
    }
    logE << "Failed to open file " << path << "\n";
    return false;
}

//-----------------------------------------------------------------------------

bool save_parameters_override(std::string path, Fit_Parameters fit_params, int detector_num)
{

    std::ifstream in_stream(path);
    std::string save_path = path+std::to_string(detector_num);
    std::ofstream out_stream(save_path);

    logI<<save_path<<"\n";

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
                    int idx = line.rfind(":mca1");
                    if(idx > -1)
                    {
                        line[idx+4] = str_det_num[0];
                    }
                    idx = line.rfind(":mca2");
                    if (idx > -1)
                    {
                        line[idx + 4] = str_det_num[0];
                    }
                    idx = line.rfind(":mca3");
                    if (idx > -1)
                    {
                        line[idx + 4] = str_det_num[0];
                    }
                    idx = line.rfind(":mca4");
                    if (idx > -1)
                    {
                        line[idx + 4] = str_det_num[0];
                    }


                    idx = line.rfind(":dxp1");
                    if(idx > -1)
                    {

                        line[idx+4] = str_det_num[0];
                    }
                    idx = line.rfind(":dxp2");
                    if (idx > -1)
                    {

                        line[idx + 4] = str_det_num[0];
                    }
                    idx = line.rfind(":dxp3");
                    if (idx > -1)
                    {

                        line[idx + 4] = str_det_num[0];
                    }
                    idx = line.rfind(":dxp4");
                    if (idx > -1)
                    {

                        line[idx + 4] = str_det_num[0];
                    }

                    out_stream<<line<<"\n";
                }
                else
                {
                    out_stream<<line<<"\n";
                }
            }
        }

        catch(std::exception& e)
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

    logE<<"Couldn't opening file "<<path<<"\n";
    return false;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} //end namespace aps
} //end namespace file
}// end namespace io
