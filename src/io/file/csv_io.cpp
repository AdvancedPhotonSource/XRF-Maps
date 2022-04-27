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



#include "csv_io.h"

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <iostream>

namespace io
{
namespace file
{
namespace csv
{

// ----------------------------------------------------------------------------

template<typename T_real>
bool load_raw_spectra(std::string filename, unordered_map<string, ArrayTr<T_real>> &data)
{
    std::ifstream file_stream(filename);
    try
    {
        std::string tmp_line;
        int num_lines;
        // First get number of lines in file so we know spectra size
        for (num_lines = 0; std::getline(file_stream, tmp_line); num_lines++ )
        {
            
        }
        // rewind file and start reading data
        file_stream.clear();
        file_stream.seekg(0);

        // subtract 1 from num lines for header
        num_lines--;
        bool first_line = true;
        int line_num = 0;
        std::vector<string> names;
        for (std::string line; std::getline(file_stream, line); )
        {
            std::stringstream strstream(line);
            if (first_line)
            {
                for (std::string value; std::getline(strstream, value, ','); )
                {
                    names.push_back(value);
                    data[value] = ArrayTr<T_real>();
                    data[value].resize(num_lines);
                    data[value].setZero(num_lines);
                }
                first_line = false;
            }
            else
            {
                int idx = 0;
                for (std::string value; std::getline(strstream, value, ','); )
                {
                    data[names[idx]](line_num) = stof(value);
                    idx++;
                }
                line_num++;
            }
        }
    }
    catch (std::exception& e)
    {
        if (file_stream.eof() == 0 && (file_stream.bad() || file_stream.fail()))
        {
            std::cerr << "ios Exception happened: " << e.what() << "\n"
                << "Error bits are: "
                << "\nfailbit: " << file_stream.fail()
                << "\neofbit: " << file_stream.eof()
                << "\nbadbit: " << file_stream.bad() << "\n";
        }
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

template<typename T_real>
bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<T_real>* energy, const data_struct::ArrayTr<T_real>* spectra, const data_struct::ArrayTr<T_real>* spectra_model, const data_struct::ArrayTr<T_real>* background)
{
    return save_fit_and_int_spectra<T_real>(fullpath, energy, spectra, spectra_model, background, nullptr);
}

// ----------------------------------------------------------------------------

template<typename T_real>
bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<T_real>* energy, const data_struct::ArrayTr<T_real>* spectra, const data_struct::ArrayTr<T_real>* spectra_model, const data_struct::ArrayTr<T_real>* background, unordered_map<string, data_struct::ArrayTr<T_real>>* labeled_spectras)
{
    if (energy == nullptr || spectra == nullptr || spectra_model == nullptr || background == nullptr)
    {
        return false;
    }

    data_struct::ArrayTr<T_real> temp_zero(energy->size());
    temp_zero.setZero(energy->size());
    // set detailed lines to zero
    data_struct::ArrayTr<T_real>* k_alpha = &temp_zero;
    data_struct::ArrayTr<T_real>* k_beta = &temp_zero;
    data_struct::ArrayTr<T_real>* l_line = &temp_zero;
    data_struct::ArrayTr<T_real>* m_line = &temp_zero;
    data_struct::ArrayTr<T_real>* step = &temp_zero;
    data_struct::ArrayTr<T_real>* tail = &temp_zero;
    data_struct::ArrayTr<T_real>* elastic = &temp_zero;
    data_struct::ArrayTr<T_real>* compton = &temp_zero;
    data_struct::ArrayTr<T_real>* pileup = &temp_zero;
    data_struct::ArrayTr<T_real>* escape = &temp_zero;


    if (labeled_spectras != nullptr)
    {      
        if (labeled_spectras->count(STR_K_A_LINES) > 0)
        {
            k_alpha = &(labeled_spectras->at(STR_K_A_LINES));
        }
        if (labeled_spectras->count(STR_K_B_LINES) > 0)
        {
            k_beta = &(labeled_spectras->at(STR_K_B_LINES));
        }
        if (labeled_spectras->count(STR_L_LINES) > 0)
        {
            l_line = &(labeled_spectras->at(STR_L_LINES));
        }
        if (labeled_spectras->count(STR_M_LINES) > 0 )
        {
            m_line = &(labeled_spectras->at(STR_M_LINES));
        }
        if (labeled_spectras->count(STR_STEP_LINES) > 0)
        {
            step = &(labeled_spectras->at(STR_STEP_LINES));
        }
        if (labeled_spectras->count(STR_TAIL_LINES) > 0)
        {
            tail = &(labeled_spectras->at(STR_TAIL_LINES));
        }
        if (labeled_spectras->count(STR_ELASTIC_LINES) > 0)
        {
            elastic = &(labeled_spectras->at(STR_ELASTIC_LINES));
        }
        if (labeled_spectras->count(STR_COMPTON_LINES) > 0)
        {
            compton = &(labeled_spectras->at(STR_COMPTON_LINES));
        }
        if (labeled_spectras->count(STR_PILEUP_LINES) > 0)
        {
            pileup = &(labeled_spectras->at(STR_PILEUP_LINES));
        }
        if (labeled_spectras->count(STR_ESCAPE_LINES) > 0)
        {
            escape = &(labeled_spectras->at(STR_ESCAPE_LINES));
        }
    }

    std::ofstream file_stream(fullpath);
    if (file_stream.is_open())
    {
        file_stream << "Energy,Spectrum,Fitted,Background,K alpha, K beta, L Lines, M Lines, step, tail, elastic, compton, pileip, escape" << "\n";

        for (int i = 0; i < energy->size(); i++)
        {
            file_stream << (*energy)(i) << "," << (*spectra)(i) << "," << (*spectra_model)(i) << "," << (*background)(i) << "," << (*k_alpha)(i) << "," << (*k_beta)(i) << "," << (*l_line)(i) << "," << (*m_line)(i) << "," << (*step)(i) << "," << (*tail)(i) << "," << (*elastic)(i) << "," << (*compton)(i) << "," << (*pileup)(i) << "," << (*escape)(i) << "\n";
        }

        file_stream.close();
    }
    else
    {
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------------

template<typename T_real>
void save_quantification(std::string path, Detector<T_real>* detector)
{
    if (detector == nullptr)
    {
        logW << "Detector == nullptr, can't save quantification\n";
    }

    //iterate through proc_type {roi, nnls, fitted}
    for (auto& itr1 : detector->fitting_quant_map)
    {
        //iterate through quantifier {sr_current, us_ic, ds_ic}
        for (auto& itr2 : itr1.second.quant_scaler_map)
        {
            std::string str_path_full = path + "calib_" + Fitting_Routine_To_Str.at(itr1.first) + "_" + itr2.first + "_K_det";
            if (detector->number() != -1)
            {
                str_path_full += std::to_string(detector->number()) + ".csv";
            }
            else
            {
                str_path_full += ".csv";
            }
            save_calibration_curve(str_path_full, detector, &(detector->quantification_standards), itr1.first, itr2.first, &(itr2.second));
        }
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
bool save_calibration_curve(std::string path, Detector<T_real>* detector, std::map<string, Quantification_Standard<T_real>>* standards, Fitting_Routines routine, string quantifier_scaler_name, Quantification_Scaler_Struct<T_real>* quants_map)
{
    if (standards == nullptr || quants_map == nullptr || detector == nullptr)
    {
        logW << "standards or quants_map or detector are null. Cannot save csv " << path << ". \n";
        return false;
    }

    std::ofstream file_stream(path);
    if (file_stream.is_open())
    {

        for (const auto& itr : detector->quantification_standards)
        {
            file_stream << "Standard Filename: " << itr.first << "\n";
            file_stream << " SR_Current: " << itr.second.sr_current << "\n";
            file_stream << " US_IC: " << itr.second.US_IC << "\n";
            file_stream << " DS_IC: " << itr.second.DS_IC << "\n";
            file_stream << "\n\n";
        }
        file_stream << "beryllium_window_thickness : " << detector->beryllium_window_thickness << "\n";
        file_stream << "germanium_dead_layer : " << detector->germanium_dead_layer << "\n";
        file_stream << "detector_chip_thickness : " << detector->detector_chip_thickness << "\n";
        file_stream << "incident_energy : " << detector->incident_energy << "\n";
        file_stream << "airpath : " << detector->airpath << "\n";
        file_stream << "detector_element : " << detector->detector_element->name << "\n";

        if (detector->avg_quantification_scaler_map.count(quantifier_scaler_name) > 0)
        {
            file_stream << quantifier_scaler_name << ": " << detector->avg_quantification_scaler_map.at(quantifier_scaler_name) << "\n";
        }

        file_stream << "\n\n";

        for (const auto& shell_itr : Shells_Quant_List)
        {
            file_stream << "\n\n";
            file_stream << "Element,Z,Counts,e_cal_ratio,absorption,transmission_Be,transmission_Ge,yield,transmission_through_Si_detector,transmission_through_air,weight  \n";

            for (const auto& itr : quants_map->curve_quant_map[shell_itr])
            {
                string name = itr.name;
                T_real counts = 0.0;
                if (shell_itr == Electron_Shell::L_SHELL)
                {
                    name += "_L";
                }
                else if (shell_itr == Electron_Shell::M_SHELL)
                {
                    name += "_M";
                }

                for (const auto& s_itr : *standards)
                {
                    if (s_itr.second.element_counts.at(routine).count(name) > 0)
                    {
                        counts = s_itr.second.element_counts.at(routine).at(name);
                        break;
                    }
                }
                
                file_stream << name << "," <<
                    itr.Z << "," <<
                    counts << "," <<
                    itr.e_cal_ratio << "," <<
                    itr.absorption << "," <<
                    itr.transmission_Be << "," <<
                    itr.transmission_Ge << "," <<
                    itr.yield << "," <<
                    itr.transmission_through_Si_detector << "," <<
                    itr.transmission_through_air << "," <<
                    itr.weight << "\n";
            }
        }
        file_stream << "\n\n";
        file_stream << "\n\n";
            
        file_stream << "Element,Z,K Shell,L Shell,M Shell\n";
        for (int i=0; i < quants_map->curve_quant_map[Electron_Shell::K_SHELL].size() ; i++)
        {
            file_stream << quants_map->curve_quant_map[Electron_Shell::K_SHELL][i].name << ","
                        << i+1 << "," 
                        << quants_map->curve_quant_map[Electron_Shell::K_SHELL][i].calib_curve_val << ","
                        << quants_map->curve_quant_map[Electron_Shell::L_SHELL][i].calib_curve_val << "," 
                        << quants_map->curve_quant_map[Electron_Shell::M_SHELL][i].calib_curve_val << "\n";
        }
        file_stream << "\n\n";
        
        file_stream.close();
    }
    else
    {
        logE << "Could not open file " << path << "\n";
        return false;
    }
    return true;
}

// ----------------------------------------------------------------------------

template<typename T_real>
bool load_element_info(std::string filename)
{
    data_struct::Element_Info_Map<T_real>* element_map = data_struct::Element_Info_Map<T_real>::inst();

    std::ifstream file_stream(filename);
    try
    {
        std::string value;
        for (std::string line; std::getline(file_stream, line); )
        {
            std::stringstream strstream(line);
            std::getline(strstream, value, ',');
            //if( std::stoi(value) > 0)
            if (value[0] >= 48 && value[0] <= 57) // 0 - 9
            {

                //logD<< "value = "<< value<<"\n";
                Element_Info<T_real>* element = nullptr;
                int element_number = std::stoi(value);
                element = element_map->get_element(element_number);
                std::string el_name;
                std::getline(strstream, el_name, ',');
                if (element == nullptr)
                {
                    element = new Element_Info<T_real>();
                    element->number = element_number;
                    element->name = el_name;
                    element_map->add_element(element);
                }
                else
                {
                    element->number = element_number;
                    element->name = el_name;
                }
                //logD<< element->number << " : "<< element->name <<"\n";
                std::getline(strstream, value, ',');
                element->xrf["ka1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ka2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["kb1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["kb2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["la1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["la2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lb1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lb2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lb3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lb4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lg1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lg2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lg3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["lg4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ll"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ln"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ma1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["ma2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["mb"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf["mg"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["k"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["l1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["l2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["l3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->yieldD["m"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ka1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ka2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["kb1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["kb2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["la1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["la2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lb1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lb2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lb3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lb4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lg1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lg2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lg3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["lg4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ll"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ln"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ma1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["ma2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["mb"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->xrf_abs_yield["mg"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->density = std::stof(value);
                std::getline(strstream, value, ',');
                element->mass = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["K"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["L1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["L2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["L3"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["M1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["M2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["M3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["M4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["M5"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["N1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N5"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N6"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["N7"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["O1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["O2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["O3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["O4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["O5"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->bindingE["P1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["P2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->bindingE["P3"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["K"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["L1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["L2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["L3"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["M1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["M2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["M3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["M4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["M5"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["N1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["N2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["N3"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["N4"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["N5"] = std::stof(value);

                std::getline(strstream, value, ',');
                element->jump["O1"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["O2"] = std::stof(value);
                std::getline(strstream, value, ',');
                element->jump["O3"] = std::stof(value);

                //element_information.emplace(std::make_pair(element->name, element));
            }
        }
    }
    catch (std::exception & e)
    {
        if (file_stream.eof() == 0 && (file_stream.bad() || file_stream.fail()))
        {
            std::cerr << "ios Exception happened: " << e.what() << "\n"
                << "Error bits are: "
                << "\nfailbit: " << file_stream.fail()
                << "\neofbit: " << file_stream.eof()
                << "\nbadbit: " << file_stream.bad() << "\n";
        }
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

}// end namespace csv
}// end namespace file
}// end namespace io
