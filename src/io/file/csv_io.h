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



#ifndef CSV_IO_H
#define CSV_IO_H

#include <fstream>

#include "data_struct/spectra.h"
#include "data_struct/element_info.h"
#include "data_struct/fit_parameters.h"
#include "data_struct/detector.h"

using namespace data_struct;

namespace io
{
namespace file
{
namespace csv
{
   
    template<typename T_real>
    DLL_EXPORT bool load_raw_spectra(std::string filename, unordered_map<string, ArrayTr<T_real>> &data)
    {
        std::ifstream file_stream(filename);
        try
        {
            std::string tmp_line;
            int num_lines;
            // First get number of lines in file so we know spectra size
            for (num_lines = 0; std::getline(file_stream, tmp_line); num_lines++)
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
                        data[names[idx]](line_num) = parse_input_real<T_real>(value);
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
    DLL_EXPORT bool save_calibration_curve(std::string path,
                                            Detector<T_real>* detector,
                                            std::map<string, Quantification_Standard<T_real>>* standards,
                                            Fitting_Routines routine,
                                            string quantifier_scaler_name, 
                                            Quantification_Scaler_Struct<T_real>* quants_map)
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
            for (int i = 0; i < quants_map->curve_quant_map[Electron_Shell::K_SHELL].size(); i++)
            {
                file_stream << quants_map->curve_quant_map[Electron_Shell::K_SHELL][i].name << ","
                    << i + 1 << ","
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
    DLL_EXPORT bool save_v9_specfit(std::string path,
        std::map<std::string, data_struct::Fit_Parameters<T_real>>& roi_files_fits_map)
    {
        const std::vector<std::string> e_list = { "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "dummy", "dummy", "Mo_L", "Tc_L", "Ru_L", "Rh_L", "Pd_L", "Ag_L", "Cd_L", "In_L", "Sn_L", "Sb_L", "Te_L", "I_L", "Xe_L", " ", "Cs_L", "Ba_L", "La_L", "Ce_L", "Pr_L", "Nd_L", "Pm_L", "Sm_L", "Eu_L", "Gd_L", "Tb_L", "Dy_L", "Br_L", "Er_L", "Tm_L", "Yb_L", "Lu_L", "Hf_L", "Ta_L", "W_L", "Re_L", "Os_L", "Ir_L", "Pt_L", "Au_L", "Hg_L", "Tl_L", "Pb_L", "Bi_L", "Po_L", "At_L", "Rn_L", "Ho_L", "Ac_L", "Th_L", "Pa_L", "U_L", "Np_L", "Pu_L", "Zr_L", "Au_M", "Pb_M", "U_M", "Hg_M", "Pt_M", "Os_M", "Bi_M", "dummy", "dummy", "real_time", "live_time", "SRcurrent", "us_IC", "ds_IC", "total_counts", "status", "niter", "total_perror", "abs_error", "relative_error", "roi_areas", "roi_pixels", "US_num", "US_unit", "US_sensfactor", "DS_num", "DS_unit", "DS_sensfactor" };

        const std::vector<std::string> p_list = { "perror_Na", "perror_Mg", "perror_Al", "perror_Si", "perror_P", "perror_S", "perror_Cl", "perror_Ar", "perror_K", "perror_Ca", "perror_Sc", "perror_Ti", "perror_V", "perror_Cr", "perror_Mn", "perror_Fe", "perror_Co", "perror_Ni", "perror_Cu", "perror_Zn", "perror_Ga", "perror_Ge", "perror_As", "perror_Se", "perror_Br", "perror_Kr", "perror_Rb", "perror_Sr", "perror_Y", "perror_Zr", "perror_Nb", "perror_Mo", "perror_Tc", "perror_Ru", "perror_Rh", "perror_Pd", "perror_Ag", "perror_Cd", "perror_In", "perror_Sn", "perror_Sb", "perror_Te", "perror_I", "perror_dummy", "perror_dummy", "perror_Mo_L", "perror_Tc_L", "perror_Ru_L", "perror_Rh_L", "perror_Pd_L", "perror_Ag_L", "perror_Cd_L", "perror_In_L", "perror_Sn_L", "perror_Sb_L", "perror_Te_L", "perror_I_L", "perror_Xe_L", "", "perror_Cs_L", "perror_Ba_L", "perror_La_L", "perror_Ce_L", "perror_Pr_L", "perror_Nd_L", "perror_Pm_L", "perror_Sm_L", "perror_Eu_L", "perror_Gd_L", "perror_Tb_L", "perror_Dy_L", "perror_Br_L", "perror_Er_L", "perror_Tm_L", "perror_Yb_L", "perror_Lu_L", "perror_Hf_L", "perror_Ta_L", "perror_W_L", "perror_Re_L", "perror_Os_L", "perror_Ir_L", "perror_Pt_L", "perror_Au_L", "perror_Hg_L", "perror_Tl_L", "perror_Pb_L", "perror_Bi_L", "perror_Po_L", "perror_At_L", "perror_Rn_L", "perror_Ho_L", "perror_Ac_L", "perror_Th_L", "perror_Pa_L", "perror_U_L", "perror_Np_L", "perror_Pu_L", "perror_Zr_L", "perror_Au_M", "perror_Pb_M", "perror_U_M", "perror_Hg_M", "perror_Pt_M", "perror_Os_M", "perror_Bi_M", "perror_dummy", "perror_dummy" };

        const std::vector<std::string> l_list = { "chisquare", "chisqred", "gen_pars_at_bndry", "ele_pars_at_bndry", "free_pars" };

        logI << "Exporting roi results to " << path << "\n";

        std::ofstream file_stream(path);
        if (file_stream.is_open())
        {

            file_stream << "spec_name, e_offset, e_linear, e_quadratic, fwhm_offset, fwhm_fanoprime, coherent_sct_energy, coherent_sct_amplitude, compton_angle, compton_fwhm_corr, compton_amplitude, compton_f_step, compton_f_tail, compton_gamma, compton_hi_f_tail, compton_hi_gamma, snip_width, si_escape, ge_escape, linear, pileup, pileup, pileup, pileup, pileup, pileup, pileup, pileup, pileup, f_step_offset, f_step_linear, f_step_quadratic, f_tail_offset, f_tail_linear, f_tail_quadratic, gamma_offset, gamma_linear, gamma_quadratic, kb_f_tail_offset, kb_f_tail_linear, kb_f_tail_quadratic, Na, Mg, Al, Si, P, S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, dummy, dummy, Mo_L, Tc_L, Ru_L, Rh_L, Pd_L, Ag_L, Cd_L, In_L, Sn_L, Sb_L, Te_L, I_L, Xe_L, , Cs_L, Ba_L, La_L, Ce_L, Pr_L, Nd_L, Pm_L, Sm_L, Eu_L, Gd_L, Tb_L, Dy_L, Br_L, Er_L, Tm_L, Yb_L, Lu_L, Hf_L, Ta_L, W_L, Re_L, Os_L, Ir_L, Pt_L, Au_L, Hg_L, Tl_L, Pb_L, Bi_L, Po_L, At_L, Rn_L, Ho_L, Ac_L, Th_L, Pa_L, U_L, Np_L, Pu_L, Zr_L, Au_M, Pb_M, U_M, Hg_M, Pt_M, Os_M, Bi_M, dummy, dummy, real_time, live_time, SRcurrent, us_IC, ds_IC, total_counts, status, niter, total_perror, abs_error, relative_error, roi_areas, roi_pixels, US_num, US_unit, US_sensfactor, DS_num, DS_unit, DS_sensfactor, perror_Na, perror_Mg, perror_Al, perror_Si, perror_P, perror_S, perror_Cl, perror_Ar, perror_K, perror_Ca, perror_Sc, perror_Ti, perror_V, perror_Cr, perror_Mn, perror_Fe, perror_Co, perror_Ni, perror_Cu, perror_Zn, perror_Ga, perror_Ge, perror_As, perror_Se, perror_Br, perror_Kr, perror_Rb, perror_Sr, perror_Y, perror_Zr, perror_Nb, perror_Mo, perror_Tc, perror_Ru, perror_Rh, perror_Pd, perror_Ag, perror_Cd, perror_In, perror_Sn, perror_Sb, perror_Te, perror_I, perror_dummy, perror_dummy, perror_Mo_L, perror_Tc_L, perror_Ru_L, perror_Rh_L, perror_Pd_L, perror_Ag_L, perror_Cd_L, perror_In_L, perror_Sn_L, perror_Sb_L, perror_Te_L, perror_I_L, perror_Xe_L, , perror_Cs_L, perror_Ba_L, perror_La_L, perror_Ce_L, perror_Pr_L, perror_Nd_L, perror_Pm_L, perror_Sm_L, perror_Eu_L, perror_Gd_L, perror_Tb_L, perror_Dy_L, perror_Br_L, perror_Er_L, perror_Tm_L, perror_Yb_L, perror_Lu_L, perror_Hf_L, perror_Ta_L, perror_W_L, perror_Re_L, perror_Os_L, perror_Ir_L, perror_Pt_L, perror_Au_L, perror_Hg_L, perror_Tl_L, perror_Pb_L, perror_Bi_L, perror_Po_L, perror_At_L, perror_Rn_L, perror_Ho_L, perror_Ac_L, perror_Th_L, perror_Pa_L, perror_U_L, perror_Np_L, perror_Pu_L, perror_Zr_L, perror_Au_M, perror_Pb_M, perror_U_M, perror_Hg_M, perror_Pt_M, perror_Os_M, perror_Bi_M, perror_dummy, perror_dummy, chisquare, chisqred, gen_pars_at_bndry, ele_pars_at_bndry, free_pars,\n";
            for (const auto& itr : roi_files_fits_map)
            {
                file_stream << itr.first << "," << itr.second.at(STR_ENERGY_OFFSET).value << "," << itr.second.at(STR_ENERGY_SLOPE).value << "," << itr.second.at(STR_ENERGY_QUADRATIC).value << "," << itr.second.at(STR_FWHM_OFFSET).value
                    << "," << itr.second.at(STR_FWHM_FANOPRIME).value << "," << itr.second.at(STR_COHERENT_SCT_ENERGY).value << "," << itr.second.at(STR_COHERENT_SCT_AMPLITUDE).value << "," << itr.second.at(STR_COMPTON_ANGLE).value
                    << "," << itr.second.at(STR_COMPTON_FWHM_CORR).value << "," << itr.second.at(STR_COMPTON_AMPLITUDE).value << "," << itr.second.at(STR_COMPTON_F_STEP).value << "," << itr.second.at(STR_COMPTON_F_TAIL).value
                    << "," << itr.second.at(STR_COMPTON_GAMMA).value << "," << itr.second.at(STR_COMPTON_HI_F_TAIL).value << "," << itr.second.at(STR_COMPTON_HI_GAMMA).value << "," << itr.second.at(STR_SNIP_WIDTH).value
                    << ",0,0,0,0,0,0,0,0,0,0,0,0," << itr.second.at(STR_F_STEP_OFFSET).value << "," << itr.second.at(STR_F_STEP_LINEAR).value << "," << itr.second.at(STR_F_STEP_QUADRATIC).value << "," << itr.second.at(STR_F_TAIL_OFFSET).value
                    << "," << itr.second.at(STR_F_TAIL_LINEAR).value << "," << itr.second.at(STR_F_TAIL_QUADRATIC).value << "," << itr.second.at(STR_GAMMA_OFFSET).value << "," << itr.second.at(STR_GAMMA_LINEAR).value
                    << "," << itr.second.at(STR_GAMMA_QUADRATIC).value << "," << itr.second.at(STR_KB_F_TAIL_OFFSET).value << "," << itr.second.at(STR_KB_F_TAIL_LINEAR).value << "," << itr.second.at(STR_KB_F_TAIL_QUADRATIC).value << ",";
                for (auto& e_itr : e_list)
                {
                    if (itr.second.contains(e_itr))
                    {
                        file_stream << itr.second.at(e_itr).value << ",";
                    }
                    else
                    {
                        file_stream << "1.0e-10, ";
                    }
                }

                for (auto& p_itr : p_list)
                {
                    file_stream << "1.0e-10, ";
                }
                
                for (auto& l_itr : l_list)
                {
                    if (itr.second.contains(l_itr))
                    {
                        file_stream << itr.second.at(l_itr).value << ",";
                    }
                    else
                    {
                        file_stream << "1.0e-10, ";
                    }
                }

                file_stream << "\n";
            }

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
    DLL_EXPORT bool save_v9_specfit_quantified(std::string path,
        Detector<T_real>* detector,
        std::map<std::string, data_struct::Fit_Parameters<double>>& roi_files_fits_map)
    {
        if (detector == nullptr)
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
            /*
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
            for (int i = 0; i < quants_map->curve_quant_map[Electron_Shell::K_SHELL].size(); i++)
            {
                file_stream << quants_map->curve_quant_map[Electron_Shell::K_SHELL][i].name << ","
                    << i + 1 << ","
                    << quants_map->curve_quant_map[Electron_Shell::K_SHELL][i].calib_curve_val << ","
                    << quants_map->curve_quant_map[Electron_Shell::L_SHELL][i].calib_curve_val << ","
                    << quants_map->curve_quant_map[Electron_Shell::M_SHELL][i].calib_curve_val << "\n";
            }
            file_stream << "\n\n";
            */
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
    DLL_EXPORT void save_quantification(std::string path, Detector<T_real>* detector)
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
    DLL_EXPORT bool save_fit_and_int_spectra_labeled(const std::string fullpath, const data_struct::ArrayTr<T_real>* energy, const data_struct::ArrayTr<T_real>* spectra, const data_struct::ArrayTr<T_real>* spectra_model, const data_struct::ArrayTr<T_real>* background, unordered_map<string, data_struct::ArrayTr<T_real>>* labeled_spectras)
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
            if (labeled_spectras->count(STR_M_LINES) > 0)
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
    DLL_EXPORT bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<T_real>* energy, const data_struct::ArrayTr<T_real>* spectra, const data_struct::ArrayTr<T_real>* spectra_model, const data_struct::ArrayTr<T_real>* background, unordered_map<string, ArrayTr<T_real>>* labeled_spectras = nullptr)
    {
        return save_fit_and_int_spectra_labeled(fullpath, energy, spectra, spectra_model, background, labeled_spectras);
    }

    // ----------------------------------------------------------------------------

    template<typename T_real>
    DLL_EXPORT bool load_element_info(std::string filename)
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

                    std::getline(strstream, value, ',');
                    element->xrf["ka1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["ka2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["kb1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["kb2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["la1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["la2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["lb1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["lb2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["lb3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["lb4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["lg1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["lg2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["lg3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["lg4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["ll"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["ln"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["ma1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["ma2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["mb"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf["mg"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->yieldD["k"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->yieldD["l1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->yieldD["l2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->yieldD["l3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->yieldD["m"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["ka1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["ka2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["kb1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["kb2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["la1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["la2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["lb1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["lb2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["lb3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["lb4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["lg1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["lg2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["lg3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["lg4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["ll"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["ln"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["ma1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["ma2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["mb"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->xrf_abs_yield["mg"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->density = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->mass = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->bindingE["K"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->bindingE["L1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["L2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["L3"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->bindingE["M1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["M2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["M3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["M4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["M5"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->bindingE["N1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["N2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["N3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["N4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["N5"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["N6"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["N7"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->bindingE["O1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["O2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["O3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["O4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["O5"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->bindingE["P1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["P2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->bindingE["P3"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->jump["K"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->jump["L1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["L2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["L3"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->jump["M1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["M2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["M3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["M4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["M5"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->jump["N1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["N2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["N3"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["N4"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["N5"] = parse_input_real<T_real>(value);

                    std::getline(strstream, value, ',');
                    element->jump["O1"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["O2"] = parse_input_real<T_real>(value);
                    std::getline(strstream, value, ',');
                    element->jump["O3"] = parse_input_real<T_real>(value);

                    //element_information.emplace(std::make_pair(element->name, element));
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

}// end namespace CSV
}// end namespace file
}// end namespace io

#endif // CSV_IO_H
