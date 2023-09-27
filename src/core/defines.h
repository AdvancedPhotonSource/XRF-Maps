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
#ifndef __XRF_DEFINES__
#define __XRF_DEFINES__

#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>

static std::time_t now_c;
#define logit now_c=std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());std::cout<<std::put_time(std::localtime(&now_c),"[%F_%T]\t")<<__FILE__<<"::"<<__FUNCTION__<<"():"<<__LINE__<<"\t"
#define logit_s std::cout
#define logE logit<<"Error: "
#define logW logit<<"Warning: "
#define logI logit<<"Info: "

#if defined _WIN32 || defined __CYGWIN__
  #pragma warning( disable : 4251 4127 4996 4505 4244 )
  #define DIR_END_CHAR '\\'
  #define DIR_END_CHAR_OPPOSITE '/'
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define DLL_EXPORT __attribute__ ((dllexport))
    #else
      #define DLL_EXPORT __declspec(dllexport) 
    #endif
  #else
    #ifdef __GNUC__
      #define DLL_EXPORT __attribute__ ((dllimport))
    #else
      #define DLL_EXPORT __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
    #endif
  #endif
  #define DLL_LOCAL
#else
  #define DIR_END_CHAR '/'
  #define DIR_END_CHAR_OPPOSITE '\\'
  #if __GNUC__ >= 4
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define DLL_EXPORT
    #define DLL_LOCAL
  #endif
#endif


#if defined _WIN32 || defined __CYGWIN__
#define TEMPLATE_CLASS_DLL_EXPORT template DLL_EXPORT class
#define TEMPLATE_STRUCT_DLL_EXPORT template DLL_EXPORT struct
#else
#define TEMPLATE_CLASS_DLL_EXPORT template class DLL_EXPORT
#define TEMPLATE_STRUCT_DLL_EXPORT template struct DLL_EXPORT
#endif

template<typename T_real>
T_real parse_input_real(std::string value)
{
    if (std::is_same<T_real, float>::value)
    {
        return std::stof(value);
    }
    else if (std::is_same<T_real, double>::value)
    {
        return std::stod(value);
    }
}

#define CALIBRATION_CURVE_SIZE 92

//namespace keys
//{

#define AVOGADRO 6.02204531e23
#define HC_ANGSTROMS 12398.52
#define RE 2.817938070e-13		// in cm
#define ENERGY_RES_OFFSET 150.0
#define ENERGY_RES_SQRT 12.0

/**
* @brief String defines for fit parameters string value pair.
*/
const std::string STR_MAX_ENERGY_TO_FIT = "MAX_ENERGY_TO_FIT";
const std::string STR_MIN_ENERGY_TO_FIT = "MIN_ENERGY_TO_FIT";

const std::string STR_SI_ESCAPE = "SI_ESCAPE";
const std::string STR_GE_ESCAPE = "GE_ESCAPE";
const std::string STR_ESCAPE_LINEAR = "ESCAPE_LINEAR";
const std::string STR_PILEUP0 = "PILEUP0";
const std::string STR_PILEUP1 = "PILEUP1";
const std::string STR_PILEUP2 = "PILEUP2";
const std::string STR_PILEUP3 = "PILEUP3";
const std::string STR_PILEUP4 = "PILEUP4";
const std::string STR_PILEUP5 = "PILEUP5";
const std::string STR_PILEUP6 = "PILEUP6";
const std::string STR_PILEUP7 = "PILEUP7";
const std::string STR_PILEUP8 = "PILEUP8";

const std::string STR_NUM_ITR = "Num_Iter";
const std::string STR_RESIDUAL = "Fit_Residual";
const std::string STR_OUTCOME = "Fit_Outcome";
const std::string STR_CHISQUARE = "chisquare";
const std::string STR_CHISQRED = "chisqred";
const std::string STR_FREE_PARS = "free_pars";
const std::string STR_DETECTOR_ELEMENT = "DETECTOR_ELEMENT";
const std::string STR_BE_WINDOW_THICKNESS = "BE_WINDOW_THICKNESS";
const std::string STR_DET_CHIP_THICKNESS = "DET_CHIP_THICKNESS";
const std::string STR_GE_DEAD_LAYER = "GE_DEAD_LAYER";

//namespace gaussian_model
//{
const std::string STR_ENERGY_OFFSET = "ENERGY_OFFSET";
const std::string STR_ENERGY_SLOPE = "ENERGY_SLOPE";
const std::string STR_ENERGY_QUADRATIC = "ENERGY_QUADRATIC";

const std::string STR_FWHM_OFFSET = "FWHM_OFFSET";
const std::string STR_FWHM_FANOPRIME = "FWHM_FANOPRIME";

const std::string STR_COHERENT_SCT_ENERGY = "COHERENT_SCT_ENERGY";
const std::string STR_COHERENT_SCT_AMPLITUDE = "COHERENT_SCT_AMPLITUDE";

const std::string STR_COMPTON_ANGLE = "COMPTON_ANGLE";
const std::string STR_COMPTON_FWHM_CORR = "COMPTON_FWHM_CORR";
const std::string STR_COMPTON_AMPLITUDE = "COMPTON_AMPLITUDE";
const std::string STR_COMPTON_F_STEP = "COMPTON_F_STEP";
const std::string STR_COMPTON_F_TAIL = "COMPTON_F_TAIL";
const std::string STR_COMPTON_GAMMA = "COMPTON_GAMMA";
const std::string STR_COMPTON_HI_F_TAIL = "COMPTON_HI_F_TAIL";
const std::string STR_COMPTON_HI_GAMMA = "COMPTON_HI_GAMMA";

const std::string STR_SNIP_WIDTH = "SNIP_WIDTH";
const std::string STR_FIT_SNIP_WIDTH = "FIT_SNIP_WIDTH";

const std::string STR_F_STEP_OFFSET = "F_STEP_OFFSET";
const std::string STR_F_STEP_LINEAR = "F_STEP_LINEAR";
const std::string STR_F_STEP_QUADRATIC = "F_STEP_QUADRATIC";

const std::string STR_F_TAIL_OFFSET = "F_TAIL_OFFSET";
const std::string STR_F_TAIL_LINEAR = "F_TAIL_LINEAR";
const std::string STR_F_TAIL_QUADRATIC = "F_TAIL_QUADRATIC";

const std::string STR_GAMMA_OFFSET = "GAMMA_OFFSET";
const std::string STR_GAMMA_LINEAR = "GAMMA_LINEAR";
const std::string STR_GAMMA_QUADRATIC = "GAMMA_QUADRATIC";

const std::string STR_KB_F_TAIL_OFFSET = "KB_F_TAIL_OFFSET";
const std::string STR_KB_F_TAIL_LINEAR = "KB_F_TAIL_LINEAR";
const std::string STR_KB_F_TAIL_QUADRATIC = "KB_F_TAIL_QUADRATIC";

const std::string STR_SUM_ELASTIC_INELASTIC_AMP = "Sum_Elastic_Inelastic";

const std::string STR_TOTAL_FLUORESCENCE_YIELD = "Total_Fluorescence_Yield";

const std::string STR_FIT_ROI = "ROI";
const std::string STR_FIT_SVD = "SVD";
const std::string STR_FIT_NNLS = "NNLS";
const std::string STR_FIT_GAUSS_MATRIX = "Fitted";
const std::string STR_FIT_GAUSS_TAILS = "gaussian_parameter";
const std::string STR_FIT_GAUSS_NNLS_TAILS = "Hybrid_NNLS";

const std::string STR_SR_CURRENT = "SR_Current";
const std::string STR_US_IC = "US_IC";
const std::string STR_DS_IC = "DS_IC";
const std::string STR_CFG_2 = "CFG_2";
const std::string STR_CFG_3 = "CFG_3";
const std::string STR_CFG_4 = "CFG_4";
const std::string STR_CFG_5 = "CFG_5";


const std::string STR_K_A_LINES = "K Alpha";
const std::string STR_K_B_LINES = "K Beta";
const std::string STR_L_LINES = "L Lines";
const std::string STR_M_LINES = "M Lines";
const std::string STR_STEP_LINES = "Step";
const std::string STR_TAIL_LINES = "Tail";
const std::string STR_ELASTIC_LINES = "Elastic";
const std::string STR_COMPTON_LINES = "Compton";
const std::string STR_PILEUP_LINES = "Pile Up";
const std::string STR_ESCAPE_LINES = "Escape";

const std::string STR_INT_SPEC = "Integrated_Spectra";
const std::string STR_FIT_INT_SPEC = "Fitted_Integrated_Spectra";
const std::string STR_MAX_CHANNELS_INT_SPEC = "Max_Channels_Integrated_Spectra";
const std::string STR_MAX10_INT_SPEC = "Max_10_Channels_Integrated_Spectra";
const std::string STR_FIT_INT_BACKGROUND = "FIT_Integrated_Background";

const std::string STR_CALIB_CURVE_SR_CUR = "Calibration_Curve_SR_Current";
const std::string STR_CALIB_CURVE_US_IC = "Calibration_Curve_US_IC";
const std::string STR_CALIB_CURVE_DS_IC = "Calibration_Curve_DS_IC";
const std::string STR_CALIB_LABELS = "Calibration_Curve_Labels";

const std::string STR_DS_IC_ELEMENT_INFO_VALUES = "DS_IC_Element_Info_Values";
const std::string STR_US_IC_ELEMENT_INFO_VALUES = "US_IC_Element_Info_Values";
const std::string STR_SR_CURRENT_ELEMENT_INFO_VALUES = "SR_Current_Element_Info_Values";

const std::string STR_BEAMLINES = "BeamLines";
const std::string STR_SCALERS = "Scalers";
const std::string STR_TIME_NORMALIZED_SCALERS = "TimeNormalizedScalers";
const std::string STR_SUMMED_SCALERS = "Summed_Scalers";
const std::string STR_TIMING = "Timing";
const std::string STR_ELT = "ELT";
const std::string STR_ERT = "ERT";
const std::string STR_ICR = "ICR";
const std::string STR_OCR = "OCR";
const std::string STR_DEAD_TIME = "Dead_Time";

const std::string STR_GENERAL_BEAMLINE = "General";

const std::string STR_ENERGY = "Energy";
const std::string STR_ENERGY_CALIB = "Energy_Calibration";

const std::string STR_ELAPSED_REAL_TIME = "Elapsed_Realtime";
const std::string STR_ELAPSED_LIVE_TIME = "Elapsed_Livetime";
const std::string STR_INPUT_COUNTS = "Input_Counts";
const std::string STR_OUTPUT_COUNTS = "Output_Counts";

const std::string STR_MAPS = "MAPS";
const std::string STR_SCAN = "Scan";
const std::string STR_SPECTRA = "Spectra";
const std::string STR_VERSION = "version";
const std::string STR_XRF_ANALYZED = "XRF_Analyzed";
const std::string STR_COUNTS_PER_SEC = "Counts_Per_Sec";
const std::string STR_CHANNEL_NAMES = "Channel_Names";
const std::string STR_CHANNEL_UNITS = "Channel_Units";
const std::string STR_FIT_PARAMETERS_OVERRIDE = "Fit_Parameters_Override";
const std::string STR_FIT_PARAMETERS = "Fit_Parameters";

const std::string STR_QUANTIFICATION = "Quantification";
const std::string STR_CALIBRATION = "Calibration";
const std::string STR_NUMBER_OF_STANDARDS = "Number_Of_Standards";
const std::string STR_ELEMENT_WEIGHTS = "Element_Weights";
const std::string STR_ELEMENT_WEIGHTS_NAMES = "Element_Weights_Names";
const std::string STR_STANDARD_NAME = "Standard_Name";

const std::string STR_Y_AXIS = "y_axis";
const std::string STR_X_AXIS = "x_axis";
const std::string STR_REQUESTED_ROWS = "requested_rows";
const std::string STR_REQUESTED_COLS = "requested_cols";
const std::string STR_THETA = "theta";
const std::string STR_SCAN_TIME_STAMP = "scan_time_stamp";
const std::string STR_NAME = "name";


const std::string STR_EXTRA_PVS = "Extra_PVs";
const std::string STR_NAMES = "Names";
const std::string STR_VALUES = "Values";
const std::string STR_DESCRIPTION = "Description";
const std::string STR_UNIT = "Unit";
const std::string STR_UNITS = "Units";

const std::string STR_US_AMP = "us_amp";
const std::string STR_DS_AMP = "ds_amp";
const std::string STR_US_AMP_NUM = "us_amp_num";
const std::string STR_DS_AMP_NUM = "ds_amp_num";
const std::string STR_US_AMP_UNIT = "us_amp_unit";
const std::string STR_DS_AMP_UNIT = "ds_amp_unit";
const std::string STR_US_AMP_NUM_UPPR = "US_AMP_NUM";
const std::string STR_DS_AMP_NUM_UPPR = "DS_AMP_NUM";
const std::string STR_US_AMP_UNIT_UPPR = "US_AMP_UNIT";
const std::string STR_DS_AMP_UNIT_UPPR = "DS_AMP_UNIT";

const std::string STR_NBS_1832 = "nbs1832";

// ROI
const std::string STR_MAPS_ROIS = "MAPS_ROIS";
const std::string STR_MAP_ROI_NAME = "Name";
const std::string STR_MAP_ROI_COLOR = "Color";
const std::string STR_MAP_ROI_COLOR_ALPHA = "Color_Alpha";
const std::string STR_MAP_ROI_PIXEL_LOC = "Pixel_Loc";
const std::string STR_MAP_ROI_INT_SPEC = "Integrated_Spectras";
const std::string STR_MAP_ROI_INT_SPEC_FILENAME = "File_Name";

//}

//}// namespace keys




#endif
