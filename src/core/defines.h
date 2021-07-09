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

#if defined _REAL_FLOAT
  #define real_t float
  #define H5T_NATIVE_REAL H5T_NATIVE_FLOAT
  #define H5T_INTEL_R H5T_INTEL_F32
#elif defined _REAL_DOUBLE
  #define real_t double
  #define H5T_NATIVE_REAL H5T_NATIVE_DOUBLE
  #define H5T_INTEL_R H5T_INTEL_F64
#endif

#if defined _WIN32 || defined __CYGWIN__
  #pragma warning( disable : 4251 4127 4996 4505 4244 )
  #define DIR_END_CHAR '\\'
  #define DIR_END_CHAR_OPPOSITE '/'
  #ifdef BUILDING_DLL
    #ifdef __GNUC__
      #define DLL_EXPORT __attribute__ ((dllexport))
    #else
      #define DLL_EXPORT __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
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

// STRING KEYS
using namespace std;

#define CALIBRATION_CURVE_SIZE 92

//namespace keys
//{

/**
* @brief String defines for fit parameters string value pair.
*/
const string STR_MAX_ENERGY_TO_FIT = "MAX_ENERGY_TO_FIT";
const string STR_MIN_ENERGY_TO_FIT = "MIN_ENERGY_TO_FIT";

const string STR_SI_ESCAPE = "SI_ESCAPE";
const string STR_GE_ESCAPE = "GE_ESCAPE";
const string STR_ESCAPE_LINEAR = "ESCAPE_LINEAR";
const string STR_PILEUP0 = "PILEUP0";
const string STR_PILEUP1 = "PILEUP1";
const string STR_PILEUP2 = "PILEUP2";
const string STR_PILEUP3 = "PILEUP3";
const string STR_PILEUP4 = "PILEUP4";
const string STR_PILEUP5 = "PILEUP5";
const string STR_PILEUP6 = "PILEUP6";
const string STR_PILEUP7 = "PILEUP7";
const string STR_PILEUP8 = "PILEUP8";

const string STR_NUM_ITR = "Num_Iter";
const string STR_RESIDUAL = "Fit_Residual";
const string STR_OUTCOME = "Fit_Outcome";
const string STR_DETECTOR_ELEMENT = "DETECTOR_ELEMENT";
const string STR_BE_WINDOW_THICKNESS = "BE_WINDOW_THICKNESS";
const string STR_DET_CHIP_THICKNESS = "DET_CHIP_THICKNESS";
const string STR_GE_DEAD_LAYER = "GE_DEAD_LAYER";

//namespace gaussian_model
//{
const string STR_ENERGY_OFFSET = "ENERGY_OFFSET";
const string STR_ENERGY_SLOPE = "ENERGY_SLOPE";
const string STR_ENERGY_QUADRATIC = "ENERGY_QUADRATIC";

const string STR_FWHM_OFFSET = "FWHM_OFFSET";
const string STR_FWHM_FANOPRIME = "FWHM_FANOPRIME";

const string STR_COHERENT_SCT_ENERGY = "COHERENT_SCT_ENERGY";
const string STR_COHERENT_SCT_AMPLITUDE = "COHERENT_SCT_AMPLITUDE";

const string STR_COMPTON_ANGLE = "COMPTON_ANGLE";
const string STR_COMPTON_FWHM_CORR = "COMPTON_FWHM_CORR";
const string STR_COMPTON_AMPLITUDE = "COMPTON_AMPLITUDE";
const string STR_COMPTON_F_STEP = "COMPTON_F_STEP";
const string STR_COMPTON_F_TAIL = "COMPTON_F_TAIL";
const string STR_COMPTON_GAMMA = "COMPTON_GAMMA";
const string STR_COMPTON_HI_F_TAIL = "COMPTON_HI_F_TAIL";
const string STR_COMPTON_HI_GAMMA = "COMPTON_HI_GAMMA";

const string STR_SNIP_WIDTH = "SNIP_WIDTH";
const string STR_FIT_SNIP_WIDTH = "FIT_SNIP_WIDTH";

const string STR_F_STEP_OFFSET = "F_STEP_OFFSET";
const string STR_F_STEP_LINEAR = "F_STEP_LINEAR";
const string STR_F_STEP_QUADRATIC = "F_STEP_QUADRATIC";

const string STR_F_TAIL_OFFSET = "F_TAIL_OFFSET";
const string STR_F_TAIL_LINEAR = "F_TAIL_LINEAR";
const string STR_F_TAIL_QUADRATIC = "F_TAIL_QUADRATIC";

const string STR_GAMMA_OFFSET = "GAMMA_OFFSET";
const string STR_GAMMA_LINEAR = "GAMMA_LINEAR";
const string STR_GAMMA_QUADRATIC = "GAMMA_QUADRATIC";

const string STR_KB_F_TAIL_OFFSET = "KB_F_TAIL_OFFSET";
const string STR_KB_F_TAIL_LINEAR = "KB_F_TAIL_LINEAR";
const string STR_KB_F_TAIL_QUADRATIC = "KB_F_TAIL_QUADRATIC";

const string STR_SUM_ELASTIC_INELASTIC_AMP = "Sum_Elastic_Inelastic";

const string STR_TOTAL_FLUORESCENCE_YIELD = "Total_Fluorescence_Yield";

const string STR_FIT_ROI = "ROI";
const string STR_FIT_SVD = "SVD";
const string STR_FIT_NNLS = "NNLS";
const string STR_FIT_GAUSS_MATRIX = "Fitted";
const string STR_FIT_GAUSS_TAILS = "gaussian_parameter";

const string STR_SR_CURRENT = "SR_Current";
const string STR_US_IC = "US_IC";
const string STR_DS_IC = "DS_IC";

const string STR_K_A_LINES = "K Alpha";
const string STR_K_B_LINES = "K Beta";
const string STR_L_LINES = "L Lines";
const string STR_M_LINES = "M Lines";
const string STR_STEP_LINES = "Step";
const string STR_TAIL_LINES = "Tail";
const string STR_ELASTIC_LINES = "Elastic";
const string STR_COMPTON_LINES = "Compton";
const string STR_PILEUP_LINES = "Pile Up";
const string STR_ESCAPE_LINES = "Escape";

const string STR_INT_SPEC = "Integrated_Spectra";
const string STR_FIT_INT_SPEC = "Fitted_Integrated_Spectra";
const string STR_MAX_CHANNELS_INT_SPEC = "Max_Channels_Integrated_Spectra";
const string STR_MAX10_INT_SPEC = "Max_10_Channels_Integrated_Spectra";
const string STR_FIT_INT_BACKGROUND = "FIT_Integrated_Background";

const string STR_CALIB_CURVE_SR_CUR = "Calibration_Curve_SR_Current";
const string STR_CALIB_CURVE_US_IC = "Calibration_Curve_US_IC";
const string STR_CALIB_CURVE_DS_IC = "Calibration_Curve_DS_IC";
const string STR_CALIB_LABELS = "Calibration_Curve_Labels";
//}

//}// namespace keys




#endif
