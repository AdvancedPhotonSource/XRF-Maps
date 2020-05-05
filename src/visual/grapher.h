/***

Copyright (c) 2016 Arthur Glowacki

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

***/


#include <valarray>
#include "core/defines.h"
#include "data_struct/spectra.h"
#include "data_struct/detector.h"

using namespace data_struct;

namespace visual
{

void SavePlotSpectras(std::string path,
                      data_struct::ArrayXr *energy,
                      data_struct::ArrayXr *spectra,
                      data_struct::ArrayXr *model,
                      data_struct::ArrayXr *background,
                      bool log_them);

void find_shell_Z_offset(quantification::models::Electron_Shell shell_idx, unordered_map<string, Element_Quant*>* all_elements_with_weights, int& zstart, int& zstop);

bool contains_shell(quantification::models::Electron_Shell shell_idx, unordered_map<string, Element_Quant*>* element_quants);

void SavePlotQuantification(std::string path, Detector* detector);

void SavePlotCalibrationCurve(std::string path,
                              Detector* detector,
                              string quantifier_scaler_name,
                              unordered_map<string, Element_Quant*>* all_elements_with_weights,
                              vector<Element_Quant>* calibration_curve,
                              int zstart,
                              int zstop);

//void SavePlotCalibrationCurve(std::string path, std::string standard_name, data_struct::Calibration_Curve *calib_curve, map<string, real_t> fitted_e_cal_raitos,  int shell_idx, int zstart, int zstop);

}
