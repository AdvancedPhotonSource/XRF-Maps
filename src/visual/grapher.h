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
#include "data_struct/quantification_standard.h"

namespace visual
{

void SavePlotSpectras(std::string path,
                      data_struct::ArrayXr *energy,
                      data_struct::ArrayXr *spectra,
                      data_struct::ArrayXr *model,
                      data_struct::ArrayXr *background,
                      bool log_them);

void SavePlotQuantification(std::string path, data_struct::Quantification_Standard *standard, int detector_num, int zstart=13, int zstop=30);

void SavePlotCalibrationCurve(std::string path, data_struct::Calibration_Curve *calib_curve, map<int, real_t> fitted_e_cal_raitos, int zstart, int zstop);

}
