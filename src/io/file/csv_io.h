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
    DLL_EXPORT bool load_element_info(std::string filename);

    template<typename T_real>
    DLL_EXPORT bool load_raw_spectra(std::string filename, unordered_map<string, ArrayTr<T_real>> &data);

    template DLL_EXPORT bool load_raw_spectra(std::string filename, unordered_map<string, ArrayTr<float>>& data);
    template DLL_EXPORT bool load_raw_spectra(std::string filename, unordered_map<string, ArrayTr<double>>& data);

    template<typename T_real>
    DLL_EXPORT bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<T_real>* energy, const data_struct::ArrayTr<T_real>* spectra, const data_struct::ArrayTr<T_real>* spectra_model, const data_struct::ArrayTr<T_real>* background);

    template DLL_EXPORT bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<float>* energy, const data_struct::ArrayTr<float>* spectra, const data_struct::ArrayTr<float>* spectra_model, const data_struct::ArrayTr<float>* background);
    template DLL_EXPORT bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<double>* energy, const data_struct::ArrayTr<double>* spectra, const data_struct::ArrayTr<double>* spectra_model, const data_struct::ArrayTr<double>* background);

    template<typename T_real>
    DLL_EXPORT bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<T_real>* energy, const data_struct::ArrayTr<T_real>* spectra, const data_struct::ArrayTr<T_real>* spectra_model, const data_struct::ArrayTr<T_real>* background, unordered_map<string, data_struct::ArrayTr<T_real>>* labeled_spectras);

    template DLL_EXPORT bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<float>* energy, const data_struct::ArrayTr<float>* spectra, const data_struct::ArrayTr<float>* spectra_model, const data_struct::ArrayTr<float>* background, unordered_map<string, data_struct::ArrayTr<float>>* labeled_spectras);
    template DLL_EXPORT bool save_fit_and_int_spectra(const std::string fullpath, const data_struct::ArrayTr<double>* energy, const data_struct::ArrayTr<double>* spectra, const data_struct::ArrayTr<double>* spectra_model, const data_struct::ArrayTr<double>* background, unordered_map<string, data_struct::ArrayTr<double>>* labeled_spectras);

    template<typename T_real>
    DLL_EXPORT void save_quantification(std::string path, Detector<T_real>* detector);

    template DLL_EXPORT void save_quantification(std::string path, Detector<float>* detector);
    template DLL_EXPORT void save_quantification(std::string path, Detector<double>* detector);

    template<typename T_real>
    DLL_EXPORT bool save_calibration_curve(std::string path,
                                            Detector<T_real>* detector,
                                            std::map<string, Quantification_Standard<T_real>>* standards,
                                            Fitting_Routines routine,
                                            string quantifier_scaler_name, 
                                            Quantification_Scaler_Struct<T_real>* quants_map);

    template DLL_EXPORT bool save_calibration_curve(std::string path, Detector<float>* detector, std::map<string, Quantification_Standard<float>>* standards, Fitting_Routines routine, string quantifier_scaler_name, Quantification_Scaler_Struct<float>* quants_map);
    template DLL_EXPORT bool save_calibration_curve(std::string path, Detector<double>* detector, std::map<string, Quantification_Standard<double>>* standards, Fitting_Routines routine, string quantifier_scaler_name, Quantification_Scaler_Struct<double>* quants_map);

}// end namespace CSV
}// end namespace file
}// end namespace io

#endif // CSV_IO_H
