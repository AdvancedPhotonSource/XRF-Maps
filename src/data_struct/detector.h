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



#ifndef DETECTOR_H
#define DETECTOR_H

#include "core/defines.h"
#include "data_struct/element_info.h"
#include "fitting/routines/base_fit_routine.h"
#include <vector>
#include <string>
#include "data_struct/quantification_standard.h"
#include "data_struct/params_override.h"
#include "fitting/optimizers/lmfit_optimizer.h"
#include "fitting/optimizers/mpfit_optimizer.h"

namespace data_struct
{


//-----------------------------------------------------------------------------

///
/// \brief The Detector class
///
template<typename T_real>
class DLL_EXPORT Detector
{
public:
    Detector(unsigned int detector_num = -1);

    ~Detector();

    void append_element(Fitting_Routines routine, string name, string quant_scaler, T_real weight);

    void update_element_quants(Fitting_Routines routine,
                                string quantifier_scaler,
                                Quantification_Standard<T_real>* standard,
                                Quantification_Model<T_real>* quantification_model,
                                T_real ic_quantifier);

    void update_calibration_curve(Fitting_Routines routine,
                                string quantifier_scaler,
                                Quantification_Model<T_real>* quantification_model,
                                T_real val);

    void update_from_fit_paramseters();

    void generage_avg_quantification_scalers();

    unsigned int number() { return _number; }

    // Fitting routines map
    std::unordered_map<Fitting_Routines, fitting::routines::Base_Fit_Routine<T_real> *> fit_routines;

    // Fitting model
    fitting::models::Base_Model<T_real>* model;

    // Quantification
    std::map<string, Quantification_Standard<T_real>> quantification_standards;

    unordered_map <Fitting_Routines, struct Fitting_Quantification_Struct<T_real>> fitting_quant_map;

    //  proc_type          quantifier            element    quant_prop
    map<Fitting_Routines, map<string, unordered_map<string, Element_Quant<T_real>*>>> all_element_quants;

    // Fit Parameters Override for model
    Params_Override<T_real> fit_params_override_dict;

    T_real beryllium_window_thickness;
    T_real germanium_dead_layer;
    T_real detector_chip_thickness;
    T_real incident_energy;
    T_real airpath;
    data_struct::Element_Info<T_real>* detector_element;

    // SR_CURRENT, US_IC, DS_IC  : average if we have multiple standards
    unordered_map<string, T_real> avg_quantification_scaler_map;

private:
    unsigned int _number;

};

#if defined _WIN32 || defined __CYGWIN__
template DLL_EXPORT class Detector<float>;
template DLL_EXPORT class Detector<double>;
#else
template class DLL_EXPORT Detector<float>;
template class DLL_EXPORT Detector<double>;
#endif

} //namespace data_struct

#endif // DETECTOR_H
