/***
-Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.
-
-Copyright 2016. UChicago Argonne, LLC. This software was produced
-under U.S. Government contract DE-AC02-06CH11357 for Argonne National
-Laboratory (ANL), which is operated by UChicago Argonne, LLC for the
-U.S. Department of Energy. The U.S. Government has rights to use,
-reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
-UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
-ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
-modified to produce derivative works, such modified software should
-be clearly marked, so as not to confuse it with the version available
-from ANL.
-
-Additionally, redistribution and use in source and binary forms, with
-or without modification, are permitted provided that the following
-conditions are met:
-
-    * Redistributions of source code must retain the above copyright
-      notice, this list of conditions and the following disclaimer.
-
-    * Redistributions in binary form must reproduce the above copyright
-      notice, this list of conditions and the following disclaimer in
-      the documentation and/or other materials provided with the
-      distribution.
-
-    * Neither the name of UChicago Argonne, LLC, Argonne National
-      Laboratory, ANL, the U.S. Government, nor the names of its
-      contributors may be used to endorse or promote products derived
-      from this software without specific prior written permission.
-
-THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS
-"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
-LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
-FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago
-Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
-INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
-BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
-LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
-CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
-LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
-ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
-POSSIBILITY OF SUCH DAMAGE.
-***/

#ifndef PARAMS_OVERRIDE_H
#define PARAMS_OVERRIDE_H


#include <map>
#include <string>
#include <list>
#include "core/defines.h"

#include "data_struct/fit_parameters.h"
#include "data_struct/fit_element_map.h"

namespace data_struct
{

using namespace std;

struct Summed_Scaler
{
  string scaler_name;
  list<string> scalers_to_sum;
};

//-----------------------------------------------------------------------------
/**
 * @brief The Params_Override struct
 */
class DLL_EXPORT Params_Override
{

public:

    Params_Override()
    {
        detector_num = 0;
        si_escape_factor = 0.0;
        ge_escape_factor = 0.0;
        si_escape_enabled = false;
        ge_escape_enabled = false;
        fit_snip_width = 0.0;
        us_amp_sens_num = 0.0;
        us_amp_sens_unit = 0.0;
        ds_amp_sens_num = 0.0;
        ds_amp_sens_unit = 0.0;
        theta_pv = "";
    }
    Params_Override(string dir, int detector)
    {

        si_escape_factor = 0.0;
        ge_escape_factor = 0.0;
        si_escape_enabled = false;
        ge_escape_enabled = false;
        fit_snip_width = 0.0;
        us_amp_sens_num = 0.0;
        us_amp_sens_unit = 0.0;
        ds_amp_sens_num = 0.0;
        ds_amp_sens_unit = 0.0;
        theta_pv = "";

        dataset_directory = dir;
        detector_num = detector;
    }
    ~Params_Override()
    {
        time_normalized_scalers.clear();
        scaler_pvs.clear();
    }

    string dataset_directory;
    int detector_num;
    Fit_Parameters fit_params;
    Fit_Element_Map_Dict elements_to_fit;
    string time_scaler;
    string time_scaler_clock;
    map< string, string > time_normalized_scalers;
    string detector_element;
    map< string, string > scaler_pvs;
    real_t si_escape_factor;
    real_t ge_escape_factor;
    bool si_escape_enabled;
    bool ge_escape_enabled;
    real_t fit_snip_width;

    string be_window_thickness;
    string det_chip_thickness;
    string ge_dead_layer;

    string elt_pv;
    string ert_pv;
    string in_cnt_pv;
    string out_cnt_pv;

    string us_amp_sens_num_pv;
    string us_amp_sens_unit_pv;
    string ds_amp_sens_num_pv;
    string ds_amp_sens_unit_pv;

    string theta_pv;

    list<struct Summed_Scaler> summed_scalers;

    real_t us_amp_sens_num;
    real_t us_amp_sens_unit;
    real_t ds_amp_sens_num;
    real_t ds_amp_sens_unit;

};

//-----------------------------------------------------------------------------

} //namespace data_struct

#endif // PARAMS_OVERRIDE_H
