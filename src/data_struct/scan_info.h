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

#ifndef SCAN_INFO_H
#define SCAN_INFO_H


#include <map>
#include <string>
#include <vector>
#include "core/defines.h"

#include <Eigen/core>

namespace data_struct
{

using namespace std;

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

struct Extra_PV
{
    string name;
    string description;
    string value;
    string unit;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

struct Scaler_Map
{
    string name;
    string unit;
    Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> values;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

struct Scan_Meta_Info
{
    string name;
    string scan_time_stamp;
    vector<real_t> x_axis;
    vector<real_t> y_axis;
    int requested_cols;
    int requested_rows;
    real_t theta;
    
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
/**
 * @brief The Params_Override struct
 */
class DLL_EXPORT Scan_Info
{

public:

    Scan_Info()
    {
        
    }
    
    ~Scan_Info()
    {
        
    }

    Scan_Meta_Info meta_info;
    vector<Scaler_Map> scaler_maps;
    vector<Extra_PV> extra_pvs;
    bool has_netcdf; 
};

//-----------------------------------------------------------------------------

} //namespace data_struct

#endif // Scan_Info_H
