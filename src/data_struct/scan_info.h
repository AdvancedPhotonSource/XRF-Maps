/***
-Copyright (c) 2020, UChicago Argonne, LLC. All rights reserved.
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

#include "fit_parameters.h" // for ArrayXXr 

namespace data_struct
{

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

struct Extra_PV
{
    std::string name;
    std::string description;
    std::string value;
    std::string unit;
};

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
template<typename T_real>
struct Scaler_Map
{
    std::string name;
    std::string unit;
    bool time_normalized;
    //bool is_timer;
    ArrayXXr<T_real> values;
};

TEMPLATE_STRUCT_DLL_EXPORT Scaler_Map<float>;
TEMPLATE_STRUCT_DLL_EXPORT Scaler_Map<double>;

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<typename T_real>
struct Scan_Meta_Info
{
    std::string name;
    std::string scan_time_stamp;
    std::string scan_type;
    ArrayTr<T_real> x_axis;
    ArrayTr<T_real> y_axis;
    int requested_cols;
    int requested_rows;
    std::vector<int> detectors;
    float theta;
    
};

TEMPLATE_STRUCT_DLL_EXPORT Scan_Meta_Info<float>;
TEMPLATE_STRUCT_DLL_EXPORT Scan_Meta_Info<double>;

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
/**
 * @brief The Params_Override struct
 */
template<typename T_real>
class DLL_EXPORT Scan_Info
{

public:

    Scan_Info()
    {
        has_netcdf = false;
    }
    
    ~Scan_Info()
    {
        
    }

    const ArrayXXr<T_real>* scaler_values(const std::string& scaler_name) const
    {
        for (const auto& itr : scaler_maps)
        {
            if (itr.name == scaler_name)
            {
                return &(itr.values);
            }
        }
        return nullptr;
    }

    T_real scaler_avg_value(const std::string& scaler_name)
    {
        for (auto& itr : scaler_maps)
        {
            if (itr.name == scaler_name)
            {
                return itr.values.mean();
            }
        }
        return (T_real)0.0;
    }

    Scan_Meta_Info<T_real> meta_info;
    std::vector<Scaler_Map<T_real>> scaler_maps;
    std::vector<Extra_PV> extra_pvs;
    bool has_netcdf; 
};

TEMPLATE_CLASS_DLL_EXPORT Scan_Info<float>;
TEMPLATE_CLASS_DLL_EXPORT Scan_Info<double>;

//-----------------------------------------------------------------------------

} //namespace data_struct

#endif // Scan_Info_H
