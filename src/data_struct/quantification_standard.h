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



#ifndef Quantification_Standard_H
#define Quantification_Standard_H

#include "core/defines.h"

#include <string>
#include <unordered_map>

#include "data_struct/element_info.h"
#include "data_struct/element_quant.h"
#include "fitting/optimizers/optimizer.h"

namespace data_struct
{

using namespace quantification::models;

const static std::vector<Electron_Shell> Shells_Quant_List({ Electron_Shell::K_SHELL, Electron_Shell::L_SHELL, Electron_Shell::M_SHELL });

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<typename T_real>
struct DLL_EXPORT Quantification_Scaler_Struct
{
    Quantification_Scaler_Struct(unsigned int max_z= CALIBRATION_CURVE_SIZE)
    {
        std::vector<Element_Quant<T_real>> e_quants;
        for (int i = 0; i < max_z; i++)
        {
            e_quants.emplace_back(Element_Quant<T_real>(i + 1));
        }
        for (const auto& itr : Shells_Quant_List)
        {
            curve_quant_map[itr] = e_quants;
        }
    }
    std::unordered_map<Electron_Shell, std::vector<Element_Quant<T_real>> > curve_quant_map;
};

TEMPLATE_STRUCT_DLL_EXPORT Quantification_Scaler_Struct<float>;
TEMPLATE_STRUCT_DLL_EXPORT Quantification_Scaler_Struct<double>;


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<typename T_real>
struct DLL_EXPORT Fitting_Quantification_Struct
{
    Fitting_Quantification_Struct(unsigned int max_z = CALIBRATION_CURVE_SIZE) // 1 (H) - 92 (U)
    {
        std::vector<std::string> quant_scalers = { STR_SR_CURRENT, STR_US_IC, STR_DS_IC };
        init(quant_scalers, max_z);
    }

    Fitting_Quantification_Struct(std::vector<std::string> quantifier_scalers, unsigned int max_z= CALIBRATION_CURVE_SIZE)
    {
        init(quantifier_scalers, max_z);
    }

    void init(std::vector<std::string> quantifier_scalers, unsigned int max_z)
    {
        for (const auto& itr : quantifier_scalers)
        {
            quant_scaler_map[itr] = Quantification_Scaler_Struct<T_real>(max_z);
        }
    }

    void update_weight(Electron_Shell shell, unsigned int Z, T_real weight)
    {
        for (auto& itr : quant_scaler_map)
        {
            itr.second.curve_quant_map.at(shell).at(Z - 1).weight = weight;
        }
    }

    //            Quantifier {SR_Current, US_IC, DS_IC}
    std::unordered_map<std::string, Quantification_Scaler_Struct<T_real>> quant_scaler_map;
};


TEMPLATE_STRUCT_DLL_EXPORT Fitting_Quantification_Struct<float>;
TEMPLATE_STRUCT_DLL_EXPORT Fitting_Quantification_Struct<double>;


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

///
/// \brief The Quantification_Standard class:
///
template<typename T_real>
class DLL_EXPORT Quantification_Standard
{

public:
    Quantification_Standard();

    Quantification_Standard(std::string standard_file, std::vector<std::string> element_names, std::vector<T_real> element_weights);

    Quantification_Standard(std::string standard_file, std::unordered_map<std::string, T_real> e_standard_weights);

    Quantification_Standard(std::string standard_file, std::vector<std::string> element_names, std::vector<T_real> element_weights, bool disable_Ka, bool disable_La);

    ~Quantification_Standard();

    void init_defaults();

    void init_weights_struct(std::string standard_file, std::vector<std::string> element_names, std::vector<T_real> element_weights);

    void normalize_counts_by_time(Fitting_Routines routine);

    //per standard
    std::string standard_filename;
    std::unordered_map<std::string, T_real> element_standard_weights;

    // element name       cts
    std::unordered_map<Fitting_Routines, std::unordered_map<std::string, T_real> > element_counts;

    Spectra<T_real> integrated_spectra;

    T_real sr_current;

    T_real US_IC;

    T_real DS_IC;

    bool disable_Ka_for_quantification;

    bool disable_La_for_quantification;

};


//-----------------------------------------------------------------------------



//-----------------------------------------------------------------------------



} //namespace data_struct

#endif // Quantification_Standard_H
