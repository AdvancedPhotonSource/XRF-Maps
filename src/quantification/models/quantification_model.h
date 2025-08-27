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



#ifndef Quantification_Model_H
#define Quantification_Model_H

#include "core/defines.h"

#include <string>
#include <unordered_map>

#include "data_struct/element_info.h"
#include "data_struct/element_quant.h"

namespace quantification
{
namespace models
{

using namespace data_struct;




//-----------------------------------------------------------------------------

///
/// \brief The Quantification_Model class:
///
template<typename T_real>
class DLL_EXPORT Quantification_Model
{

public:
    Quantification_Model();

    ~Quantification_Model();

    void init_element_quant(Element_Quant<T_real>& out_quant,
                            T_real incident_energy,
                            Element_Info<T_real>* detector_element,
                            data_struct::Electron_Shell shell,
                            T_real airpath,
                            T_real detector_chip_thickness,
                            T_real beryllium_window_thickness,
                            T_real germanium_dead_layer,
                            int z_number);

    T_real transmission(T_real thickness, T_real beta, T_real llambda) const;

    T_real absorption(T_real thickness, T_real beta, T_real llambda, T_real shell_factor=1) const;

    std::unordered_map<std::string, T_real> model_calibrationcurve(std::unordered_map<std::string, Element_Quant<T_real>> quant_map, T_real p);

    void model_calibrationcurve(std::vector<Element_Quant<T_real>>* quant_vec, T_real p);

protected:

};

//-----------------------------------------------------------------------------



} //namespace models

} //namespace quantification

#endif // Quantification_Model_H
