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


enum class Electron_Shell {K_SHELL, L_SHELL, M_SHELL, N_SHELL, O_SHELL, P_SHELL, Q_SHELL};

const static map<Electron_Shell, string> Shell_To_String = { {Electron_Shell::K_SHELL, "K"},
                                                            {Electron_Shell::L_SHELL, "L"} ,
                                                            {Electron_Shell::M_SHELL, "M"} ,
                                                            {Electron_Shell::N_SHELL, "N"} ,
                                                            {Electron_Shell::O_SHELL, "O"} ,
                                                            {Electron_Shell::P_SHELL, "P"},
                                                            {Electron_Shell::Q_SHELL, "Q"} };

//-----------------------------------------------------------------------------

///
/// \brief The Quantification_Model class:
///
class DLL_EXPORT Quantification_Model
{

public:
    Quantification_Model();

    ~Quantification_Model();

    void init_element_quant(Element_Quant& out_quant,
                            real_t incident_energy,
                            Element_Info* detector_element,
                            Electron_Shell shell,
                            real_t airpath,
                            real_t detector_chip_thickness,
                            real_t beryllium_window_thickness,
                            real_t germanium_dead_layer,
                            size_t z_number);

    real_t transmission(real_t thickness, real_t beta, real_t llambda) const;

    real_t absorption(real_t thickness, real_t beta, real_t llambda, real_t shell_factor=1) const;

    std::unordered_map<std::string, real_t> model_calibrationcurve(std::unordered_map<std::string, Element_Quant> quant_map, real_t p);

    void model_calibrationcurve(std::vector<Element_Quant>* quant_vec, real_t p);

protected:

};

//-----------------------------------------------------------------------------

DLL_EXPORT Electron_Shell get_shell_by_name(std::string element_name);


} //namespace models

} //namespace quantification

#endif // Quantification_Model_H
