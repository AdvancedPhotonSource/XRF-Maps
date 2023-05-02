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



#include "optimizer.h"


namespace fitting
{
namespace optimizers
{

    // ----------------------------------------------------------------------------
    
    std::string optimizer_outcome_to_str(OPTIMIZER_OUTCOME outcome)
    {
        switch (outcome)
        {
        case fitting::optimizers::OPTIMIZER_OUTCOME::CONVERGED:
            return "Converged";
        case fitting::optimizers::OPTIMIZER_OUTCOME::F_TOL_LT_TOL:
            return "Ftol less than TOL";
        case fitting::optimizers::OPTIMIZER_OUTCOME::X_TOL_LT_TOL:
            return "XTol less than TOL";
        case fitting::optimizers::OPTIMIZER_OUTCOME::G_TOL_LT_TOL:
            return "GTol less than TOL";
        case fitting::optimizers::OPTIMIZER_OUTCOME::EXHAUSTED:
            return "EXHAUSTED";
        case fitting::optimizers::OPTIMIZER_OUTCOME::CRASHED:
            return "CRASHED";
        case fitting::optimizers::OPTIMIZER_OUTCOME::EXPLODED:
            return "Exploded";
        case fitting::optimizers::OPTIMIZER_OUTCOME::FAILED:
            return "FAILED";
        case fitting::optimizers::OPTIMIZER_OUTCOME::FOUND_NAN:
            return "Found NAN";
        case fitting::optimizers::OPTIMIZER_OUTCOME::FOUND_ZERO:
            return "Found Zero";
        case fitting::optimizers::OPTIMIZER_OUTCOME::STOPPED:
            return "Stopped";
        case fitting::optimizers::OPTIMIZER_OUTCOME::TRAPPED:
            return "Trapped";
        }
        return "";
    }

    // ----------------------------------------------------------------------------


TEMPLATE_STRUCT_DLL_EXPORT User_Data<float>;
TEMPLATE_STRUCT_DLL_EXPORT User_Data<double>;
TEMPLATE_STRUCT_DLL_EXPORT Gen_User_Data<float>;
TEMPLATE_STRUCT_DLL_EXPORT Gen_User_Data<double>;
TEMPLATE_STRUCT_DLL_EXPORT Quant_User_Data<float>;
TEMPLATE_STRUCT_DLL_EXPORT Quant_User_Data<double>;
TEMPLATE_CLASS_DLL_EXPORT Optimizer<float>;
TEMPLATE_CLASS_DLL_EXPORT Optimizer<double>;

} //namespace optimizers
} //namespace fitting
