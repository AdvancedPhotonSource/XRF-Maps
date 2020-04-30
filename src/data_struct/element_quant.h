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



#ifndef Element_Quant_H
#define Element_Quant_H

#include "core/defines.h"

namespace data_struct
{

//-----------------------------------------------------------------------------

struct DLL_EXPORT Element_Quant
{
    Element_Quant()
    {
        zero();
    }
    Element_Quant(real_t weight_)
    {
        zero();
        weight = weight_;
    }
    Element_Quant(const Element_Quant& e)
    {
        weight = e.weight;
        absorption = e.absorption;
        transmission_Be = e.transmission_Be;
        transmission_Ge = e.transmission_Ge;
        yield = e.yield;
        transmission_through_Si_detector = e.transmission_through_Si_detector;
        transmission_through_air = e.transmission_through_air;
        e_cal_ratio = e.e_cal_ratio;
        Z = e.Z;
        name = e.name;
    }
    void zero()
    {
        weight = 0.0;
        absorption = 0.0;
        transmission_Be = 0.0;
        transmission_Ge = 0.0; // or Si dead layer
        yield = 0.0;
        transmission_through_Si_detector = 0.0;
        transmission_through_air = 0.0;// (N2)
        e_cal_ratio = 0.0;
        Z = 0;
        name = "";
    }

    real_t weight;  // in ug/cm2
    real_t absorption;
    real_t transmission_Be;
    real_t transmission_Ge; // or Si dead layer
    real_t yield;
    real_t transmission_through_Si_detector;
    real_t transmission_through_air;// (N2)
    real_t Z;
    real_t e_cal_ratio;
    std::string name;

};

//-----------------------------------------------------------------------------

} //namespace data_struct

#endif // Element_Quant_H
