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



#ifndef Calibration_Standard_H
#define Calibration_Standard_H

#include "defines.h"
#include <valarray>
#include <vector>
#include <string>
#include <unordered_map>
#include "spectra.h"

namespace data_struct
{
namespace xrf
{

//-----------------------------------------------------------------------------

///
/// \brief The Calibration_Standard class: Class of resulting fit per element
///
class DLL_EXPORT Calibration_Standard
{

public:
    Calibration_Standard();

    ~Calibration_Standard();

    void append_element_weight(std::string name, real_t weight) {_element_weights[name] = weight;}

    void offset(real_t val) { _e_offset = val; }

    const real_t& offset() const { return _e_offset; }

    void slope(real_t val) { _e_slope = val; }

    const real_t& slope() const { return _e_slope; }

    void quad(real_t val) { _e_quad = val; }

    const real_t& quad() const { return _e_quad; }

    void live_time(real_t val) { _live_time = val; }

    const real_t& live_time() const { return _live_time; }

    void real_time(real_t val) { _real_time = val; }

    const real_t& real_time() const { return _real_time; }

    void current(real_t val) { _current = val; }

    const real_t& current() const { return _current; }

    data_struct::xrf::Spectra* spectra() { return &_spectra; }

    std::unordered_map<std::string, real_t>* element_weights() { return &_element_weights; }

    real_t element_weight(std::string element_symb) { return _element_weights[element_symb]; }

protected:

    real_t _e_offset;

    real_t _e_quad;

    real_t _e_slope;

    //date,
    real_t _live_time;
    real_t _real_time;
    real_t _current;
    real_t _IC_US;
    real_t _IC_DS;

    std::unordered_map<std::string, real_t> _element_weights; // in ug/cm2

    std::valarray<real_t> _us_amp;
    std::valarray<real_t> _ds_amp;

    data_struct::xrf::Spectra _spectra;

};

//-----------------------------------------------------------------------------


} //namespace xrf

} //namespace data_struct

#endif // Calibration_H
