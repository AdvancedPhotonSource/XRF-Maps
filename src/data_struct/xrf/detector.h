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



#ifndef Detector_H
#define Detector_H

#include "defines.h"
#include "element_info.h"
#include <vector>
#include <string>
#include <math.h>

namespace data_struct
{
namespace xrf
{



/**
 * @brief The Range struct to determine size of spectra we want to fit or model
 */
struct Range
{
    Range() {min = 0; max = 0;}
    Range(int rmin, int rmax) {min = rmin; max = rmax;}
    size_t count() const  {return (max - min) + 1;}
    int min;
    int max;
};


struct detector_meta_struct
{
    short       number;
    std::string name;
    std::string description;
    std::string unit;
    real_t      value;
};

//-----------------------------------------------------------------------------

///
/// \brief The Detector class:
///
class DLL_EXPORT Detector
{

public:
    Detector();

    ~Detector();

    void energy_offset(real_t val) { _calib_energy_offset = val; }

    const real_t& energy_offset() const { return _calib_energy_offset; }

    void energy_slope(real_t val) { _calib_energy_slope = val; }

    const real_t& energy_slope() const { return _calib_energy_slope; }

    void energy_quadratic(real_t val) { _calib_energy_quad = val; }

    const real_t& energy_quadratic() const { return _calib_energy_quad; }

    void set_element(Element_Info* detector_element)  { _detector_element = detector_element; }

    const Element_Info * const get_element() const { return _detector_element; }

    std::vector<detector_meta_struct> meta_array;

protected:

    real_t _calib_energy_offset;

    real_t _calib_energy_slope;

    real_t _calib_energy_quad;

    real_t _chip_thickness;

    Element_Info* _detector_element;

};

//-----------------------------------------------------------------------------

/**
 * @brief get_energy_range: genereates a range which consists of min and max. This represents the min energy and max enegry of the spectra to fit.
 * @param min_energy
 * @param max_energy
 * @param spectra_size
 * @param calibration: energy calibration
 * @return Range structure with the min energy and max enegry of the spectra to fit.
 */
DLL_EXPORT Range get_energy_range(real_t min_energy, real_t max_energy, size_t spectra_size, real_t energy_offset, real_t energy_slope);


DLL_EXPORT void gen_energy_vector(real_t number_channels, real_t energy_offset, real_t energy_slope, std::vector<real_t> *out_vec);

} //namespace xrf

} //namespace data_struct

#endif // Detector_H
