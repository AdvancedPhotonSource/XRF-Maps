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



#include "detector.h"

namespace data_struct
{
namespace xrf
{


Range get_energy_range(real_t min_energy, real_t max_energy, size_t spectra_size, real_t energy_offset, real_t energy_slope)
{

    real_t MIN_ENERGY_TO_FIT = min_energy;
    real_t MAX_ENERGY_TO_FIT = max_energy;


    struct Range energy_range;
    energy_range.min = (int)ceil( (MIN_ENERGY_TO_FIT - energy_offset) / energy_slope );
    energy_range.max = (int)ceil( (MAX_ENERGY_TO_FIT - energy_offset) / energy_slope );
    //if (xmax > used_chan - 1) or (xmax <= np.amin([xmin, used_chan / 20.])):
    if ( (energy_range.max > spectra_size - 1) || (energy_range.max <= energy_range.min) )
    {
        energy_range.max = spectra_size - 1;
    }
    if (energy_range.min < 0 || energy_range.min > energy_range.max)
    {
        energy_range.min = 0;
    }
    return energy_range;

}

void gen_energy_vector(real_t number_channels, real_t energy_offset, real_t energy_slope, std::vector<real_t> *out_vec)
{

    out_vec->resize(number_channels);
    for(int i=0; i<number_channels; i++)
    {
        (*out_vec)[i] = (i * energy_slope) + energy_offset;
    }

}

//===============================================================================================

Detector::Detector()
{

    _calib_energy_offset = 0.0;

    _calib_energy_slope = 1.0;

    _calib_energy_quad = 0.0;

    _chip_thickness = 0.0;

    _detector_element = nullptr;

}

Detector::~Detector()
{

}


} //namespace xrf
} //namespace data_struct
