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



#include "quantification_standard.h"

#include <string>
#include <cmath>

#include "quantification/models/quantification_model.h"

namespace data_struct
{

//-----------------------------------------------------------------------------

Quantification_Standard::Quantification_Standard()
{
    init_defaults();
}

//-----------------------------------------------------------------------------

Quantification_Standard::Quantification_Standard(std::string standard_file, std::vector<std::string> element_names, std::vector<real_t> element_weights)
{
    init_defaults();
    init_weights_struct(standard_file, element_names, element_weights);
}

//-----------------------------------------------------------------------------

Quantification_Standard::Quantification_Standard(std::string standard_file, std::unordered_map<std::string, real_t> e_standard_weights)
{
    init_defaults();
    this->standard_filename = standard_file;
    for (const auto& itr : e_standard_weights)
    {
        element_standard_weights[itr.first] = itr.second;
    }
}

//-----------------------------------------------------------------------------

Quantification_Standard::~Quantification_Standard()
{

}

//-----------------------------------------------------------------------------

void Quantification_Standard::init_defaults()
{

    sr_current = 1.0;
    US_IC = 1.0;
    DS_IC = 1.0;
}

//-----------------------------------------------------------------------------

void Quantification_Standard::init_weights_struct(std::string standard_file, std::vector<std::string> element_names, std::vector<real_t> element_weights)
{
    standard_filename = standard_file;
    for (size_t i = 0; i < element_names.size(); i++)
    {
        element_standard_weights[element_names[i]] = element_weights[i];
    }
}

//-----------------------------------------------------------------------------

void Quantification_Standard::normalize_counts_by_time(Fitting_Routines routine)
{
   
    if (element_counts.count(routine) > 0)
    {
        for (auto& itr : element_counts.at(routine))
        {
            itr.second /= integrated_spectra.elapsed_livetime();
        }
    }
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

} //namespace data_struct
