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


#ifndef Fit_Element_Map_H
#define Fit_Element_Map_H

#include "defines.h"
#include "fit_parameters.h"
#include "element_info.h"

namespace data_struct
{
namespace xrf
{

const real_t ENERGY_RES_OFFSET = 150.0;
const real_t ENERGY_RES_SQRT = 12.0;

const real_t AVOGADRO = (real_t)6.02204531e23;
const real_t HC_ANGSTROMS = (real_t)12398.52;
const real_t RE = (real_t)2.817938070e-13;		// in cm

//-----------------------------------------------------------------------------

enum Element_Param_Type{ None = 0, Ka_Line = 1, Kb_Line = 2, L_Line = 3, M_Line = 7 };

//-----------------------------------------------------------------------------

struct Element_Energy_Ratio
{

    Element_Energy_Ratio(real_t e, real_t r, real_t m, Element_Param_Type et)
    {
        energy = e;
        ratio = r;
        mu_fraction = m;
        ptype = et;
    }

    real_t energy;
    real_t ratio;
    real_t mu_fraction;
    Element_Param_Type ptype;
};

//-----------------------------------------------------------------------------

///
/// \brief The Fit_Element class: Class that hold element information and the results of the fitting to a spectra.
///                                Able to store 2d image of the fit an element
///
class DLL_EXPORT Fit_Element_Map
{

public:
    Fit_Element_Map(std::string name, Element_Info* element_info);

    ~Fit_Element_Map();

    const real_t center() const { return _center; }

    const real_t width() const { return _width; }

    void set_custom_multiply_ratio(unsigned int idx, real_t multi);

    void init_energy_ratio_for_detector_element(Element_Info* detector_element);

    const std::string full_name() const { return _full_name; }

    const std::vector<Element_Energy_Ratio>& energy_ratios() const { return _energy_ratios; }

    const real_t width_multi() const { return _width_multi; }

protected:

    void generate_energy_ratio(real_t energy, real_t ratio, Element_Param_Type et, Element_Info* detector_element);

    // reference to element information from Database
    Element_Info* _element_info;

    std::string _full_name;
    std::vector<Element_Energy_Ratio> _energy_ratios;
    std::vector<real_t> _energy_ratio_custom_multipliers;

    real_t _center;
    real_t _width;
    real_t _width_multi;

};

//-----------------------------------------------------------------------------

typedef std::unordered_map<std::string, Fit_Element_Map*> Fit_Element_Map_Dict;


} //namespace xrf

} //namespace data_struct

#endif // Fit_Element_Map_H
