/***

Copyright (c) 2016 Arthur Glowacki

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

***/

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
