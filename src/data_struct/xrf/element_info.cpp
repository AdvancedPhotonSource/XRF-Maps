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



#include "element_info.h"

namespace data_struct
{
namespace xrf
{

Element_Info::Element_Info()
{

    number = 0;
    name = " ";
    density = 1.0;
    mass = 1.0;

    xrf.emplace(std::make_pair("Ka1", 0.0));
    xrf.emplace(std::make_pair("Ka2", 0.0));
    xrf.emplace(std::make_pair("Kb1", 0.0));
    xrf.emplace(std::make_pair("Kb2", 0.0));

    xrf.emplace(std::make_pair("La1", 0.0));
    xrf.emplace(std::make_pair("La2", 0.0));
    xrf.emplace(std::make_pair("Lb1", 0.0));
    xrf.emplace(std::make_pair("Lb2", 0.0));
    xrf.emplace(std::make_pair("Lb3", 0.0));
    xrf.emplace(std::make_pair("Lb4", 0.0));
    xrf.emplace(std::make_pair("Lb5", 0.0));

    xrf.emplace(std::make_pair("Lg1", 0.0));
    xrf.emplace(std::make_pair("Lg2", 0.0));
    xrf.emplace(std::make_pair("Lg3", 0.0));
    xrf.emplace(std::make_pair("Lg4", 0.0));
    xrf.emplace(std::make_pair("Ll", 0.0));
    xrf.emplace(std::make_pair("Ln", 0.0));

    xrf.emplace(std::make_pair("Ma1", 0.0));
    xrf.emplace(std::make_pair("Ma2", 0.0));
    xrf.emplace(std::make_pair("Mb", 0.0));
    xrf.emplace(std::make_pair("Mg", 0.0));


    xrf_abs_yield.emplace(std::make_pair("Ka1", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Ka2", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Kb1", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Kb2", 0.0));

    xrf_abs_yield.emplace(std::make_pair("La1", 0.0));
    xrf_abs_yield.emplace(std::make_pair("La2", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Lb1", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Lb2", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Lb3", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Lb4", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Lb5", 0.0));

    xrf_abs_yield.emplace(std::make_pair("Lg1", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Lg2", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Lg3", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Lg4", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Ll", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Ln", 0.0));

    xrf_abs_yield.emplace(std::make_pair("Ma1", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Ma2", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Mb", 0.0));
    xrf_abs_yield.emplace(std::make_pair("Mg", 0.0));

    yieldD.emplace(std::make_pair("k", 0.0));
    yieldD.emplace(std::make_pair("l1", 0.0));
    yieldD.emplace(std::make_pair("l2", 0.0));
    yieldD.emplace(std::make_pair("l3", 0.0));
    yieldD.emplace(std::make_pair("m", 0.0));

    bindingE.emplace(std::make_pair("K", 0.0));

    bindingE.emplace(std::make_pair("L1", 0.0));
    bindingE.emplace(std::make_pair("L2", 0.0));
    bindingE.emplace(std::make_pair("L3", 0.0));

    bindingE.emplace(std::make_pair("M1", 0.0));
    bindingE.emplace(std::make_pair("M2", 0.0));
    bindingE.emplace(std::make_pair("M3", 0.0));
    bindingE.emplace(std::make_pair("M4", 0.0));
    bindingE.emplace(std::make_pair("M5", 0.0));

    bindingE.emplace(std::make_pair("N1", 0.0));
    bindingE.emplace(std::make_pair("N2", 0.0));
    bindingE.emplace(std::make_pair("N3", 0.0));
    bindingE.emplace(std::make_pair("N4", 0.0));
    bindingE.emplace(std::make_pair("N5", 0.0));
    bindingE.emplace(std::make_pair("N6", 0.0));
    bindingE.emplace(std::make_pair("N7", 0.0));

    bindingE.emplace(std::make_pair("O1", 0.0));
    bindingE.emplace(std::make_pair("O2", 0.0));
    bindingE.emplace(std::make_pair("O3", 0.0));
    bindingE.emplace(std::make_pair("O4", 0.0));
    bindingE.emplace(std::make_pair("O5", 0.0));

    bindingE.emplace(std::make_pair("P1", 0.0));
    bindingE.emplace(std::make_pair("P2", 0.0));
    bindingE.emplace(std::make_pair("P3", 0.0));


    jump.emplace(std::make_pair("K", 0.0));

    jump.emplace(std::make_pair("L1", 0.0));
    jump.emplace(std::make_pair("L2", 0.0));
    jump.emplace(std::make_pair("L3", 0.0));

    jump.emplace(std::make_pair("M1", 0.0));
    jump.emplace(std::make_pair("M2", 0.0));
    jump.emplace(std::make_pair("M3", 0.0));
    jump.emplace(std::make_pair("M4", 0.0));
    jump.emplace(std::make_pair("M5", 0.0));

    jump.emplace(std::make_pair("N1", 0.0));
    jump.emplace(std::make_pair("N2", 0.0));
    jump.emplace(std::make_pair("N3", 0.0));
    jump.emplace(std::make_pair("N4", 0.0));
    jump.emplace(std::make_pair("N5", 0.0));

    jump.emplace(std::make_pair("O1", 0.0));
    jump.emplace(std::make_pair("O2", 0.0));
    jump.emplace(std::make_pair("O3", 0.0));

    /*
    f1_atomic_scattering_real = nullptr;
    f2_atomic_scattering_imaginary = nullptr; 
    atomic_scattering_len = -1;
    extra_energies = nullptr;
    extra_f1 = nullptr;
    extra_f2 = nullptr;
    extra_energies_len = -1;
    */
}

Element_Info::~Element_Info()
{
    xrf.clear();
    xrf_abs_yield.clear();
    yieldD.clear();
    bindingE.clear();
    jump.clear();

    f1_atomic_scattering_real.clear();
    f2_atomic_scattering_imaginary.clear();
    extra_energies.clear();
    extra_f1.clear();
    extra_f2.clear();

    /*
    if(f1_atomic_scattering_real != nullptr)
        delete [] f1_atomic_scattering_real;
    if(f2_atomic_scattering_imaginary != nullptr)
         delete [] f2_atomic_scattering_imaginary;
    f1_atomic_scattering_real = nullptr;
    f2_atomic_scattering_imaginary = nullptr;

    if(extra_energies != nullptr)
        delete [] extra_energies;
    if(extra_f1 != nullptr)
        delete [] extra_f1;
    if(extra_f2 != nullptr)
        delete [] extra_f2;
    extra_energies = nullptr;
    extra_f1 = nullptr;
    extra_f2 = nullptr;
    */
}

void Element_Info::init_f_energies(int len)
{
    f1_atomic_scattering_real.resize(len);
    f2_atomic_scattering_imaginary.resize(len);
    /*
    atomic_scattering_len = len;
    if(f1_atomic_scattering_real != nullptr)
        delete [] f1_atomic_scattering_real;
    f1_atomic_scattering_real = new real_t[atomic_scattering_len];

    if(f2_atomic_scattering_imaginary != nullptr)
        delete [] f2_atomic_scattering_imaginary;
    f2_atomic_scattering_imaginary = new real_t[atomic_scattering_len];
    */
}

void Element_Info::init_extra_energies(int len)
{
    extra_energies.resize(len);
    extra_f1.resize(len);
    extra_f2.resize(len);
    /*
    extra_energies_len = len;
    if(extra_energies != nullptr)
        delete [] extra_energies;
    extra_energies = new real_t[extra_energies_len];

    if(extra_f1 != nullptr)
        delete [] extra_f1;
    extra_f1 = new real_t[extra_energies_len];

    if(extra_f2 != nullptr)
        delete [] extra_f2;
    extra_f2 = new real_t[extra_energies_len];
    */
}

void Element_Info::get_energies_between(real_t energy, real_t* out_low, real_t* out_high, size_t* out_low_idx, size_t* out_high_idx)
{
    *out_low_idx = 0;
    *out_high_idx = energies->size()-1;
    for(size_t l = 0; l < energies->size()-1; l++)
    {
        if(energy < (*energies)[l+1])
        {
            *out_low_idx = l;
            *out_low = (*energies)[l];
            break;
        }
    }
    for(size_t h = energies->size()-1; h > 1; h--)
    {
        if(energy > (*energies)[h-1])
        {
            *out_high_idx = h;
            *out_high = (*energies)[h];
            break;
        }
    }
}


// -----------------------------Element_Info_Map-------------------------------

Element_Info_Map* Element_Info_Map::_this_inst(0);

Element_Info_Map* Element_Info_Map::inst()
{
    if (_this_inst == nullptr)
    {
        _this_inst = new Element_Info_Map();
    }
    return _this_inst;
}

Element_Info_Map::Element_Info_Map()
{

}

Element_Info_Map::~Element_Info_Map()
{
    clear();
}

void Element_Info_Map::clear()
{
    /*
    for(auto itr : _number_element_info_map)
    {
        Element_Info* element = _number_element_info_map[itr.first];
        if (element != nullptr)
            delete element;
        element = nullptr;
        _number_element_info_map[itr.first] = nullptr;
    }
    */
    _number_element_info_map.clear();
    _name_element_info_map.clear();
}

void Element_Info_Map::add_element(Element_Info* element)
{
    // set global energies pointer
    element->energies = &_energies;
    // TODO: check if it exists
    _number_element_info_map.insert(std::pair<int, Element_Info*>(element->number, element));
    _name_element_info_map.insert(std::pair<std::string, Element_Info*>(element->name, element));
}

Element_Info* Element_Info_Map::get_element(int element_number)
{
    if (_number_element_info_map.count(element_number) > 0)
        return _number_element_info_map[element_number];
    return nullptr;
}

Element_Info* Element_Info_Map::get_element(std::string element_name)
{
    if (_name_element_info_map.count(element_name) > 0)
        return _name_element_info_map[element_name];
    return nullptr;
}

void Element_Info_Map::generate_default_elements(int start_element, int end_element)
{
    for (int i = start_element; i <end_element; i++)
    {
        Element_Info* element = new Element_Info();
        element->number = i;
        element->name = Element_Symbols[i];
       _number_element_info_map.insert(std::pair<int, Element_Info*>(element->number, element));
       _name_element_info_map.insert(std::pair<std::string, Element_Info*>(element->name, element));
    }
}
/*
void Element_Info_Map::set_energies(float* energy_arr, int num_energies)
{

    _energies.resize(num_energies);
    for(int i=0; i<num_energies; i++)
    {
        _energies[i] = energy_arr[i];
    }

}
*/
} //namespace data_struct
} //namespace xrf
