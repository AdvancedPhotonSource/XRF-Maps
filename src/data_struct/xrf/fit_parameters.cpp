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



#include "fit_parameters.h"
#include <algorithm>
#include <vector>
#include <iostream>

namespace data_struct
{
namespace xrf
{

std::string Fit_Param::bound_type_str()
{
    std::string str_bound_type = " ";
    switch (bound_type)
    {
        case data_struct::xrf::NOT_INIT:
            str_bound_type = "Not Initialized";
            break;
        case data_struct::xrf::FIXED:
            str_bound_type = "Fixed";
            break;
        case data_struct::xrf::LIMITED_LO_HI:
            str_bound_type = "LIMITED LO HI";
            break;
        case data_struct::xrf::LIMITED_LO:
            str_bound_type = "LIMITED LO";
            break;
        case data_struct::xrf::LIMITED_HI:
            str_bound_type = "LIMITED HI";
            break;
        case data_struct::xrf::FIT:
            str_bound_type = "FIT";
            break;
    }
    return str_bound_type;
}

Fit_Parameters::Fit_Parameters(const Fit_Parameters& fit_pars)
{
    _params.clear();
    for(const auto& itr : fit_pars._params)
    {
        _params[itr.first] = itr.second;
    }
}
void Fit_Parameters::add_parameter(std::string name, Fit_Param param)
{
    _params[name] = param;
}

std::vector<real_t> Fit_Parameters::to_array()
{
    std::vector<real_t> arr;
    for(const auto& itr : _params)
    {
        if (itr.second.bound_type > E_Bound_Type::FIXED)
        {
            _params[itr.first].opt_array_index = arr.size();
            arr.push_back(itr.second.value);
        }
    }
    return arr;
}

void Fit_Parameters::moving_average_with(Fit_Parameters fit_params)
{
    for(const auto &itr : _params)
    {
        if(fit_params.contains(itr.first))
        {
            _params[itr.first].value += fit_params[itr.first].value;
            _params[itr.first].value *= 0.5;;
        }
    }
}

void Fit_Parameters::from_array(std::vector<real_t> arr)
{
    from_array(&arr[0], arr.size());
}

void Fit_Parameters::from_array(const real_t* arr, size_t arr_size)
{
    //std::cout<<std::endl;
    for(auto& itr : _params)
    {
        if (itr.second.opt_array_index > -1 && itr.second.opt_array_index < arr_size)
        {
            //std::cout<<" [ "<<_params[itr.first].value<<" : "<<arr[itr.second.opt_array_index]<< " ]  ";
            //std::cout<<" [ "<<_params[itr.first].value - arr[itr.second.opt_array_index]<< " ]  ";
            itr.second.value = arr[itr.second.opt_array_index];
            //_params[itr.first].value = arr[itr.second.opt_array_index];

        }
    }
    //std::cout<<std::endl;
}

void Fit_Parameters::update_value_by_idx(real_t * val, int idx)
{
    std::cout<<"idx "<<idx<<" val "<<*val<<std::endl;
    for(auto& itr : _params)
    {
        if (itr.second.opt_array_index == idx)
        {
            //_params[itr.first].value = *val;
            itr.second.value = *val;
            break;
        }
    }
}

void Fit_Parameters::set_all_value(real_t value, data_struct::xrf::E_Bound_Type btype)
{
    for(auto& itr : _params)
    {
        if (itr.second.bound_type == btype)
            //_params[itr.first].value = value;
            itr.second.value = value;
    }
}

void Fit_Parameters::set_all(data_struct::xrf::E_Bound_Type btype)
{
    for(auto& itr : _params)
    {
        //_params[itr.first].bound_type = btype;
        itr.second.bound_type = btype;
    }
}

void Fit_Parameters::update_values(Fit_Parameters override_fit_params)
{
    for(auto& itr : _params)
    {
        if(override_fit_params.contains(itr.first))
        {
            if( false == std::isnan(override_fit_params[itr.first].value) )
            {
                //_params[itr.first].value = override_fit_params[itr.first].value;
                itr.second.value = override_fit_params[itr.first].value;
            }
            if( false == std::isnan(override_fit_params[itr.first].min_val) )
            {
                //_params[itr.first].min_val = override_fit_params[itr.first].min_val;
                itr.second.min_val = override_fit_params[itr.first].min_val;
            }
            if( false == std::isnan(override_fit_params[itr.first].max_val) )
            {
                //_params[itr.first].max_val = override_fit_params[itr.first].max_val;
                itr.second.max_val = override_fit_params[itr.first].max_val;
            }
            if( override_fit_params[itr.first].bound_type != E_Bound_Type::NOT_INIT)
            {
                //_params[itr.first].bound_type = override_fit_params[itr.first].bound_type;
                itr.second.bound_type = override_fit_params[itr.first].bound_type;
            }
        }
    }
}

void Fit_Parameters::print()
{
    for(const auto& itr : _params)
    {
        std::cout<<" [ "<<itr.first<<" ] = "<<itr.second.value<<"  "<<itr.second.bound_type << std::endl;
    }
    std::cout<<std::endl;

}

void Fit_Parameters::print_non_fixed()
{
    for(const auto& itr : _params)
    {
        if(itr.second.bound_type != E_Bound_Type::FIXED)
        {
            std::cout<<" [ "<<itr.first<<" ] = "<<itr.second.value<<"  "<<itr.second.bound_type << std::endl;
        }
    }
    std::cout<<std::endl;

}

void Fit_Parameters::pow10values()
{
    /* this messes with the fit parameters for next fitting
    for (auto itr : _params)
    {
        if(_params[itr.first].bound_type != E_Bound_Type::FIXED)
        {
            real_t val = _params[itr.first].value;
            _params[itr.first].value = std::pow(10.0, val);
        }
    }
    */
}

} //namespace xrf
} //namespace data_struct
