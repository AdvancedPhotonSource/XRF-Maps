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

//-----------------------------------------------------------------------------

template<typename T_real>
const std::string Fit_Param<T_real>::bound_type_str() const
{
    switch (bound_type)
    {
		case E_Bound_Type::NOT_INIT:
            return "Not Initialized";
            break;
        case E_Bound_Type::FIXED:
			return "Fixed";
            break;
        case E_Bound_Type::LIMITED_LO_HI:
			return "LIMITED LO HI";
            break;
        case E_Bound_Type::LIMITED_LO:
			return "LIMITED LO";
            break;
        case E_Bound_Type::LIMITED_HI:
			return "LIMITED HI";
            break;
        case E_Bound_Type::FIT:
			return "FIT";
            break;
    }
    return "N/A";
}

//-----------------------------------------------------------------------------

template<typename T_real>
Fit_Parameters<T_real>::Fit_Parameters(const Fit_Parameters<T_real>& fit_pars)
{
    _params.clear();
    for(const auto& itr : fit_pars._params)
    {
        _params[itr.first] = itr.second;
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::add_parameter(Fit_Param<T_real> param)
{
    _params[param.name] = param;
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::append_and_update(const Fit_Parameters& fit_params)
{
	for (typename std::unordered_map<std::string, Fit_Param<T_real>>::const_iterator itr = fit_params.begin(); itr != fit_params.end(); itr++ )
	{
		_params[itr->first] = itr->second;	
	}
}

//-----------------------------------------------------------------------------

template<typename T_real>
std::vector<T_real> Fit_Parameters<T_real>::to_array()
{
    std::vector<T_real> arr;
    for(const auto& itr : _params)
    {
        if (itr.second.bound_type != E_Bound_Type::FIXED)
        {
            _params[itr.first].opt_array_index = arr.size();
            arr.push_back(itr.second.value);
        }
    }
    return arr;
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::to_array_with_bounds(std::vector<double>& fitp, std::vector<double>&lb, std::vector<double>& ub, std::vector<double>& step )
{
    for(const auto& itr : _params)
    {
        if (itr.second.bound_type != E_Bound_Type::FIXED)
        {
            _params[itr.first].opt_array_index = fitp.size();
            fitp.push_back(static_cast<double>(itr.second.value));
            lb.push_back(static_cast<double>(itr.second.min_val));
            ub.push_back(static_cast<double>(itr.second.max_val));
            step.push_back(static_cast<double>(itr.second.step_size));
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
std::vector<std::string> Fit_Parameters<T_real>::names_to_array()
{
    std::vector<std::string> arr;
    for(const auto& itr : _params)
    {
        _params[itr.first].opt_array_index = arr.size();
        arr.push_back(itr.first);
    }
    return arr;
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::sum_values(Fit_Parameters<T_real>  fit_params)
{
    for(const auto &itr : _params)
    {
        if(fit_params.contains(itr.first) && itr.second.bound_type != E_Bound_Type::FIXED)
        {
            _params[itr.first].value += fit_params[itr.first].value;
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::divide_fit_values_by(T_real divisor)
{
    for(const auto &itr : _params)
    {
        if (itr.second.bound_type != E_Bound_Type::FIXED)
        {
            _params[itr.first].value /= divisor;
        }
    }

}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::from_array(std::vector<T_real> &arr)
{
    from_array(&arr[0], arr.size());
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::from_array_d(const std::vector<double> &arr)
{
    for(auto& itr : _params)
    {
        if (itr.second.opt_array_index > -1 && itr.second.opt_array_index < (int)arr.size())
        {
            itr.second.value = arr[itr.second.opt_array_index];
        }
    }
}


//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::from_array(const T_real* arr, size_t arr_size)
{
    //logit_s<<"\n";
    for(auto& itr : _params)
    {
        if (itr.second.opt_array_index > -1 && itr.second.opt_array_index < (int)arr_size)
        {
            itr.second.value = arr[itr.second.opt_array_index];
        }
    }
    //logit_s<<"\n";
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::set_all_value(T_real value, E_Bound_Type btype)
{
    for(auto& itr : _params)
    {
        if (itr.second.bound_type == btype)
            //_params[itr.first].value = value;
            itr.second.value = value;
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::set_all(E_Bound_Type btype)
{
    for(auto& itr : _params)
    {
        //_params[itr.first].bound_type = btype;
        itr.second.bound_type = btype;
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::set_all_except(E_Bound_Type btype,const std::vector<std::string>& exception_list)
{
    for (auto& itr : _params)
    {
        bool found = false;
        for (const auto& eitr : exception_list)
        {
            if (itr.first == eitr)
            {
                found = true;
                break;
            }
        }
        if (false == found)
        {
            itr.second.bound_type = btype;
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::update_values(const Fit_Parameters<T_real>  * const override_fit_params)
{
    for(auto& itr : _params)
    {
        if(override_fit_params->contains(itr.first))
        {
            if( std::isfinite(override_fit_params->at(itr.first).value) )
            {
                itr.second.value = override_fit_params->at(itr.first).value;
            }
            else
            {
                logW << itr.first << "new value = " << override_fit_params->at(itr.first).value << " Reverting to original value "<< itr.second.value << " .\n";
            }
            if( std::isfinite(override_fit_params->at(itr.first).min_val) )
            {
                itr.second.min_val = override_fit_params->at(itr.first).min_val;
            }
            else
            {
                //logW << itr.first << "new min_val = " << override_fit_params->at(itr.first).min_val << " Reverting to original min_val "<< itr.second.min_val << " . Not updating min_val\n";
            }
            if( std::isfinite(override_fit_params->at(itr.first).max_val) )
            {
                itr.second.max_val = override_fit_params->at(itr.first).max_val;
            }
            else
            {
                //logW << itr.first << " max_val = " << itr.second.max_val << " . Not updating max_val\n";
            }
            if( override_fit_params->at(itr.first).bound_type != E_Bound_Type::NOT_INIT)
            {
                itr.second.bound_type = override_fit_params->at(itr.first).bound_type;
            }
            else
            {
                //logW << itr.first << " bound_type == E_Bound_Type::NOT_INIT. Not updating bound_type\n";
            }
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::update_and_add_values(Fit_Parameters<T_real>  *override_fit_params)
{
    for(auto& itr : *override_fit_params)
    {
        _params[itr.first] = itr.second;
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::update_and_add_values_gt_zero(Fit_Parameters<T_real>  *override_fit_params)
{
    for(auto& itr : *override_fit_params)
    {
        if(itr.second.value > 0.0)
        {
            _params[itr.first] = itr.second;
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::update_follow_constraints(Fit_Parameters<T_real>  *override_fit_params)
{
    for(auto& itr : *override_fit_params)
    {
        if(itr.second.value != _params.at(itr.first).value)
        {
            if(itr.second.value <= _params.at(itr.first).max_val && itr.second.value >= _params.at(itr.first).min_val)
            {
                _params[itr.first] = itr.second;
            }
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::update_value_to_constraints()
{
    for(auto& itr : _params)
    {
        if(itr.second.bound_type != E_Bound_Type::FIXED)
        {
            if(itr.second.value > itr.second.max_val)
            {
                itr.second.value = itr.second.max_val;
            }
            if(itr.second.value < itr.second.min_val)
            {
                itr.second.value = itr.second.min_val;
            }
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::print()
{
    logit_s << "     Name  \t value  \t min  \t max  \t step size \t fitting\n\n";
    for(const auto& itr : _params)
    {
        if(itr.second.value > itr.second.max_val || itr.second.value < itr.second.min_val)
        {
            logit_s<<"\033[1;31m "<<" "<<itr.first<<" \t "<<itr.second.value<<" \t " << itr.second.min_val << " \t " << itr.second.max_val << " \t " << itr.second.step_size << " \t " <<itr.second.bound_type_str() << "\033[0;m \n";    
        }
        else
        {
            logit_s<<" "<<itr.first<<" \t "<<itr.second.value<<" \t " << itr.second.min_val << " \t " << itr.second.max_val << " \t " << itr.second.step_size << " \t " <<itr.second.bound_type_str() << "\n";
        }
    }
    logit_s<<"\n";

}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::print_non_fixed()
{
    for(const auto& itr : _params)
    {
        if(itr.second.bound_type != E_Bound_Type::FIXED)
        {
            if(itr.second.value > itr.second.max_val || itr.second.value < itr.second.min_val)
            {
                logit_s<<"\033[1;31m [ "<<itr.first<<" ] = "<<itr.second.value<<"  "<<itr.second.bound_type_str() << "\033[0;m \n";    
            }
            else
            {
                logit_s<<" [ "<<itr.first<<" ] = "<<itr.second.value<<"  "<<itr.second.bound_type_str() << "\n";
            }
        }
    }
    logit_s<<"\n";

}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::remove(Fit_Parameters* override_fit_params)
{
    for (auto& itr : *override_fit_params)
    {
        if (_params.count(itr.first) > 0)
        {
            _params.erase(itr.first);
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Fit_Parameters<T_real>::remove(std::string key)
{
    if (_params.count(key) > 0)
    {
        _params.erase(key);
    }

}

//-----------------------------------------------------------------------------

template<typename T_real>
size_t Fit_Parameters<T_real>::size_non_fixed()
{
    size_t sum = 0;
    for (auto& itr : _params)
    {
        bool is_fit = (itr.second.bound_type == E_Bound_Type::LIMITED_LO_HI) || (itr.second.bound_type == E_Bound_Type::LIMITED_LO) || (itr.second.bound_type == E_Bound_Type::LIMITED_HI) || (itr.second.bound_type == E_Bound_Type::FIT);
        if (is_fit)
        {
            sum++;
        }
    }
    return sum;
}

//-----------------------------------------------------------------------------


TEMPLATE_STRUCT_DLL_EXPORT Fit_Param<float>;
TEMPLATE_STRUCT_DLL_EXPORT Fit_Param<double>;

TEMPLATE_CLASS_DLL_EXPORT Fit_Parameters<float>;
TEMPLATE_CLASS_DLL_EXPORT Fit_Parameters<double>;

} //namespace data_struct
