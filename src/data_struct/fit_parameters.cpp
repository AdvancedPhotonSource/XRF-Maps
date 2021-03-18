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

const std::string Fit_Param::bound_type_str() const
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

Fit_Parameters::Fit_Parameters(const Fit_Parameters& fit_pars)
{
    _params.clear();
    for(const auto& itr : fit_pars._params)
    {
        _params[itr.first] = itr.second;
    }
}

void Fit_Parameters::add_parameter(Fit_Param param)
{
    _params[param.name] = param;
}

void Fit_Parameters::append_and_update(Fit_Parameters* fit_params)
{
	for (auto& itr : *(fit_params->Params()))
	{
		_params[itr.first] = itr.second;	
	}
}

std::vector<float> Fit_Parameters::to_array_f()
{
    std::vector<float> arr;
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

std::vector<double> Fit_Parameters::to_array_d()
{
    std::vector<double> arr;
    for (const auto& itr : _params)
    {
        if (itr.second.bound_type > E_Bound_Type::FIXED)
        {
            _params[itr.first].opt_array_index = arr.size();
            arr.push_back(itr.second.value);
        }
    }
    return arr;
}


std::vector<std::string> Fit_Parameters::names_to_array()
{
    std::vector<std::string> arr;
    for(const auto& itr : _params)
    {
        _params[itr.first].opt_array_index = arr.size();
        arr.push_back(itr.first);
    }
    return arr;
}

void Fit_Parameters::sum_values(Fit_Parameters fit_params)
{
    for(const auto &itr : _params)
    {
        if(fit_params.contains(itr.first) && itr.second.bound_type > E_Bound_Type::FIXED)
        {
            _params[itr.first].value += fit_params[itr.first].value;
        }
    }
}

void Fit_Parameters::divide_fit_values_by(double divisor)
{
    for(const auto &itr : _params)
    {
        if (itr.second.bound_type > E_Bound_Type::FIXED)
        {
            _params[itr.first].value /= divisor;
        }
    }

}

void Fit_Parameters::from_array(std::vector<float> &arr)
{
    from_array(&arr[0], arr.size());
}

void Fit_Parameters::from_array(std::vector<double>& arr)
{
    from_array(&arr[0], arr.size());
}

void Fit_Parameters::from_array(const float* arr, size_t arr_size)
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

void Fit_Parameters::from_array(const double* arr, size_t arr_size)
{
    //logit_s<<"\n";
    for (auto& itr : _params)
    {
        if (itr.second.opt_array_index > -1 && itr.second.opt_array_index < (int)arr_size)
        {
            itr.second.value = arr[itr.second.opt_array_index];
        }
    }
    //logit_s<<"\n";
}

void Fit_Parameters::set_all_value(double value, E_Bound_Type btype)
{
    for(auto& itr : _params)
    {
        if (itr.second.bound_type == btype)
            //_params[itr.first].value = value;
            itr.second.value = value;
    }
}

void Fit_Parameters::set_all(E_Bound_Type btype)
{
    for(auto& itr : _params)
    {
        //_params[itr.first].bound_type = btype;
        itr.second.bound_type = btype;
    }
}

void Fit_Parameters::update_values(Fit_Parameters *override_fit_params)
{
    for(auto& itr : _params)
    {
        if(override_fit_params->contains(itr.first))
        {
            if( std::isfinite(override_fit_params->at(itr.first).value) )
            {
                itr.second.value = override_fit_params->at(itr.first).value;
            }
            if( std::isfinite(override_fit_params->at(itr.first).min_val) )
            {
                itr.second.min_val = override_fit_params->at(itr.first).min_val;
            }
            if( std::isfinite(override_fit_params->at(itr.first).max_val) )
            {
                itr.second.max_val = override_fit_params->at(itr.first).max_val;
            }
            if( override_fit_params->at(itr.first).bound_type != E_Bound_Type::NOT_INIT)
            {
                itr.second.bound_type = override_fit_params->at(itr.first).bound_type;
            }
        }
    }
}

void Fit_Parameters::update_and_add_values(Fit_Parameters *override_fit_params)
{
    for(auto& itr : *override_fit_params)
    {
        _params[itr.first] = itr.second;
    }
}

void Fit_Parameters::update_and_add_values_gt_zero(Fit_Parameters *override_fit_params)
{
    for(auto& itr : *override_fit_params)
    {
        if(itr.second.value > 0.0)
        {
            _params[itr.first] = itr.second;
        }
    }
}

void Fit_Parameters::print()
{
    for(const auto& itr : _params)
    {
        logit_s<<" [ "<<itr.first<<" ] = "<<itr.second.value<<"  "<<itr.second.bound_type_str() << "\n";
    }
    logit_s<<"\n";

}

void Fit_Parameters::print_non_fixed()
{
    for(const auto& itr : _params)
    {
        if(itr.second.bound_type != E_Bound_Type::FIXED)
        {
            logit_s<<" [ "<<itr.first<<" ] = "<<itr.second.value<<"  "<<itr.second.bound_type_str() << "\n";
        }
    }
    logit_s<<"\n";

}

void Fit_Parameters::remove(Fit_Parameters* override_fit_params)
{
    for (auto& itr : *override_fit_params)
    {
        if (_params.count(itr.first) > 0)
        {
            _params.erase(itr.first);
        }
    }
}

void Fit_Parameters::remove(std::string key)
{
    if (_params.count(key) > 0)
    {
        _params.erase(key);
    }

}

Range get_energy_range(size_t spectra_size, Fit_Parameters* params)
{
	return get_energy_range(params->value(STR_MIN_ENERGY_TO_FIT),
							params->value(STR_MAX_ENERGY_TO_FIT),
							spectra_size,
							params->value(STR_ENERGY_OFFSET),
							params->value(STR_ENERGY_SLOPE));
}


Range get_energy_range(real_t min_energy, real_t max_energy, size_t spectra_size, real_t energy_offset, real_t energy_slope)
{

	struct Range energy_range;
	energy_range.min = static_cast<size_t>(round((min_energy - energy_offset) / energy_slope));
	energy_range.max = static_cast<size_t>(round((max_energy - energy_offset) / energy_slope));
	//if (xmax > used_chan - 1) or (xmax <= np.amin([xmin, used_chan / 20.])):
	if ((energy_range.max > spectra_size - 1) || (energy_range.max <= energy_range.min))
	{
		energy_range.max = spectra_size - 1;
	}
    if (energy_range.min > energy_range.max)
	{
		energy_range.min = 0;
	}
	return energy_range;

}

} //namespace data_struct
