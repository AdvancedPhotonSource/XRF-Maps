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

/// Initial Author <2021>: Arthur Glowacki

#include "scaler_lookup.h"

namespace data_struct
{


// ----------------------------------------------------------------------------
// -----------------------------Scaler_Lookup-------------------------------
// ----------------------------------------------------------------------------

Scaler_Lookup* Scaler_Lookup::_this_inst(0);

// ----------------------------------------------------------------------------

Scaler_Lookup* Scaler_Lookup::inst()
{
    if (_this_inst == nullptr)
    {
        _this_inst = new Scaler_Lookup();
    }
    return _this_inst;
}

// ----------------------------------------------------------------------------

Scaler_Lookup::Scaler_Lookup()
{

}

// ----------------------------------------------------------------------------

Scaler_Lookup::~Scaler_Lookup()
{
    clear();
}

// ----------------------------------------------------------------------------

void Scaler_Lookup::clear()
{
	for (auto &itr : _beamline_map)
	{
		itr.second.time_normalized_scaler_pv_label_map.clear();
		itr.second.scaler_pv_label_map.clear();
		itr.second.summed_scalers.clear();
		itr.second.timing_info.clear();
	}
}

// ----------------------------------------------------------------------------

void Scaler_Lookup::add_beamline_scaler(const string& beamline, const string& scaler_label, const string& scaler_pv, bool is_time_normalized)
{
    if (is_time_normalized)
    {
		_beamline_map[beamline].time_normalized_scaler_pv_label_map[scaler_pv] = scaler_label;
    }
    else
    {
		_beamline_map[beamline].scaler_pv_label_map[scaler_pv] = scaler_label;
    }
}

// ----------------------------------------------------------------------------

void Scaler_Lookup::add_timing_info(const string& beamline, const string& time_pv, double clock)
{
	_beamline_map[beamline].timing_info[time_pv] = clock;
}

// ----------------------------------------------------------------------------

void Scaler_Lookup::add_summed_scaler(const string& beamline, const string& scaler_label, const vector<string>& scaler_list)
{
    data_struct::Summed_Scaler s_scal;
    s_scal.scaler_name = scaler_label;
    for (const auto& itr : scaler_list)
    {
        s_scal.scalers_to_sum.push_back(itr);
    }
	_beamline_map[beamline].summed_scalers.push_back(s_scal);
}

// ----------------------------------------------------------------------------

bool Scaler_Lookup::search_for_timing_info(const vector<string>& pv_list, string& out_pv, double& out_clock, string& out_beamline)
{
	
	for (const auto& itr : pv_list)
	{
		for (const auto& bItr : _beamline_map)
		{
			if (bItr.second.timing_info.count(itr) > 0)
			{
				out_pv = itr;
				out_clock = bItr.second.timing_info.at(itr);
				out_beamline = bItr.first;
				return true;
			}
		}
	}
    return false;
}

// ----------------------------------------------------------------------------

bool Scaler_Lookup::search_for_timing_info(const unordered_map<string, double>& pv_map, string& out_pv, double& out_clock, string& out_beamline)
{
    for (const auto& itr : pv_map)
    {
		for (const auto& bItr : _beamline_map)
		{
			if (bItr.second.timing_info.count(itr.first) > 0)
			{
				out_pv = itr.first;
				out_clock = bItr.second.timing_info.at(itr.first);
				out_beamline = bItr.first;
				return true;
			}
		}
    }
    return false;
}

// ----------------------------------------------------------------------------

bool Scaler_Lookup::search_pv(const string& pv, string& out_label, bool& out_is_time_normalized, string& out_beamline)
{
	for (const auto& bItr : _beamline_map)
	{
		if (bItr.second.time_normalized_scaler_pv_label_map.count(pv) > 0)
		{
			out_is_time_normalized = true;
			out_label = bItr.second.time_normalized_scaler_pv_label_map.at(pv);
			if (bItr.first != STR_GENERAL_BEAMLINE)
			{
				out_beamline = bItr.first;
			}
			return true;
		}
		if (bItr.second.scaler_pv_label_map.count(pv) > 0)
		{
			out_is_time_normalized = false;
			out_label = bItr.second.scaler_pv_label_map.at(pv);
			if (bItr.first != STR_GENERAL_BEAMLINE)
			{
				out_beamline = bItr.first;
			}
			return true;
		}
	}
    return false;
}

// ----------------------------------------------------------------------------

const vector<struct Summed_Scaler>* Scaler_Lookup::get_summed_scaler_list(string beamline) const
{
	if(_beamline_map.count(beamline) > 0)
	{
        return &(_beamline_map.at(beamline).summed_scalers);
    }
    return nullptr;
}

// ----------------------------------------------------------------------------

} //namespace data_struct
