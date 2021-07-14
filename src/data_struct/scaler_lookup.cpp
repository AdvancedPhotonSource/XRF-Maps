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
    _time_normalized_scaler_pv_label_map.clear();
    _scaler_pv_label_map.clear();
    _summed_scalers.clear();
}

// ----------------------------------------------------------------------------

void Scaler_Lookup::add_beamline_scaler(const string& beamline, const string& scaler_label, const string& scaler_pv, bool is_time_normalized)
{
    if (is_time_normalized)
    {
        _time_normalized_scaler_pv_label_map[scaler_pv] = scaler_label;
    }
    else
    {
        _scaler_pv_label_map[scaler_pv] = scaler_label;
    }
}

// ----------------------------------------------------------------------------

void Scaler_Lookup::add_timing_info(const string& time_pv, double clock)
{
    _timing_info[time_pv] = clock;
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
    _summed_scalers[beamline].push_back(s_scal);
}

// ----------------------------------------------------------------------------

bool Scaler_Lookup::search_for_timing_info(const vector<string>& pv_list, string& out_pv, double& out_clock)
{
    for (const auto& itr : pv_list)
    {
        if (_timing_info.count(itr) > 0)
        {
            out_pv = itr;
            out_clock = _timing_info.at(itr);
            return true;
        }
    }
    return false;
}

// ----------------------------------------------------------------------------

bool Scaler_Lookup::search_for_timing_info(const unordered_map<string, real_t>& pv_map, string& out_pv, double& out_clock)
{
    for (const auto& itr : pv_map)
    {
        if (_timing_info.count(itr.first) > 0)
        {
            out_pv = itr.first;
            out_clock = _timing_info.at(itr.first);
            return true;
        }
    }
    return false;
}

// ----------------------------------------------------------------------------

bool Scaler_Lookup::search_pv(const string& pv, string& out_label, bool& out_is_time_normalized)
{
    if (_time_normalized_scaler_pv_label_map.count(pv) > 0)
    {
        out_is_time_normalized = true;
        out_label = _time_normalized_scaler_pv_label_map.at(pv);
        return true;
    }
    if (_scaler_pv_label_map.count(pv) > 0)
    {
        out_is_time_normalized = false;
        out_label = _scaler_pv_label_map.at(pv);
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------

const vector<struct Summed_Scaler>* Scaler_Lookup::get_summed_scaler_list(string beamline) const
{
    if (_summed_scalers.count(beamline) > 0)
    {
        return &_summed_scalers.at(beamline);
    }
    return nullptr;
}

// ----------------------------------------------------------------------------

} //namespace data_struct
