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



#ifndef Scaler_Lookup_H
#define Scaler_Lookup_H

#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include "core/defines.h"

namespace data_struct
{
    struct Summed_Scaler
    {
        string scaler_name;
        std::vector<string> scalers_to_sum;
    };

//singleton
class DLL_EXPORT Scaler_Lookup
{
public:

    static Scaler_Lookup* inst();

	~Scaler_Lookup();

	void clear();

	void add_beamline_scaler(const string& beamline, const string& scaler_label, const string& scaler_pv, bool is_time_normalized);

    void add_timing_info(const string& beamline, const string& time_pv, double clock);

    void add_summed_scaler(const string& beamline, const string& scaler_label, const vector<string>& scaler_list);

    bool search_for_timing_info(const vector<string>& pv_list, string& out_pv, double& out_clock, string& out_beamline);

    bool search_for_timing_info(const unordered_map<string, double>& pv_map, string& out_pv, double& out_clock, string& out_beamline);

    bool search_pv(const string& pv, string& out_label, bool& out_is_time_normalized, string& out_beamline);

    const vector<struct Summed_Scaler>* get_summed_scaler_list(string beamline) const;

private:

    Scaler_Lookup();

    static Scaler_Lookup *_this_inst;

	struct BeamLine
	{
		//     PV    Label
		map< string, string > scaler_pv_label_map;
		//     PV    Label
		map< string, string > time_normalized_scaler_pv_label_map;
		//    Time_PV  Clock
		map< string, double > timing_info;
		//   beamline      summed scalers
		vector<struct Summed_Scaler> summed_scalers;
	};
	
    //   beamline     Label
    map< string, struct BeamLine > _beamline_map;
 
};
     
} //namespace data_struct

#endif // Scaler_Lookup_H
