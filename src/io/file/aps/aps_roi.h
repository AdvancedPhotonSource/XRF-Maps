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

/// Initial Author <2022>: Arthur Glowacki



#ifndef _APS_ROI_
#define _APS_ROI_

#include "core/defines.h"
#include <string>
#include <vector>
#include <map>
#include<unordered_map>
#include <json/json.h>
#include "data_struct/spectra.h"


namespace io
{
namespace file
{
namespace aps
{

typedef std::pair<int, int> int_point;

template<class _real>
DLL_EXPORT bool load_v10_rois(std::string path, std::map<std::string, std::vector<int_point>>& rois, std::unordered_map<std::string, data_struct::Spectra<_real> > &int_specs)
{

    bool parsed = false;

    std::ifstream file(path);
    // json reader
    Json::Reader reader;
    // this will contain complete JSON data
    Json::Value completeJsonData;
    // reader reads the data and stores it in completeJsonData
    parsed = reader.parse(file, completeJsonData);
    
    if (completeJsonData.isMember(STR_MAPS_ROIS))
    {
        const Json::Value json_rois = completeJsonData[STR_MAPS_ROIS];
        for (int i = 0; i < json_rois.size(); ++i)
        {
            const Json::Value json_map_roi = json_rois[i];
            
            std::string roi_name = std::to_string(i);
            // get the value associated with grade key
            if (json_map_roi.isMember(STR_MAP_ROI_NAME))
            {
                roi_name = json_map_roi[STR_MAP_ROI_NAME].asString();
            }
            logI << "Loading roi name :" << roi_name << " from file " << path << "\n";
            // load pixel locations
            const Json::Value json_pixel_locs = json_map_roi[STR_MAP_ROI_PIXEL_LOC];
            for (int index = 0; index < json_pixel_locs.size(); ++index)
            {
                const Json::Value json_ppair = json_pixel_locs[index];

                rois[roi_name].push_back({ json_ppair[0].asInt() ,json_ppair[1].asInt() });
            }
            // load integrated spectras
            if (json_map_roi.isMember(STR_MAP_ROI_INT_SPEC))
            {
                const Json::Value json_int_specs = json_map_roi[STR_MAP_ROI_INT_SPEC];
                for (int j = 0; j < json_int_specs.size(); ++j)
                {
                    const Json::Value json_spectra = json_int_specs[j];
                    if (json_spectra.isMember(STR_MAP_ROI_INT_SPEC_FILENAME)
                        && json_spectra.isMember(STR_SPECTRA)
                        && json_spectra.isMember(STR_ELT)
                        && json_spectra.isMember(STR_ERT)
                        && json_spectra.isMember(STR_ICR)
                        && json_spectra.isMember(STR_OCR))
                    {
                        std::string filename = json_spectra[STR_MAP_ROI_INT_SPEC_FILENAME].asString();
                        const Json::Value json_spectra_values = json_spectra[STR_SPECTRA];
                        int_specs[filename].resize(json_spectra_values.size());
                        int_specs[filename].elapsed_livetime(json_spectra[STR_ELT].asDouble());
                        int_specs[filename].elapsed_realtime(json_spectra[STR_ERT].asDouble());
                        int_specs[filename].input_counts(json_spectra[STR_ICR].asDouble());
                        int_specs[filename].output_counts(json_spectra[STR_OCR].asDouble());

                        for (int k = 0; k < json_spectra_values.size(); ++k)
                        {
                            int_specs[filename](k) = json_spectra_values[k].asDouble();
                        }
                    }
                }
            }
        }
    }
    else
    {
        logE << "Could not find key: " << STR_MAPS_ROIS << " . Skipping file: " << path << "\n";
        parsed = false;
    }
    return parsed;
}

DLL_EXPORT bool load_v9_rois(std::string path, std::map<std::string, std::vector<int_point>> &rois);

}// end namespace aps
}// end namespace file
}// end namespace io

#endif // APS_ROI
