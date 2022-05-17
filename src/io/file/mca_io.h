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



#ifndef MCA_IO_H
#define MCA_IO_H

#include "data_struct/spectra.h"
#include "data_struct/fit_parameters.h"

using namespace data_struct;

namespace io
{
namespace file
{

namespace mca
{

//-----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT bool load_integrated_spectra(std::string path, data_struct::Spectra<T_real>* spectra, unordered_map<string, T_real>& pv_map)
{
    std::ifstream paramFileStream(path);

    logI << "Loading:  " << path << "\n";

    if (paramFileStream.is_open())
    {
        //paramFileStream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        paramFileStream.exceptions(std::ifstream::failbit);
        //std::string line;
        std::string tag;
        try
        {
            for (std::string line; std::getline(paramFileStream, line, '\n'); )
                //while(std::getline(paramFileStream, line))
            {
                std::istringstream strstream(line);
                std::getline(strstream, tag, ':');
                //logD<<"tag : "<<tag<<"\n";
                if (tag == "VERSION")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    if (value != "3.1")
                    {
                        logW << "MCA_IO supports version 3.1. MAy not load this version " << value << " properly.\n";
                    }
                    try
                    {
                        pv_map[tag] = parse_input_real<T_real>(value);

                    }
                    catch (std::exception& e)
                    {
                        pv_map[tag] = 0.;
                    }
                }
                else if (tag == "CHANNELS")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    int ivalue = std::stoi(value);
                    try
                    {
                        pv_map[tag] = parse_input_real<T_real>(value);

                    }
                    catch (std::exception& e)
                    {
                        pv_map[tag] = 0.;
                    }
                    spectra->resize(ivalue);
                    spectra->Zero(ivalue);
                }
                else if (tag == "DATE")
                {

                }
                else if (tag == "ELEMENTS")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    int ivalue = std::stoi(value);
                    if (ivalue > 1)
                    {
                        logW << "MCA_IO only supports loading 1 channel. This file has " << value << " channels.\n";
                    }
                }
                else if (tag == "T_realIME")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    float fvalue = parse_input_real<T_real>(value);
                    spectra->elapsed_realtime(fvalue);
                }
                else if (tag == "LIVE_TIME")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    float fvalue = parse_input_real<T_real>(value);
                    spectra->elapsed_livetime(fvalue);
                }
                else if (tag == "CAL_OFFSET")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    try
                    {
                        pv_map[tag] = parse_input_real<T_real>(value);

                    }
                    catch (std::exception& e)
                    {
                        pv_map[tag] = 0.;
                    }
                }
                else if (tag == "CAL_SLOPE")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    try
                    {
                        pv_map[tag] = parse_input_real<T_real>(value);

                    }
                    catch (std::exception& e)
                    {
                        pv_map[tag] = 0.;
                    }
                }
                else if (tag == "CAL_QUAD")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    try
                    {
                        pv_map[tag] = parse_input_real<T_real>(value);

                    }
                    catch (std::exception& e)
                    {
                        pv_map[tag] = 0.;
                    }
                }
                else if (tag == "ENVIRONMENT")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());

                    int f_idx = value.find('=');
                    if (f_idx > -1)
                    {
                        string pv_name = value.substr(0, f_idx);
                        value = value.substr(f_idx + 1);
                        if (value[0] == '"')
                        {
                            f_idx = value.find('"', 1);
                            if (f_idx > -1)
                            {
                                string pv_value = value.substr(1, f_idx - 1);
                                try
                                {
                                    pv_map[pv_name] = parse_input_real<T_real>(pv_value);

                                }
                                catch (std::exception& e)
                                {
                                    pv_map[pv_name] = 0.;
                                }
                            }
                        }
                    }
                }
                else if (tag == "DATA")
                {
                    for (int i = 0; i < spectra->size(); i++)
                    {
                        std::string value;
                        std::getline(paramFileStream, value, '\n');
                        value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                        value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                        value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                        float fvalue = parse_input_real<T_real>(value);
                        (*spectra)[i] = fvalue;
                    }
                }
            }
        }
        catch (std::exception& e)
        {
            if (paramFileStream.eof() == 0 && (paramFileStream.bad() || paramFileStream.fail()))
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << paramFileStream.fail()
                    << "\neofbit: " << paramFileStream.eof()
                    << "\nbadbit: " << paramFileStream.bad() << "\n";
                return false;
            }
        }

        paramFileStream.close();
        return true;
    }
    return false;

}

//-----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT bool save_integrated_spectra(std::string path, data_struct::Spectra<T_real>* spectra, unordered_map<string, T_real>& pv_map)
{
    std::ofstream paramFileStream(path);

    T_real sr_current = 0.0;
    T_real us_ic = 0.0;
    T_real ds_ic = 0.0;
    T_real offset = 0.0;
    T_real slope = 0.0;
    T_real quad = 0.0;

    if (pv_map.count(STR_SR_CURRENT) > 0)
    {
        sr_current = pv_map.at(STR_SR_CURRENT);
    }
    if (pv_map.count(STR_US_IC) > 0)
    {
        us_ic = pv_map.at(STR_US_IC);
    }
    if (pv_map.count(STR_DS_IC) > 0)
    {
        ds_ic = pv_map.at(STR_DS_IC);
    }
    if (pv_map.count(STR_ENERGY_OFFSET) > 0)
    {
        offset = pv_map.at(STR_ENERGY_OFFSET);
    }
    if (pv_map.count(STR_ENERGY_SLOPE) > 0)
    {
        slope = pv_map.at(STR_ENERGY_SLOPE);
    }
    if (pv_map.count(STR_ENERGY_QUADRATIC) > 0)
    {
        quad = pv_map.at(STR_ENERGY_QUADRATIC);
    }


    logI << "Saving:  " << path << "\n";

    if (spectra == nullptr)
    {
        return false;
    }
    if (paramFileStream.is_open())
    {


        paramFileStream << "VERSION: 3.1\n";
        paramFileStream << "ELEMENTS: 1\n";
        paramFileStream << "DATE: \n";
        paramFileStream << "CHANNELS: " << spectra->size() << "\n";
        paramFileStream << "T_realIME: " << spectra->elapsed_realtime() << "\n";
        paramFileStream << "LIVE_TIME: " << spectra->elapsed_livetime() << "\n";
        paramFileStream << "CAL_OFFSET: " << offset << "\n";
        paramFileStream << "CAL_SLOPE: " << slope << "\n";
        paramFileStream << "CAL_QUAD: " << quad << "\n";
        paramFileStream << "ENVIRONMENT: SRCURRENT=\"" << sr_current << "\"\n";
        paramFileStream << "ENVIRONMENT: UPSTREAM_IONCHAMBER=\"" << us_ic << "\"\n";
        paramFileStream << "ENVIRONMENT: DOWNSTREAM_IONCHAMBER=\"" << ds_ic << "\"\n";
        paramFileStream << "DATA: \n";
        for (int i = 0; i < spectra->size(); i++)
        {
            paramFileStream << (*spectra)(i) << "\n";
        }

        paramFileStream.close();
        return true;
    }
    return false;

}

//-----------------------------------------------------------------------------

}// end namespace mca


//DLL_EXPORT bool load_element_info_from_csv(std::string filename);

}// end namespace file
}// end namespace io

#endif // MCA_IO_H
