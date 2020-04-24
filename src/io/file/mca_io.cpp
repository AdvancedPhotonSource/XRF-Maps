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



#include "mca_io.h"

#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <iostream>



namespace io
{
namespace file
{
namespace mca
{

bool load_integrated_spectra(std::string path, data_struct::Spectra* spectra, unordered_map<string, string>& pv_map)
{
    std::ifstream paramFileStream(path);

    logI<<"Loading:  "<<path<<"\n";

    if (paramFileStream.is_open() )
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
                    if(value != "3.1")
                    {
                        logW<<"MCA_IO supports version 3.1. MAy not load this version "<<value<<" properly.\n";
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
                    if(ivalue > 1)
                    {
                        logW<<"MCA_IO only supports loading 1 channel. This file has "<<value<<" channels.\n";
                    }
                }
                else if (tag == "REAL_TIME")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    float fvalue = std::stof(value);
                    spectra->elapsed_realtime(fvalue);
                }
                else if (tag == "LIVE_TIME")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    float fvalue = std::stof(value);
                    spectra->elapsed_livetime(fvalue);
                }
                else if (tag == "CAL_OFFSET")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    ///float fvalue = std::stof(value);
                }
                else if (tag == "CAL_SLOPE")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    ///float fvalue = std::stof(value);
                }
                else if (tag == "CAL_QUAD")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                    ///float fvalue = std::stof(value);
                }
                else if (tag == "ENVIRONMENT")
                {
                    std::string value;
                    std::getline(strstream, value);
                    value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                    value.erase(std::remove(value.begin(), value.end(), ' '), value.end());

                    int f_idx = value.find('=');
                    if(f_idx > -1)
                    {
                        string pv_name = value.substr(0,f_idx);
                        value = value.substr(f_idx+1);
                        if(value[0] == '"')
                        {
                            f_idx = value.find('"', 1);
                            if(f_idx > -1)
                            {
                                string pv_value = value.substr(1, f_idx-1);
                                pv_map[pv_name] = pv_value;
                            }
                        }
                    }
                }
                else if (tag == "DATA")
                {
                    for (int i=0; i<spectra->size(); i++)
                    {
                        std::string value;
                        std::getline(paramFileStream, value, '\n');
                        value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                        value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                        value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                        float fvalue = std::stof(value);
                        (*spectra)[i] = fvalue;
                    }
                }
            }
        }
        catch(std::exception& e)
        {
            if (paramFileStream.eof() == 0 && (paramFileStream.bad() || paramFileStream.fail()) )
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << paramFileStream.fail()
                    << "\neofbit: " << paramFileStream.eof()
                    << "\nbadbit: " << paramFileStream.bad() << "\n";
            }
        }

        paramFileStream.close();
        return true;
    }
    return false;

}

} //end namespace mca
} //end namespace file
}// end namespace io
