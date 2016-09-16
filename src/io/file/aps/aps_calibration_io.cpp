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



#include "aps_calibration_io.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "mda_io.h"
#include "spectra_volume.h"

namespace io
{
namespace file
{
namespace aps
{

APS_Calibration_IO::APS_Calibration_IO()
{



}

APS_Calibration_IO::~APS_Calibration_IO()
{


}

bool APS_Calibration_IO::load(std::string path, std::string standard_filename, data_struct::xrf::Calibration_Standard* calibration, size_t detector_num)
{

    bool status = _load_standards_txt(path+standard_filename, calibration);

    if (!status)
    {
        //TODO: log error
        std::cout<<"error";
        return status;
    }

    // check if standard_file is mca ext
    _standard_filename.erase(std::remove_if(_standard_filename.begin(), _standard_filename.end(), ::isspace), _standard_filename.end());
    if(_standard_filename.substr(_standard_filename.find_last_of(".") + 1) == "mca")
    {

        std::ifstream standardFileStream(path+_standard_filename);
        if (standardFileStream.good())
        {
            status = _load_mca(path+_standard_filename, 0, calibration);
        }
        else
        {
            std::cout<<"Could not find standard file "<<_standard_filename<<std::endl;

            //check to see if 0 is on the end of it ( detector number )
            std::ifstream standardFileStream(path+_standard_filename+"0");
            if (standardFileStream.good())
            {
                status = _load_mca(path+_standard_filename+"0", 0, calibration);
                //for number of detectors
                //status = _load_mca(_standard_filename+"1"); ...

            }
            else
            {
                std::cout<<"Could not find standard file "<<_standard_filename+"0"<<std::endl;
                // check if mda file exists and convert to mca
                //TODO:

                std::string filename = _standard_filename.substr(0, _standard_filename.find_last_of(".") );
                std::string mda_filename = path + "mda/" + filename + ".mda";
                status = _load_mda(mda_filename, calibration);
                if (status)
                {
                    //_convert_mda_to_mca()
                }
                else
                {
                    std::cout<<"Could not find standard file "<<mda_filename<<std::endl;
                }
            }
        }

    }
    else if(_standard_filename.substr(_standard_filename.find_last_of(".") + 1) == "mda")
    {
        status = _load_mda(_standard_filename, calibration);
    }
    else
    {
        std::cout<<"Unknown standard file format "<<_standard_filename<<std::endl;
    }

    return status;

}


bool APS_Calibration_IO::_load_standards_txt(std::string path, data_struct::xrf::Calibration_Standard* calibration)
{

    std::ifstream paramFileStream(path);

    if (paramFileStream.is_open() )
    {
        paramFileStream.exceptions(std::ifstream::failbit);
        bool has_filename = false;
        bool has_elements = false;
        bool has_weights = false;
        //std::string line;
        std::string tag;

        std::vector<std::string> element_names;
        std::vector<real_t> element_weights;

        try
        {

            for (std::string line; std::getline(paramFileStream, line); )
            //while(std::getline(paramFileStream, line))
            {
                std::istringstream strstream(line);
                std::getline(strstream, tag, ':');
                //std::cout<<"tag : "<<tag<<std::endl;
                if (tag == "FILENAME")
                {
                    std::cout << line << std::endl;
                    std::getline(strstream, _standard_filename, ':');
                    std::cout << "Standard file name = "<< _standard_filename << std::endl;
                    has_filename = true;
                }
                else if (tag == "ELEMENTS_IN_STANDARD")
                {
                    std::string element_symb;
                    while(std::getline(strstream, element_symb, ','))
                    {
                        element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());
                        std::cout<<"Element : "<<element_symb<<std::endl;
                        element_names.push_back(element_symb);
                    }
                    has_elements = true;
                }
                else if (tag == "WEIGHT")
                {
                    std::string element_weight_str;
                    while(std::getline(strstream, element_weight_str, ','))
                    {
                        element_weight_str.erase(std::remove_if(element_weight_str.begin(), element_weight_str.end(), ::isspace), element_weight_str.end());
                        std::cout<<"Element weight: "<<element_weight_str<<std::endl;
                        real_t weight = std::stof(element_weight_str);
                        element_weights.push_back(weight);
                    }
                    has_weights = true;
                }

            }
        }
        catch(std::exception e)
        {
            if (paramFileStream.eof() == 0 && (paramFileStream.bad() || paramFileStream.fail()) )
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << paramFileStream.fail()
                    << "\neofbit: " << paramFileStream.eof()
                    << "\nbadbit: " << paramFileStream.bad() << std::endl;
            }
        }


        paramFileStream.close();
        if(has_filename && has_elements && has_weights)
        {
            if(element_names.size() == element_weights.size())
            {
                for(size_t i=0; i<element_names.size(); i++)
                {
                    calibration->append_element_weight(element_names[i], element_weights[i]);
                }
            }
            else
            {
                std::cout<<"Error: number of element names ["<<element_names.size()<<"] does not match number of element weights ["<<element_weights.size()<<"]!"<<std::endl;
            }

            return true;
        }

    }
    else
    {
        std::cout<<"Failed to open file "<<path<<std::endl;
    }
    return false;

}

bool APS_Calibration_IO::_load_mca(std::string path, int detector_num, data_struct::xrf::Calibration_Standard* calibration)
{

    std::ifstream paramFileStream(path);

    if (paramFileStream.is_open() )
    {
        paramFileStream.exceptions(std::ifstream::failbit);
        //std::string line;
        std::string tag;
        size_t num_channels = 0;
        data_struct::xrf::Spectra* spectra = calibration->spectra();
        std::string str_value;

        try
        {
            for (std::string line; std::getline(paramFileStream, line); )
            //while(std::getline(paramFileStream, line))
            {
                std::istringstream strstream(line);
                std::getline(strstream, tag, ':');

                //std::cout<<"tag : "<<tag<<std::endl;
                if (tag == "VERSION" || tag == "DATE")
                {
                    std::cout << line << std::endl;
                }
                else if (tag == "ELEMENTS")
                {

                }
                else if (tag == "CHANNELS")
                {
                    std::getline(strstream, str_value, ':');
                    num_channels = std::stoi(str_value);
                    spectra->resize(num_channels);

                }
                else if (tag == "REAL_TIME")
                {
                    std::getline(strstream, str_value, ':');
                    calibration->real_time(std::stoi(str_value));
                }
                else if (tag == "LIVE_TIME")
                {
                    std::getline(strstream, str_value, ':');
                    calibration->live_time(std::stoi(str_value));
                }
                else if (tag == "CAL_OFFSET")
                {
                    std::getline(strstream, str_value, ':');
                    calibration->offset(std::stoi(str_value));
                }
                else if (tag == "CAL_SLOPE")
                {
                    std::getline(strstream, str_value, ':');
                    calibration->slope(std::stoi(str_value));
                }
                else if (tag == "CAL_QUAD")
                {
                    std::getline(strstream, str_value, ':');
                    calibration->quad(std::stoi(str_value));
                }
                else if (tag == "ENVIRONMENT")
                {
                    /*
                    if (etag == "S:SRcurrentAI")
                        calibration->current(std::stoi(str_value));
                    elif etag == '2xfm:scaler1_cts1.B':
                        if IC_US == 0 :
                            IC_US = real_t(temp)
                    elif etag == '2xfm:scaler1_cts1.C':
                        if IC_DS == 0 :
                            IC_DS = real_t(temp)
                    elif etag == '2xfm:scaler3_cts1.B':
                        IC_US = real_t(temp)
                    elif etag == '2xfm:scaler3_cts1.C':
                        IC_DS = real_t(temp)
                    elif etag == '2idd:scaler1_cts1.C':
                        IC_US = real_t(temp)
                    elif etag == '2idd:scaler1_cts1.B':
                        IC_DS = real_t(temp)
                    elif etag == '8bmb:3820:scaler1_cts1.B':
                        IC_US = real_t(temp)
                    elif etag == '8bmb:3820:scaler1_cts1.C':
                        IC_DS = real_t(temp)
                    elif etag[5:] == 'A1sens_num.VAL':
                        amp[0, 0] = real_t(temp)
                    elif etag[5:] == 'A2sens_num.VAL':
                        amp[1, 0] = real_t(temp)
                    elif etag[5:] == 'A3sens_num.VAL':
                        amp[2, 0] = real_t(temp)
                    elif etag[5:] == 'A4sens_num.VAL':
                        amp[3, 0] = real_t(temp)
                    elif etag[5:] == 'A1sens_unit.VAL':
                        if (temp == "nA/V") or	(temp == "pA/V") or (temp == "uA/V") or (temp == "mA/V"):
                            if (temp == "pA/V") : amp[0, 1] = 0
                            if (temp == "nA/V") : amp[0, 1] = 1
                            if (temp == "uA/V") : amp[0, 1] = 2
                            if (temp == "mA/V") : amp[0, 1] = 3
                        else:
                            amp[0, 1] = real_t(temp)
                    elif etag[5:] == 'A2sens_unit.VAL':
                        if (temp == "nA/V") or	(temp == "pA/V") or (temp == "uA/V") or (temp == "mA/V"):
                            if (temp == "pA/V") : amp[1, 1] = 0
                            if (temp == "nA/V") : amp[1, 1] = 1
                            if (temp == "uA/V") : amp[1, 1] = 2
                            if (temp == "mA/V") : amp[1, 1] = 3
                        else:
                            amp[1, 1] = real_t(temp)
                    elif etag[5:] == 'A3sens_unit.VAL':
                        if (temp == "nA/V") or	(temp == "pA/V") or (temp == "uA/V") or (temp == "mA/V"):
                            if (temp == "pA/V") : amp[2, 1] = 0
                            if (temp == "nA/V") : amp[2, 1] = 1
                            if (temp == "uA/V") : amp[2, 1] = 2
                            if (temp == "mA/V") : amp[2, 1] = 3
                        else:
                            amp[2, 1] = real_t(temp)
                    elif etag[5:] == 'A4sens_unit.VAL':
                        if (temp == "nA/V") or (temp == "pA/V") or (temp == "uA/V") or (temp == "mA/V") :
                            if (temp == "pA/V") : amp[3, 1] = 0
                            if (temp == "nA/V") : amp[3, 1] = 1
                            if (temp == "uA/V") : amp[3, 1] = 2
                            if (temp == "mA/V") : amp[3, 1] = 3
                        else:
                            amp[3, 1] = real_t(temp)
                    */
                }
                else if (tag == "DATA")
                {
                    for(size_t i=0; i<num_channels; i++)
                    {
                        std::getline(paramFileStream, line);
                        (*spectra)[i] = std::stof(line);
                    }
                }
            }
        }
        catch(std::exception e)
        {
            if (paramFileStream.eof() == 0 && (paramFileStream.bad() || paramFileStream.fail()) )
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << paramFileStream.fail()
                    << "\neofbit: " << paramFileStream.eof()
                    << "\nbadbit: " << paramFileStream.bad() << std::endl;
            }
        }

        paramFileStream.close();

        return true;

    }
    return false;

}

bool APS_Calibration_IO::_load_mda(std::string path, data_struct::xrf::Calibration_Standard* calibration)
{

    io::file::MDA_IO mda_io;
    data_struct::xrf::Spectra_Volume dset;
    //bool status = mda_io.load_dataset(path, &dset);
    //return status;
	return false;
}



} //end namespace aps
} //end namespace file
}// end namespace io
