/***
Copyright (c) 2023, UChicago Argonne, LLC. All rights reserved.

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

/// Initial Author <2023>: Arthur Glowacki



#ifndef EDF_IO_H
#define EDF_IO_H

#include <fstream>

#include "data_struct/spectra_line.h"
#include "data_struct/element_info.h"
#include "data_struct/fit_parameters.h"
#include "data_struct/detector.h"

using namespace data_struct;

namespace io
{
namespace file
{
namespace edf
{
   
    // ----------------------------------------------------------------------------
    template<typename T_real>
    bool check_header_and_load(std::string filename, char** buffer, int& length, int& idx, std::map<std::string, std::string> &header_map)
    {
        char *buf = nullptr;
        std::ifstream file_stream(filename, std::ios_base::binary);
        try
        {
            file_stream.seekg(0, std::ios_base::end);
            length = file_stream.tellg();

            if (length < 256)
            {
                logW << "File length smaller than 256\n";
                return false;
            }

            buf = new char[length];
            // rewind file and start reading data
            file_stream.clear();
            file_stream.seekg(0);

            file_stream.read(buf, length);

            file_stream.close();
            // first byte should be '{'
            if (buf[0] != '{')
            {
                logW << "File header does not start with '{'\n";

                delete[] buf;
                buffer = nullptr;
                return false;
            }

            std::string header = "";
            idx = 0;
            int last_i = 0;
            // find end of header '}'
            for (int i = 1; i < length; i++)
            {
                if (buf[i] == '}')
                {
                    header = std::string(buf, i + 1);
                    idx = i + 2; // 1 char for '}' and 1 for '\n'
                    break;
                }
                else if (buf[i] == '\n')
                {
                    std::string entry(buf+last_i, (i - last_i) + 1);
                    int e_idx = entry.find('=');
                    if (e_idx >= 0)
                    {
                        std::string key = entry.substr(0, e_idx);
                        key.erase(std::remove(key.begin(), key.end(), '\n'), key.end());
                        key.erase(std::remove(key.begin(), key.end(), '\r'), key.end());
                        key.erase(std::remove(key.begin(), key.end(), ' '), key.end());
                        std::string value = entry.substr(e_idx + 1);
                        value.erase(std::remove(value.begin(), value.end(), '\n'), value.end());
                        value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
                        value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
                        value.erase(std::remove(value.begin(), value.end(), ';'), value.end());
                        header_map[key] = value;
                    }
                    last_i = i+1;
                }
            }
        }
        catch (std::exception& e)
        {
            if (buf != nullptr)
            {
                delete[] buf;
                buffer = nullptr;
            }

            if (file_stream.eof() == 0 && (file_stream.bad() || file_stream.fail()))
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << file_stream.fail()
                    << "\neofbit: " << file_stream.eof()
                    << "\nbadbit: " << file_stream.bad() << "\n";
            }
            return false;
        }
        *buffer = buf;
        return true;
    }

    // ----------------------------------------------------------------------------

    template<typename T_real>
    DLL_EXPORT bool load_spectra_line(std::string filename, data_struct::Spectra_Line<T_real> *spec_line)
    {
        char* buffer = nullptr;
        std::map<std::string, std::string> header;
        try
        {
            int length = 0;
            int idx = 0;
            if (false == check_header_and_load<T_real>(filename, &buffer, length, idx, header))
            {
                return false;
            }

            int* val = (int*)(buffer + idx);
            for (int col = 0; col < spec_line->size(); col++)
            {
                for (int samp = 0; samp < 2048; samp++)
                {
                    (*spec_line)[col][samp] = static_cast<T_real>(*val);
                    val++;
                }
            }
            if (buffer != nullptr)
            {
                delete[] buffer;
            }
        }
        catch (std::exception& e)
        {
            if (buffer != nullptr)
            {
                delete[] buffer;
            }
            logE << e.what() << "\n";
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------
   
    template<typename T_real>
    DLL_EXPORT bool load_spectra_line_meta(std::string filename, int detector_num, int row, data_struct::Spectra_Line<T_real>* spec_line, data_struct::Scaler_Map<T_real> &dead_time_map, data_struct::Scaler_Map<T_real>& output_map)
    {
        char* buffer = nullptr;
        std::map<std::string, std::string> header;
        try
        {
            int length = 0;
            int idx = 0;
            if (false == check_header_and_load<T_real>(filename, &buffer, length, idx, header))
            {
                return false;
            }

            int num_vals = (length - idx) / sizeof(double);
            int col = 0;
            double* val = (double*)(buffer + idx);
            for (int i = 0; i < num_vals; i+=6)
            {
                //  chnum output icr ocr livetime deadtime
                double det = static_cast<double>(*val);
                val++;
                double output = static_cast<double>(*val);
                val++;
                double icr = static_cast<double>(*val);
                val++;
                double ocr = static_cast<double>(*val);
                val++;
                double livetime = static_cast<double>(*val);
                val++;
                double deadtime = static_cast<double>(*val);
                val++;

                if (detector_num == (int)det)
                {
                    output_map.values(row, col) = static_cast<T_real>(output);
                    dead_time_map.values(row, col) = static_cast<T_real>(deadtime);
                    (*spec_line)[col].elapsed_livetime(static_cast<T_real>(livetime));
                    (*spec_line)[col].elapsed_realtime(static_cast<T_real>(livetime));
                    (*spec_line)[col].input_counts(static_cast<T_real>(icr));
                    (*spec_line)[col].output_counts(static_cast<T_real>(ocr));
                    col++;
                }
            }

            if (buffer != nullptr)
            {
                delete[] buffer;
            }
        }
        catch (std::exception& e)
        {
            if (buffer != nullptr)
            {
                delete[] buffer;
            }
            logE << e.what() << "\n";
            return false;
        }

        return true;
    }

    // ----------------------------------------------------------------------------

    template<typename T_real>
    DLL_EXPORT bool load_scaler(std::string filename, data_struct::Scaler_Map<T_real>& scaler_map)
    {
        char* buffer = nullptr;
        std::map<std::string, std::string> header;
        try
        {
            int length = 0;
            int idx = 0;
            if (false == check_header_and_load<T_real>(filename, &buffer, length, idx, header))
            {
                return false;
            }

            float* val = (float*)(buffer + idx);
            for (int col = 0; col < scaler_map.values.cols(); col++)
            {
                for (int row = 0; row < scaler_map.values.rows(); row++)
                {   
                    scaler_map.values(row, col) = static_cast<T_real>(*val);
                    val++;
                }
            }

            if (buffer != nullptr)
            {
                delete[] buffer;
            }
        }
        catch (std::exception& e)
        {
            if (buffer != nullptr)
            {
                delete[] buffer;
            }
            logE << e.what() << "\n";
            return false;
        }

        return true;
    }
    
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

}// end namespace edf
}// end namespace file
}// end namespace io

#endif // EDF_IO_H
