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
    DLL_EXPORT bool load_spectra_line(std::string filename, data_struct::Spectra_Line<T_real> *spec_line)
    {
        std::ifstream file_stream(filename, std::ios_base::binary);
        try
        {
            file_stream.seekg(0, std::ios_base::end);
            int length = file_stream.tellg();

            if (length < 256)
            {
                logW << "File length smaller than 256\n";
                return false;
            }

            char* buffer = new char[length];
            // rewind file and start reading data
            file_stream.clear();
            file_stream.seekg(0);

            file_stream.read(buffer, length);

            file_stream.close();
            // first byte should be '{'
            if (buffer[0] != '{')
            {
                logW << "File header does not start with '{'\n";

                delete[] buffer;
                return false;
            }

            std::string header = "";
            int idx = 0;
            // find end of header '}'
            for (int i = 1; i < length; i++)
            {
                if (buffer[i] == '}')
                {
                    header = std::string(buffer, i+1);
                    idx = i + 2; // 1 char for '}' and 1 for '\n'
                    break;
                }
            }
     
            logI << header << "\n";
            
            float* val = (float*)(buffer + idx);
            for (int col = 0; col < spec_line->size(); col++)
            {
                for (int samp = 0; samp < 2048; samp++)
                {
                    (*spec_line)[col][samp] = *val;
                    val++;
                }
            }

            delete[] buffer;
        }
        catch (std::exception& e)
        {
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

        return true;
    }

    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------
    // ----------------------------------------------------------------------------

}// end namespace edf
}// end namespace file
}// end namespace io

#endif // EDF_IO_H
