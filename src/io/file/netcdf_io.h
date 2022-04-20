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



#ifndef NetCDF_IO_H
#define NetCDF_IO_H

#include "data_struct/spectra_volume.h"
#include <netcdf.h>
#include <mutex>

namespace io
{
namespace file
{

#ifdef _REAL_DOUBLE
	#define nc_get_vars_real nc_get_vars_double 
#else
	#define nc_get_vars_real nc_get_vars_float
#endif

enum class E_load_type { LINE = 0, CALLBACKF = 1, INTEGRATED = 2 };

template<typename T_real>
class DLL_EXPORT NetCDF_IO
{
public:

    static NetCDF_IO* inst();

    ~NetCDF_IO();

    /**
     * @brief load_spectra_line
     * @param path
     * @param detector
     * @param spec_line
     * @return the number of spectra loaded. 0 if fail.
     */
    size_t load_spectra_line(std::string path, size_t detector, data_struct::Spectra_Line<T_real>* spec_line);

    bool load_spectra_line_with_callback(std::string path,
										const std::vector<size_t>& detector_num_arr,
                                        int row,
                                        size_t max_rows,
                                        size_t max_cols,
                                        data_struct::IO_Callback_Func_Def<T_real> callback_fun,
                                        void* user_data);

    size_t load_spectra_line_integrated(std::string path, size_t detector, size_t line_size, data_struct::Spectra<T_real>* spectra);

private:
    NetCDF_IO();

    size_t _load_spectra(E_load_type ltype,
                        std::string path,
                        size_t detector,
                        data_struct::Spectra_Line<T_real>* spec_line,
                        size_t line_size, 
                        data_struct::Spectra<T_real>* spectra);

    static NetCDF_IO *_this_inst;

    static std::mutex _mutex;

};

}// end namespace file
}// end namespace io

#endif // NetCDF_IO_H
