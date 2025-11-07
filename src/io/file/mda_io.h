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


#ifndef MDA_IO_H
#define MDA_IO_H

#include "support/mdautils-1.4.1/mda-load.h"
#include "data_struct/element_info.h"
#include "data_struct/spectra_volume.h"
#include "data_struct/quantification_standard.h"
#include "data_struct/scan_info.h"
#include "data_struct/analysis_job.h"

#include <string>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "data_struct/scaler_lookup.h"


#ifdef XDR_HACK
#include "support/mdautils-1.4.1/xdr_hack.h"
#else
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif


namespace io
{
namespace file
{

template<typename T_real>
class DLL_EXPORT MDA_IO
{
public:

    /**
     * @brief MDA_IO
     * @param filename
     * @return
     */
    MDA_IO();

    /**
     * @brief ~MDA_IO
     * @return
     */
    ~MDA_IO();

    bool load_scalers(std::string path);

    bool load_spectra_volume(std::string path,
                            size_t detector_num,
                            data_struct::Spectra_Volume<T_real>* vol,
                            bool hasNetCDF,
                            bool subtract_two_cols);

    bool load_spectra_volume_with_callback(std::string path,
										                    const std::vector<size_t>& detector_num_arr,
                                        bool hasNetCDF,
                                        data_struct::Analysis_Job<T_real>* analysis_job,
                                        size_t& out_rows,
                                        size_t& out_cols,
										                    data_struct::IO_Callback_Func_Def<T_real> callback_func,
                                        void *user_data);

	bool load_integrated_spectra(std::string path,
								size_t detector_num,
								data_struct::Spectra<T_real>*out_integrated_spectra,
								bool hasNetCDF);

    bool load_quantification_scalers(std::string path,
                                     data_struct::Params_Override<T_real>*override_values,
                                     bool hasNetCDF);

    void unload();
	
    //void search_and_update_amps(std::string us_amp_pv_str, std::string ds_amp_pv_str, T_real& out_us_amp, T_real& out_ds_amp);

    data_struct::Scan_Info<T_real>* get_scan_info() { return &_scan_info; }

    size_t get_num_integreated_spectra() { return _integrated_spectra_map.size(); }

    data_struct::ArrayTr<T_real>* get_integrated_spectra(unsigned int detector);

    void unload_int_spectra() { _integrated_spectra_map.clear(); }

    bool load_henke_from_xdr(std::string filename);

private:

    void _load_scalers(bool load_int_spec, bool hasNetCDF, bool subtract_two_cols);

    void _load_extra_pvs_vector();

    void _load_meta_info(bool hasNetCDF, bool subtract_two_cols);

    bool _find_theta(std::string pv_name, float* theta_out);

    /**
     * @brief _mda_file: mda helper structure
     */
    struct mda_file* _mda_file;

    /**
     * @brief _mda_file_info: lazy load struct
     */
    mda_fileinfo *_mda_file_info;

    bool _hasNetcdf;

    data_struct::Scan_Info<T_real> _scan_info;

    std::string _theta_pv_str;

    std::unordered_map<unsigned int, data_struct::ArrayTr<T_real>> _integrated_spectra_map;
};

DLL_EXPORT int mda_get_multiplied_dims(std::string path);

DLL_EXPORT int mda_get_rank_and_dims(std::string path, size_t* dims);


}// end namespace file
}// end namespace io

#endif // MDA_IO_H
