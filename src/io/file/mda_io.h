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

#include "support/mda_utils/mda-load.h"
#include "data_struct/element_info.h"
#include "data_struct/spectra_volume.h"
#include "data_struct/quantification_standard.h"
#include "data_struct/params_override.h"
#include "data_struct/analysis_job.h"

namespace io
{
namespace file
{

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

    void unload();

    struct mda_file* get_scan_ptr() { return _mda_file; }

    bool load_spectra_volume(std::string path,
                            size_t detector_num,
                            data_struct::Spectra_Volume* vol,
                            bool hasNetCDF,
                            data_struct::Params_Override *override_values,
                            data_struct::Quantification_Standard * quantification_standard);

    bool load_spectra_volume_with_callback(std::string path,
                                        size_t detector_num_start,
                                        size_t detector_num_end,
                                        bool hasNetCDF,
                                        data_struct::Analysis_Job *analysis_job,
										data_struct::IO_Callback_Func_Def callback_func,
                                        void *user_data);

    int find_scaler_index(struct mda_file* mda_file, std::string det_name, real_t& val);

    int get_multiplied_dims(std::string path);

    int get_rank_and_dims(std::string path, int* dims);

    int rows() { return _rows; }

    int cols() { return _cols; }

    inline bool is_single_row_scan() {return _is_single_row;}

private:

    bool _find_theta(std::string pv_name, float* theta_out);

    //void _load_detector_meta_data(data_struct::Detector * detector);
    bool _is_single_row;

    /**
     * @brief _mda_file: mda helper structure
     */
    struct mda_file* _mda_file;

    /**
     * @brief _mda_file_info: lazy load struct
     */
    mda_fileinfo *_mda_file_info;

    int _rows;

    int _cols;

};

DLL_EXPORT bool load_henke_from_xdr(std::string filename, data_struct::Element_Info_Map* element_map);

}// end namespace file
}// end namespace io

#endif // MDA_IO_H
