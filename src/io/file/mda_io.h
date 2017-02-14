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

#include "base_file_io.h"
#include "mda-load.h"
#include "element_info.h"
#include "spectra_volume.h"
#include "quantification_standard.h"
#include "params_override.h"

namespace io
{
namespace file
{

class DLL_EXPORT MDA_IO : public Base_File_IO
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

    /**
     * @brief lazy_load : Only load in the meta info, not the actual datasets
     * @param filename
     */
    virtual void lazy_load();

    /**
     * @brief load : Load the full dataset
     * @param filename
     */
    virtual bool load_dataset(std::string path, Base_Dataset* dset);

    virtual bool load_spectra_volume(std::string path,
                                     size_t detector_num,
                                     data_struct::xrf::Spectra_Volume* vol,
                                     bool hasNetCDF,
                                     data_struct::xrf::Params_Override *override_values,
                                     data_struct::xrf::Quantification_Standard * quantification_standard);

    int find_2d_detector_index(struct mda_file* mda_file, std::string det_name, int detector_num, real_t& val);

    int get_multiplied_dims(std::string path);

    int get_rank_and_dims(std::string path, int* dims);

private:

    void _load_detector_meta_data(data_struct::xrf::Detector * detector);

    /**
     * @brief _mda_file: mda helper structure
     */
    struct mda_file* _mda_file;

    /**
     * @brief _mda_file_info: lazy load struct
     */
    mda_fileinfo *_mda_file_info;

};

DLL_EXPORT bool load_henke_from_xdr(std::string filename, data_struct::xrf::Element_Info_Map* element_map);

}// end namespace file
}// end namespace io

#endif // MDA_IO_H
