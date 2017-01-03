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


#ifndef HDF5_IO_H
#define HDF5_IO_H

#include <mutex>
#include <list>
#include "base_file_io.h"
#include "base_dataset.h"
#include "spectra_volume.h"
#include "fit_element_map.h"
#include "hdf5.h"

//Include mda data structures to save scalers
#include "mda_io.h"

#include "quantification_standard.h"

namespace io
{
namespace file
{

struct HDF5_Spectra_Layout
{
   HDF5_Spectra_Layout()
   {
      spectrum_dim = 0;
      row_dim = 1;
      col_dim = 2;
   }

   HDF5_Spectra_Layout(unsigned int row, unsigned int col, unsigned int spectrum)
   {
      row_dim = row;
      col_dim = col;
      spectrum_dim = spectrum;
   }

   unsigned int row_dim;
   unsigned int col_dim;
   unsigned int spectrum_dim;
};


struct HDF5_Range
{
    HDF5_Range(unsigned int rank)
    {
        this->rank = rank;
        offset = new unsigned long[rank];
        count = new unsigned long[rank];
        for (unsigned int i=0; i < rank; i++)
        {
            offset[i] = 0;
            count[i] = 0;
        }
    }

    unsigned int rank;
    unsigned long* offset;
    unsigned long* count;
};

class DLL_EXPORT HDF5_IO : public Base_File_IO
{
public:
    HDF5_IO();

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

    //void load_dataset2(std::string path, HDF5_Spectra_Layout layout, HDF5_Range range, data_struct::xrf::Spectra_Volume* spec_vol);

    bool load_spectra_volume(std::string path, size_t detector_num, data_struct::xrf::Spectra_Volume* spec_vol);

    //DLL_EXPORT void load_spectra_volume(std::string path, HDF5_Spectra_Layout layout, data_struct::xrf::Spectra_Volume* spec_vol);

    hid_t start_save_seq(const std::string filename);

    bool save_spectra_volume(const hid_t file_id,
                             const std::string path,
                             data_struct::xrf::Spectra_Volume * spectra_volume,
                             size_t row_idx_start=0,
                             int row_idx_end=-1,
                             size_t col_idx_start=0,
                             int col_idx_end=-1);

    bool save_element_fits(const hid_t file_id,
                           const std::string path,
                           const data_struct::xrf::Fit_Count_Dict * const element_counts,
                           size_t row_idx_start=0,
                           int row_idx_end=-1,
                           size_t col_idx_start=0,
                           int col_idx_end=-1);

    bool save_quantification(const hid_t file_id,
                             data_struct::xrf::Quantification_Standard * quantification_standard,
                             size_t row_idx_start=0,
                             int row_idx_end=-1,
                             size_t col_idx_start=0,
                             int col_idx_end=-1);

    bool save_scalers(const std::string filename,
                      size_t detector_num,
                      struct mda_file *mda_scalers,
                      std::unordered_map< std::string, std::string > *extra_override_values,
                      size_t row_idx_start=0,
                      int row_idx_end=-1,
                      size_t col_idx_start=0,
                      int col_idx_end=-1);

    bool end_save_seq(const hid_t);

private:


    //bool save_scalar(const hid_t group_id,  mda_scan *mda_scalers)

    struct scaler_struct
    {
        scaler_struct(std::string name, int mda_idx_, int hdf_idx_)
        {
             hdf_name = name;
             mda_idx = mda_idx_;
             hdf_idx = hdf_idx_;
        }
        int mda_idx;
        int hdf_idx;
        std::string hdf_name;
    };

    //void _parse_group_info(hid_t h5file, hid_t id);
    //void _parse_dataset_info(hid_t h5file, hid_t id);
    //void _parse_attr_info(hid_t h5file, hid_t id);

    std::mutex _mutex;

};

}// end namespace file
}// end namespace io

#endif // HDF5_IO_H
