/***

Copyright (c) 2016 Arthur Glowacki

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

***/

#ifndef HDF5_IO_H
#define HDF5_IO_H

#include <mutex>
#include "base_file_io.h"
#include "base_dataset.h"
#include "spectra_volume.h"
#include "fit_element_map.h"
#include "hdf5.h"

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

    virtual void load_spectra_volume(std::string path, HDF5_Spectra_Layout layout, data_struct::xrf::Spectra_Volume* spec_vol);

    //DLL_EXPORT void load_spectra_volume(std::string path, HDF5_Spectra_Layout layout, data_struct::xrf::Spectra_Volume* spec_vol);

    bool save_spectra_volume(const std::string filename,
                             const std::string path,
                             const data_struct::xrf::Spectra_Volume * const spectra_volume,
                             size_t row_idx_start=0,
                             int row_idx_end=-1,
                             size_t col_idx_start=0,
                             int col_idx_end=-1);

    bool save_element_fits(const std::string filename,
                           const std::string path,
                           const data_struct::xrf::Fit_Count_Dict * const element_counts,
                           size_t row_idx_start=0,
                           int row_idx_end=-1,
                           size_t col_idx_start=0,
                           int col_idx_end=-1);

private:
    void parse_group_info(hid_t h5file, hid_t id);
    void parse_dataset_info(hid_t h5file, hid_t id);
    void parse_attr_info(hid_t h5file, hid_t id);

    std::mutex _mutex;

};

}// end namespace file
}// end namespace io

#endif // HDF5_IO_H
