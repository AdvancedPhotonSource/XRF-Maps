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

#ifndef MDA_IO_H
#define MDA_IO_H

#include "base_file_io.h"
#include "mda-load.h"
#include "element_info.h"
#include "spectra_volume.h"
#include "detector.h"

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
                                     data_struct::xrf::Detector *detector,
                                     data_struct::xrf::Spectra_Volume* vol,
                                     std::unordered_map< std::string, std::string > *extra_override_values);

private:

    int _find_2d_detector_index(std::string det_name);

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
