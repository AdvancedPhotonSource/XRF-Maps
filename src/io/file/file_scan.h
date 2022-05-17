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

/// Initial Author <2022>: Arthur Glowacki

#ifndef _FILE_SCAN_H
#define _FILE_SCAN_H

#if defined _WIN32
#include "support/direct/dirent.h"
#else
#include <dirent.h>
#endif

#include "io/file/netcdf_io.h"
#include "io/file/mda_io.h"
#include "io/file/mca_io.h"
#include "io/file/hdf5_io.h"
#include "io/file/csv_io.h"

#include "data_struct/spectra_volume.h"



namespace io
{
    namespace file
    {

        struct file_name_size
        {
            file_name_size(std::string name, long size) { filename = name; total_rank_size = size; }
            std::string filename;
            long total_rank_size;
        };

        bool compare_file_size(const file_name_size& first, const file_name_size& second);

        class DLL_EXPORT File_Scan
        {

        public:
            static File_Scan* inst();

            ~File_Scan();

            void populate_netcdf_hdf5_files(std::string dataset_dir);

            //void check_and_create_dirs(std::string dataset_directory);

            std::vector<std::string> find_all_dataset_files(std::string dataset_directory, std::string search_str);

            void sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string>* dataset_files);

            const std::vector<std::string>& netcdf_files() {  return _netcdf_files; }

            const std::vector<std::string>& bnp_netcdf_files() { return _bnp_netcdf_files; }

            const std::vector<std::string>& hdf_files() { return _hdf_files; }

            const std::vector<std::string>& hdf_xspress_files() { return _hdf_xspress_files; }

            const std::vector<std::string>& hdf_emd_files() { return _hdf_emd_files; }

        private:

            File_Scan();

            static File_Scan* _this_inst;

            std::vector<std::string> _netcdf_files;
            std::vector<std::string> _bnp_netcdf_files;
            std::vector<std::string> _hdf_files;
            std::vector<std::string> _hdf_xspress_files;
            //std::vector<std::string> _hdf_confocal_files;
            std::vector<std::string> _hdf_emd_files;
        };


        // ----------------------------------------------------------------------------
    }
}// end namespace io

#endif // HL_FILE_IO_H
