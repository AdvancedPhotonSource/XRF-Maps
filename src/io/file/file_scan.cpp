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

/// Initial Author <2017>: Arthur Glowacki

#include "file_scan.h"

namespace io
{
    namespace file
    {
        //static std::vector<std::string> netcdf_files;
        //static std::vector<std::string> bnp_netcdf_files;
        //static std::vector<std::string> hdf_files;
        //static std::vector<std::string> hdf_xspress_files;
        //std::vector<std::string> hdf_confocal_files;
        //static std::vector<std::string> hdf_emd_files;

        // ----------------------------------------------------------------------------

        bool compare_file_size(const file_name_size& first, const file_name_size& second)
        {
            return (first.total_rank_size > second.total_rank_size);
        }

        //-----------------------------------------------------------------------------

        File_Scan* File_Scan::_this_inst(nullptr);

        //-----------------------------------------------------------------------------

        File_Scan::File_Scan()
        {

        }

        //-----------------------------------------------------------------------------

        File_Scan* File_Scan::inst()
        {
            //std::lock_guard<std::mutex> lock(_mutex);

            if (_this_inst == nullptr)
            {
                _this_inst = new File_Scan();
            }
            return _this_inst;
        }

        //-----------------------------------------------------------------------------

        File_Scan::~File_Scan()
        {
            _netcdf_files.clear();
            _bnp_netcdf_files.clear();
            _hdf_files.clear();
            _hdf_xspress_files.clear();
            _hdf_emd_files.clear();
        }


        // ----------------------------------------------------------------------------

        void File_Scan::populate_netcdf_hdf5_files(std::string dataset_dir)
        {
            _netcdf_files.clear();
            _bnp_netcdf_files.clear();
            _hdf_files.clear();
            _hdf_xspress_files.clear();
            _hdf_emd_files.clear();

            std::replace(dataset_dir.begin(), dataset_dir.end(), '/', DIR_END_CHAR);
            if (dataset_dir[dataset_dir.length() - 1] != DIR_END_CHAR)
            {
                dataset_dir += DIR_END_CHAR;
            }

            //populate netcdf and hdf5 files for fly scans
            _netcdf_files = find_all_dataset_files(dataset_dir + "flyXRF" + DIR_END_CHAR, "_0.nc");
            _bnp_netcdf_files = find_all_dataset_files(dataset_dir + "flyXRF" + DIR_END_CHAR, "_001.nc");
            _hdf_files = find_all_dataset_files(dataset_dir + "flyXRF.h5" + DIR_END_CHAR, "_0.h5");
            _hdf_xspress_files = find_all_dataset_files(dataset_dir + "flyXspress" + DIR_END_CHAR, "_0.h5");
            //_hdf_confocal_files = find_all_dataset_files(dataset_dir , ".hdf5");
            _hdf_emd_files = find_all_dataset_files(dataset_dir, ".emd");
        }

        // ----------------------------------------------------------------------------

        std::vector<std::string> File_Scan::find_all_dataset_files(std::string dataset_directory, std::string search_str)
        {
            std::vector<std::string> dataset_files;
            logI << dataset_directory << " searching for " << search_str << "\n";
            DIR* dir;
            struct dirent* ent;
            size_t search_str_size = search_str.length();
            if ((dir = opendir(dataset_directory.c_str())) != NULL)
            {
                /* print all the files and directories within directory */
                while ((ent = readdir(dir)) != NULL)
                {
                    std::string fname(ent->d_name);
                    // check if extension is .mda
                    if (fname.size() > 4)
                    {
                        if (fname.rfind(search_str) == fname.size() - search_str_size)
                        {
                            dataset_files.push_back(fname);
                        }
                    }
                }
                closedir(dir);
            }
            else
            {
                /* could not open directory */
                logW << "Could not open directory " << dataset_directory << " using search string " << search_str << "\n";
            }

            logI << "found " << dataset_files.size() << "\n";
            return dataset_files;
        }

        // ----------------------------------------------------------------------------
        
        void File_Scan::find_all_dataset_files_by_list(std::string dataset_directory, std::vector<std::string>& search_strs, std::vector<std::string>& out_dataset_files)
        {
            DIR* dir;
            struct dirent* ent;
            if ((dir = opendir(dataset_directory.c_str())) != NULL)
            {
                /* print all the files and directories within directory */
                while ((ent = readdir(dir)) != NULL)
                {
                    std::string fname(ent->d_name);
                    // check if extension is .mda
                    if (fname.size() > 4)
                    {
                        for (auto& itr : search_strs)
                        {
                            size_t search_str_size = itr.length();
                            if (fname.rfind(itr) == fname.size() - search_str_size)
                            {
                                out_dataset_files.push_back(fname);
                            }
                        }
                    }
                }
                closedir(dir);
            }
            else
            {
                /* could not open directory */
                logW << "Could not open directory " << dataset_directory << "\n";
            }
            
        }

        // ----------------------------------------------------------------------------

        /*
        void File_Scan::check_and_create_dirs(std::string dataset_directory)
        {

            bool found_img_dat = false;
            bool found_output = false;
            logI << dataset_directory << "\n";
            DIR* dir;
            struct dirent* ent;
            if ((dir = opendir(dataset_directory.c_str())) != NULL)
            {
                / print all the files and directories within directory 
                while ((ent = readdir(dir)) != NULL)
                {
                    if (strcmp(ent->d_name, "img.dat") == 0)
                    {
                        found_img_dat = true;
                    }
                    if (strcmp(ent->d_name, "output") == 0)
                    {
                        found_output = true;
                    }
                    if (found_img_dat && found_output)
                    {
                        break;
                    }
                }
                closedir(dir);
            }
            else
            {
                // could not open directory 
                logW << "Could not open directory " << dataset_directory << "\n";
            }

            if (false == found_img_dat)
            {
                int retval = system(nullptr);
                std::string cmd = "mkdir " + dataset_directory + "img.dat";
                logI << cmd << "\n";
                retval = system(cmd.c_str());
                if (retval != 0)
                {
                    logE << "Could not create directory: " << cmd << " . May not be able to save results!\n";
                }
            }
            if (false == found_output)
            {
                int retval = system(nullptr);
                std::string cmd = "mkdir " + dataset_directory + "output";
                logI << cmd << "\n";
                retval = system(cmd.c_str());
                if (retval != 0)
                {
                    logE << "Could not create directory: " << cmd << " . May not be able to save results!\n";
                }
            }
            logI << "done" << "\n";

        }
        */
        // ----------------------------------------------------------------------------

        void File_Scan::sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string>* dataset_files)
        {
            // only supports soring mda files
            std::string ending = ".mda";
            io::file::MDA_IO<float> mda_io;
            logI << dataset_directory << " " << dataset_files->size() << " files" << "\n";
            std::list<file_name_size> f_list;

            for (const auto& itr : *dataset_files)
            {
                //check if file ends with .mda
                if (itr.compare(itr.length() - ending.length(), ending.length(), ending) == 0)
                {
                    std::string full_path = dataset_directory + DIR_END_CHAR + "mda" + DIR_END_CHAR + itr;
                    long fsize = file::mda_get_multiplied_dims(full_path);
                    f_list.push_back(file_name_size(itr, fsize));
                }
                else
                {
                    std::string full_path = dataset_directory + DIR_END_CHAR + itr;
                    std::ifstream in(full_path.c_str(), std::ifstream::ate | std::ifstream::binary);
                    long fsize = in.tellg();
                    if (fsize > -1)
                    {
                        f_list.push_back(file_name_size(itr, fsize));
                    }
                }
            }

            f_list.sort(compare_file_size);
            if (f_list.size() > 0)
            {
                dataset_files->clear();

                for (auto& itr : f_list)
                {
                    dataset_files->push_back(itr.filename);
                }
            }
            logI << "done" << "\n";
        }

        // ----------------------------------------------------------------------------
        // ----------------------------------------------------------------------------
    }// end namespace file
}// end namespace io
