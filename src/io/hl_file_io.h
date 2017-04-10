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

#ifndef HL_FILE_IO_H
#define HL_FILE_IO_H

#include <dirent.h>

#include "threadpool.h"

#include "netcdf_io.h"
#include "mda_io.h"
#include "hdf5_io.h"
#include "csv_io.h"

#include "spectra_volume.h"

#include "gaussian_model.h"

#include "element_info.h"

#include "aps_fit_params_import.h"

#include "base_fit_routine.h"

#include "fit_element_map.h"
#include "params_override.h"

#include "quantification_standard.h"

#include "stream_block.h"

#ifdef _BUILD_WITH_VTK
  #include "vtk_graph.h"
#endif


namespace io
{

// ----------------------------------------------------------------------------

struct file_name_size
{
    file_name_size(std::string name, int size) { filename = name; total_rank_size = size;}
    std::string filename;
    int total_rank_size;
};

struct file_name_fit_params
{
    std::string dataset_dir;
    std::string dataset_filename;
    int detector_num;
    data_struct::xrf::Fit_Parameters fit_params;
    data_struct::xrf::Spectra spectra;
    data_struct::xrf::Fit_Element_Map_Dict elements_to_fit;
    bool success;
};

// ----------------------------------------------------------------------------

void check_and_create_dirs(std::string dataset_directory);

bool compare_file_size (const file_name_size& first, const file_name_size& second);

std::vector<std::string> find_all_dataset_files(std::string dataset_directory, std::string search_str);

void generate_h5_averages(std::string dataset_directory,
                          std::string dataset_file,
                          ThreadPool* tp,
                          size_t detector_num_start,
                          size_t detector_num_end);

bool load_element_info(std::string element_henke_filename,
                       std::string element_csv_filename,
                       data_struct::xrf::Element_Info_Map *element_info_map);

bool load_and_integrate_spectra_volume(std::string dataset_directory,
                                       std::string dataset_file,
                                       data_struct::xrf::Spectra *integrated_spectra,
                                       size_t detector_num,
                                       data_struct::xrf::Params_Override * params_override);

bool load_override_params(std::string dataset_directory,
                          int detector_num,
                          data_struct::xrf::Params_Override *params_override);

bool load_quantification_standard(std::string dataset_directory,
                                  std::string quantification_info_file,
                                  std::string *standard_file_name,
                                  std::unordered_map<std::string, real_t> *element_standard_weights);

bool load_spectra_volume(std::string dataset_directory,
                         std::string dataset_file,
                         data_struct::xrf::Spectra_Volume *spectra_volume,
                         size_t detector_num,
                         data_struct::xrf::Params_Override * params_override,
                         data_struct::xrf::Quantification_Standard * quantification_standard,
                         bool save_scalers);

bool load_spectra_volume_with_callback(std::string dataset_directory,
                                       std::string dataset_file,
                                       size_t detector_num,
                                       data_struct::xrf::Params_Override * params_override,
                                       file::IO_Callback_Func_Def callback_fun,
                                       void* user_data);

void populate_netcdf_hdf5_files(std::string dataset_dir);

void save_averaged_fit_params(std::string dataset_dir,
                              std::vector<data_struct::xrf::Fit_Parameters> fit_params_avgs,
                              int detector_num_start,
                              int detector_num_end);

bool save_results(std::string save_loc,
                  const data_struct::xrf::Fit_Count_Dict * const element_counts,
                  fitting::routines::Base_Fit_Routine* fit_routine,
                  std::queue<std::future<bool> >* job_queue,
                  std::chrono::time_point<std::chrono::system_clock> start);

void save_optimized_fit_params(struct file_name_fit_params file_and_fit_params);

bool save_volume(data_struct::xrf::Quantification_Standard * quantification_standard,
                 data_struct::xrf::Spectra_Volume *spectra_volume,
                 real_t energy_offset,
                 real_t energy_slope,
                 real_t energy_quad);

void sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string> *dataset_files);

}// end namespace io

#endif // HL_FILE_IO_H
