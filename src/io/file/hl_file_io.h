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

#if defined _WIN32
#include "support/direct/dirent.h"
#else
#include <dirent.h>
#endif

#include "io/file/netcdf_io.h"
#include "io/file/mda_io.h"
#include "io/file/hdf5_io.h"
#include "io/file/csv_io.h"

#include "data_struct/spectra_volume.h"

#include "fitting/models/gaussian_model.h"

#include "data_struct/element_info.h"

#include "io/file/aps/aps_fit_params_import.h"

#include "fitting/routines/base_fit_routine.h"

#include "data_struct/fit_element_map.h"
#include "data_struct/params_override.h"

#include "data_struct/quantification_standard.h"

#include "data_struct/stream_block.h"

#include "fitting/routines/roi_fit_routine.h"
#include "fitting/routines/svd_fit_routine.h"
#include "fitting/routines/nnls_fit_routine.h"

#include "fitting/models/gaussian_model.h"

#ifdef _BUILD_WITH_QT
  #include "visual/grapher.h"
#endif


namespace io
{

// ----------------------------------------------------------------------------

struct DLL_EXPORT file_name_size
{
    file_name_size(std::string name, int size) { filename = name; total_rank_size = size;}
    std::string filename;
    int total_rank_size;
};

struct DLL_EXPORT file_name_fit_params
{
    std::string dataset_dir;
    std::string dataset_filename;
    int detector_num;
    data_struct::Fit_Parameters fit_params;
    data_struct::Spectra spectra;
    data_struct::Fit_Element_Map_Dict elements_to_fit;
    bool success;
};

// ----------------------------------------------------------------------------

DLL_EXPORT void check_and_create_dirs(std::string dataset_directory);

DLL_EXPORT bool compare_file_size (const file_name_size& first, const file_name_size& second);

DLL_EXPORT std::vector<std::string> find_all_dataset_files(std::string dataset_directory, std::string search_str);

DLL_EXPORT void generate_h5_averages(std::string dataset_directory,
                          std::string dataset_file,
                          size_t detector_num_start,
                          size_t detector_num_end);

DLL_EXPORT fitting::routines::Base_Fit_Routine* generate_fit_routine(data_struct::Fitting_Routines proc_type,
                                                                     fitting::optimizers::Optimizer* optimizer);

/**
 * @brief init_analysis_job_detectors : Read in maps_fit_parameters_override.txt[0-3] and initialize data structres
 *                                      Override file contains information about which element to fit, updated branching
 *                                      conditions and other custom properties of the dataset.
 * @param analysis_job : data structure that holds information about the analysis to be perfomred.
 * @return True if successful
 */
DLL_EXPORT bool init_analysis_job_detectors(data_struct::Analysis_Job* analysis_job);

DLL_EXPORT bool load_element_info(std::string element_henke_filename,
                       std::string element_csv_filename);

DLL_EXPORT bool load_and_integrate_spectra_volume(std::string dataset_directory,
                                       std::string dataset_file,
                                       data_struct::Spectra *integrated_spectra,
                                       size_t detector_num,
                                       data_struct::Params_Override * params_override);

DLL_EXPORT bool load_override_params(std::string dataset_directory,
                          int detector_num,
                          data_struct::Params_Override *params_override);

DLL_EXPORT bool load_quantification_standard(std::string dataset_directory,
                                  std::string quantification_info_file,
                                  std::string *standard_file_name,
                                  std::unordered_map<std::string, real_t> *element_standard_weights);

DLL_EXPORT bool load_spectra_volume(std::string dataset_directory,
                         std::string dataset_file,
                         data_struct::Spectra_Volume *spectra_volume,
                         size_t detector_num,
                         data_struct::Params_Override * params_override,
                         data_struct::Quantification_Standard * quantification_standard,
                         bool *is_loaded_from_analyazed_h5,
                         bool save_scalers);

DLL_EXPORT void populate_netcdf_hdf5_files(std::string dataset_dir);

DLL_EXPORT void save_averaged_fit_params(std::string dataset_dir,
                              std::unordered_map<int, data_struct::Fit_Parameters> fit_params_avgs,
                              size_t detector_num_start,
                              size_t detector_num_end);

DLL_EXPORT bool save_results(std::string save_loc,
                  const data_struct::Fit_Count_Dict * const element_counts,
                  std::queue<std::future<bool> >* job_queue,
                  std::chrono::time_point<std::chrono::system_clock> start);

DLL_EXPORT void save_optimized_fit_params(struct file_name_fit_params* file_and_fit_params);

DLL_EXPORT bool save_volume(data_struct::Quantification_Standard * quantification_standard,
                 data_struct::Spectra_Volume *spectra_volume,
                 real_t energy_offset,
                 real_t energy_slope,
                 real_t energy_quad);

DLL_EXPORT void sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string> *dataset_files);

}// end namespace io

#endif // HL_FILE_IO_H
