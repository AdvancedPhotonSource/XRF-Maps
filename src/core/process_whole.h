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


#ifndef PROCESS_WHOLE
#define PROCESS_WHOLE

#include <iostream>
#include <queue>
#include <string>
#include <array>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <ctime>
#include <limits>
#include <sstream>
#include <fstream>

#include <stdlib.h>

#if defined _WIN32
#include "support/direct/dirent.h"
#else
#include <dirent.h>
#endif

#include "core/defines.h"

#include "workflow/threadpool.h"

#include "io/file/hl_file_io.h"
#include "io/file/mca_io.h"

#include "data_struct/spectra_volume.h"

#include "fitting/models/gaussian_model.h"

#include "data_struct/element_info.h"

#include "io/file/aps/aps_fit_params_import.h"

#include "fitting/routines/roi_fit_routine.h"
#include "fitting/routines/svd_fit_routine.h"
#include "fitting/routines/nnls_fit_routine.h"

#include "fitting/optimizers/lmfit_optimizer.h"
#include "fitting/optimizers/mpfit_optimizer.h"

#include "data_struct/fit_element_map.h"
#include "data_struct/params_override.h"

#include "data_struct/quantification_standard.h"


using namespace std::placeholders; //for _1, _2,

// ----------------------------------------------------------------------------

template<typename T>
data_struct::Fit_Count_Dict* generate_fit_count_dict(std::unordered_map<std::string, T> *elements_to_fit,
                                                          size_t height,
                                                          size_t width,
                                                          bool alloc_iter_count);

// ----------------------------------------------------------------------------

DLL_EXPORT bool fit_single_spectra(fitting::routines::Base_Fit_Routine * fit_routine,
                        const fitting::models::Base_Model * const model,
                        const data_struct::Spectra * const spectra,
                        const data_struct::Fit_Element_Map_Dict * const elements_to_fit,
                        data_struct::Fit_Count_Dict * out_fit_counts,
                        size_t i,
                        size_t j);

// ----------------------------------------------------------------------------

DLL_EXPORT struct io::file_name_fit_params* optimize_integrated_fit_params(std::string dataset_directory,
                                                                std::string  dataset_filename,
                                                                size_t detector_num,
                                                                fitting::models::Fit_Params_Preset,
                                                                fitting::optimizers::Optimizer *optimizer);

// ----------------------------------------------------------------------------

//DLL_EXPORT struct io::file_name_fit_params* optimize_integrated_fit_params_preloaded(std::string dataset_directory,
//                                                                                     std::string  dataset_filename,
//                                                                                     size_t detector_num,
//                                                                                     data_struct::Spectra * integrated_spectra,
//                                                                                     data_struct::Params_Override *params_override,
//                                                                                     fitting::models::Fit_Params_Preset optimize_fit_params_preset,
//                                                                                     fitting::optimizers::Optimizer *optimizer);

// ----------------------------------------------------------------------------

DLL_EXPORT void generate_optimal_params(data_struct::Analysis_Job* analysis_job);

// ----------------------------------------------------------------------------

//DLL_EXPORT void generate_optimal_params_mp(data_struct::Analysis_Job* analysis_job);

// ----------------------------------------------------------------------------

DLL_EXPORT void proc_spectra(data_struct::Spectra_Volume* spectra_volume,
                             data_struct::Detector* detector_struct,
                             ThreadPool* tp,
                             bool save_spec_vol);

// ----------------------------------------------------------------------------

DLL_EXPORT void process_dataset_files(data_struct::Analysis_Job* analysis_job);

// ----------------------------------------------------------------------------

void process_dataset_files_quick_and_dirty(std::string dataset_file, data_struct::Analysis_Job* analysis_job, ThreadPool& tp);

// ----------------------------------------------------------------------------

DLL_EXPORT bool perform_quantification(data_struct::Analysis_Job* analysis_job);

// ----------------------------------------------------------------------------

DLL_EXPORT void average_quantification(std::vector<data_struct::Quantification_Standard>* quant_stand_list,
										const std::vector<size_t>& detector_num_arr);

// ----------------------------------------------------------------------------

DLL_EXPORT void interate_datasets_and_update(data_struct::Analysis_Job& analysis_job);

// ----------------------------------------------------------------------------

#endif
