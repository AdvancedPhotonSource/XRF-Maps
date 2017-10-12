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

#ifndef PROCESS_STREAMING_H
#define PROCESS_STREAMING_H

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

#include "core/defines.h"

#include "io/file/hl_file_io.h"
#include "data_struct/xrf/element_info.h"
#include "data_struct/xrf/fit_element_map.h"
#include "data_struct/xrf/analysis_job.h"
#include "core/command_line_parser.h"
#include "data_struct/xrf/stream_block.h"
#include "workflow/pipeline.h"
#include "workflow/xrf/spectra_file_source.h"
#include "workflow/xrf/integrated_spectra_source.h"
#include "workflow/xrf/detector_sum_spectra_source.h"
#include "workflow/xrf/spectra_net_streamer.h"
#include "workflow/xrf/spectra_stream_saver.h"

#include "fitting/routines/param_optimized_fit_routine.h"
#include "fitting/models/gaussian_model.h"

DLL_EXPORT data_struct::xrf::Stream_Block* proc_spectra_block( data_struct::xrf::Stream_Block* stream_block );

DLL_EXPORT void run_stream_pipeline(data_struct::xrf::Analysis_Job* job);

DLL_EXPORT struct io::file_name_fit_params* optimize_integrated_fit_params( data_struct::xrf::Stream_Block* stream_block );

DLL_EXPORT void save_optimal_params(struct io::file_name_fit_params* f_struct);

DLL_EXPORT void run_optimization_stream_pipeline(data_struct::xrf::Analysis_Job* job);

DLL_EXPORT void run_quick_n_dirty_pipeline(data_struct::xrf::Analysis_Job* job);

DLL_EXPORT bool perform_quantification_streaming(std::string dataset_directory,
                            std::string quantification_info_file,
                            std::vector<data_struct::xrf::Fitting_Routines> proc_types,
                            std::vector<data_struct::xrf::Quantification_Standard>* quant_stand_list,
                            std::unordered_map<int, data_struct::xrf::Params_Override> *fit_params_override_dict,
                            size_t detector_num_start,
                            size_t detector_num_end);

#endif
