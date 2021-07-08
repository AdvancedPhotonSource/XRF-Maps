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



#ifndef Analysis_Job_H
#define Analysis_Job_H

#include "core/defines.h"
#include "data_struct/detector.h"
#include "data_struct/element_info.h"
#include "fitting/routines/base_fit_routine.h"
#include <vector>
#include <string>
#include <thread>
#include "data_struct/quantification_standard.h"
#include "data_struct/params_override.h"
#include "fitting/optimizers/lmfit_optimizer.h"
#include "fitting/optimizers/mpfit_optimizer.h"
#include <iostream>

namespace data_struct
{

// ----------------------------------------------------------------------------

///
/// \brief The Analysis_Job class
///
class DLL_EXPORT Analysis_Job
{
public:

    Analysis_Job();

    ~Analysis_Job();

    Detector* get_first_detector();

    Detector* get_detector(int detector_num);

    void set_optimizer(std::string optimizer);

    fitting::optimizers::Optimizer *optimizer(){return _optimizer;}

    void init_fit_routines(size_t spectra_samples, bool force=false);

    std::string command_line;

    std::string dataset_directory;

    std::string quantification_standard_filename;

    std::string theta_pv;

    std::string network_source_ip;

    std::string network_source_port;

    std::string network_stream_port;

    float theta;

    std::vector<std::string> dataset_files;

    std::vector<std::string> optimize_dataset_files;

    std::vector<Fitting_Routines> fitting_routines;

    std::map<int, Detector> detectors_meta_data;

    fitting::models::Fit_Params_Preset optimize_fit_params_preset;

	std::string update_theta_str;

	std::vector<size_t> detector_num_arr;

    size_t num_threads;

    bool update_scalers;

    bool quick_and_dirty;

    bool generate_average_h5;

    bool add_v9_layout;

    bool add_exchange_layout;

    bool is_network_source;

    bool stream_over_network;

    bool export_int_fitted_to_csv;

	long long mem_limit;

	std::string update_us_amps_str;

	std::string update_ds_amps_str;

	std::string update_quant_us_amps_str;

	std::string update_quant_ds_amps_str;

    //list of quantification standards to use
    vector<Quantification_Standard> standard_element_weights;

protected:

    //Optimizers for fitting models
    fitting::optimizers::LMFit_Optimizer _lmfit_optimizer;
    fitting::optimizers::MPFit_Optimizer _mpfit_optimizer;
    fitting::optimizers::Optimizer *_optimizer;

    size_t _last_init_sample_size;

private:
	bool _first_init;

};

} //namespace data_struct

#endif // Analysis_Job_H
