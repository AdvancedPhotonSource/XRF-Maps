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
#include "data_struct/xrf/element_info.h"
#include "fitting/routines/base_fit_routine.h"
#include <vector>
#include <string>
#include "data_struct/xrf/quantification_standard.h"
#include "data_struct/xrf/params_override.h"
#include "fitting/optimizers/lmfit_optimizer.h"
#include "fitting/optimizers/mpfit_optimizer.h"
#include <iostream>

namespace data_struct
{
namespace xrf
{

//-----------------------------------------------------------------------------

///
/// \brief The Analysis_Sub_Struct class
///
struct DLL_EXPORT Analysis_Sub_Struct
{

    std::unordered_map<int, fitting::routines::Base_Fit_Routine *> fit_routines;

    fitting::models::Base_Model * model;

    data_struct::xrf::Quantification_Standard quant_standard;

    data_struct::xrf::Params_Override fit_params_override_dict;

};

///
/// \brief The Analysis_Job class
///
class DLL_EXPORT Analysis_Job
{
public:

    Analysis_Job(size_t num_threads=1);

    ~Analysis_Job();

    const std::string& dataset_directory() { return _dataset_directory; }

    const std::vector<std::string>& dataset_files() { return _dataset_files; }

    struct Analysis_Sub_Struct* get_sub_struct(int detector_num);

    //truct Analysis_Sub_Struct* get_default_sub_struct() { return &_default_sub_struct; }
/*
    inline auto get_detector_begin() { return _detectors_meta_data.begin(); }

    inline auto get_detector_end() { return _detectors_meta_data.begin(); }
*/
    bool load(std::string dataset_directory,
              std::vector<std::string> dataset_files,
              std::vector<Fitting_Routines> fitting_routines,
              size_t detector_num_start,
              size_t detector_num_end);

    const size_t& num_threads() { return _num_threads; }

    size_t get_num_detectors() { return _detectors_meta_data.size(); }

    void set_optimizer(std::string optimizer);

    const size_t& detector_num_start() { return _detector_num_start; }

    const size_t& detector_num_end() { return _detector_num_end; }

    void fit_params_preset(fitting::models::Fit_Params_Preset preset) { _optimize_fit_params_preset = preset;}

    const fitting::models::Fit_Params_Preset& fit_params_preset() {return _optimize_fit_params_preset;}

protected:

    fitting::routines::Base_Fit_Routine* _generate_fit_routine(Fitting_Routines proc_type);

    std::string _dataset_directory;

    std::vector<std::string> _dataset_files;

    std::vector<Fitting_Routines> _fitting_routines;

    //struct Analysis_Sub_Struct _default_sub_struct;
    std::map<int, struct Analysis_Sub_Struct> _detectors_meta_data;

    fitting::models::Fit_Params_Preset _optimize_fit_params_preset;

    //Optimizers for fitting models
    fitting::optimizers::LMFit_Optimizer _lmfit_optimizer;
    fitting::optimizers::MPFit_Optimizer _mpfit_optimizer;
    fitting::optimizers::Optimizer *_optimizer;

    size_t _detector_num_start;
    size_t _detector_num_end;

    size_t _num_threads;
};

} //namespace xrf

} //namespace data_struct

#endif // Analysis_Job_H
