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



#include "analysis_job.h"

#include "fitting/routines/roi_fit_routine.h"
#include "fitting/routines/svd_fit_routine.h"
#include "fitting/routines/nnls_fit_routine.h"

#include "fitting/models/gaussian_model.h"

#include "io/hl_file_io.h"

namespace data_struct
{
namespace xrf
{

//-----------------------------------------------------------------------------

Analysis_Job::Analysis_Job(size_t num_threads)
{
    _optimizer = &_lmfit_optimizer;
    _num_threads = num_threads;
    _detector_num_start = 0;
    _detector_num_end = 0;
    //default mode for which parameters to fit when optimizing fit parameters
    _optimize_fit_params_preset = fitting::models::BATCH_FIT_NO_TAILS;
    _quick_and_dirty = false;
    _optimize_fit_override_params = false;
    _generate_average_h5 = false;
    _command_line = "";
}

//-----------------------------------------------------------------------------

Analysis_Job::~Analysis_Job()
{

}

//-----------------------------------------------------------------------------

fitting::routines::Base_Fit_Routine* Analysis_Job::_generate_fit_routine(Fitting_Routines proc_type)
{
    //Fitting routines
    fitting::routines::Base_Fit_Routine *fit_routine = nullptr;
    switch(proc_type)
    {
        case GAUSS_TAILS:
            fit_routine = new fitting::routines::Param_Optimized_Fit_Routine();
            ((fitting::routines::Param_Optimized_Fit_Routine*)fit_routine)->set_optimizer(_optimizer);
            break;
        case GAUSS_MATRIX:
            fit_routine = new fitting::routines::Matrix_Optimized_Fit_Routine();
            ((fitting::routines::Matrix_Optimized_Fit_Routine*)fit_routine)->set_optimizer(_optimizer);
            break;
        case ROI:
            fit_routine = new fitting::routines::ROI_Fit_Routine();
            break;
        case SVD:
            fit_routine = new fitting::routines::SVD_Fit_Routine();
            break;
        case NNLS:
            fit_routine = new fitting::routines::NNLS_Fit_Routine();
            break;
        default:
            break;
    }
    return fit_routine;
}

//-----------------------------------------------------------------------------

struct Analysis_Sub_Struct* Analysis_Job::get_sub_struct(int detector_num)
{
       struct Analysis_Sub_Struct* sub_struct = nullptr;
       if(_detectors_meta_data.count(detector_num) > 0)
       {
            sub_struct = &(_detectors_meta_data.at(detector_num));
       }
       return sub_struct;
}

//-----------------------------------------------------------------------------

bool Analysis_Job::init(size_t detector_num_start, size_t detector_num_end)
{

    _detector_num_start = detector_num_start;
    _detector_num_end = detector_num_end;

    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = 2000;

    /*
    _default_sub_struct.
    if( false == io::load_override_params(_dataset_directory, -1, override_params) )
    {
        return false;
    }
    */

    _detectors_meta_data.clear();

    //initialize models and fit routines for all detectors
    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        _detectors_meta_data[detector_num] = Analysis_Sub_Struct();
        Analysis_Sub_Struct *sub_struct = &_detectors_meta_data[detector_num];

        sub_struct->model = new fitting::models::Gaussian_Model();
        data_struct::xrf::Params_Override * override_params = &(sub_struct->fit_params_override_dict);

        override_params->dataset_directory = _dataset_directory;
        override_params->detector_num = detector_num;

        if( false == io::load_override_params(_dataset_directory, detector_num, override_params) )
        {
            if( false == io::load_override_params(_dataset_directory, -1, override_params) )
            {
                return false;
            }
        }

        if (override_params->elements_to_fit.size() < 1)
        {
            logit<<"Error, no elements to fit. Check  maps_fit_parameters_override.txt0 - 3 exist"<<std::endl;
            return false;
        }

        for(auto proc_type : _fitting_routines)
        {
            //Fitting models
            fitting::routines::Base_Fit_Routine *fit_routine = _generate_fit_routine(proc_type);
            sub_struct->fit_routines[proc_type] = fit_routine;

            logit << "Generating model for "<< fit_routine->get_name() <<" detector "<<detector_num<<std::endl;


            //reset model fit parameters to defaults
            sub_struct->model->reset_to_default_fit_params();
            //Update fit parameters by override values
            sub_struct->model->update_fit_params_values(override_params->fit_params);

            //Initialize model
            fit_routine->initialize(sub_struct->model, &override_params->elements_to_fit, energy_range);
        }
    }

    return true;
}

//-----------------------------------------------------------------------------

void Analysis_Job::set_optimizer(std::string optimizer)
{
    if(optimizer == "mpfit")
    {
        _optimizer = &_mpfit_optimizer;
    }
    else
    {
        _optimizer = &_lmfit_optimizer;
    }
}

//-----------------------------------------------------------------------------

} //namespace xrf
} //namespace data_struct
