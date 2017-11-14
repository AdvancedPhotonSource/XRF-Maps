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

namespace data_struct
{
namespace xrf
{

//-----------------------------------------------------------------------------

Analysis_Job::Analysis_Job()
{
    _optimizer = &_lmfit_optimizer;
    _last_init_sample_size = 0;
	_first_init = true;
    num_threads = 1;
    detector_num_start = 0;
    detector_num_end = 0;
    //default mode for which parameters to fit when optimizing fit parameters
    optimize_fit_params_preset = fitting::models::BATCH_FIT_NO_TAILS;
    quick_and_dirty = false;
    optimize_fit_override_params = false;
    generate_average_h5 = false;
    is_network_source = false;
    stream_over_network = false;
    command_line = "";
    theta_pv = "";
    theta = 0.f;
}

//-----------------------------------------------------------------------------

Analysis_Job::~Analysis_Job()
{

    for(auto& itr : detectors_meta_data)
    {
        Analysis_Sub_Struct *sub = &itr.second;
        if (sub->model != nullptr)
        {
            delete sub->model;
            sub->model = nullptr;
        }
        for(auto &itr2 : sub->fit_routines)
        {
            fitting::routines::Base_Fit_Routine *fit_routine = itr2.second;
            if(fit_routine != nullptr)
            {
                delete fit_routine;
                itr2.second = nullptr;
            }
        }
        sub->fit_routines.clear();
    }
    detectors_meta_data.clear();

}

//-----------------------------------------------------------------------------

struct Analysis_Sub_Struct* Analysis_Job::get_sub_struct(int detector_num)
{
       struct Analysis_Sub_Struct* sub_struct = nullptr;
       if(detectors_meta_data.count(detector_num) > 0)
       {
            sub_struct = &(detectors_meta_data.at(detector_num));
       }
       return sub_struct;
}

//-----------------------------------------------------------------------------

void Analysis_Job::init_fit_routines(size_t spectra_samples)
{
    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = spectra_samples - 1;

    if(_first_init && _last_init_sample_size != spectra_samples)
    {
		_first_init = false;
        _last_init_sample_size = spectra_samples;
        for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
        {
            data_struct::xrf::Analysis_Sub_Struct *sub_struct = &detectors_meta_data[detector_num];

            for(auto &proc_type : fitting_routines)
            {
                //Fitting models
                fitting::routines::Base_Fit_Routine *fit_routine = sub_struct->fit_routines[proc_type];
                logit << "Updating fit routine "<< fit_routine->get_name() <<" detector "<<detector_num<<std::endl;

                Fit_Element_Map_Dict *elements_to_fit = &(sub_struct->fit_params_override_dict.elements_to_fit);
                //Initialize model
                fit_routine->initialize(sub_struct->model, elements_to_fit, energy_range);
            }
        }
    }
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
