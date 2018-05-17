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

#include "core/process_streaming.h"

// ----------------------------------------------------------------------------

data_struct::Stream_Block* proc_spectra_block( data_struct::Stream_Block* stream_block )
{

    for(auto &itr : stream_block->fitting_blocks)
    {
        int i = itr.first;
        std::unordered_map<std::string, real_t> counts_dict = stream_block->fitting_blocks[i].fit_routine->fit_spectra(stream_block->model, stream_block->spectra, stream_block->elements_to_fit);
        //make count / sec
        for (auto& el_itr : *(stream_block->elements_to_fit))
        {
            stream_block->fitting_blocks[i].fit_counts[el_itr.first] = counts_dict[el_itr.first] / stream_block->spectra->elapsed_livetime();
        }
        stream_block->fitting_blocks[i].fit_counts[STR_NUM_ITR] = counts_dict[STR_NUM_ITR];
    }
    return stream_block;
}

// ----------------------------------------------------------------------------

void run_stream_pipeline(data_struct::Analysis_Job* job)
{
    workflow::Source<data_struct::Stream_Block*> *source;
    workflow::Distributor<data_struct::Stream_Block*, data_struct::Stream_Block*> distributor(job->num_threads);
    workflow::Sink<data_struct::Stream_Block*> *sink;

    //setup input
    if(job->quick_and_dirty)
    {
        source = new workflow::xrf::Detector_Sum_Spectra_Source(job);
    }
    else if(job->is_network_source)
    {
        source = new workflow::xrf::Spectra_Net_Source(job);
    }
    else
    {
        source = new workflow::xrf::Spectra_File_Source(job);
    }

    //setup output
    if(job->stream_over_network)
    {
        sink = new workflow::xrf::Spectra_Net_Streamer();
    }
    else
    {
        sink = new workflow::xrf::Spectra_Stream_Saver();
    }

    distributor.set_function(proc_spectra_block);
    source->connect(&distributor);
    sink->connect(&distributor);

    sink->start();
    source->run();
    sink->wait_and_stop();

    delete source;
    delete sink;
}

// ----------------------------------------------------------------------------

struct io::file_name_fit_params* optimize_integrated_fit_params( data_struct::Stream_Block* stream_block )
{
    struct io::file_name_fit_params* ret_struct = new struct io::file_name_fit_params();

    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = ret_struct->spectra.size() -1;

//    fitting::models::Range energy_range = get_energy_range(sub_struct->fit_params_override_dict.min_energy,
//                                                           sub_struct->fit_params_override_dict.min_energy,
//                                                           spectra_samples,
//                                                           sub_struct->fit_params_override_dict.fit_params[STR_ENERGY_OFFSET].value,
//                                                           sub_struct->fit_params_override_dict.fit_params[STR_ENERGY_SLOPE].value);


    fitting::models::Gaussian_Model model;
    //Fitting routines
    fitting::routines::Param_Optimized_Fit_Routine fit_routine;
    //TODO:
////    fit_routine.set_optimizer(optimizer);

    //reset model fit parameters to defaults
    model.reset_to_default_fit_params();
    //Update fit parameters by override values
    //TODO:
////    model.update_fit_params_values(params_override.fit_params);
    model.set_fit_params_preset(stream_block->optimize_fit_params_preset);
    //Initialize the fit routine
    fit_routine.initialize(&model, &ret_struct->elements_to_fit, energy_range);
    //Fit the spectra saving the element counts in element_fit_count_dict
    ret_struct->fit_params = fit_routine.fit_spectra_parameters(stream_block->model, stream_block->spectra, stream_block->elements_to_fit);

    ret_struct->success = true;

    delete stream_block;

    return ret_struct;
}

// ----------------------------------------------------------------------------

void save_optimal_params(struct io::file_name_fit_params* f_struct)
{
    static bool first = false;
    if(f_struct->success)
    {
        if(first)
        {
            //TODO:
///            fit_params_avgs[f_struct->detector_num] = f_struct->fit_params;
            first = false;
        }
        else
        {
            //TODO:
///            fit_params_avgs[f_struct->detector_num].moving_average_with(f_struct->fit_params);
        }
        io::save_optimized_fit_params(f_struct);
    }
    delete f_struct;
}

// ----------------------------------------------------------------------------

void run_optimization_stream_pipeline(data_struct::Analysis_Job* job)
{
    //TODO: Run only on 8 largest files if no files are specified
    workflow::xrf::Integrated_Spectra_Source spectra_stream_producer(job);
    workflow::Distributor<data_struct::Stream_Block*, struct io::file_name_fit_params*> distributor(job->num_threads);
    workflow::Sink<struct io::file_name_fit_params*> sink;
    sink.set_function(save_optimal_params);

    distributor.set_function(optimize_integrated_fit_params);
    sink.connect(&distributor);
    spectra_stream_producer.connect(&distributor);


    sink.start();
    spectra_stream_producer.run();
    sink.wait_and_stop();
    /*
    io::save_averaged_fit_params(dataset_directory, fit_params_avgs, detector_num_start, detector_num_end);
    */
}

// ----------------------------------------------------------------------------
