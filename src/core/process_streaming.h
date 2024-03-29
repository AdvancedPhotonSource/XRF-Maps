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

#include "core/defines.h"
#include "data_struct/analysis_job.h"
#include "data_struct/stream_block.h"
#include "workflow/xrf/detector_sum_spectra_source.h"
#include "workflow/xrf/integrated_spectra_source.h"
#include "workflow/xrf/spectra_file_source.h"
#include "workflow/xrf/spectra_net_source.h"
#include "workflow/xrf/spectra_net_streamer.h"
#include "workflow/xrf/spectra_stream_saver.h"

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT data_struct::Stream_Block<T_real>* proc_spectra_block( data_struct::Stream_Block<T_real>* stream_block )
{

    for (auto& itr : stream_block->fitting_blocks)
    {
        std::unordered_map<std::string, T_real> counts_dict;
        stream_block->fitting_blocks[itr.first].fit_routine->fit_spectra(stream_block->model, stream_block->spectra, stream_block->elements_to_fit, counts_dict);
        //make count / sec
        for (auto& el_itr : *(stream_block->elements_to_fit))
        {
            stream_block->fitting_blocks[itr.first].fit_counts[el_itr.first] = counts_dict[el_itr.first] / stream_block->spectra->elapsed_livetime();
        }
        stream_block->fitting_blocks[itr.first].fit_counts[STR_NUM_ITR] = counts_dict[STR_NUM_ITR];
    }
    return stream_block;
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT void run_stream_pipeline(data_struct::Analysis_Job<T_real>* job)
{
    workflow::Source<data_struct::Stream_Block<T_real>*>* source;
    workflow::Distributor<data_struct::Stream_Block<T_real>*, data_struct::Stream_Block<T_real>*> distributor(job->num_threads);
    workflow::Sink<data_struct::Stream_Block<T_real>*>* sink;

    //setup input
    if (job->quick_and_dirty)
    {
        source = new workflow::xrf::Detector_Sum_Spectra_Source<T_real>(job);
    }
    else if (job->is_network_source)
    {
        if (job->network_source_ip.length() > 0)
        {
            source = new workflow::xrf::Spectra_Net_Source<T_real>(job, job->network_source_ip, job->network_source_port);
        }
        else
        {
            source = new workflow::xrf::Spectra_Net_Source<T_real>(job);
        }
    }
    else
    {
        source = new workflow::xrf::Spectra_File_Source<T_real>(job);
    }

    //setup output
    if (job->stream_over_network)
    {
        sink = new workflow::xrf::Spectra_Net_Streamer<T_real>(job->network_stream_port);
    }
    else
    {
        sink = new workflow::xrf::Spectra_Stream_Saver<T_real>();
    }

    distributor.set_function(proc_spectra_block<T_real>);
    source->connect(&distributor);
    sink->connect(&distributor);


    sink->start();
    source->run();
    sink->wait_and_stop();

    delete source;
    delete sink;
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT void stream_spectra(data_struct::Analysis_Job<T_real>* job)
{
    workflow::Source<data_struct::Stream_Block<T_real>*>* source;
    workflow::xrf::Spectra_Net_Streamer<T_real> sink(job->network_stream_port);
    sink.set_send_counts(false);
    sink.set_send_spectra(true);

    //setup input
    if (job->quick_and_dirty)
    {
        source = new workflow::xrf::Detector_Sum_Spectra_Source<T_real>(job);
    }
    else
    {
        source = new workflow::xrf::Spectra_File_Source<T_real>(job);
    }

    source->connect(&sink);
    source->run();

    delete source;
}

// ----------------------------------------------------------------------------

#endif
