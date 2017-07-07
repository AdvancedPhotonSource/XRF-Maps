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



#ifndef Spectra_Net_Source_H
#define Spectra_Net_Source_H

#include "defines.h"

#include "source.h"
#include "stream_block.h"
#include "analysis_job.h"
#include "zmq.hpp"

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

class DLL_EXPORT Spectra_Net_Source : public Source<data_struct::xrf::Stream_Block*>
{

public:

    Spectra_Net_Source(data_struct::xrf::Analysis_Job* analysis_job);

    ~Spectra_Net_Source();

    //virtual void cb_load_spectra_data(size_t row, size_t col, size_t height, size_t width, size_t detector_num, data_struct::xrf::Spectra* spectra, void* user_data);

    virtual void run();

protected:
/*
    virtual bool _load_spectra_volume_with_callback(std::string dataset_directory,
                                                    std::string dataset_file,
                                                    size_t detector_num_start,
                                                    size_t detector_num_end,
                                                    io::file::IO_Callback_Func_Def callback_fun);
*/
    zmq::context_t _context;

    zmq::socket_t *_receiver;

    std::string *_current_dataset_directory;
    std::string *_current_dataset_name;

    data_struct::xrf::Analysis_Job* _analysis_job;

    bool _running;
//    std::function <void (size_t, size_t, size_t, size_t, size_t, data_struct::xrf::Spectra*, void*)> _cb_function;

};

} //namespace xrf
} //namespace workflow

#endif // Spectra_Net_Source_H
