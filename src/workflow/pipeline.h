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



#ifndef Pipeline_H
#define Pipeline_H

#include "core/defines.h"
#include "workflow/source.h"
#include "workflow/distributor.h"
#include "workflow/sink.h"

namespace workflow
{

//-----------------------------------------------------------------------------
template <typename _T>
class DLL_EXPORT Simple_Pipeline
{

public:

    Simple_Pipeline(int num_threads)
    {
        _distributor = new Distributor<_T, _T>(num_threads);

        _producer.connect(_distributor->get_callback_func());
        _distributor->connect(_sink.get_job_queue());
    }

    ~Simple_Pipeline()
    {
        delete _distributor;
    }

    void set_producer_func(std::function<_T (void*)> func)
    {
        _producer.set_function(func);
    }

    void set_distributor_func(std::function<_T (_T)> func)
    {
        _distributor->set_function(func);
    }

    void set_sink_func(std::function<void (_T)> func)
    {
        _sink.set_function(func);
    }

    void run()
    {
        _sink.start();
        _producer.run();
        _sink.wait_and_stop();
    }

protected:

    Source<_T> _producer;
    Distributor<_T, _T> *_distributor;
    Sink<_T> _sink;

};

} //namespace workflow

#endif // Pipeline_H
