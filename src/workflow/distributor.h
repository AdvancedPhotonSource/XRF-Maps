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



#ifndef Distributor_H
#define Distributor_H

#include "defines.h"
#include "threadpool.h"
#include <functional>

namespace workflow
{

//-----------------------------------------------------------------------------
template <typename T_IN, typename T_OUT>
class DLL_EXPORT Distributor
{

public:

    Distributor(size_t num_threads)
    {
        _thread_pool = new ThreadPool(num_threads);
        _job_queue = nullptr;
        _callback_func = std::bind(&Distributor::distribute, this, std::placeholders::_1);
    }

    ~Distributor()
    {
        delete _thread_pool;
    }

    void connect(std::queue<std::future<T_OUT> > *job_queue)
    {
        _job_queue = job_queue;
    }

    void distribute(T_IN input)
    {
        if(_job_queue != nullptr)
        {
            _job_queue->emplace( _thread_pool->enqueue(_dist_func, input) );
        }
    }

    std::function<void (T_IN)> get_callback_func()
    {
        return _callback_func;
    }

    void set_function(std::function<T_OUT (T_IN)> dist_func)
    {
        _dist_func = dist_func;
    }

protected:

    std::function<void (T_IN)> _callback_func;

    std::function<T_OUT (T_IN)> _dist_func;

    ThreadPool *_thread_pool;

    std::queue<std::future<T_OUT> > *_job_queue;

};

} //namespace workflow

#endif // Distributor_H
