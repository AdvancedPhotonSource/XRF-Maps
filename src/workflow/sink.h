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



#ifndef Sink_H
#define Sink_H

#include "core/defines.h"
#include <functional>
#include <future>
#include <memory>
#include <thread>
#include <atomic>
#include <mutex>
#include "workflow/distributor.h"

namespace workflow
{

//-----------------------------------------------------------------------------
template<typename T_IN>
class DLL_EXPORT Sink
{

public:

	
    Sink() : _running(false) {}

    Sink(const Sink&) = delete;

    Sink& operator=(const Sink&) = delete;

    virtual ~Sink()
    {
        stop();
    }

    template<typename _T>
    void connect(Distributor<_T, T_IN> *distributor)
    {
        _check_func = std::bind(&Distributor<_T, T_IN>::is_queue_empty, distributor);
        _get_func = std::bind(&Distributor<_T, T_IN>::front_chunk, distributor, std::placeholders::_1);
        //_get_func = std::bind(&Distributor<_T, T_IN>::front_pop, distributor);
    }

    virtual void set_function(std::function<void (T_IN)> func)
    {
        _callback_func = func;
    }

    void start()
    {
        if(_thread != nullptr)
        {
            stop();
        }
        std::packaged_task<void(void)> task([this](){ this->_execute(); });
        _running = true;
        _thread = std::make_unique<std::thread>(std::move(task));
    }

    void stop()
    {
        _running.store(false);
        // stop() is reachable without a prior start() (directly, and via
        // wait_and_stop()), so the worker may never have been created.
        if(_thread != nullptr && _thread->joinable())
        {
            _thread->join();
        }
        _thread.reset();
    }

    void wait_and_stop()
    {
        for (;;)
        {
            bool idle = false;
            {
                // Holding the local queue lock across the distributor check
                // prevents observing the distributor as empty in the middle of
                // a worker transfer into this queue.
                std::lock_guard<std::mutex> lock(_queue_mutex);
                idle = _job_queue.empty() && (!_check_func || _check_func());
            }
            if (idle)
            {
                break;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(400));
        }
        stop();
    }

    void sink_function(T_IN val)
    {
        _callback_func(val);
        // No manual cleanup: T_IN owns itself (shared_ptr), so the block is
        // released once this frame and any downstream copies go out of scope.
    }

protected:

    void _execute()
    {
        while(_running.load())
        {
            std::future<T_IN> next;
            bool have_job = false;
            {
                std::lock_guard<std::mutex> lock(_queue_mutex);
                if (_check_func && _get_func && !_check_func())
                {
                    _get_func(&_job_queue);
                }
                if (!_job_queue.empty())
                {
                    next = std::move(_job_queue.front());
                    _job_queue.pop();
                    have_job = true;
                }
            }

            if (have_job)
            {
                try
                {
                    T_IN input_block = next.get();
                    if (_callback_func)
                    {
                        _callback_func(input_block);
                    }
                    else
                    {
                        logE << "Sink callback is not configured; dropping input block\n";
                    }
                }
                catch (const std::exception& e)
                {
                    // Do not strand later futures in the queue if one fitting
                    // job or sink callback fails.
                    logE << "Sink worker failed to process input: " << e.what() << "\n";
                }
                catch (...)
                {
                    logE << "Sink worker failed to process input with an unknown error\n";
                }
            }
            else
            {
                std::this_thread::sleep_for(std::chrono::milliseconds(10));
            }
        }
    }


    std::function<bool (void)> _check_func;

    std::function<void (std::queue<std::future<T_IN> > *)> _get_func;
    //std::function<T_IN (void)> _get_func;

    std::function<void (T_IN)> _callback_func;

    std::queue<std::future<T_IN> > _job_queue;

    std::mutex _queue_mutex;

    std::atomic<bool> _running;

    // Owned: the worker thread started by start() and joined by stop()/~Sink().
    std::unique_ptr<std::thread> _thread;

};

} //namespace workflow

#endif // Sink_H
