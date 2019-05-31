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
#include <thread>
#include "workflow/distributor.h"

namespace workflow
{

//-----------------------------------------------------------------------------
template<typename T_IN>
class DLL_EXPORT Sink
{

public:

	
    Sink()
    {
        _thread = nullptr;
        _running = false;
        _delete_block = true;
    }

	Sink(const Sink &)
	{

	}

	Sink& operator=(const Sink&)
	{
		return *this;
	}

    virtual ~Sink()
    {
        if(_thread != nullptr)
        {
            _thread->join();
            delete _thread;
        }
        _thread = nullptr;
    }

    void set_delete_block(bool val) { _delete_block = val; }

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
        _thread = new std::thread(std::move(task));
    }

    void stop()
    {
        _running = false;
        _thread->join();
        delete _thread;
        _thread = nullptr;
    }

    void wait_and_stop()
    {
        while( _check_func() == false || !_job_queue.empty())
        {
            std::this_thread::sleep_for(std::chrono::milliseconds(400));
        }
        stop();
    }

    void sink_function(T_IN val)
    {
        _callback_func(val);
		// if sink thread is not running we have to delete the stream_block
		if (_delete_block && _running == false)
		{
			delete val;
		}
    }

protected:

    void _execute()
    {
        while(_running)
        {
            if( _check_func() == false)
            {
                _get_func(&_job_queue);
                while(! _job_queue.empty())
                {
                    auto ret = std::move(_job_queue.front());
                    _job_queue.pop();
                    T_IN input_block = ret.get();

                    _callback_func(input_block);

                    if(_delete_block && input_block != nullptr)
                    {
                        delete input_block;
						input_block = nullptr;
                    }
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

    bool _running;

    std::thread *_thread;

    bool _delete_block;

};

} //namespace workflow

#endif // Sink_H
