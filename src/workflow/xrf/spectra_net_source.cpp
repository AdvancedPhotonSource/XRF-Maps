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



#include "spectra_net_source.h"

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

template<typename T_real>
Spectra_Net_Source<T_real>::Spectra_Net_Source(data_struct::Analysis_Job<T_real>* analysis_job, std::string ip_addr, std::string port) : Source<data_struct::Stream_Block<T_real>*>()
{
    _analysis_job = analysis_job;
#ifdef _BUILD_WITH_ZMQ
    _conn_str = "tcp://"+ip_addr+":"+port;
    logI<<"Connecting to "<<_conn_str<<"\n";
	_context = new zmq::context_t(1);
	_zmq_socket = new zmq::socket_t(*_context, ZMQ_SUB);
    _zmq_socket->connect(_conn_str);
    _zmq_socket->set(zmq::sockopt::subscribe, "XRF-Spectra");
    //_zmq_socket->setsockopt(ZMQ_RCVTIMEO, 1000); //set timeout to 1000ms
#else
    logE<<"Spectra_Net_Source needs ZeroMQ to work. Recompile with option -DBUILD_WITH_ZMQ\n";
#endif

}

//-----------------------------------------------------------------------------

template<typename T_real>
Spectra_Net_Source<T_real>::~Spectra_Net_Source()
{
#ifdef _BUILD_WITH_ZMQ
	if (_zmq_socket != nullptr)
	{
		_zmq_socket->close();
		delete _zmq_socket;
	}
	if (_context != nullptr)
	{
		_context->close();
		delete _context;
	}
	_zmq_socket = nullptr;
	_context = nullptr;
#endif
}

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Net_Source<T_real>::run()
{
#ifdef _BUILD_WITH_ZMQ
    _running = true;
    zmq::message_t token, message;
    data_struct::Stream_Block<T_real> *stream_block;
    while (_running)
    {
        zmq::recv_result_t result = _zmq_socket->recv(token, zmq::recv_flags::none);
        if (result.has_value())
        {
            std::string s1 (static_cast<char*>(token.data()), token.size()); 
            if(s1 == "XRF-Spectra")
            {
                zmq::recv_result_t result2 = _zmq_socket->recv(message, zmq::recv_flags::none);
                if (result2.has_value())
                {
                    if(this->_output_callback_func != nullptr && _analysis_job != nullptr)
                    {
                        stream_block = _serializer.decode_spectra((char*)message.data(), message.size());
                        _analysis_job->init_fit_routines(stream_block->spectra->size());
                        data_struct::Detector<T_real>* cp = _analysis_job->get_detector(stream_block->detector_number());

                        if(cp == nullptr)
                        {
                            cp = _analysis_job->get_first_detector();
                        }
                        if(cp != nullptr)
                        {
                            stream_block->init_fitting_blocks(&(cp->fit_routines), &(cp->fit_params_override_dict.elements_to_fit));
                            stream_block->model = cp->model;
                            logI<<"Detector number "<<stream_block->detector_number()<<" is "<<cp->get_name()<<"\n";
                        }

                        this->_output_callback_func(stream_block);
                    }
                }
            }
        }
    }
    _zmq_socket->close();
    logI<<"ZMQ socket closed\n";
#endif
}

// ----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT Spectra_Net_Source<float>;
TEMPLATE_CLASS_DLL_EXPORT Spectra_Net_Source<double>;

} //namespace xrf
} //namespace workflow
