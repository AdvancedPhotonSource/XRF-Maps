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



#include "spectra_net_streamer.h"

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

template<typename T_real>
Spectra_Net_Streamer<T_real>::Spectra_Net_Streamer(std::string port) : Sink<data_struct::Stream_Block<T_real>*>()
{
#ifdef _BUILD_WITH_ZMQ
    _send_counts = true;
    _send_spectra = true;
    _verbose = false;
    this->_callback_func = std::bind(&Spectra_Net_Streamer<T_real>::stream, this, std::placeholders::_1);

    std::string conn_str = "tcp://*:" + port;
	_context = new zmq::context_t(1);
	_zmq_socket = new zmq::socket_t(*_context, ZMQ_PUB);
	_zmq_socket->bind(conn_str);
#else
    logE<<"Spectra_Net_Streamer needs ZeroMQ to work. Recompile with option -DBUILD_WITH_ZMQ\n";
#endif
}

//-----------------------------------------------------------------------------

template<typename T_real>
Spectra_Net_Streamer<T_real>::~Spectra_Net_Streamer()
{
#ifdef _BUILD_WITH_ZMQ
    if(_zmq_socket != nullptr)
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
void Spectra_Net_Streamer<T_real>::stream(data_struct::Stream_Block<T_real>* stream_block)
{
#ifdef _BUILD_WITH_ZMQ
	std::string data;

    if(_send_counts && _send_spectra)
    {
        zmq::message_t topic("XRF-Counts-and-Spectra", 22);
        _zmq_socket->send(topic, zmq::send_flags::sndmore);
        data = _serializer.encode_counts_and_spectra(stream_block);
        zmq::message_t message(data.c_str(), data.length());
        if (false == _zmq_socket->send(message, zmq::send_flags::none))
        {
            logE << "sending ZMQ counts and spectra message"<<"\n";
        }
    }
    else
    {
        if(_send_counts)
        {
            if(_verbose)
            {
                logI<<"Sending counts "<< stream_block->dataset_name <<" "<<stream_block->row()<<" " <<stream_block->col()<<"\n";
            }
            zmq::message_t topic("XRF-Counts", 10);
            _zmq_socket->send(topic, zmq::send_flags::sndmore);
            data = _serializer.encode_counts(stream_block);
            zmq::message_t message(data.c_str(), data.length());
            if (false == _zmq_socket->send(message, zmq::send_flags::none))
            {
                logE << "sending ZMQ counts message"<<"\n";
            }
        }
        if(_send_spectra)
        {
            if(_verbose)
            {
                logI<<"Sending spectra "<< stream_block->dataset_name <<" "<<stream_block->row()<<" "<< stream_block->col()<<"\n";
            }
            zmq::message_t topic("XRF-Spectra", 11);
            _zmq_socket->send(topic, zmq::send_flags::sndmore);
            data = _serializer.encode_spectra(stream_block);
            zmq::message_t message(data.c_str(), data.length());
            if (false == _zmq_socket->send(message, zmq::send_flags::none))
            {
                logE << "sending ZMQ spectra message"<<"\n";
            }
        }
    }
#else
    logE<<"Spectra_Net_Streamer needs ZeroMQ to work. Recompile with option -DBUILD_WITH_ZMQ\n";
#endif
}

// ----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT Spectra_Net_Streamer<float>;
TEMPLATE_CLASS_DLL_EXPORT Spectra_Net_Streamer<double>;

} //namespace xrf
} //namespace workflow
