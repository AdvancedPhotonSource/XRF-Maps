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



#ifndef BASIC_SERIALIZER_H
#define BASIC_SERIALIZER_H

#include "core/defines.h"
#include "data_struct/stream_block.h"

namespace io
{
namespace net
{

class DLL_EXPORT Basic_Serializer
{
public:

    Basic_Serializer();

    ~Basic_Serializer();

    std::string encode_counts(data_struct::Stream_Block* in_stream_block);

    data_struct::Stream_Block* decode_counts(char* message, size_t message_len);

    std::string encode_spectra(data_struct::Stream_Block* in_stream_block);

    data_struct::Stream_Block* decode_spectra(char* message, size_t message_len);

    std::string encode_counts_and_spectra(data_struct::Stream_Block* in_stream_block);

    data_struct::Stream_Block* decode_counts_and_spectra(char* message, size_t message_len);

protected:
	template <typename T>
	void _convert_var_to_bytes(std::string * str, char* bytes_temp, T variable, size_t size)
	{
		memcpy(bytes_temp, (char*)(&variable), size);
        for (size_t i = 0; i < size; i++)
			*str += bytes_temp[i];
	}

    void _encode_meta(data_struct::Stream_Block* stream_block, std::string& raw_msg);

    void _encode_counts(data_struct::Stream_Block* stream_block, std::string& raw_msg);

    void _encode_spectra(data_struct::Stream_Block* stream_block, std::string& raw_msg);

    data_struct::Stream_Block* _decode_meta(char* message, size_t message_len, size_t& idx);

    void _decode_counts(char* message, size_t message_len, size_t& idx, data_struct::Stream_Block* out_stream_block);

    void _decode_spectra(char* message, size_t message_len, size_t& idx, data_struct::Stream_Block* out_stream_block);
};

}// end namespace net
}// end namespace io

#endif // Basic_Serializer_H
