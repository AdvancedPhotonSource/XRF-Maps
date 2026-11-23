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


#include "io/net/basic_serializer.h"

#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>


namespace io
{
namespace net
{

namespace
{

bool can_read(size_t message_len, size_t idx, size_t byte_count)
{
    return idx <= message_len && byte_count <= (message_len - idx);
}

template<typename T>
bool read_wire_value(const char* message, size_t message_len, size_t& idx, T& value)
{
    if (message == nullptr || !can_read(message_len, idx, sizeof(T)))
    {
        return false;
    }

    std::memcpy(&value, message + idx, sizeof(T));
    idx += sizeof(T);
    return true;
}

bool read_wire_string(const char* message, size_t message_len, size_t& idx, std::string& value)
{
    if (message == nullptr || idx >= message_len)
    {
        return false;
    }

    const void* terminator = std::memchr(message + idx, '\0', message_len - idx);
    if (terminator == nullptr)
    {
        return false;
    }

    const char* end = static_cast<const char*>(terminator);
    value.assign(message + idx, static_cast<size_t>(end - (message + idx)));
    idx += value.size() + 1;
    return true;
}

} // namespace

//-----------------------------------------------------------------------------
    
template<typename T_real>
Basic_Serializer<T_real>::Basic_Serializer()
{
    // _tmp_* are fixed-size member arrays; nothing to allocate.
}

template<typename T_real>
Basic_Serializer<T_real>::~Basic_Serializer()
{
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool Basic_Serializer<T_real>::_encode_meta(const data_struct::Stream_Block_Ptr<T_real>& stream_block, std::string& raw_msg)
{
    if (stream_block == nullptr)
    {
        return false;
    }

    //TODO:
    // add something to tell if 4 or 8 byte real
    _convert_var_to_bytes(&raw_msg, _tmp_uint, stream_block->detector_number(), sizeof(unsigned int));
    _convert_var_to_bytes(&raw_msg, _tmp_size_t, stream_block->row(), sizeof(size_t));
    _convert_var_to_bytes(&raw_msg, _tmp_size_t, stream_block->col(), sizeof(size_t));
    _convert_var_to_bytes(&raw_msg, _tmp_size_t, stream_block->height(), sizeof(size_t));
    _convert_var_to_bytes(&raw_msg, _tmp_size_t, stream_block->width(), sizeof(size_t));
    _convert_var_to_bytes(&raw_msg, _tmp_real, stream_block->theta, sizeof(T_real));
    // dataset name
    if (stream_block->dataset_name != nullptr)
    {
        raw_msg += *stream_block->dataset_name;
    }
    raw_msg += '\0';
    // dataset directory
    if (stream_block->dataset_directory != nullptr)
    {
        raw_msg += *stream_block->dataset_directory;
    }
    raw_msg += '\0';

    return true;
}

//-----------------------------------------------------------------------------

template<typename T_real>
data_struct::Stream_Block_Ptr<T_real> Basic_Serializer<T_real>::_decode_meta(char* message, size_t message_len, size_t& idx)
{
    unsigned int detector_number = 0;
    size_t row = 0;
    size_t col = 0;
    size_t height = 0;
    size_t width = 0;
    T_real theta = 0;
    std::string dataset_name;
    std::string dataset_directory;

    idx = 0;
    if (!read_wire_value(message, message_len, idx, detector_number) ||
        !read_wire_value(message, message_len, idx, row) ||
        !read_wire_value(message, message_len, idx, col) ||
        !read_wire_value(message, message_len, idx, height) ||
        !read_wire_value(message, message_len, idx, width) ||
        !read_wire_value(message, message_len, idx, theta) ||
        !read_wire_string(message, message_len, idx, dataset_name) ||
        !read_wire_string(message, message_len, idx, dataset_directory))
    {
        return nullptr;
    }

    data_struct::Stream_Block_Ptr<T_real> out_stream_block =
        std::make_shared<data_struct::Stream_Block<T_real> >(
            static_cast<int>(detector_number), row, col, height, width,
            dataset_directory, dataset_name);
    out_stream_block->theta = theta;

    return out_stream_block;
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Basic_Serializer<T_real>::_encode_counts(const data_struct::Stream_Block_Ptr<T_real>& stream_block, std::string& raw_msg)
{
    _convert_var_to_bytes(&raw_msg, _tmp_uint, stream_block->fitting_blocks.size(), 4);

    // iterate through fitting routine
    for( auto& itr : stream_block->fitting_blocks)
    {
        _convert_var_to_bytes(&raw_msg, _tmp_uint, itr.first, 4);
        _convert_var_to_bytes(&raw_msg, _tmp_uint, itr.second.fit_counts.size(), 4);
        // iterate through elements counts
        for(auto &itr2 : itr.second.fit_counts)
        {
            raw_msg += itr2.first;
            raw_msg += '\0';
            _convert_var_to_bytes(&raw_msg, _tmp_real, itr2.second, sizeof(T_real));
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
std::string Basic_Serializer<T_real>::encode_counts(const data_struct::Stream_Block_Ptr<T_real>& stream_block)
{
    std::string raw_msg = "";
    if (!_encode_meta(stream_block, raw_msg))
    {
        return {};
    }
    _encode_counts(stream_block, raw_msg);
	return raw_msg;
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool Basic_Serializer<T_real>::_decode_counts(char* message, size_t message_len, size_t& idx, const data_struct::Stream_Block_Ptr<T_real>& out_stream_block)
{
    std::uint32_t proc_type_count = 0;
    if (out_stream_block == nullptr || !read_wire_value(message, message_len, idx, proc_type_count))
    {
        return false;
    }

    // Each process block needs at least a 4-byte type and a 4-byte count.
    if (proc_type_count > (message_len - idx) / (sizeof(std::uint32_t) * 2))
    {
        return false;
    }

    for (std::uint32_t proc_type_itr = 0; proc_type_itr < proc_type_count; proc_type_itr++)
    {
        std::uint32_t proc_type_value = 0;
        std::uint32_t fit_block_size = 0;
        if (!read_wire_value(message, message_len, idx, proc_type_value) ||
            !read_wire_value(message, message_len, idx, fit_block_size))
        {
            return false;
        }

        // Even an empty element name needs a terminator followed by one value.
        const size_t minimum_entry_size = 1 + sizeof(T_real);
        if (fit_block_size > (message_len - idx) / minimum_entry_size)
        {
            return false;
        }

        const auto proc_type = static_cast<data_struct::Fitting_Routines>(proc_type_value);
        if (!data_struct::Fitting_Routine_To_Str.contains(proc_type))
        {
            return false;
        }
        auto& fit_counts = out_stream_block->fitting_blocks[proc_type].fit_counts;
        for (std::uint32_t i = 0; i < fit_block_size; i++)
        {
            std::string name;
            T_real val = 0.0;
            if (!read_wire_string(message, message_len, idx, name) ||
                !read_wire_value(message, message_len, idx, val))
            {
                return false;
            }
            fit_counts[std::move(name)] = val;
        }
    }
    return true;
}

//-----------------------------------------------------------------------------

template<typename T_real>
data_struct::Stream_Block_Ptr<T_real> Basic_Serializer<T_real>::decode_counts(char* message, size_t message_len)
{
    size_t idx = 0;
    data_struct::Stream_Block_Ptr<T_real> out_stream_block = _decode_meta(message, message_len, idx);
	if (out_stream_block == nullptr || !_decode_counts(message, message_len, idx, out_stream_block) || idx != message_len)
	{
		return nullptr;
	}
	return out_stream_block;
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool Basic_Serializer<T_real>::_encode_spectra(const data_struct::Stream_Block_Ptr<T_real>& stream_block, std::string& raw_msg)
{
	if (stream_block == nullptr || stream_block->spectra == nullptr)
	{
		return false;
	}

    const auto spectra_size = static_cast<size_t>(stream_block->spectra->size());
    constexpr size_t max_spectrum_size = static_cast<size_t>(std::numeric_limits<unsigned short>::max()) + 1;
    if (spectra_size == 0 || spectra_size > max_spectrum_size)
    {
        logE << "Cannot encode spectrum with " << spectra_size << " channels.\n";
        return false;
    }

    std::map<unsigned short, T_real> spec_indx_values;

    // Collect first so count overflow cannot leave a partially encoded packet.
    for (size_t i = 0; i < spectra_size; i++)
    {
        T_real spec_val = (*(stream_block->spectra))[static_cast<Eigen::Index>(i)];
        if(spec_val > (T_real)0.0)
        {
            if (spec_indx_values.size() == std::numeric_limits<unsigned short>::max())
            {
                logE << "Cannot encode more than " << std::numeric_limits<unsigned short>::max()
                     << " nonzero spectrum channels.\n";
                return false;
            }
            spec_indx_values.emplace(static_cast<unsigned short>(i), spec_val);
        }
    }

    const auto send_cnt = static_cast<unsigned short>(spec_indx_values.size());

    _convert_var_to_bytes(&raw_msg, _tmp_real, stream_block->spectra->elapsed_livetime(), sizeof(T_real));
    _convert_var_to_bytes(&raw_msg, _tmp_real, stream_block->spectra->elapsed_realtime(), sizeof(T_real));
    _convert_var_to_bytes(&raw_msg, _tmp_real, stream_block->spectra->input_counts(), sizeof(T_real));
    _convert_var_to_bytes(&raw_msg, _tmp_real, stream_block->spectra->output_counts(), sizeof(T_real));
    _convert_var_to_bytes(&raw_msg, _tmp_uint, stream_block->spectra->size(), 4);

    //if we have more than half the spec it doesn't help to to send index/values pairs
//    if(send_cnt < (stream_block->spectra->size() / 2))
//    {

    _convert_var_to_bytes(&raw_msg, _tmp_ushort, send_cnt, sizeof(unsigned short));

    //encode values
    for(auto const & itr : spec_indx_values)
    {
        _convert_var_to_bytes(&raw_msg, _tmp_ushort, itr.first, sizeof(unsigned short));
        _convert_var_to_bytes(&raw_msg, _tmp_real, itr.second, sizeof(T_real));
    }

//    }

    return true;
}

//-----------------------------------------------------------------------------

template<typename T_real>
std::string Basic_Serializer<T_real>::encode_spectra(const data_struct::Stream_Block_Ptr<T_real>& stream_block)
{
    std::string raw_msg = "";
    if (!_encode_meta(stream_block, raw_msg) || !_encode_spectra(stream_block, raw_msg))
    {
        return {};
    }
    return raw_msg;
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool Basic_Serializer<T_real>::_decode_spectra(char* message, size_t message_len, size_t& idx, const data_struct::Stream_Block_Ptr<T_real>& out_stream_block)
{
    T_real elt = 0;
    T_real ert = 0;
    T_real incnt = 0;
    T_real outcnt = 0;
    std::uint32_t spectra_size = 0;
    unsigned short recv_cnt = 0;
    unsigned short spec_index = 0;
    T_real spec_value = 0.0f;

    if (out_stream_block == nullptr ||
        !read_wire_value(message, message_len, idx, elt) ||
        !read_wire_value(message, message_len, idx, ert) ||
        !read_wire_value(message, message_len, idx, incnt) ||
        !read_wire_value(message, message_len, idx, outcnt) ||
        !read_wire_value(message, message_len, idx, spectra_size) ||
        !read_wire_value(message, message_len, idx, recv_cnt))
    {
        return false;
    }

    // Spectrum indexes use an unsigned-short wire representation.
    constexpr size_t max_spectrum_size = static_cast<size_t>(std::numeric_limits<unsigned short>::max()) + 1;
    if (spectra_size < 1 || static_cast<size_t>(spectra_size) > max_spectrum_size)
    {
        logE<<"Invalid spectra_size " << spectra_size << "\n";
        return false;
    }

    const size_t entry_size = sizeof(unsigned short) + sizeof(T_real);
    if (static_cast<size_t>(recv_cnt) > (message_len - idx) / entry_size)
    {
        return false;
    }

    auto spectra = std::make_unique<data_struct::Spectra<T_real> >(spectra_size);
    spectra->elapsed_livetime(elt);
    spectra->elapsed_realtime(ert);
    spectra->input_counts(incnt);
    spectra->output_counts(outcnt);

    for (unsigned short i = 0; i < recv_cnt; i++)
    {
        if (!read_wire_value(message, message_len, idx, spec_index) ||
            !read_wire_value(message, message_len, idx, spec_value) ||
            static_cast<size_t>(spec_index) >= static_cast<size_t>(spectra_size))
        {
            return false;
        }
        (*spectra)[spec_index] = spec_value;
    }

    delete out_stream_block->spectra;
    out_stream_block->spectra = spectra.release();
    return true;
}

//-----------------------------------------------------------------------------

template<typename T_real>
data_struct::Stream_Block_Ptr<T_real> Basic_Serializer<T_real>::decode_spectra(char* message, size_t message_len)
{
    size_t idx = 0;
    data_struct::Stream_Block_Ptr<T_real> out_stream_block = _decode_meta(message, message_len, idx);
	if (out_stream_block == nullptr || !_decode_spectra(message, message_len, idx, out_stream_block) || idx != message_len)
	{
		return nullptr;
	}
	return out_stream_block;
}

//-----------------------------------------------------------------------------

template<typename T_real>
std::string Basic_Serializer<T_real>::encode_counts_and_spectra(const data_struct::Stream_Block_Ptr<T_real>& in_stream_block)
{
    std::string raw_msg = "";
    if (!_encode_meta(in_stream_block, raw_msg))
    {
        return {};
    }
    _encode_counts(in_stream_block, raw_msg);
    if (!_encode_spectra(in_stream_block, raw_msg))
    {
        return {};
    }
    return raw_msg;
}

//-----------------------------------------------------------------------------

template<typename T_real>
data_struct::Stream_Block_Ptr<T_real> Basic_Serializer<T_real>::decode_counts_and_spectra(char* message, size_t message_len)
{
    size_t idx = 0;
    data_struct::Stream_Block_Ptr<T_real> out_stream_block = _decode_meta(message, message_len, idx);
	if (out_stream_block != nullptr)
	{
		if (idx >= message_len || !_decode_counts(message, message_len, idx, out_stream_block))
		{
			return nullptr;
		}
		if (idx >= message_len || !_decode_spectra(message, message_len, idx, out_stream_block))
		{
			return nullptr;
		}
		if (idx != message_len)
		{
			return nullptr;
		}
	}
    return out_stream_block;
}

//-----------------------------------------------------------------------------


TEMPLATE_CLASS_DLL_EXPORT Basic_Serializer<float>;
TEMPLATE_CLASS_DLL_EXPORT Basic_Serializer<double>;

} //end namespace net
}// end namespace io
