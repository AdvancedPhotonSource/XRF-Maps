/***

Copyright (c) 2016 Arthur Glowacki

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

***/

#ifndef SPECTRA_H
#define SPECTRA_H

#include "defines.h"
#include <valarray>

namespace data_struct
{
namespace xrf
{

template<typename _T>
class Spectra_T : public std::valarray< _T >
{
public:
    /**
     * @brief Spectra : Constructor
     */
	Spectra_T() : std::valarray<_T>()
	{
		_elapsed_lifetime = 1.0;
		_elapsed_realtime = 1.0;
		_input_counts = 1.0;
		_output_counts = 1.0;
	}

	Spectra_T(size_t sample_size) : std::valarray<_T>(0.0, sample_size)
	{
		_elapsed_lifetime = 1.0;
		_elapsed_realtime = 1.0;
		_input_counts = 1.0;
		_output_counts = 1.0;
	}

    void recalc_elapsed_lifetime()
    {
        if(_input_counts == 0 || _output_counts == 0)
        {
            _elapsed_lifetime = _elapsed_realtime;
        }
        else
        {
            _elapsed_lifetime = _elapsed_realtime * _output_counts / _input_counts;
        }
    }

    void elapsed_lifetime(_T val) { _elapsed_lifetime = val; }

    const _T elapsed_lifetime() const { return _elapsed_lifetime; }

    void elapsed_realtime(_T val) { _elapsed_realtime = val; }

    const _T elapsed_realtime() const { return _elapsed_realtime; }

    void input_counts(_T val) { _input_counts = val; }

    const _T input_counts() const { return _input_counts; }

    void output_counts(_T val) { _output_counts = val; }

    const _T output_counts() const { return _output_counts; }

private:

    _T _elapsed_lifetime;
    _T _elapsed_realtime;
    _T _input_counts;
    _T _output_counts;

};


template DLL_EXPORT class Spectra_T<real_t>;
typedef Spectra_T<real_t> Spectra;

std::valarray<real_t> convolve1d(std::valarray<real_t> arr, size_t boxcar_size);
std::valarray<real_t> convolve1d(std::valarray<real_t> arr, std::valarray<real_t> boxcar);
std::valarray<real_t> snip_background(const Spectra * const spectra, real_t energy_offset, real_t energy_linear, real_t energy_quadratic, real_t spectral_binning, real_t width, real_t xmin, real_t xmax);

} //namespace xrf
} //namespace data_struct

#endif // SPECTRA_H
