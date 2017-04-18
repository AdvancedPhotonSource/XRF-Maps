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

/// Initial Author <2016>: Arthur Glowacki


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

    Spectra_T(size_t sample_size, _T elt, _T ert, _T incnt, _T outcnt) : std::valarray<_T>(0.0, sample_size)
    {
        _elapsed_lifetime = elt;
        _elapsed_realtime = ert;
        _input_counts = incnt;
        _output_counts = outcnt;
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

    void add(const Spectra_T& spectra)
    {
        *this += (std::valarray<_T>)spectra;
        _elapsed_lifetime += spectra.elapsed_lifetime();
        _elapsed_realtime += spectra.elapsed_realtime();
        _input_counts += spectra.input_counts();
        _output_counts += spectra.output_counts();
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
