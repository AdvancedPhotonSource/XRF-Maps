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

#include "core/defines.h"
#include <Eigen/Core>
#include <vector>
#include <functional>

namespace data_struct
{

using namespace std;

template<typename _T>
using ArrayTr = Eigen::Array<_T, Eigen::Dynamic, Eigen::RowMajor>;

template<typename _T>
using VectorTr = Eigen::Vector<_T, Eigen::Dynamic>;

template<typename _T>
class Spectra : public ArrayTr<_T>
{
public:
    /**
     * @brief Spectra : Constructor
     */
    Spectra() : ArrayTr<_T>()
	{
        _elapsed_livetime = 1.0;
		_elapsed_realtime = 1.0;
		_input_counts = 1.0;
		_output_counts = 1.0;
	}

    Spectra(const Spectra &spectra) : ArrayTr<_T>(spectra)
	{
        _elapsed_livetime = spectra._elapsed_livetime;
		_elapsed_realtime = spectra._elapsed_realtime;
		_input_counts = spectra._input_counts;
		_output_counts = spectra._output_counts;
	}

    Spectra(size_t sample_size) : ArrayTr<_T>(sample_size)
	{
		this->setZero();
        _elapsed_livetime = 1.0;
		_elapsed_realtime = 1.0;
		_input_counts = 1.0;
		_output_counts = 1.0;
	}

    Spectra(size_t sample_size, _T elt, _T ert, _T incnt, _T outcnt) : ArrayTr<_T>(sample_size)
    {
        this->setZero();
        _elapsed_livetime = elt;
        _elapsed_realtime = ert;
        _input_counts = incnt;
        _output_counts = outcnt;
    }

    Spectra(const ArrayTr<_T>& arr) : ArrayTr<_T>(arr)
    {
        _elapsed_livetime = 1.0;
        _elapsed_realtime = 1.0;
        _input_counts = 1.0;
        _output_counts = 1.0;
    }

    Spectra(const ArrayTr<_T>& arr, _T livetime, _T realtime, _T incnt, _T outnt) : ArrayTr<_T>(arr)
    {
        _elapsed_livetime = livetime;
        _elapsed_realtime = realtime;
        _input_counts = incnt;
        _output_counts = outnt;
    }

    Spectra(const ArrayTr<_T>&& arr) : ArrayTr<_T>(arr)
    {
        _elapsed_livetime = 1.0;
        _elapsed_realtime = 1.0;
        _input_counts = 1.0;
        _output_counts = 1.0;
    }

    Spectra(const ArrayTr<_T>&& arr, _T livetime, _T realtime, _T incnt, _T outnt) : ArrayTr<_T>(arr)
    {
        _elapsed_livetime = livetime;
        _elapsed_realtime = realtime;
        _input_counts = incnt;
        _output_counts = outnt;
    }

    Spectra(Eigen::Index& rows, Eigen::Index& cols) : ArrayTr<_T>(rows, cols)
	{
        _elapsed_livetime = 1.0;
        _elapsed_realtime = 1.0;
        _input_counts = 1.0;
        _output_counts = 1.0;
	}

    virtual ~Spectra()
    {

    }

    void recalc_elapsed_livetime()
    {
        if(_input_counts == 0 || _output_counts == 0)
        {
            _elapsed_livetime = _elapsed_realtime;
        }
        else
        {
            _elapsed_livetime = _elapsed_realtime * _output_counts / _input_counts;
        }
    }

    void add(const Spectra& spectra)
    {
        *this += (ArrayTr<_T>)spectra;
        _T val = spectra.elapsed_livetime();
        if(std::isfinite(val))
        {
            _elapsed_livetime += val;
        }
        val = spectra.elapsed_realtime();
        if(std::isfinite(val))
        {
            _elapsed_realtime += val;
        }
        val = spectra.input_counts();
        if(std::isfinite(val))
        {
            _input_counts += val;
        }
        val = spectra.output_counts();
        if(std::isfinite(val))
        {
            _output_counts += val;
        }
    }

    void elapsed_livetime(_T val) { _elapsed_livetime = val; }

    const _T elapsed_livetime() const { return _elapsed_livetime; }

    void elapsed_realtime(_T val) { _elapsed_realtime = val; }

    const _T elapsed_realtime() const { return _elapsed_realtime; }

    void input_counts(_T val) { _input_counts = val; }

    const _T input_counts() const { return _input_counts; }

    void output_counts(_T val) { _output_counts = val; }

    const _T output_counts() const { return _output_counts; }

    Spectra sub_spectra(size_t start, size_t count) const
	{
        return Spectra(this->segment(start, count), _elapsed_livetime, _elapsed_realtime, _input_counts, _output_counts);
	}

private:

    _T _elapsed_livetime;
    _T _elapsed_realtime;
    _T _input_counts;
    _T _output_counts;

};

TEMPLATE_CLASS_DLL_EXPORT Spectra<float>;
TEMPLATE_CLASS_DLL_EXPORT Spectra<double>;

template<typename T_real>
DLL_EXPORT ArrayTr<T_real> convolve1d(const ArrayTr<T_real>& arr, size_t boxcar_size);

template DLL_EXPORT ArrayTr<float> convolve1d<float>(const ArrayTr<float>& arr, size_t boxcar_size);
template DLL_EXPORT ArrayTr<double> convolve1d<double>(const ArrayTr<double>& arr, size_t boxcar_size);

template<typename T_real>
DLL_EXPORT ArrayTr<T_real> convolve1d(const ArrayTr<T_real>& arr, const ArrayTr<T_real>& boxcar);

template DLL_EXPORT ArrayTr<float> convolve1d<float>(const ArrayTr<float>& arr, const ArrayTr<float>& boxcar);
template DLL_EXPORT ArrayTr<double> convolve1d<double>(const ArrayTr<double>& arr, const ArrayTr<double>& boxcar);

template<typename T_real>
DLL_EXPORT ArrayTr<T_real> snip_background(const Spectra<T_real> * const spectra, T_real energy_offset, T_real energy_linear, T_real energy_quadratic, T_real width, T_real xmin, T_real xmax);

template DLL_EXPORT ArrayTr<float> snip_background<float>(const Spectra<float>* const spectra, float energy_offset, float energy_linear, float energy_quadratic, float width, float xmin, float xmax);
template DLL_EXPORT ArrayTr<double> snip_background<double>(const Spectra<double>* const spectra, double energy_offset, double energy_linear, double energy_quadratic, double width, double xmin, double xmax);

template<typename T_real>
DLL_EXPORT void gen_energy_vector(T_real number_channels, T_real energy_offset, T_real energy_slope, std::vector<T_real> *out_vec);

template DLL_EXPORT void gen_energy_vector<float>(float number_channels, float energy_offset, float energy_slope, std::vector<float>* out_vec);
template DLL_EXPORT void gen_energy_vector<double>(double number_channels, double energy_offset, double energy_slope, std::vector<double>* out_vec);

template<typename T_real>
using IO_Callback_Func_Def = std::function<void(size_t, size_t, size_t, size_t, size_t, Spectra<T_real>*, void*)>;

} //namespace data_struct

#endif // SPECTRA_H
