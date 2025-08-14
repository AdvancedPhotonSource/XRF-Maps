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
#include <algorithm>
#include <functional>

namespace data_struct
{

#define default_time_and_io_counts 0.000000001

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
        _elapsed_livetime = (_T)default_time_and_io_counts;
		_elapsed_realtime = (_T)default_time_and_io_counts;
		_input_counts = 0.0;
		_output_counts = 0.0;
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
        _elapsed_livetime = (_T)default_time_and_io_counts;
		_elapsed_realtime = (_T)default_time_and_io_counts;
		_input_counts = 0.0;
		_output_counts = 0.0;
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
        _elapsed_livetime = (_T)default_time_and_io_counts;
        _elapsed_realtime = (_T)default_time_and_io_counts;
        _input_counts = 0.0;
        _output_counts = 0.0;
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
        _elapsed_livetime = (_T)default_time_and_io_counts;
        _elapsed_realtime = (_T)default_time_and_io_counts;
        _input_counts = 0.0;
        _output_counts = 0.0;
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
        _elapsed_livetime = (_T)default_time_and_io_counts;
        _elapsed_realtime = (_T)default_time_and_io_counts;
        _input_counts = 0.0;
        _output_counts = 0.0;
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

    void add(const Spectra<_T>& spectra)
    {
        *this += (ArrayTr<_T>)spectra;
        _T val = spectra.elapsed_livetime();
        if(std::isfinite(val))
        {
            _elapsed_livetime += val;
        }
        else
        {
            logW << " Spectra _elapsed_livetime = " << val << " . Not updating!\n";
        }
        val = spectra.elapsed_realtime();
        if(std::isfinite(val))
        {
            _elapsed_realtime += val;
        }
        else
        {
            logW << " Spectra _elapsed_realtime = " << val << " . Not updating!\n";
        }
        val = spectra.input_counts();
        if(std::isfinite(val))
        {
            _input_counts += val;
        }
        else
        {
            logW << " Spectra _input_counts = " << val << " . Not updating!\n";
        }
        val = spectra.output_counts();
        if(std::isfinite(val))
        {
            _output_counts += val;
        }
        else
        {
            logW << " Spectra _output_counts = " << val << " . Not updating!\n";
        }
    }

    // Note: T_real is different template typename, this is to convert double to float or the other way around.
    template<typename T_real>
    void add(const Spectra<T_real>* spectra)
    {
        if (spectra != nullptr)
        {
            for (int i = 0; i < this->size(); i++)
            {
                (ArrayTr<_T>(*this))(i) += static_cast<_T>( (*spectra)(i) );
            }
            _T val = spectra->elapsed_livetime();
            if (std::isfinite(val))
            {
                _elapsed_livetime += val;
            }
            else
            {
                logW << " Spectra _elapsed_livetime = " << val << " . Not updating!\n";
            }
            val = spectra->elapsed_realtime();
            if (std::isfinite(val))
            {
                _elapsed_realtime += val;
            }
            else
            {
                logW << " Spectra _elapsed_realtime = " << val << " . Not updating!\n";
            }
            val = spectra->input_counts();
            if (std::isfinite(val))
            {
                _input_counts += val;
            }
            else
            {
                logW << " Spectra _input_counts = " << val << " . Not updating!\n";
            }
            val = spectra->output_counts();
            if (std::isfinite(val))
            {
                _output_counts += val;
            }
            else
            {
                logW << " Spectra _output_counts = " << val << " . Not updating!\n";
            }
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

//-----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT ArrayTr<T_real> convolve1d(const ArrayTr<T_real>& arr, const ArrayTr<T_real>& boxcar)
{
    ArrayTr<T_real> new_background(arr.size());
    new_background.setZero(arr.size());
    //convolve 1d

    size_t const nf = arr.size();
    size_t const ng = boxcar.size();
    ArrayTr<T_real> const& min_v = (nf < ng) ? arr : boxcar;
    ArrayTr<T_real> const& max_v = (nf < ng) ? boxcar : arr;
    size_t const n = std::max(nf, ng) - std::min(nf, ng) + 1;
    ArrayTr<T_real> out(n);
    out.setZero(n);
    for (size_t i = 0; i < n; ++i)
    {
        for (int j(min_v.size() - 1), k(i); j >= 0; --j)
        {
            out[i] += min_v[j] * max_v[k];
            ++k;
        }
    }
    T_real norm = 1 / T_real(boxcar.size());
    size_t j = min_v.size() / 2;
    for (size_t i = 0; i < n; i++)
    {
        new_background[j] = out[i] * norm;
        j++;
    }

    return new_background;
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT ArrayTr<T_real> convolve1d(const ArrayTr<T_real>& arr, size_t boxcar_size)
{
    ArrayTr<T_real> boxcar(boxcar_size);
    boxcar.setConstant(boxcar_size, 1.0);
    return convolve1d(arr, boxcar);
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT ArrayTr<T_real> snip_background(const Spectra<T_real> * const spectra, T_real energy_offset, T_real energy_linear, T_real energy_quadratic, T_real width, T_real xmin, T_real xmax)
{
    ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(spectra->size(), 0, spectra->size() - 1);

    ArrayTr<T_real> background;

    energy = energy_offset + (energy * energy_linear) + (Eigen::pow(energy, (T_real)2.0) * energy_quadratic);

    ArrayTr<T_real> tmp = std::pow((energy_offset / (T_real)2.3548), (T_real)2.0) + energy * (T_real)3.58 * energy_linear;
    tmp = tmp.unaryExpr([](T_real r) { return r < 0.0 ? (T_real)0.0 : r;  });

    //ArrayTr<T_real> fwhm = 2.35 * std::sqrt(tmp);
    ArrayTr<T_real> current_width = (T_real)2.35 * Eigen::sqrt(tmp);

    ArrayTr<T_real> boxcar;
    // smooth the background
    boxcar.resize(5);
    boxcar.setConstant(5, 1.0);


    if (spectra != nullptr)
    {
        if (spectra->size() > 0)
        {
            //convolve 1d
            background = convolve1d(*spectra, boxcar);
        }
        else
        {
            return background;
        }
    }
    else
    {
        return background;
    }
    //fwhm
    current_width = width * current_width / energy_linear;  // in channels

    background = Eigen::log(Eigen::log(background + (T_real)1.0) + (T_real)1.0);

    // FIRST SNIPPING
    int no_iterations = 2;
    
    int lo_limit = std::max<int>(xmin, 0);
    int up_limit = std::min<int>(xmax, (spectra->size() - 1));
    for (int j = 0; j < no_iterations; j++)
    {
        for (long int k = 0; k < background.size(); k++)
        {
            int lo_index = std::max<int>(k - current_width[k], lo_limit);
            lo_index = std::min<int>(lo_index, up_limit);
            int hi_index = std::max<int>(k + current_width[k], lo_limit);
            hi_index = std::min<int>(k + current_width[k], up_limit);
            T_real temp = (background[lo_index] + background[hi_index]) / (T_real)2.0;
            if (background[k] > temp)
            {
                background[k] = temp;
            }
        }
    }

    while (current_width.maxCoeff() >= 0.5)
    {
        for (long int k = 0; k < background.size(); k++)
        {
            int lo_index = std::max<int>(k - current_width[k], lo_limit);
            lo_index = std::min<int>(lo_index, up_limit);
            int hi_index = std::max<int>(k + current_width[k], lo_limit);
            hi_index = std::min<int>(k + current_width[k], up_limit);
            T_real temp = (background[lo_index] + background[hi_index]) / (T_real)2.0;
            if (background[k] > temp)
            {
                background[k] = temp;
            }
        }

        current_width = current_width / T_real(M_SQRT2); // window_rf
    }

    background = Eigen::exp(Eigen::exp(background) - (T_real)1.0) - (T_real)1.0;
    background = background.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });

    return background;

}

// ----------------------------------------------------------------------------
/*
template<typename T_real>
DLL_EXPORT ArrayTr<T_real> snip_background2(const Spectra<T_real> * const spectra, T_real width)
{
    width = width * 100.0;
    assert(spectra != nullptr);
    ArrayTr<T_real> boxcar;
    boxcar.resize(5);
    boxcar.setConstant(5, 1.0);
    ArrayTr<T_real> background = convolve1d(*spectra, boxcar);
    background += 1.0;
    background = Eigen::sqrt(background);
    background = Eigen::log(Eigen::log(background + (T_real)1.0) + (T_real)1.0);

    while (width >= 0.5)
    {
        for (int k = 0; k < background.size(); k++)
        {
            int lo = std::max<int>( k - (int)width, 0);
            int hi = std::min<int>( k + (int)width, background.size()-1);
            T_real temp = (background[lo] + background[hi]) / (T_real)2.0;
            if (background[k] > temp)
            {
                background[k] = temp;
            }
        }
        width = width / T_real(M_SQRT2); // window_rf
    }

    background = Eigen::exp(Eigen::exp(background) - (T_real)1.0) - (T_real)1.0;
    background = Eigen::pow(background, 2.0);
    background -= 1.0;
    background = background.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });

    return background;

}
*/
// ----------------------------------------------------------------------------


template<typename T_real>
using IO_Callback_Func_Def = std::function<void(size_t, size_t, size_t, size_t, size_t, Spectra<T_real>*, void*)>;

} //namespace data_struct

#endif // SPECTRA_H
