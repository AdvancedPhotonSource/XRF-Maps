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



#include "spectra.h"
#include <algorithm>
#include <math.h>
#include <iostream>
#include <vector>

namespace data_struct
{

/*
template<typename T>
std::valarray<T> conv_valid(std::valarray<T> const &f, std::valarray<T> const &g)
{
	size_t const nf = f.size();
	size_t const ng = g.size();
	std::valarray<T> const &min_v = (nf < ng)? f : g;
	std::valarray<T> const &max_v = (nf < ng)? g : f;
	size_t const n  = std::max(nf, ng) - std::min(nf, ng) + 1;
	std::valarray<T> out(0.0, n);
	for(auto i(0); i < n; ++i)
	{
        for(int j(min_v.size() - 1), k(i); j >= 0; --j)
		{
			out[i] += min_v[j] * max_v[k];
			++k;
		}
	}
	return out;
}
*/

// ----------------------------------------------------------------------------

template<typename T_real>
ArrayTr<T_real> convolve1d(const ArrayTr<T_real>& arr, size_t boxcar_size)
{
    ArrayTr<T_real> boxcar(boxcar_size);
    boxcar.setConstant(boxcar_size, 1.0);
    return convolve1d(arr, boxcar);
}

// ----------------------------------------------------------------------------

template<typename T_real>
ArrayTr<T_real> convolve1d(const ArrayTr<T_real>& arr, const ArrayTr<T_real> &boxcar)
{
    ArrayTr<T_real> new_background(arr.size());
    new_background.setZero(arr.size());
    //convolve 1d

    size_t const nf = arr.size();
	size_t const ng = boxcar.size();
    ArrayTr<T_real> const &min_v = (nf < ng)? arr : boxcar;
    ArrayTr<T_real> const &max_v = (nf < ng)? boxcar : arr;
	size_t const n  = std::max(nf, ng) - std::min(nf, ng) + 1;
    ArrayTr<T_real> out(n);
    out.setZero(n);
    for(size_t i = 0; i < n; ++i)
    {
        for(int j(min_v.size() - 1), k(i); j >= 0; --j)
        {
            out[i] += min_v[j] * max_v[k];
            ++k;
        }
    }
	T_real norm = 1 / T_real(boxcar.size());
	int j = min_v.size() / 2;
    for(size_t i=0; i< n; i++)
    {
        new_background[j] = out[i] * norm;
		j++;
    }

    return new_background;
}

// ----------------------------------------------------------------------------

template<typename T_real>
ArrayTr<T_real> snip_background(const Spectra<T_real>* const spectra,
									  T_real energy_offset,
									  T_real energy_linear,
									  T_real energy_quadratic,
									  T_real width,
									  T_real xmin,
									  T_real xmax)
{
	ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(spectra->size(), 0, spectra->size() - 1);
    
	ArrayTr<T_real> background;

    energy = energy_offset + (energy * energy_linear) + (Eigen::pow(energy, (T_real)2.0) * energy_quadratic);

        ArrayTr<T_real> tmp = std::pow((energy_offset / (T_real)2.3548), (T_real)2.0) + energy * (T_real)2.96 * energy_linear;
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

    int max_of_xmin = (std::max)(xmin, (T_real)0.0);
    int min_of_xmax = (std::min)(xmax, T_real(spectra->size() - 1));
	for (int j = 0; j<no_iterations; j++)
	{
        for (long int k = 0; k<background.size(); k++)
		{
            long int lo_index = k - current_width[k];
            long int hi_index = k + current_width[k];
			if (lo_index < max_of_xmin)
			{
                lo_index = max_of_xmin;
			}
            if(lo_index > min_of_xmax)
            {
                lo_index = min_of_xmax;
            }
			if (hi_index > min_of_xmax)
			{
                hi_index = min_of_xmax;
			}
            if(hi_index < max_of_xmin)
            {
                hi_index = max_of_xmin;
            }
			T_real temp = (background[lo_index] + background[hi_index]) / (T_real)2.0;
			if (background[k] > temp)
			{
				background[k] = temp;
			}
		}
	}

	while (current_width.maxCoeff() >= 0.5)
	{
        for (long int k = 0; k<background.size(); k++)
		{
            long int lo_index = k - current_width[k];
            long int hi_index = k + current_width[k];
            if (lo_index < max_of_xmin)
            {
                lo_index = max_of_xmin;
            }
            if(lo_index > min_of_xmax)
            {
                lo_index = min_of_xmax;
            }
            if (hi_index > min_of_xmax)
            {
                hi_index = min_of_xmax;
            }
            if(hi_index < max_of_xmin)
            {
                hi_index = max_of_xmin;
            }
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


} //namespace data_struct

