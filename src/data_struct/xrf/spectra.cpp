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
namespace xrf
{

template<typename T>
std::vector<T> conv_valid(std::vector<T> const &f, std::vector<T> const &g)
{
	size_t const nf = f.size();
	size_t const ng = g.size();
	std::vector<T> const &min_v = (nf < ng)? f : g;
	std::vector<T> const &max_v = (nf < ng)? g : f;
	size_t const n  = std::max(nf, ng) - std::min(nf, ng) + 1;
	std::vector<T> out(n, T());
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

std::valarray<real_t> convolve1d( std::valarray<real_t> arr, size_t boxcar_size)
{
    std::valarray<real_t> boxcar(1.0, boxcar_size);
    return convolve1d(arr, boxcar);
}

std::valarray<real_t> convolve1d( std::valarray<real_t> arr, std::valarray<real_t> boxcar)
{
    std::valarray<real_t> new_background((real_t)0.0, arr.size());
    //convolve 1d

    size_t const nf = arr.size();
	size_t const ng = boxcar.size();
    std::valarray<real_t> const &min_v = (nf < ng)? arr : boxcar;
    std::valarray<real_t> const &max_v = (nf < ng)? boxcar : arr;
	size_t const n  = std::max(nf, ng) - std::min(nf, ng) + 1;
    std::valarray<real_t> out((real_t)0.0, n);
    for(auto i(0); i < n; ++i)
    {
        for(int j(min_v.size() - 1), k(i); j >= 0; --j)
        {
            out[i] += min_v[j] * max_v[k];
            ++k;
        }
    }
    for(size_t i=0; i< arr.size(); i++)
    {
        if( out[i] != (real_t)0.0)
        {
            new_background[i] = out[i] / real_t(boxcar.size());
        }
    }


    /*
    for(size_t i=0; i< arr.size(); i++)
    {
        new_background[i] = 0.0;
        for (size_t j=0; j<boxcar.size(); j++)
        {
            if( (i-j) >= 0 )
                new_background[i] += arr[i-j] * boxcar[j];
        }
    }
    new_background /= real_t(boxcar.size());
    */
    return new_background;
}

std::valarray<real_t> snip_background(const Spectra* const spectra,
									  real_t energy_offset,
									  real_t energy_linear,
									  real_t energy_quadratic,
								      real_t spectral_binning,
									  real_t width,
									  real_t xmin,
									  real_t xmax)
{
	std::valarray<real_t> energy;
	std::valarray<real_t> index;
	std::valarray<real_t> background;
	size_t buffer_size = spectra->size();
	energy.resize(buffer_size);
	index.resize(buffer_size);
	background.resize(buffer_size);
	for (size_t i = 0; i < buffer_size; i++)
	{
		background[i] = (*spectra)[i];
		energy[i] = real_t(i);
		index = real_t(i);
	}

	if (spectral_binning > 0)
	{
		energy = energy * spectral_binning;
	}

	energy = energy_offset + energy * energy_linear + std::pow(energy, (real_t)2.0) * energy_quadratic;

	std::valarray<real_t> tmp = std::pow((energy_offset / (real_t)2.3548), (real_t)2.0) + energy * (real_t)2.96 * energy_linear;
	for (size_t i = 0; i<tmp.size(); i++)
	{
		if (tmp[i] < 0.0)
		{
			tmp[i] = 0.0;
		}
	}
	//std::valarray<real_t> fwhm = 2.35 * std::sqrt(tmp);
	std::valarray<real_t> current_width = (real_t)2.35 * std::sqrt(tmp);


	std::valarray<real_t> boxcar;
	std::valarray<real_t> new_background;
	new_background.resize(background.size());
	// smooth the background
	if (spectral_binning > 0)
	{
		boxcar.resize(3);
	}
	else
	{
		boxcar.resize(5);
	}

	for (size_t i = 0; i<boxcar.size(); i++)
	{
		boxcar[i] = 1.0;
	}

	//convolve 1d
    background = convolve1d(background, boxcar);
	//clear out
	new_background.resize(1);
	boxcar.resize(1);
	//fwhm
	current_width = width * current_width / energy_linear;  // in channels
	if (spectral_binning > 0)
	{
		current_width = current_width / (real_t)2.0;
	}

	background = std::log(std::log(background + (real_t)1.0) + (real_t)1.0);

	// FIRST SNIPPING
	int no_iterations = 2;
	if (spectral_binning > 0)
	{
		no_iterations = 3;
	}

	real_t max_of_xmin = (std::max)(xmin, (real_t)0.0);
	real_t min_of_xmax = (std::min)(xmax, real_t(buffer_size - 1));
	for (int j = 0; j<no_iterations; j++)
	{
		for (size_t k = 0; k<background.size(); k++)
		{
			real_t lo_index = k - current_width[k];
			real_t hi_index = k + current_width[k];
			if (lo_index < max_of_xmin)
			{
				lo_index = max_of_xmin;
			}
			if (hi_index > min_of_xmax)
			{
				hi_index = min_of_xmax;
			}
			real_t temp = (background[lo_index] + background[hi_index]) / (real_t)2.0;
			if (background[k] > temp)
			{
				background[k] = temp;
			}
		}
	}

	while (current_width.max() >= 0.5)
	{
		for (size_t k = 0; k<background.size(); k++)
		{
			real_t lo_index = k - current_width[k];
			real_t hi_index = k + current_width[k];
			if (lo_index < max_of_xmin)
			{
				lo_index = max_of_xmin;
			}
			if (hi_index > min_of_xmax)
			{
				hi_index = min_of_xmax;
			}
			real_t temp = (background[lo_index] + background[hi_index]) / (real_t)2.0;
			if (background[k] > temp)
			{
				background[k] = temp;
			}
		}

		current_width = current_width / real_t(M_SQRT2); // window_rf
	}

	background = std::exp(std::exp(background) - (real_t)1.0) - (real_t)1.0;

	for (size_t i = 0; i<background.size(); i++)
	{
		if (std::isnan(background[i]))
		{
			background[i] = 0.0;
		}
	}

	return background;

}

Range get_energy_range(real_t min_energy, real_t max_energy, size_t spectra_size, real_t energy_offset, real_t energy_slope)
{

   struct Range energy_range;
    energy_range.min = (int)ceil( (min_energy - energy_offset) / energy_slope );
    energy_range.max = (int)ceil( (max_energy - energy_offset) / energy_slope );
    //if (xmax > used_chan - 1) or (xmax <= np.amin([xmin, used_chan / 20.])):
    if ( (energy_range.max > spectra_size - 1) || (energy_range.max <= energy_range.min) )
    {
        energy_range.max = spectra_size - 1;
    }
    if (energy_range.min < 0 || energy_range.min > energy_range.max)
    {
        energy_range.min = 0;
    }
    return energy_range;

}

void gen_energy_vector(real_t number_channels, real_t energy_offset, real_t energy_slope, std::vector<real_t> *out_vec)
{

    out_vec->resize(number_channels);
    for(int i=0; i<number_channels; i++)
    {
        (*out_vec)[i] = (i * energy_slope) + energy_offset;
    }

}

} //namespace data_struct
} //namespace xrf
