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



#include "gaussian_model.h"

#include <iostream>
#include <algorithm>
#include <math.h>

#include "Faddeeva.hh"
#include <string.h>

#define SQRT_2xPI (real_t)2.506628275 // sqrt ( 2.0 * M_PI )

using namespace data_struct::xrf;


namespace fitting
{
namespace models
{

// ----------------------------------------------------------------------------

Gaussian_Model::Gaussian_Model() : Base_Model()
{

    _optimizer = nullptr;
    _snip_background = true;

    _fit_parameters = _generate_default_fit_parameters();


}

// ----------------------------------------------------------------------------

Gaussian_Model::~Gaussian_Model()
{

}

// ----------------------------------------------------------------------------

Fit_Parameters Gaussian_Model::_generate_default_fit_parameters()
{

    Fit_Parameters fit_params;
    //                                                        name                     min               max             val              step              use
    fit_params.add_parameter(STR_ENERGY_OFFSET,    Fit_Param(STR_ENERGY_OFFSET,        (real_t)-0.2,    (real_t)0.2,    (real_t)0.0,    (real_t)0.00001,   E_Bound_Type::LIMITED_LO_HI)); //LIMITED_LO_HI  fixed
    fit_params.add_parameter(STR_ENERGY_SLOPE,     Fit_Param(STR_ENERGY_OFFSET,        (real_t)0.001,   (real_t)0.1,    (real_t)1.0,    (real_t)0.00001,   E_Bound_Type::LIMITED_LO_HI)); //LIMITED_LO_HI  fixed
    fit_params.add_parameter(STR_ENERGY_QUADRATIC, Fit_Param(STR_ENERGY_OFFSET,        (real_t)-0.0001, (real_t)0.0001, (real_t)0.0,    (real_t)0.00001,   E_Bound_Type::LIMITED_LO_HI)); //LIMITED_LO_HI  fixed

    fit_params.add_parameter(STR_FWHM_OFFSET,    Fit_Param(STR_FWHM_OFFSET,        (real_t)0.005,    (real_t)0.5,  (real_t)0.12,    (real_t)0.00001,   E_Bound_Type::LIMITED_LO_HI)); //LIMITED_LO_HI  fixed
    fit_params.add_parameter(STR_FWHM_FANOPRIME, Fit_Param(STR_FWHM_FANOPRIME,     (real_t)0.000001, (real_t)0.05, (real_t)0.00012, (real_t)0.000001,  E_Bound_Type::LIMITED_LO_HI)); //LIMITED_LO_HI  fixed

    fit_params.add_parameter(STR_COHERENT_SCT_ENERGY, Fit_Param(STR_COHERENT_SCT_ENERGY,	   (real_t)9.4, (real_t)10.4, (real_t)9.99, (real_t)0.001,  E_Bound_Type::LIMITED_LO_HI)); // LIMITED_LO_HI fixed  0
    fit_params.add_parameter(STR_COHERENT_SCT_AMPLITUDE, Fit_Param(STR_COHERENT_SCT_AMPLITUDE, (real_t)0.000001, (real_t)10.0, (real_t)10.0,  (real_t)0.00001, E_Bound_Type::FIT));

    fit_params.add_parameter(STR_COMPTON_ANGLE, Fit_Param(STR_COMPTON_ANGLE,		 (real_t)-0.0001, (real_t)0.0001, (real_t)90.0, (real_t)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_COMPTON_FWHM_CORR, Fit_Param(STR_COMPTON_FWHM_CORR, (real_t)-0.0001, (real_t)0.0001, (real_t)1.0,  (real_t)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_COMPTON_AMPLITUDE, Fit_Param(STR_COMPTON_AMPLITUDE, (real_t)-0.0001, (real_t)0.0001, (real_t)10.0, (real_t)0.000001,       E_Bound_Type::FIT));
    fit_params.add_parameter(STR_COMPTON_F_STEP, Fit_Param(STR_COMPTON_F_STEP,		 (real_t)-0.0001, (real_t)0.0001, (real_t)0.0,  (real_t)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_COMPTON_F_TAIL, Fit_Param(STR_COMPTON_F_TAIL,		 (real_t)-0.0001, (real_t)0.0001, (real_t)0.0,  (real_t)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_COMPTON_GAMMA, Fit_Param(STR_COMPTON_GAMMA,		 (real_t)-0.0001, (real_t)0.0001, (real_t)1.0,  (real_t)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_COMPTON_HI_F_TAIL, Fit_Param(STR_COMPTON_HI_F_TAIL, (real_t)-0.0001, (real_t)0.0001, (real_t)0.0,  (real_t)0.0000001, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_COMPTON_HI_GAMMA, Fit_Param(STR_COMPTON_HI_GAMMA,	 (real_t)-0.0001, (real_t)0.0001, (real_t)1.0,  (real_t)0.01,      E_Bound_Type::FIXED));

    fit_params.add_parameter(STR_SNIP_WIDTH, Fit_Param(STR_SNIP_WIDTH,			     (real_t)0.1,   (real_t)2.828427, (real_t)0.15, (real_t)0.01,      E_Bound_Type::FIXED)); //max = 2* sqrt(2)

    fit_params.add_parameter(STR_F_STEP_OFFSET, Fit_Param(STR_F_STEP_OFFSET,		(real_t)0.0, (real_t)1.0, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_F_STEP_LINEAR, Fit_Param(STR_F_STEP_LINEAR,        (real_t)0.0, (real_t)1.0, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_F_STEP_QUADRATIC, Fit_Param(STR_F_STEP_QUADRATIC,  (real_t)0.0, (real_t)0.0, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));

    fit_params.add_parameter(STR_F_TAIL_OFFSET, Fit_Param(STR_F_TAIL_OFFSET,		(real_t)0.0, (real_t)0.1, (real_t)0.04, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_F_TAIL_LINEAR, Fit_Param(STR_F_TAIL_LINEAR,		(real_t)0.0, (real_t)1.0, (real_t)0.01, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_F_TAIL_QUADRATIC, Fit_Param(STR_F_TAIL_QUADRATIC,  (real_t)0.0, (real_t)0.01, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));

    fit_params.add_parameter(STR_GAMMA_OFFSET, Fit_Param(STR_GAMMA_OFFSET,			(real_t)0.1, (real_t)10.0, (real_t)0.0, (real_t)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_GAMMA_LINEAR, Fit_Param(STR_GAMMA_LINEAR,			(real_t)0.0, (real_t)3.0,  (real_t)0.2, (real_t)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_GAMMA_QUADRATIC, Fit_Param(STR_GAMMA_QUADRATIC,	(real_t)0.0, (real_t)0.0,  (real_t)0.0, (real_t)0.1, E_Bound_Type::FIXED));

    fit_params.add_parameter(STR_KB_F_TAIL_OFFSET, Fit_Param(STR_KB_F_TAIL_OFFSET,		 (real_t)0.0, (real_t)0.2,  (real_t)0.0, (real_t)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_KB_F_TAIL_LINEAR, Fit_Param(STR_KB_F_TAIL_LINEAR,		 (real_t)0.0, (real_t)0.02, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(STR_KB_F_TAIL_QUADRATIC, Fit_Param(STR_KB_F_TAIL_QUADRATIC, (real_t)0.0, (real_t)0.0,  (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));

    return fit_params;

}

// ----------------------------------------------------------------------------

const Spectra Gaussian_Model::model_spectrum(const Fit_Parameters * const fit_params,
                                             const unordered_map<string, Fit_Element_Map*> * const elements_to_fit,
                                             const struct Range energy_range)
{

    Spectra agr_spectra(energy_range.count());
    Spectra spectra_model;


    std::valarray<real_t> energy((real_t)0.0, energy_range.count());
    //std::valarray<real_t> background_counts((real_t)0.0, energy_range.count());
    real_t e_val = energy_range.min;
    for(int i=0; i < (energy_range.max - energy_range.min )+1; i++)
    {
        energy[i] = e_val;
        e_val += 1.0;
    }

    //real_t gain = detector->energy_slope();
    std::valarray<real_t> ev = fit_params->at(STR_ENERGY_OFFSET).value + energy * fit_params->at(STR_ENERGY_SLOPE).value + std::pow(energy, (real_t)2.0) * fit_params->at(STR_ENERGY_QUADRATIC).value;

/* move outside
    if( _snip_background )
    {
        real_t spectral_binning = 0.0;
        background_counts = snip_background(spectra, detector->energy_offset(), detector->energy_slope(), detector->energy_quadratic(), spectral_binning, fit_params->at(STR_SNIP_WIDTH).value, energy_range.min, energy_range.max);
    }
*/
    for(const auto& itr : (*elements_to_fit))
    {
        //Fit_Element_Map* element = itr.second;
        spectra_model = model_spectrum_element(fit_params, itr.second, ev, energy);
        agr_spectra += spectra_model;
    }

    agr_spectra += elastic_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);
    agr_spectra += compton_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);

    //agr_spectra += background_counts;

/*
    // pileup
    //temp_element_pos = np.array(keywords.added_params[4:13]);
    for (int ii=0; ii<9; ii++)
    {
        // skip calculation for peaks (element + scatter) that are fixed AND
        // close to zero (default 'zero' at 1e-10)
        if (pall[temp_element_pos[ii]] <= -10)
        {
            continue;
        }
        int j = 0;
        int i = np.amax(keywords.mele_pos) - np.amin(keywords.kele_pos) + 1 + ii;
        if (er_struct.energy <= 0.0)
        {
            continue;
        }
        //repetition delta_energy = ev.copy() - (er_struct.energy);
        faktor = er_struct.ratio* (10.**pall[temp_element_pos[ii]]);
        // peak, gauss
        value = faktor * this->peak(gain, sigma[i, j], delta_energy);
        //fit_counts.pileup = fit_counts.pileup + value;
        *(agr_spectra.buffer()) += value;
    }

    //(*counts) = fit_counts.ka + fit_counts.kb + fit_counts.l + fit_counts.m + fit_counts.elastic + fit_counts.compton + fit_counts.step + fit_counts.tail + fit_counts.pileup;
    //escape
    //if (np.sum(np.abs(pall[keywords.added_params[1:4]])) >= 0.0)
    {
        // si escape
        if (fit_params->at(STR_SI_ESCAPE).value > 0.0)
        //if (pall[keywords.added_params[1]] > 0.0)
        {
            real_t escape_E = (real_t)1.73998;
            for (size_t i=0; i<ev.size(); i++)
            {
                if( ev[i] > (escape_E+ev[0]) )
                {
                    escape_factor = pall[keywords.added_params[1]] + pall[keywords.added_params[3]] * ev;

                    for ( int ii=0; ii<len(wo[0]); ii++ )
                    {
                        fit_counts.escape[ii] = counts[wo[0][ii]] * std::max(np.append(escape_factor[wo[0][ii]], 0.0));
                    }
                }
            }
        }
        *(agr_spectra.buffer()) += fit_counts.escape;
    }

*/

    return agr_spectra;
}

// ----------------------------------------------------------------------------

const Spectra Gaussian_Model::model_spectrum_element(const Fit_Parameters * const fitp, const Fit_Element_Map * const element_to_fit, const std::valarray<real_t> &ev, std::valarray<real_t> energy)
{
    Spectra spectra_model(energy.size());
    std::valarray<real_t> peak_counts(0.0, energy.size());

    //real_t gain = detector->energy_slope();
    //std::valarray<real_t> ev = detector->energy_offset() + energy * detector->energy_slope() + std::pow(energy, (real_t)2.0) * detector->energy_quadratic();

    if(false == fitp->contains(element_to_fit->full_name()))
    {
        return spectra_model;
    }

    real_t pre_faktor = std::pow((real_t)10.0 , fitp->at(element_to_fit->full_name()).value);

    if(std::isnan(pre_faktor))
        return spectra_model;

    for (const Element_Energy_Ratio& er_struct : element_to_fit->energy_ratios())
    {

        real_t sigma = std::sqrt( std::pow((fitp->at(STR_FWHM_OFFSET).value / (real_t)2.3548), (real_t)2.0) + (er_struct.energy) * (real_t)2.96 * fitp->at(STR_FWHM_FANOPRIME).value );
        real_t f_step =  std::abs( er_struct.mu_fraction * ( fitp->at(STR_F_STEP_OFFSET).value + (fitp->at(STR_F_STEP_LINEAR).value * er_struct.energy)));
        real_t f_tail = std::abs( fitp->at(STR_F_TAIL_OFFSET).value + (fitp->at(STR_F_TAIL_LINEAR).value * er_struct.mu_fraction));
        real_t kb_f_tail = std::abs(  fitp->at(STR_KB_F_TAIL_OFFSET).value + (fitp->at(STR_KB_F_TAIL_LINEAR).value * er_struct.mu_fraction));
        real_t value = 1.0;

        //don't process if energy is 0
        if (er_struct.ratio == 0.0)
            continue;
        if (er_struct.energy <= 0.0)
            continue;

        // gaussian peak shape
        std::valarray<real_t> delta_energy = ev - er_struct.energy;

        real_t faktor = real_t(er_struct.ratio * pre_faktor);
        if (er_struct.ptype == Element_Param_Type::Kb_Line)
        {
            faktor = faktor / ((real_t)1.0 + kb_f_tail + f_step);
        }
        if ((er_struct.ptype == Element_Param_Type::Ka_Line) || (er_struct.ptype == Element_Param_Type::L_Line))
        {
            faktor = faktor / ((real_t)1.0 + f_tail + f_step);
        }
        // peak, gauss
        //value = faktor * this->peak(fit_params->at(STR_ENERGY_SLOPE).value, sigma, delta_energy);
        spectra_model += faktor * this->peak(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy);

/*      //separates counts into which lines they are
        peak_counts = faktor * this->peak(fit_params->at(STR_ENERGY_SLOPE).value, sigma, delta_energy);
        switch(er_struct.ptype)
        {
            case Element_Param_Type::Ka_Line:
                counts_arr->ka += peak_counts;
            break;
            case Element_Param_Type::Kb_Line:
                counts_arr->kb += peak_counts;
            break;
            case Element_Param_Type::L_Line:
                counts_arr->l += peak_counts;
            break;
            case Element_Param_Type::M_Line:
                counts_arr->m += peak_counts;
            break;
            default:
            break;
        }
        (*counts) += peak_counts;
*/


        //  peak, step
        if (f_step > 0.0)
        {
            value = faktor * f_step;
            //value = value * this->step(gain, sigma, delta_energy, er_struct.energy);
            spectra_model += value * this->step(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy, er_struct.energy);
            //counts_arr->step = fit_counts.step + value;
        }
        //  peak, tail;; use different tail for K beta vs K alpha lines
        if (er_struct.ptype == Element_Param_Type::Kb_Line)
        {
            real_t gamma = std::abs( fitp->at(STR_GAMMA_OFFSET).value + fitp->at(STR_GAMMA_LINEAR).value * (er_struct.energy) ) * element_to_fit->width_multi();
            value = faktor * kb_f_tail;
            spectra_model += value * this->tail(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy, gamma);
            //fit_counts.tail = fit_counts.tail + value;
        }
    }
    return spectra_model;
}

// ----------------------------------------------------------------------------

const std::valarray<real_t> Gaussian_Model::peak(real_t gain, real_t sigma, std::valarray<real_t>& delta_energy) const
{
    // gain / (sigma * sqrt( 2.0 * M_PI) ) * exp( -0.5 * ( (delta_energy / sigma) ** 2 )
    return gain / ( sigma * SQRT_2xPI ) *  std::exp((real_t)-0.5 * std::pow((delta_energy / sigma), (real_t)2.0) );
}

// ----------------------------------------------------------------------------

const std::valarray<real_t> Gaussian_Model::step(real_t gain, real_t sigma, std::valarray<real_t>& delta_energy, real_t peak_E) const
{
    std::valarray<real_t> counts((real_t)0.0, delta_energy.size());
    for (unsigned int i=0; i<delta_energy.size(); i++)
    {
        counts[i] = gain / (real_t)2.0 /  peak_E * Faddeeva::erfc(delta_energy[i]/(M_SQRT2 * sigma));
    }
    return counts;
}

// ----------------------------------------------------------------------------

const std::valarray<real_t> Gaussian_Model::tail(real_t gain, real_t sigma, std::valarray<real_t>& delta_energy, real_t gamma) const
{
    std::valarray<real_t> counts((real_t)0.0, delta_energy.size());

    for (unsigned int i=0; i<delta_energy.size(); i++)
    {
        real_t temp_a = 1.0;
        if (delta_energy[i] < 0.0)
        {
            temp_a = exp(delta_energy[i]/ (gamma * sigma));
        }
        counts[i] = gain / 2. / gamma / sigma / exp(-0.5/pow(gamma, 2.0)) * temp_a * Faddeeva::erfc( delta_energy[i]  /( M_SQRT2*sigma) + (1.0/(gamma*M_SQRT2) )  );
    }
    return counts;
}

// ----------------------------------------------------------------------------

const std::valarray<real_t> Gaussian_Model::elastic_peak(const Fit_Parameters * const fitp, std::valarray<real_t> ev, real_t gain) const
{
    std::valarray<real_t> counts((real_t)0.0, ev.size());
    real_t sigma = std::sqrt( std::pow( (fitp->at(STR_FWHM_OFFSET).value / 2.3548), (real_t)2.0 ) + fitp->at(STR_COHERENT_SCT_ENERGY).value * 2.96 * fitp->at(STR_FWHM_FANOPRIME).value  );
    if(std::isnan(sigma))
    {
        return counts;
    }
    std::valarray<real_t> delta_energy = ev - fitp->at(STR_COHERENT_SCT_ENERGY).value;


    // elastic peak, gaussian
    real_t fvalue = 1.0;

    fvalue = fvalue * std::pow(10.0, fitp->at(STR_COHERENT_SCT_AMPLITUDE).value);

    //std::valarray<real_t> value = fvalue * this->peak(gain, *sigma, delta_energy);
    //counts = counts + value;
    counts += ( fvalue * this->peak(gain, sigma, delta_energy) );

    return counts;
}

// ----------------------------------------------------------------------------

const std::valarray<real_t> Gaussian_Model::compton_peak(const Fit_Parameters * const fitp, std::valarray<real_t> ev, real_t  gain) const
{
    std::valarray<real_t> counts((real_t)0.0, ev.size());

    real_t compton_E = fitp->at(STR_COHERENT_SCT_ENERGY).value/(1. +(fitp->at(STR_COHERENT_SCT_ENERGY).value / 511.0 ) * (1.0 -std::cos( fitp->at(STR_COMPTON_ANGLE).value * 2.0 * M_PI / 360.0 )));

    real_t sigma = std::sqrt( std::pow( (fitp->at(STR_FWHM_OFFSET).value/2.3548), 62.0) + compton_E * 2.96 * fitp->at(STR_FWHM_FANOPRIME).value );
    if(std::isnan(sigma))
    {
        return counts;
    }
    //real_t local_sigma = (*sigma) * p[14];

    std::valarray<real_t> delta_energy = ev - compton_E;

    // compton peak, gaussian
    real_t faktor = 1.0 / (1.0 + fitp->at(STR_COMPTON_F_STEP).value + fitp->at(STR_COMPTON_F_TAIL).value + fitp->at(STR_COMPTON_HI_F_TAIL).value);

    faktor = faktor * std::pow( 10.0, fitp->at(STR_COMPTON_AMPLITUDE).value) ;

    std::valarray<real_t> value = faktor * this->peak(gain, sigma * fitp->at(STR_COMPTON_FWHM_CORR).value, delta_energy);
    counts = counts + value;

    // compton peak, step
    if ( fitp->at(STR_COMPTON_F_STEP).value > 0.0 )
    {
        real_t fvalue = faktor * fitp->at(STR_COMPTON_F_STEP).value;
        value = fvalue * this->step(gain, sigma, delta_energy, compton_E);
        counts = counts + value;
    }
    // compton peak, tail on the low side
    real_t fvalue = faktor * fitp->at(STR_COMPTON_F_TAIL).value;
    value = fvalue * this->tail(gain, sigma, delta_energy, fitp->at(STR_COMPTON_GAMMA).value);
    counts = counts + value;

    // compton peak, tail on the high side
    fvalue = faktor * fitp->at(STR_COMPTON_HI_F_TAIL).value;
    delta_energy *= -1.0;
    value = ( fvalue * this->tail(gain, sigma, delta_energy, fitp->at(STR_COMPTON_HI_GAMMA).value) );
    counts = counts + value;

    return counts;
}

// ----------------------------------------------------------------------------

} //namespace models
} //namespace fitting
