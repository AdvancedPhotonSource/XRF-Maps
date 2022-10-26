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

#include <string.h>

#include <omp.h>

#define SQRT_2xPI (T_real)2.506628275 // sqrt ( 2.0 * M_PI )

using namespace data_struct;


namespace fitting
{
namespace models
{

// ----------------------------------------------------------------------------

template<typename T_real>
Gaussian_Model<T_real>::Gaussian_Model() : Base_Model<T_real>()
{
    _fit_parameters = _generate_default_fit_parameters();
}

// ----------------------------------------------------------------------------

template<typename T_real>
Gaussian_Model<T_real>::~Gaussian_Model()
{

}

// ----------------------------------------------------------------------------

template<typename T_real>
Fit_Parameters<T_real> Gaussian_Model<T_real>::_generate_default_fit_parameters()
{

    Fit_Parameters<T_real> fit_params;
    //                                      name                     min         max             val       step              use
    fit_params.add_parameter(Fit_Param<T_real>(STR_MIN_ENERGY_TO_FIT, (T_real)0.1, (T_real)100.00, (T_real)1.0, (T_real)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_MAX_ENERGY_TO_FIT, (T_real)0.1, (T_real)100.00, (T_real)11.0, (T_real)0.1, E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param<T_real>(STR_ENERGY_OFFSET,       (T_real)-0.2,    (T_real)0.2,    (T_real)0.0,    (T_real)0.00001,   E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param<T_real>(STR_ENERGY_SLOPE,        (T_real)0.001,   (T_real)10.0,    (T_real)1.0,    (T_real)0.00001,   E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param<T_real>(STR_ENERGY_QUADRATIC,    (T_real)-0.000000001, (T_real)0.00001, (T_real)0.0,    (T_real)0.000000001,   E_Bound_Type::LIMITED_LO_HI));

    fit_params.add_parameter(Fit_Param<T_real>(STR_FWHM_OFFSET,         (T_real)0.005,    (T_real)0.5,  (T_real)0.12,    (T_real)0.00001,   E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param<T_real>(STR_FWHM_FANOPRIME,      (T_real)0.000001, (T_real)0.05, (T_real)0.00012, (T_real)0.000001,  E_Bound_Type::LIMITED_LO_HI));

    fit_params.add_parameter(Fit_Param<T_real>(STR_COHERENT_SCT_ENERGY,	   (T_real)9.4, (T_real)10.4, (T_real)9.99, (T_real)0.001,  E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param<T_real>(STR_COHERENT_SCT_AMPLITUDE, (T_real)0.000001, (T_real)10.0, (T_real)5.0,  (T_real)0.00001, E_Bound_Type::FIT));

    fit_params.add_parameter(Fit_Param<T_real>(STR_COMPTON_ANGLE,		 (T_real)-0.0001, (T_real)0.0001, (T_real)90.0, (T_real)0.1,       E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param<T_real>(STR_COMPTON_FWHM_CORR,    (T_real)-0.0001, (T_real)0.0001, (T_real)1.0,  (T_real)0.1,       E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param<T_real>(STR_COMPTON_AMPLITUDE,    (T_real)-0.0001, (T_real)0.0001, (T_real)5.0, (T_real)0.00001,  E_Bound_Type::FIT));
    fit_params.add_parameter(Fit_Param<T_real>(STR_COMPTON_F_STEP,		 (T_real)-0.0001, (T_real)0.0001, (T_real)0.0,  (T_real)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_COMPTON_F_TAIL,		 (T_real)-0.0001, (T_real)0.4, (T_real)0.1,  (T_real)0.1,       E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param<T_real>(STR_COMPTON_GAMMA,		 (T_real)-0.0001, (T_real)10., (T_real)1.0,  (T_real)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_COMPTON_HI_F_TAIL, (T_real)-0.0001, (T_real)0.1, (T_real)0.013,  (T_real)0.0000001,    E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param<T_real>(STR_COMPTON_HI_GAMMA,	 (T_real)-0.0001, (T_real)10., (T_real)1.0,  (T_real)0.01,      E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param<T_real>(STR_SNIP_WIDTH,			     (T_real)0.1,   (T_real)2.828427, (T_real)0.25, (T_real)0.01,  E_Bound_Type::FIXED)); //max = 2* sqrt(2)

    fit_params.add_parameter(Fit_Param<T_real>(STR_F_STEP_OFFSET,		(T_real)0.0, (T_real)1.0, (T_real)0.0, (T_real)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_F_STEP_LINEAR,        (T_real)0.0, (T_real)1.0, (T_real)0.0, (T_real)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_F_STEP_QUADRATIC,  (T_real)0.0, (T_real)0.0, (T_real)0.0, (T_real)0.01, E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param<T_real>(STR_F_TAIL_OFFSET,		(T_real)0.0, (T_real)0.1, (T_real)0.003, (T_real)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_F_TAIL_LINEAR,		(T_real)0.0, (T_real)1.0, (T_real)0.001, (T_real)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_F_TAIL_QUADRATIC,  (T_real)0.0, (T_real)0.01, (T_real)0.0, (T_real)0.01, E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param<T_real>(STR_GAMMA_OFFSET,			(T_real)0.1, (T_real)10.0, (T_real)2.2, (T_real)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_GAMMA_LINEAR,			(T_real)0.0, (T_real)3.0,  (T_real)0.0, (T_real)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_GAMMA_QUADRATIC,	(T_real)0.0, (T_real)0.0,  (T_real)0.0, (T_real)0.1, E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param<T_real>(STR_KB_F_TAIL_OFFSET,		 (T_real)0.0, (T_real)0.2,  (T_real)0.05, (T_real)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_KB_F_TAIL_LINEAR,		 (T_real)0.0, (T_real)0.02, (T_real)0.0, (T_real)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param<T_real>(STR_KB_F_TAIL_QUADRATIC, (T_real)0.0, (T_real)0.0,  (T_real)0.0, (T_real)0.01, E_Bound_Type::FIXED));

    return fit_params;

}

// ----------------------------------------------------------------------------

template<typename T_real>
void Gaussian_Model<T_real>::set_fit_params_preset(Fit_Params_Preset preset)
{

    switch (preset)
    {
        case Fit_Params_Preset::MATRIX_BATCH_FIT: //matrix batch fit
            _fit_parameters[STR_ENERGY_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_ENERGY_SLOPE].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_ENERGY_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_FWHM_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_FWHM_FANOPRIME].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_COHERENT_SCT_ENERGY].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COHERENT_SCT_AMPLITUDE].bound_type = E_Bound_Type::FIT;

            _fit_parameters[STR_COMPTON_ANGLE].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_FWHM_CORR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_AMPLITUDE].bound_type = E_Bound_Type::FIT;
            _fit_parameters[STR_COMPTON_F_STEP].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_F_TAIL].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_GAMMA].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_HI_F_TAIL].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_HI_GAMMA].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_SNIP_WIDTH].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_F_STEP_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_STEP_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_STEP_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_F_TAIL_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_TAIL_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_TAIL_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_GAMMA_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_GAMMA_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_GAMMA_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_KB_F_TAIL_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_KB_F_TAIL_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_KB_F_TAIL_QUADRATIC].bound_type = E_Bound_Type::FIXED;
            break;
        case Fit_Params_Preset::BATCH_FIT_NO_TAILS: // batch fit without tails
            _fit_parameters[STR_ENERGY_OFFSET].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_ENERGY_SLOPE].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_ENERGY_QUADRATIC].bound_type = E_Bound_Type::LIMITED_LO_HI;

            _fit_parameters[STR_FWHM_OFFSET].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_FWHM_FANOPRIME].bound_type = E_Bound_Type::LIMITED_LO_HI;

            _fit_parameters[STR_COHERENT_SCT_ENERGY].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COHERENT_SCT_AMPLITUDE].bound_type = E_Bound_Type::FIT;

            _fit_parameters[STR_COMPTON_ANGLE].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COMPTON_FWHM_CORR].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COMPTON_AMPLITUDE].bound_type = E_Bound_Type::FIT;
            _fit_parameters[STR_COMPTON_F_STEP].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_F_TAIL].bound_type = E_Bound_Type::LIMITED_LO;
            _fit_parameters[STR_COMPTON_GAMMA].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_HI_F_TAIL].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COMPTON_HI_GAMMA].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_SNIP_WIDTH].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_F_STEP_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_STEP_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_STEP_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_F_TAIL_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_TAIL_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_TAIL_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_GAMMA_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_GAMMA_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_GAMMA_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_KB_F_TAIL_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_KB_F_TAIL_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_KB_F_TAIL_QUADRATIC].bound_type = E_Bound_Type::FIXED;
            break;
        case Fit_Params_Preset::BATCH_FIT_WITH_TAILS: //batch fit with tails
            _fit_parameters[STR_ENERGY_OFFSET].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_ENERGY_SLOPE].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_ENERGY_QUADRATIC].bound_type = E_Bound_Type::LIMITED_LO_HI;

            _fit_parameters[STR_FWHM_OFFSET].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_FWHM_FANOPRIME].bound_type = E_Bound_Type::LIMITED_LO_HI;

            _fit_parameters[STR_COHERENT_SCT_ENERGY].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COHERENT_SCT_AMPLITUDE].bound_type = E_Bound_Type::FIT;

            _fit_parameters[STR_COMPTON_ANGLE].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COMPTON_FWHM_CORR].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COMPTON_AMPLITUDE].bound_type = E_Bound_Type::FIT;
            _fit_parameters[STR_COMPTON_F_STEP].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_F_TAIL].bound_type = E_Bound_Type::LIMITED_LO;
            _fit_parameters[STR_COMPTON_GAMMA].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COMPTON_HI_F_TAIL].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_COMPTON_HI_GAMMA].bound_type = E_Bound_Type::LIMITED_LO_HI;

            _fit_parameters[STR_SNIP_WIDTH].bound_type = E_Bound_Type::FIT;

            _fit_parameters[STR_F_STEP_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_STEP_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_STEP_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_F_TAIL_OFFSET].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_F_TAIL_LINEAR].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_F_TAIL_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_GAMMA_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_GAMMA_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_GAMMA_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_KB_F_TAIL_OFFSET].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_KB_F_TAIL_LINEAR].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_KB_F_TAIL_QUADRATIC].bound_type = E_Bound_Type::FIXED;
            break;
        case Fit_Params_Preset::BATCH_FIT_WITH_FREE_ENERGY: // batch fit with free E, everything else fixed
            _fit_parameters[STR_ENERGY_OFFSET].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_ENERGY_SLOPE].bound_type = E_Bound_Type::LIMITED_LO_HI;
            _fit_parameters[STR_ENERGY_QUADRATIC].bound_type = E_Bound_Type::LIMITED_LO_HI;

            _fit_parameters[STR_FWHM_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_FWHM_FANOPRIME].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_COHERENT_SCT_ENERGY].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COHERENT_SCT_AMPLITUDE].bound_type = E_Bound_Type::FIT;

            _fit_parameters[STR_COMPTON_ANGLE].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_FWHM_CORR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_AMPLITUDE].bound_type = E_Bound_Type::FIT;
            _fit_parameters[STR_COMPTON_F_STEP].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_F_TAIL].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_GAMMA].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_HI_F_TAIL].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_COMPTON_HI_GAMMA].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_SNIP_WIDTH].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_F_STEP_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_STEP_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_STEP_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_F_TAIL_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_TAIL_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_F_TAIL_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_GAMMA_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_GAMMA_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_GAMMA_QUADRATIC].bound_type = E_Bound_Type::FIXED;

            _fit_parameters[STR_KB_F_TAIL_OFFSET].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_KB_F_TAIL_LINEAR].bound_type = E_Bound_Type::FIXED;
            _fit_parameters[STR_KB_F_TAIL_QUADRATIC].bound_type = E_Bound_Type::FIXED;
            break;
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
const Spectra<T_real> Gaussian_Model<T_real>::model_spectrum(const Fit_Parameters<T_real> * const fit_params,
                                             const unordered_map<string, Fit_Element_Map<T_real>*> * const elements_to_fit,
                                             unordered_map<string, ArrayTr<T_real>>* labeled_spectras,
                                             const struct Range energy_range)
{

    vector<string> spectra_labels = { STR_K_A_LINES, STR_K_B_LINES, STR_L_LINES, STR_M_LINES, STR_STEP_LINES, STR_TAIL_LINES, STR_ELASTIC_LINES, STR_COMPTON_LINES, STR_PILEUP_LINES, STR_ESCAPE_LINES };

    if (labeled_spectras != nullptr) // if this stucture is not null then initialize
    {
        for (auto& itr : spectra_labels)
        {
            if (labeled_spectras->count(itr) > 1)
            {
                labeled_spectras->erase(itr);
            }
            labeled_spectras->insert({ itr, Spectra<T_real>(energy_range.count()) });
        }
    }

    Spectra<T_real> agr_spectra(energy_range.count());
    Spectra<T_real> tmp_spec(energy_range.count());

    T_real energy_offset = fit_params->value(STR_ENERGY_OFFSET);
    T_real energy_slope = fit_params->value(STR_ENERGY_SLOPE);
    T_real energy_quad = fit_params->value(STR_ENERGY_QUADRATIC);

	ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    ArrayTr<T_real> ev = energy_offset + (energy * energy_slope) + (pow(energy, (T_real)2.0) * energy_quad);

    for(const auto& itr : (*elements_to_fit))
    {
        if(itr.first == STR_COHERENT_SCT_AMPLITUDE || itr.first == STR_COMPTON_AMPLITUDE)
        {
            continue;
        }
        else
        {
            agr_spectra += model_spectrum_element(fit_params, itr.second, ev, labeled_spectras);
        }
    }

    if (labeled_spectras != nullptr)
    {
        tmp_spec = elastic_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);
        (*labeled_spectras)[STR_ELASTIC_LINES] += tmp_spec;
        agr_spectra += tmp_spec;
    }
    else
    {
        agr_spectra += elastic_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);
    }

    if (labeled_spectras != nullptr)
    {
        tmp_spec = compton_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);
        (*labeled_spectras)[STR_COMPTON_LINES] += tmp_spec;
        agr_spectra += tmp_spec;
    }
    else
    {
        agr_spectra += compton_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);
    }

 //   agr_spectra += escape_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);

    return agr_spectra;
}

// ----------------------------------------------------------------------------

template<typename T_real>
const Spectra<T_real> Gaussian_Model<T_real>::model_spectrum_mp(const Fit_Parameters<T_real> * const fit_params,
                                                const unordered_map<string, Fit_Element_Map<T_real>*> * const elements_to_fit,
                                                const struct Range energy_range)
{

    Spectra<T_real> agr_spectra(energy_range.count());

    T_real energy_offset = fit_params->value(STR_ENERGY_OFFSET);
    T_real energy_slope = fit_params->value(STR_ENERGY_SLOPE);
    T_real energy_quad = fit_params->value(STR_ENERGY_QUADRATIC);

    ArrayTr<T_real> energy = ArrayTr<T_real>::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    ArrayTr<T_real> ev = energy_offset + (energy * energy_slope) + (pow(energy, (T_real)2.0) * energy_quad);

    std::vector<std::string> keys;
    for (const auto& itr : (*elements_to_fit))
    {
        if(itr.first == STR_COHERENT_SCT_AMPLITUDE || itr.first == STR_COMPTON_AMPLITUDE)
        {
            continue;
        }
        else
        {
            keys.push_back(itr.first);
        }
    }
#pragma omp parallel for
    for (int i=0; i < (int)keys.size(); i++)
    {
        Spectra<T_real> tmp = model_spectrum_element(fit_params, elements_to_fit->at(keys[i]), ev, nullptr);
#pragma omp critical
        {
            agr_spectra += tmp;
        }
    }

    agr_spectra += elastic_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);
    agr_spectra += compton_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);

    //   agr_spectra += escape_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);

    return agr_spectra;
}

// ----------------------------------------------------------------------------

template<typename T_real>
const Spectra<T_real> Gaussian_Model<T_real>::model_spectrum_element(const Fit_Parameters<T_real> * const fitp,
                                                     const Fit_Element_Map<T_real>* const element_to_fit,
                                                     const ArrayTr<T_real> &ev,
                                                     unordered_map<string, ArrayTr<T_real> > * labeled_spectras)
{
    Spectra<T_real> spectra_model(ev.size());

    if(false == fitp->contains(element_to_fit->full_name()))
    {
        return spectra_model;
    }

    T_real pre_faktor = std::pow((T_real)10.0 , fitp->at(element_to_fit->full_name()).value);

    if(false == std::isfinite(pre_faktor))
        return spectra_model;

    //T_real fwhm_offset = fitp->value(STR_FWHM_OFFSET);
    const vector<Element_Energy_Ratio<T_real>> energy_ratios = element_to_fit->energy_ratios();

    //for (const Element_Energy_Ratio& er_struct : element_to_fit->energy_ratios())
    for (int idx = 0; idx < energy_ratios.size(); idx++)
    {
        const Element_Energy_Ratio<T_real>& er_struct = energy_ratios.at(idx);
        T_real sigma = std::sqrt(std::pow((fitp->at(STR_FWHM_OFFSET).value / (T_real)2.3548), (T_real)2.0) + (er_struct.energy) * (T_real)2.96 * fitp->at(STR_FWHM_FANOPRIME).value);
        T_real f_step =  std::abs<T_real>( er_struct.mu_fraction * ( fitp->at(STR_F_STEP_OFFSET).value + (fitp->at(STR_F_STEP_LINEAR).value * er_struct.energy)));
        T_real f_tail = std::abs<T_real>( fitp->at(STR_F_TAIL_OFFSET).value + (fitp->at(STR_F_TAIL_LINEAR).value * er_struct.mu_fraction));
        T_real kb_f_tail = std::abs<T_real>(  fitp->at(STR_KB_F_TAIL_OFFSET).value + (fitp->at(STR_KB_F_TAIL_LINEAR).value * er_struct.mu_fraction));
        T_real value = 1.0;

        //don't process if energy is 0
        if (er_struct.ratio == 0.0)
            continue;
        if (er_struct.energy <= 0.0)
            continue;

        // gaussian peak shape
		ArrayTr<T_real> delta_energy = ev - er_struct.energy;

        string label = "";

        T_real incident_energy = fitp->at(STR_COHERENT_SCT_ENERGY).value;

        T_real faktor = T_real(er_struct.ratio * pre_faktor);
		if (element_to_fit->check_binding_energy(incident_energy, idx))
		{
			switch (er_struct.ptype)
			{
			case Element_Param_Type::Kb1_Line:
				label = STR_K_A_LINES;
			case Element_Param_Type::Kb2_Line:
				label = STR_K_B_LINES;
				faktor = faktor / ((T_real)1.0 + kb_f_tail + f_step);
				break;
			case Element_Param_Type::Ka1_Line:
			case Element_Param_Type::Ka2_Line:
				label = STR_K_A_LINES;
				faktor = faktor / ((T_real)1.0 + f_tail + f_step);
				break;
			case Element_Param_Type::La1_Line:
			case Element_Param_Type::La2_Line:
			case Element_Param_Type::Lb1_Line:
			case Element_Param_Type::Lb2_Line:
			case Element_Param_Type::Lb3_Line:
			case Element_Param_Type::Lb4_Line:
			case Element_Param_Type::Lg1_Line:
			case Element_Param_Type::Lg2_Line:
			case Element_Param_Type::Lg3_Line:
			case Element_Param_Type::Lg4_Line:
			case Element_Param_Type::Ll_Line:
			case Element_Param_Type::Ln_Line:
				label = STR_L_LINES;
				faktor = faktor / ((T_real)1.0 + f_tail + f_step);
				break;
			default:
				break;
			}
		}
		else
		{
			faktor = (T_real)0.0;
		}


        if (labeled_spectras != nullptr && label.length() > 0)
        {
            Spectra<T_real> tmp_spec(ev.size());
            // peak, gauss
            tmp_spec += faktor * this->peak(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy);
            ////spectra_model += faktor * (fitp->at(STR_ENERGY_SLOPE).value / ( sigma * SQRT_2xPI ) *  Eigen::exp((T_real)-0.5 * Eigen::pow((delta_energy / sigma), (T_real)2.0) ) );

            //  peak, step
            if (f_step > 0.0)
            {
                value = faktor * f_step;
                //value = value * this->step(gain, sigma, delta_energy, er_struct.energy);
                tmp_spec += value * this->step(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy, er_struct.energy);
                //counts_arr->step = fit_counts.step + value;
            }
            //  peak, tail;; use different tail for K beta vs K alpha lines
            if (er_struct.ptype == Element_Param_Type::Kb1_Line || er_struct.ptype == Element_Param_Type::Kb2_Line)
            {
                T_real gamma = std::abs(fitp->at(STR_GAMMA_OFFSET).value + fitp->at(STR_GAMMA_LINEAR).value * (er_struct.energy)) * element_to_fit->width_multi();
                value = faktor * kb_f_tail;
                tmp_spec += value * this->tail(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy, gamma);
                //fit_counts.tail = fit_counts.tail + value;
            }

            if (element_to_fit->pileup_element() != nullptr) // check if it is pileup 
            {
                (*labeled_spectras)[STR_PILEUP_LINES] += tmp_spec;
            }
            else
            {
                (*labeled_spectras)[label] += tmp_spec;
            }
            spectra_model += tmp_spec;

        }
        else
        {
            // peak, gauss
            spectra_model += faktor * this->peak(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy);
            ////spectra_model += faktor * (fitp->at(STR_ENERGY_SLOPE).value / ( sigma * SQRT_2xPI ) *  Eigen::exp((T_real)-0.5 * Eigen::pow((delta_energy / sigma), (T_real)2.0) ) );

            //  peak, step
            if (f_step > 0.0)
            {
                value = faktor * f_step;
                //value = value * this->step(gain, sigma, delta_energy, er_struct.energy);
                spectra_model += value * this->step(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy, er_struct.energy);
                //counts_arr->step = fit_counts.step + value;
            }
            //  peak, tail;; use different tail for K beta vs K alpha lines
            if (er_struct.ptype == Element_Param_Type::Kb1_Line || er_struct.ptype == Element_Param_Type::Kb2_Line)
            {
                T_real gamma = std::abs(fitp->at(STR_GAMMA_OFFSET).value + fitp->at(STR_GAMMA_LINEAR).value * (er_struct.energy)) * element_to_fit->width_multi();
                value = faktor * kb_f_tail;
                spectra_model += value * this->tail(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy, gamma);
                //fit_counts.tail = fit_counts.tail + value;
            }
        }
    }
    return spectra_model;
}

// ----------------------------------------------------------------------------

template<typename T_real>
const ArrayTr<T_real> Gaussian_Model<T_real>::peak(T_real gain, T_real sigma, const ArrayTr<T_real>& delta_energy) const
{
    // gain / (sigma * sqrt( 2.0 * M_PI) ) * exp( -0.5 * ( (delta_energy / sigma) ** 2 )
    return gain / ( sigma * (T_real)(SQRT_2xPI) ) *  Eigen::exp((T_real)-0.5 * Eigen::pow((delta_energy / sigma), (T_real)2.0) );
}

// ----------------------------------------------------------------------------

template<typename T_real>
const ArrayTr<T_real> Gaussian_Model<T_real>::step(T_real gain, T_real sigma, const ArrayTr<T_real>& delta_energy, T_real peak_E) const
{
    //delta_energy = delta_energy.unaryExpr([sigma](T_real v) { return  (T_real)std::erfc(v  /( M_SQRT2*sigma)); } );
    //return (gain / (T_real)2.0 / peak_E * delta_energy );
    ArrayTr<T_real>counts(delta_energy.size());
	for (unsigned int i = 0; i < delta_energy.size(); i++)
	{
        counts[i] = gain / (T_real)2.0 / peak_E * std::erfc((T_real)delta_energy[i] / ((T_real)(M_SQRT2) * sigma));
	}
    return counts;

}

// ----------------------------------------------------------------------------

template<typename T_real>
const ArrayTr<T_real> Gaussian_Model<T_real>::tail(T_real gain, T_real sigma, ArrayTr<T_real> delta_energy, T_real gamma) const
{
    delta_energy = delta_energy.unaryExpr([sigma,gamma](T_real v) { return  (v < (T_real)0.0) ? std::exp(v/ (gamma * sigma)) * std::erfc(v / ((T_real)(M_SQRT2)*sigma) + ((T_real)1.0/(gamma*(T_real)(M_SQRT2)))) : std::erfc(v / ((T_real)(M_SQRT2)*sigma) + ((T_real)1.0/(gamma*(T_real)(M_SQRT2)))); } );
    return( gain / (T_real)2.0 / gamma / sigma / exp((T_real)-0.5/pow(gamma, (T_real)2.0)) * delta_energy);
}

// ----------------------------------------------------------------------------

template<typename T_real>
const ArrayTr<T_real> Gaussian_Model<T_real>::elastic_peak(const Fit_Parameters<T_real> * const fitp, const ArrayTr<T_real>& ev, T_real gain) const
{
    Spectra<T_real> counts(ev.size());
	counts.setZero();
    T_real sigma = std::sqrt( std::pow( (fitp->at(STR_FWHM_OFFSET).value / (T_real)2.3548), (T_real)2.0 ) + fitp->at(STR_COHERENT_SCT_ENERGY).value * (T_real)2.96 * fitp->at(STR_FWHM_FANOPRIME).value  );
    if(false == std::isfinite(sigma))
    {
        return counts;
    }
	ArrayTr<T_real>delta_energy = ev - fitp->at(STR_COHERENT_SCT_ENERGY).value;


    // elastic peak, gaussian
    T_real fvalue = (T_real)1.0;

    fvalue = fvalue * std::pow((T_real)10.0, fitp->at(STR_COHERENT_SCT_AMPLITUDE).value);

    //Spectra value = fvalue * this->peak(gain, *sigma, delta_energy);
    //counts = counts + value;
    counts += ( fvalue * this->peak(gain, sigma, delta_energy) );
    ////counts += fvalue * (gain / ( sigma * (T_real)(SQRT_2xPI) ) * Eigen::exp((T_real)-0.5 * Eigen::pow((delta_energy / sigma), (T_real)2.0) ) );

    return counts;
}

// ----------------------------------------------------------------------------

template<typename T_real>
const ArrayTr<T_real> Gaussian_Model<T_real>::compton_peak(const Fit_Parameters<T_real> * const fitp, const ArrayTr<T_real>& ev, T_real  gain) const
{
	ArrayTr<T_real>counts(ev.size());
	counts.setZero();

    T_real compton_E = fitp->at(STR_COHERENT_SCT_ENERGY).value/((T_real)1.0 +(fitp->at(STR_COHERENT_SCT_ENERGY).value / (T_real)511.0 ) * ((T_real)1.0 -std::cos( fitp->at(STR_COMPTON_ANGLE).value * (T_real)2.0 * (T_real)(M_PI) / (T_real)360.0 )));

    T_real sigma = std::sqrt( std::pow( (fitp->at(STR_FWHM_OFFSET).value/(T_real)2.3548), (T_real)62.0) + compton_E * (T_real)2.96 * fitp->at(STR_FWHM_FANOPRIME).value );
    if(false == std::isfinite(sigma))
    {
        return counts;
    }
    //T_real local_sigma = (*sigma) * p[14];

	ArrayTr<T_real>delta_energy = ev - compton_E;

    // compton peak, gaussian
    T_real faktor = (T_real)1.0 / ((T_real)1.0 + fitp->at(STR_COMPTON_F_STEP).value + fitp->at(STR_COMPTON_F_TAIL).value + fitp->at(STR_COMPTON_HI_F_TAIL).value);

    faktor = faktor * std::pow((T_real)10.0, fitp->at(STR_COMPTON_AMPLITUDE).value) ;

    counts += faktor * this->peak(gain, sigma * fitp->at(STR_COMPTON_FWHM_CORR).value, delta_energy);
    ////counts += faktor * (gain / ( (sigma * fitp->at(STR_COMPTON_FWHM_CORR).value) * (T_real)(SQRT_2xPI) ) *  Eigen::exp((T_real)-0.5 * Eigen::pow((delta_energy / (sigma*fitp->at(STR_COMPTON_FWHM_CORR).value)), (T_real)2.0) ) );

    // compton peak, step
    if ( fitp->at(STR_COMPTON_F_STEP).value > 0.0 )
    {
        T_real fvalue = faktor * fitp->at(STR_COMPTON_F_STEP).value;
		counts += fvalue * this->step(gain, sigma, delta_energy, compton_E);
    }
    // compton peak, tail on the low side
    T_real fvalue = faktor * fitp->at(STR_COMPTON_F_TAIL).value;
    counts += fvalue * this->tail(gain, sigma, delta_energy, fitp->at(STR_COMPTON_GAMMA).value);

    // compton peak, tail on the high side
    fvalue = faktor * fitp->at(STR_COMPTON_HI_F_TAIL).value;
    delta_energy *= (T_real)-1.0;
    counts += ( fvalue * this->tail(gain, sigma, delta_energy, fitp->at(STR_COMPTON_HI_GAMMA).value) );
    return counts;
}

// ----------------------------------------------------------------------------

template<typename T_real>
const ArrayTr<T_real> Gaussian_Model<T_real>::escape_peak(const Fit_Parameters<T_real>* const fitp, const ArrayTr<T_real>& ev, T_real  gain) const
{
    // check detector element for how much to move escape peaks.
    ArrayTr<T_real>counts(ev.size());
    counts.setZero();
    /*
    //escape
    //if (np.sum(np.abs(pall[keywords.added_params[1:4]])) >= 0.0)
    {
        // si escape
        if (fit_params->at(STR_SI_ESCAPE).value > 0.0)
            //if (pall[keywords.added_params[1]] > 0.0)
        {
            T_real escape_E = (T_real)1.73998;
            for (size_t i = 0; i < ev.size(); i++)
            {
                if (ev[i] > (escape_E + ev[0]))
                {
                    escape_factor = pall[keywords.added_params[1]] + pall[keywords.added_params[3]] * ev;

                    for (int ii = 0; ii < len(wo[0]); ii++)
                    {
                        fit_counts.escape[ii] = counts[wo[0][ii]] * std::max(np.append(escape_factor[wo[0][ii]], 0.0));
                    }
                }
            }
        }
    }
    */
    return counts;
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

template<typename T_real>
ArrayTr<T_real> generate_ev_array(Range energy_range, Fit_Parameters<T_real>& fit_params)
{
    T_real energy_offset = 0.0;
    T_real energy_slope = 0.0;
    T_real energy_quad = 0.0;
    if (fit_params.contains(STR_ENERGY_OFFSET))
    {
        energy_offset = fit_params.at(STR_ENERGY_OFFSET).value;
    }
    if (fit_params.contains(STR_ENERGY_SLOPE))
    {
        energy_slope = fit_params.at(STR_ENERGY_SLOPE).value;
    }
    if (fit_params.contains(STR_ENERGY_QUADRATIC))
    {
        energy_quad = fit_params.at(STR_ENERGY_QUADRATIC).value;
    }


    return generate_ev_array(energy_range, energy_offset, energy_slope, energy_quad);
}

// ----------------------------------------------------------------------------

template<typename T_real>
ArrayTr<T_real> generate_ev_array(Range energy_range, T_real energy_offset, T_real energy_slope, T_real energy_quad)
{
    ArrayTr<T_real>energy = data_struct::ArrayTr<T_real>::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    ArrayTr<T_real>ev = energy_offset + (energy * energy_slope) + (Eigen::pow(energy, (T_real)2.0) * energy_quad);
    return ev;
}

// ----------------------------------------------------------------------------

} //namespace models
} //namespace fitting
