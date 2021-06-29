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

#define SQRT_2xPI (real_t)2.506628275 // sqrt ( 2.0 * M_PI )

using namespace data_struct;


namespace fitting
{
namespace models
{

// ----------------------------------------------------------------------------

Gaussian_Model::Gaussian_Model() : Base_Model()
{
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
    //                                      name                     min         max             val       step              use
    fit_params.add_parameter(Fit_Param(STR_MIN_ENERGY_TO_FIT, (real_t)0.1, (real_t)100.00, (real_t)1.0, (real_t)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_MAX_ENERGY_TO_FIT, (real_t)0.1, (real_t)100.00, (real_t)11.0, (real_t)0.1, E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param(STR_ENERGY_OFFSET,       (real_t)-0.2,    (real_t)0.2,    (real_t)0.0,    (real_t)0.00001,   E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param(STR_ENERGY_SLOPE,        (real_t)0.001,   (real_t)10.0,    (real_t)1.0,    (real_t)0.00001,   E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param(STR_ENERGY_QUADRATIC,    (real_t)-0.00001, (real_t)0.00001, (real_t)0.0,    (real_t)0.000000001,   E_Bound_Type::LIMITED_LO_HI));

    fit_params.add_parameter(Fit_Param(STR_FWHM_OFFSET,         (real_t)0.005,    (real_t)0.5,  (real_t)0.12,    (real_t)0.00001,   E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param(STR_FWHM_FANOPRIME,      (real_t)0.000001, (real_t)0.05, (real_t)0.00012, (real_t)0.000001,  E_Bound_Type::LIMITED_LO_HI));

    fit_params.add_parameter(Fit_Param(STR_COHERENT_SCT_ENERGY,	   (real_t)9.4, (real_t)10.4, (real_t)9.99, (real_t)0.001,  E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param(STR_COHERENT_SCT_AMPLITUDE, (real_t)0.000001, (real_t)10.0, (real_t)5.0,  (real_t)0.00001, E_Bound_Type::FIT));

    fit_params.add_parameter(Fit_Param(STR_COMPTON_ANGLE,		 (real_t)-0.0001, (real_t)0.0001, (real_t)90.0, (real_t)0.1,       E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param(STR_COMPTON_FWHM_CORR,    (real_t)-0.0001, (real_t)0.0001, (real_t)1.0,  (real_t)0.1,       E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param(STR_COMPTON_AMPLITUDE,    (real_t)-0.0001, (real_t)0.0001, (real_t)5.0, (real_t)0.00001,  E_Bound_Type::FIT));
    fit_params.add_parameter(Fit_Param(STR_COMPTON_F_STEP,		 (real_t)-0.0001, (real_t)0.0001, (real_t)0.0,  (real_t)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_COMPTON_F_TAIL,		 (real_t)-0.0001, (real_t)0.4, (real_t)0.1,  (real_t)0.1,       E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param(STR_COMPTON_GAMMA,		 (real_t)-0.0001, (real_t)10., (real_t)1.0,  (real_t)0.1,       E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_COMPTON_HI_F_TAIL, (real_t)-0.0001, (real_t)0.1, (real_t)0.013,  (real_t)0.0000001,    E_Bound_Type::LIMITED_LO_HI));
    fit_params.add_parameter(Fit_Param(STR_COMPTON_HI_GAMMA,	 (real_t)-0.0001, (real_t)10., (real_t)1.0,  (real_t)0.01,      E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param(STR_SNIP_WIDTH,			     (real_t)0.1,   (real_t)2.828427, (real_t)0.25, (real_t)0.01,  E_Bound_Type::FIXED)); //max = 2* sqrt(2)

    fit_params.add_parameter(Fit_Param(STR_F_STEP_OFFSET,		(real_t)0.0, (real_t)1.0, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_F_STEP_LINEAR,        (real_t)0.0, (real_t)1.0, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_F_STEP_QUADRATIC,  (real_t)0.0, (real_t)0.0, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param(STR_F_TAIL_OFFSET,		(real_t)0.0, (real_t)0.1, (real_t)0.003, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_F_TAIL_LINEAR,		(real_t)0.0, (real_t)1.0, (real_t)0.001, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_F_TAIL_QUADRATIC,  (real_t)0.0, (real_t)0.01, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param(STR_GAMMA_OFFSET,			(real_t)0.1, (real_t)10.0, (real_t)2.2, (real_t)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_GAMMA_LINEAR,			(real_t)0.0, (real_t)3.0,  (real_t)0.0, (real_t)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_GAMMA_QUADRATIC,	(real_t)0.0, (real_t)0.0,  (real_t)0.0, (real_t)0.1, E_Bound_Type::FIXED));

    fit_params.add_parameter(Fit_Param(STR_KB_F_TAIL_OFFSET,		 (real_t)0.0, (real_t)0.2,  (real_t)0.05, (real_t)0.1, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_KB_F_TAIL_LINEAR,		 (real_t)0.0, (real_t)0.02, (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));
    fit_params.add_parameter(Fit_Param(STR_KB_F_TAIL_QUADRATIC, (real_t)0.0, (real_t)0.0,  (real_t)0.0, (real_t)0.01, E_Bound_Type::FIXED));

    return fit_params;

}

// ----------------------------------------------------------------------------

void Gaussian_Model::set_fit_params_preset(Fit_Params_Preset preset)
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

const Spectra Gaussian_Model::model_spectrum(const Fit_Parameters * const fit_params,
                                             const unordered_map<string, Fit_Element_Map*> * const elements_to_fit,
                                             unordered_map<string, ArrayXr>* labeled_spectras,
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
            labeled_spectras->insert({ itr, Spectra(energy_range.count()) });
        }
    }

    Spectra agr_spectra(energy_range.count());
    Spectra tmp_spec(energy_range.count());

    real_t energy_offset = fit_params->value(STR_ENERGY_OFFSET);
    real_t energy_slope = fit_params->value(STR_ENERGY_SLOPE);
    real_t energy_quad = fit_params->value(STR_ENERGY_QUADRATIC);

	ArrayXr energy = ArrayXr::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    ArrayXr ev = energy_offset + (energy * energy_slope) + (pow(energy, (real_t)2.0) * energy_quad);

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

const Spectra Gaussian_Model::model_spectrum_mp(const Fit_Parameters * const fit_params,
                                                const unordered_map<string, Fit_Element_Map*> * const elements_to_fit,
                                                const struct Range energy_range)
{

    Spectra agr_spectra(energy_range.count());

    real_t energy_offset = fit_params->value(STR_ENERGY_OFFSET);
    real_t energy_slope = fit_params->value(STR_ENERGY_SLOPE);
    real_t energy_quad = fit_params->value(STR_ENERGY_QUADRATIC);

    ArrayXr energy = ArrayXr::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    ArrayXr ev = energy_offset + (energy * energy_slope) + (pow(energy, (real_t)2.0) * energy_quad);

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
        Spectra tmp = model_spectrum_element(fit_params, elements_to_fit->at(keys[i]), ev, nullptr);
#pragma omp critical
        {
            agr_spectra += tmp;
        }
    }

    agr_spectra += elastic_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);
    agr_spectra += compton_peak(fit_params, ev, fit_params->at(STR_ENERGY_SLOPE).value);

    return agr_spectra;
}

// ----------------------------------------------------------------------------

const Spectra Gaussian_Model::model_spectrum_element(const Fit_Parameters * const fitp,
                                                     const Fit_Element_Map * const element_to_fit,
                                                     const ArrayXr &ev,
                                                     unordered_map<string, ArrayXr>* labeled_spectras)
{
    Spectra spectra_model(ev.size());

    if(false == fitp->contains(element_to_fit->full_name()))
    {
        return spectra_model;
    }

    real_t pre_faktor = std::pow((real_t)10.0 , fitp->at(element_to_fit->full_name()).value);

    if(false == std::isfinite(pre_faktor))
        return spectra_model;

    //real_t fwhm_offset = fitp->value(STR_FWHM_OFFSET);
    const vector<Element_Energy_Ratio> energy_ratios = element_to_fit->energy_ratios();

    //for (const Element_Energy_Ratio& er_struct : element_to_fit->energy_ratios())
    for (int idx = 0; idx < energy_ratios.size(); idx++)
    {
        const Element_Energy_Ratio& er_struct = energy_ratios.at(idx);
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
		ArrayXr delta_energy = ev - er_struct.energy;

        string label = "";

        real_t incident_energy = fitp->at(STR_COHERENT_SCT_ENERGY).value;

        real_t faktor = real_t(er_struct.ratio * pre_faktor);
		if (element_to_fit->check_binding_energy(incident_energy, idx))
		{
			switch (er_struct.ptype)
			{
			case Element_Param_Type::Kb1_Line:
				label = STR_K_A_LINES;
			case Element_Param_Type::Kb2_Line:
				label = STR_K_B_LINES;
				faktor = faktor / ((real_t)1.0 + kb_f_tail + f_step);
				break;
			case Element_Param_Type::Ka1_Line:
			case Element_Param_Type::Ka2_Line:
				label = STR_K_A_LINES;
				faktor = faktor / ((real_t)1.0 + f_tail + f_step);
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
				faktor = faktor / ((real_t)1.0 + f_tail + f_step);
				break;
			default:
				break;
			}
		}
		else
		{
			faktor = (real_t)0.0;
		}


        if (labeled_spectras != nullptr && label.length() > 0)
        {
            Spectra tmp_spec(ev.size());
            // peak, gauss
            tmp_spec += faktor * this->peak(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy);
            ////spectra_model += faktor * (fitp->at(STR_ENERGY_SLOPE).value / ( sigma * SQRT_2xPI ) *  Eigen::exp((real_t)-0.5 * Eigen::pow((delta_energy / sigma), (real_t)2.0) ) );

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
                real_t gamma = std::abs(fitp->at(STR_GAMMA_OFFSET).value + fitp->at(STR_GAMMA_LINEAR).value * (er_struct.energy)) * element_to_fit->width_multi();
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
            ////spectra_model += faktor * (fitp->at(STR_ENERGY_SLOPE).value / ( sigma * SQRT_2xPI ) *  Eigen::exp((real_t)-0.5 * Eigen::pow((delta_energy / sigma), (real_t)2.0) ) );

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
                real_t gamma = std::abs(fitp->at(STR_GAMMA_OFFSET).value + fitp->at(STR_GAMMA_LINEAR).value * (er_struct.energy)) * element_to_fit->width_multi();
                value = faktor * kb_f_tail;
                spectra_model += value * this->tail(fitp->at(STR_ENERGY_SLOPE).value, sigma, delta_energy, gamma);
                //fit_counts.tail = fit_counts.tail + value;
            }
        }
    }
    return spectra_model;
}

// ----------------------------------------------------------------------------

const ArrayXr Gaussian_Model::peak(real_t gain, real_t sigma, const ArrayXr& delta_energy) const
{
    // gain / (sigma * sqrt( 2.0 * M_PI) ) * exp( -0.5 * ( (delta_energy / sigma) ** 2 )
    return gain / ( sigma * (real_t)(SQRT_2xPI) ) *  Eigen::exp((real_t)-0.5 * Eigen::pow((delta_energy / sigma), (real_t)2.0) );
}

// ----------------------------------------------------------------------------

const ArrayXr Gaussian_Model::step(real_t gain, real_t sigma, const ArrayXr& delta_energy, real_t peak_E) const
{
    //delta_energy = delta_energy.unaryExpr([sigma](real_t v) { return  (real_t)std::erfc(v  /( M_SQRT2*sigma)); } );
    //return (gain / (real_t)2.0 / peak_E * delta_energy );
    ArrayXr counts(delta_energy.size());
	for (unsigned int i = 0; i < delta_energy.size(); i++)
	{
        counts[i] = gain / (real_t)2.0 / peak_E * std::erfc((real_t)delta_energy[i] / ((real_t)(M_SQRT2) * sigma));
	}
    return counts;

}

// ----------------------------------------------------------------------------

const ArrayXr Gaussian_Model::tail(real_t gain, real_t sigma, ArrayXr delta_energy, real_t gamma) const
{
    delta_energy = delta_energy.unaryExpr([sigma,gamma](real_t v) { return  (v < (real_t)0.0) ? std::exp(v/ (gamma * sigma)) * std::erfc(v / ((real_t)(M_SQRT2)*sigma) + ((real_t)1.0/(gamma*(real_t)(M_SQRT2)))) : std::erfc(v / ((real_t)(M_SQRT2)*sigma) + ((real_t)1.0/(gamma*(real_t)(M_SQRT2)))); } );
    return( gain / (real_t)2.0 / gamma / sigma / exp((real_t)-0.5/pow(gamma, (real_t)2.0)) * delta_energy);
}

// ----------------------------------------------------------------------------

const ArrayXr Gaussian_Model::elastic_peak(const Fit_Parameters * const fitp, const ArrayXr& ev, real_t gain) const
{
    Spectra counts(ev.size());
	counts.setZero();
    real_t sigma = std::sqrt( std::pow( (fitp->at(STR_FWHM_OFFSET).value / (real_t)2.3548), (real_t)2.0 ) + fitp->at(STR_COHERENT_SCT_ENERGY).value * (real_t)2.96 * fitp->at(STR_FWHM_FANOPRIME).value  );
    if(false == std::isfinite(sigma))
    {
        return counts;
    }
	ArrayXr delta_energy = ev - fitp->at(STR_COHERENT_SCT_ENERGY).value;


    // elastic peak, gaussian
    real_t fvalue = (real_t)1.0;

    fvalue = fvalue * std::pow((real_t)10.0, fitp->at(STR_COHERENT_SCT_AMPLITUDE).value);

    //Spectra value = fvalue * this->peak(gain, *sigma, delta_energy);
    //counts = counts + value;
    counts += ( fvalue * this->peak(gain, sigma, delta_energy) );
    ////counts += fvalue * (gain / ( sigma * (real_t)(SQRT_2xPI) ) * Eigen::exp((real_t)-0.5 * Eigen::pow((delta_energy / sigma), (real_t)2.0) ) );

    return counts;
}

// ----------------------------------------------------------------------------

const ArrayXr Gaussian_Model::compton_peak(const Fit_Parameters * const fitp, const ArrayXr& ev, real_t  gain) const
{
	ArrayXr counts(ev.size());
	counts.setZero();

    real_t compton_E = fitp->at(STR_COHERENT_SCT_ENERGY).value/((real_t)1.0 +(fitp->at(STR_COHERENT_SCT_ENERGY).value / (real_t)511.0 ) * ((real_t)1.0 -std::cos( fitp->at(STR_COMPTON_ANGLE).value * (real_t)2.0 * (real_t)(M_PI) / (real_t)360.0 )));

    real_t sigma = std::sqrt( std::pow( (fitp->at(STR_FWHM_OFFSET).value/(real_t)2.3548), (real_t)62.0) + compton_E * (real_t)2.96 * fitp->at(STR_FWHM_FANOPRIME).value );
    if(false == std::isfinite(sigma))
    {
        return counts;
    }
    //real_t local_sigma = (*sigma) * p[14];

	ArrayXr delta_energy = ev - compton_E;

    // compton peak, gaussian
    real_t faktor = (real_t)1.0 / ((real_t)1.0 + fitp->at(STR_COMPTON_F_STEP).value + fitp->at(STR_COMPTON_F_TAIL).value + fitp->at(STR_COMPTON_HI_F_TAIL).value);

    faktor = faktor * std::pow((real_t)10.0, fitp->at(STR_COMPTON_AMPLITUDE).value) ;

    counts += faktor * this->peak(gain, sigma * fitp->at(STR_COMPTON_FWHM_CORR).value, delta_energy);
    ////counts += faktor * (gain / ( (sigma * fitp->at(STR_COMPTON_FWHM_CORR).value) * (real_t)(SQRT_2xPI) ) *  Eigen::exp((real_t)-0.5 * Eigen::pow((delta_energy / (sigma*fitp->at(STR_COMPTON_FWHM_CORR).value)), (real_t)2.0) ) );

    // compton peak, step
    if ( fitp->at(STR_COMPTON_F_STEP).value > 0.0 )
    {
        real_t fvalue = faktor * fitp->at(STR_COMPTON_F_STEP).value;
		counts += fvalue * this->step(gain, sigma, delta_energy, compton_E);
    }
    // compton peak, tail on the low side
    real_t fvalue = faktor * fitp->at(STR_COMPTON_F_TAIL).value;
    counts += fvalue * this->tail(gain, sigma, delta_energy, fitp->at(STR_COMPTON_GAMMA).value);

    // compton peak, tail on the high side
    fvalue = faktor * fitp->at(STR_COMPTON_HI_F_TAIL).value;
    delta_energy *= (real_t)-1.0;
    counts += ( fvalue * this->tail(gain, sigma, delta_energy, fitp->at(STR_COMPTON_HI_GAMMA).value) );
    return counts;
}

// ----------------------------------------------------------------------------

const ArrayXr Gaussian_Model::escape_peak(const Fit_Parameters* const fitp, const ArrayXr& ev, real_t  gain) const
{
    ArrayXr counts(ev.size());
    counts.setZero();
    /*
    //escape
    //if (np.sum(np.abs(pall[keywords.added_params[1:4]])) >= 0.0)
    {
        // si escape
        if (fit_params->at(STR_SI_ESCAPE).value > 0.0)
            //if (pall[keywords.added_params[1]] > 0.0)
        {
            real_t escape_E = (real_t)1.73998;
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

ArrayXr generate_ev_array(Range energy_range, Fit_Parameters& fit_params)
{
    real_t energy_offset = 0.0;
    real_t energy_slope = 0.0;
    real_t energy_quad = 0.0;
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

ArrayXr generate_ev_array(Range energy_range, real_t energy_offset, real_t energy_slope, real_t energy_quad)
{
    data_struct::ArrayXr energy = data_struct::ArrayXr::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
    data_struct::ArrayXr ev = energy_offset + (energy * energy_slope) + (Eigen::pow(energy, (real_t)2.0) * energy_quad);
    return ev;
}

// ----------------------------------------------------------------------------

} //namespace models
} //namespace fitting
