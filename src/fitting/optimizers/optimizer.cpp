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



#include "optimizer.h"


namespace fitting
{
namespace optimizers
{

    std::string optimizer_outcome_to_str(OPTIMIZER_OUTCOME outcome)
    {
        switch (outcome)
        {
        case fitting::optimizers::OPTIMIZER_OUTCOME::CONVERGED:
            return "Converged";
        case fitting::optimizers::OPTIMIZER_OUTCOME::F_TOL_LT_TOL:
            return "Ftol less than TOL";
        case fitting::optimizers::OPTIMIZER_OUTCOME::X_TOL_LT_TOL:
            return "XTol less than TOL";
        case fitting::optimizers::OPTIMIZER_OUTCOME::G_TOL_LT_TOL:
            return "GTol less than TOL";
        case fitting::optimizers::OPTIMIZER_OUTCOME::EXHAUSTED:
            return "EXHAUSTED";
        case fitting::optimizers::OPTIMIZER_OUTCOME::CRASHED:
            return "CRASHED";
        case fitting::optimizers::OPTIMIZER_OUTCOME::EXPLODED:
            return "Exploded";
        case fitting::optimizers::OPTIMIZER_OUTCOME::FAILED:
            return "FAILED";
        case fitting::optimizers::OPTIMIZER_OUTCOME::FOUND_NAN:
            return "Found NAN";
        case fitting::optimizers::OPTIMIZER_OUTCOME::FOUND_ZERO:
            return "Found Zero";
        case fitting::optimizers::OPTIMIZER_OUTCOME::STOPPED:
            return "Stopped";
        case fitting::optimizers::OPTIMIZER_OUTCOME::TRAPPED:
            return "Trapped";
        }
        return "";
    }

    void fill_user_data(User_Data& ud,
                        Fit_Parameters* fit_params,
                        const Spectra* const spectra,
                        const Fit_Element_Map_Dict* const elements_to_fit,
                        const Base_Model* const model,
                        const Range energy_range,
                        Callback_Func_Status_Def* status_callback,
                        size_t total_itr,
                        bool use_weights)
	{
		ud.fit_model = (Base_Model*)model;
		// set spectra to fit
        ud.spectra = spectra->sub_spectra(energy_range.min, energy_range.count());
        //not allocating memory. see https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
        //new (&(ud.spectra)) Eigen::Map<const ArrayXr>(spectra->data() + energy_range.min, energy_range.count());
        ud.orig_spectra = spectra;
		ud.fit_parameters = fit_params;
		ud.elements = (Fit_Element_Map_Dict *)elements_to_fit;
        ud.energy_range.min = energy_range.min;
        ud.energy_range.max = energy_range.max;

        ud.status_callback = status_callback;
        ud.cur_itr = 0;
        ud.total_itr = total_itr;

        if (use_weights)
        {
            ArrayXr weights = (real_t)1.0 / ((real_t)1.0 + (*spectra));
            weights = convolve1d(weights, 5);
            weights = Eigen::abs(weights);
            weights /= weights.maxCoeff();
            ud.weights = weights.segment(energy_range.min, energy_range.count());
        }
        else
        {
            ud.weights.resize(energy_range.count());
            ud.weights.fill(1.0);
        }

        ArrayXr background(spectra->size());
        background.setZero(spectra->size());
        if(fit_params->contains(STR_SNIP_WIDTH))
        {
            background = snip_background(spectra,
                                         fit_params->value(STR_ENERGY_OFFSET),
                                         fit_params->value(STR_ENERGY_SLOPE),
                                         fit_params->value(STR_ENERGY_QUADRATIC),
                                         fit_params->value(STR_SNIP_WIDTH),
                                         energy_range.min,
                                         energy_range.max);
        }
        ud.spectra_background = background.segment(energy_range.min, energy_range.count());
        ud.spectra_background = ud.spectra_background.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });
		ud.spectra_model.resize(energy_range.count());
	}

	void fill_gen_user_data(Gen_User_Data &ud,
                            Fit_Parameters *fit_params,
                            const Spectra * const spectra,
                            const Range energy_range,
                            const ArrayXr* background,
                            Gen_Func_Def gen_func,
                            bool use_weights)
	{
		ud.func = gen_func;
		// set spectra to fit
        ud.spectra = spectra->sub_spectra(energy_range.min, energy_range.count());;
        //not allocating memory. see https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
        //new (&ud.spectra) Eigen::Map<const ArrayXr>(spectra->data() + energy_range.min, energy_range.count());
		ud.fit_parameters = fit_params;
        ud.energy_range.min = energy_range.min;
        ud.energy_range.max = energy_range.max;

        if(use_weights)
        {
		    ArrayXr weights = (real_t)1.0 / ((real_t)1.0 + (*spectra));
		    weights = convolve1d(weights, 5);
		    weights = Eigen::abs(weights);
		    weights /= weights.maxCoeff();
		    ud.weights = weights.segment(energy_range.min, energy_range.count());
        }
        else
        {
            ud.weights.resize(energy_range.count());
            ud.weights.fill(1.0);
        }
        ud.spectra_background = *background;
        ud.spectra_background = ud.spectra_background.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });
        
		ud.spectra_model.resize(energy_range.count());
	}

    void update_background_user_data(User_Data *ud)
    {
        if(ud->fit_parameters->contains(STR_SNIP_WIDTH))
        {
            Fit_Param fit_snip_width = ud->fit_parameters->at(STR_SNIP_WIDTH);
            if(fit_snip_width.bound_type != E_Bound_Type::FIXED && ud->orig_spectra != nullptr)
            {
                //ud->spectra_background = snip_background(ud->orig_spectra,
				ArrayXr background = snip_background(ud->orig_spectra,
                                             ud->fit_parameters->value(STR_ENERGY_OFFSET),
                                             ud->fit_parameters->value(STR_ENERGY_SLOPE),
                                             ud->fit_parameters->value(STR_ENERGY_QUADRATIC),
                                             fit_snip_width.value,
                                             ud->energy_range.min,
                                             ud->energy_range.max);

				ud->spectra_background = background.segment(ud->energy_range.min, ud->energy_range.count());
                ud->spectra_background = ud->spectra_background.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });
				
            }
        }


    }

} //namespace optimizers
} //namespace fitting
