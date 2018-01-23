/***
-Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.
-
-Copyright 2016. UChicago Argonne, LLC. This software was produced
-under U.S. Government contract DE-AC02-06CH11357 for Argonne National
-Laboratory (ANL), which is operated by UChicago Argonne, LLC for the
-U.S. Department of Energy. The U.S. Government has rights to use,
-reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
-UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
-ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
-modified to produce derivative works, such modified software should
-be clearly marked, so as not to confuse it with the version available
-from ANL.
-
-Additionally, redistribution and use in source and binary forms, with
-or without modification, are permitted provided that the following
-conditions are met:
-
-    * Redistributions of source code must retain the above copyright
-      notice, this list of conditions and the following disclaimer.
-
-    * Redistributions in binary form must reproduce the above copyright
-      notice, this list of conditions and the following disclaimer in
-      the documentation and/or other materials provided with the
-      distribution.
-
-    * Neither the name of UChicago Argonne, LLC, Argonne National
-      Laboratory, ANL, the U.S. Government, nor the names of its
-      contributors may be used to endorse or promote products derived
-      from this software without specific prior written permission.
-
-THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS
-"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
-LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
-FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago
-Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
-INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
-BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
-LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
-CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
-LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
-ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
-POSSIBILITY OF SUCH DAMAGE.
-***/

#include "lmfit_optimizer.h"

#include <iostream>
#include <algorithm>
#include <math.h>

#include <string.h>

#include "support/lmfit_6.1/lmmin.hpp"

using namespace data_struct::xrf;


namespace fitting
{
namespace optimizers
{


void residuals_lmfit( const real_t *par, int m_dat, const void *data, real_t *fvec, int *userbreak )
{
    User_Data* ud = (User_Data*)(data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(par, m_dat);
    //Model spectra based on new fit parameters
    update_background_user_data(ud);
    ud->spectra_model = ud->fit_model->model_spectrum(ud->fit_parameters, ud->elements, ud->energy_range);
    ud->spectra_model += ud->spectra_background;
    //Calculate residuals
    for (int i = 0; i < m_dat; i++ )
    {
		fvec[i] = (ud->spectra[i] - ud->spectra_model[i]) * ud->weights[i];
    }

}


void general_residuals_lmfit( const real_t *par, int m_dat, const void *data, real_t *fvec, int *userbreak )
{

    Gen_User_Data* ud = (Gen_User_Data*)(data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(par, m_dat);
    //Model spectra based on new fit parameters
    ud->func(ud->fit_parameters, &(ud->energy_range), &(ud->spectra_model));
    ud->spectra_model += ud->spectra_background;
    //Calculate residuals
    for (int i = 0; i < m_dat; i++ )
    {
        fvec[i] = ( ud->spectra[i] - ud->spectra_model[i] ) * ud->weights[i];
    }

}


//-----------------------------------------------------------------------------

void quantification_residuals( const real_t *par, int m_dat, const void *data, real_t *fvec, int *userbreak )
{
    ///(std::valarray<real_t> p, std::valarray<real_t> y, std::valarray<real_t> x)

    //y is array of elements standards
    //x is indexes of elements in standard
    //p is array size 2 but seems only first index is used
    ///return (y - this->fit_calibrationcurve(x, p));

    Quant_User_Data* ud = (Quant_User_Data*)(data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(par, m_dat);
    //Model spectra based on new fit parameters

    //Calculate residuals
    std::unordered_map<std::string, real_t> result_map = ud->quantification_model->model_calibrationcurve(*(ud->quant_map), par[0]);

    int idx = 0;
    for(auto& itr : *(ud->quant_map))
    {
        fvec[idx] = itr.second.e_cal_ratio - result_map[itr.first];
        idx++;
    }
}


// =====================================================================================================================


LMFit_Optimizer::LMFit_Optimizer() : Optimizer()
{

}

LMFit_Optimizer::~LMFit_Optimizer()
{


}

void LMFit_Optimizer::minimize(Fit_Parameters *fit_params,
                               const Spectra * const spectra,
                               const Fit_Element_Map_Dict * const elements_to_fit,
                               const Base_Model * const model,
                               const Range energy_range)
{

    User_Data ud;

    fill_user_data(ud, fit_params, spectra, elements_to_fit, model, energy_range);

    std::vector<real_t> fitp_arr = fit_params->to_array();
    std::vector<real_t> perror(fitp_arr.size());

    lm_status_struct<real_t> status;

    /* Predefined control parameter sets (msgfile=NULL means stdout).

        */

    lm_control_struct<real_t> control = {LM_USERTOL, LM_USERTOL, LM_USERTOL, LM_USERTOL, 100., 1000, 1, NULL, 0, -1, -1};

//    control.ftol = 1.0e-10;
//    /* Relative error desired in the sum of squares.
//                         Termination occurs when both the actual and
//                         predicted relative reductions in the sum of squares
//                         are at most ftol. */
//    control.xtol = 1.0e-10;
//    /* Relative error between last two approximations.
//                         Termination occurs when the relative error between
//                         two consecutive iterates is at most xtol. */
//    control.gtol = 1.0e-10;
//    /* Orthogonality desired between fvec and its derivs.
//                         Termination occurs when the cosine of the angle
//                         between fvec and any column of the Jacobian is at
//                         most gtol in absolute value. */
//    control.epsilon = 2.2204460e-4;
//    control.epsilon = 1.0e-5;
//    /* Step used to calculate the Jacobian, should be
//                         slightly larger than the relative error in the
//                         user-supplied functions. */
//    control.stepbound = 100.0;
//    /* Used in determining the initial step bound. This
//                         bound is set to the product of stepbound and the
//                         Euclidean norm of diag*x if nonzero, or else to
//                         stepbound itself. In most cases stepbound should lie
//                         in the interval (0.1,100.0). Generally, the value
//                         100.0 is recommended. */
//    control.patience = 1000;
    /* Used to set the maximum number of function evaluations
                         to patience*(number_of_parameters+1). */
    //control.scale_diag;
    /* If 1, the variables will be rescaled internally.
                         Recommended value is 1. */
    //FILE* msgfile;    /* Progress messages will be written to this file. */
    //int verbosity;    /* OR'ed: 1: print some messages; 2: print Jacobian. */
    //int n_maxpri;     /* -1, or max number of parameters to print. */
    //int m_maxpri;     /* -1, or max number of residuals to print. */


    //control.verbosity = 3;

    /* perform the fit */
    //printf( "Fitting:\n" );
    lmmin( fitp_arr.size(), &fitp_arr[0], energy_range.count(), (const void*) &ud, residuals_lmfit, &control, &status );
	printf(".");
    /* print results */
    //printf( "\nResults:\n" );
    //printf( "status after %d function evaluations:\n  %s\n",  status.nfev, lm_infmsg[status.outcome] );
    //if(residuals_count_lmfit < 100)
     //   logit<<"=-=-=-=-=-=-=-=-=-=--==--=-=-==--==---==--==--== bad fit =-=-=-=-=-=-=-=-=-=-=-=-=-=-= "<<residuals_count_lmfit<<"\n";

    //logit<<"residuals count = "<<residuals_count_lmfit<<"\n";
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = static_cast<real_t>(status.nfev);
    }

}


void LMFit_Optimizer::minimize_func(Fit_Parameters *fit_params,
                                    const Spectra * const spectra,
                                    const Range energy_range,
                                    Gen_Func_Def gen_func)
{

    Gen_User_Data ud;

    fill_gen_user_data(ud, fit_params, spectra, energy_range, gen_func);

    std::vector<real_t> fitp_arr = fit_params->to_array();
    std::vector<real_t> perror(fitp_arr.size());

    lm_status_struct<real_t> status;

    lm_control_struct<real_t> control = {LM_USERTOL, LM_USERTOL, LM_USERTOL, LM_USERTOL, 100., 100, 1, NULL, 0, -1, -1};


    lmmin( fitp_arr.size(), &fitp_arr[0], energy_range.count(), (const void*) &ud, general_residuals_lmfit, &control, &status );
    printf(".");
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = status.nfev;
    }

}

void LMFit_Optimizer::minimize_quantification(Fit_Parameters *fit_params,
                                              std::unordered_map<std::string, Element_Quant> * quant_map,
                                              quantification::models::Quantification_Model * quantification_model)
{
    Quant_User_Data ud;

    ud.quant_map = quant_map;
    ud.quantification_model = quantification_model;
    ud.fit_parameters = fit_params;

    std::vector<real_t> fitp_arr = fit_params->to_array();
    std::vector<real_t> perror(fitp_arr.size());

    lm_status_struct<real_t> status;

    lm_control_struct<real_t> control = {LM_USERTOL, LM_USERTOL, LM_USERTOL, LM_USERTOL, 100., 100, 1, NULL, 0, -1, -1};

    lmmin( fitp_arr.size(), &fitp_arr[0], quant_map->size(), (const void*) &ud, quantification_residuals, &control, &status );
    printf(".");
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = static_cast<real_t>(status.nfev);
    }
}

} //namespace optimizers
} //namespace fitting
