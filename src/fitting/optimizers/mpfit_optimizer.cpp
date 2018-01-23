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



#include "mpfit_optimizer.h"

#include <iostream>
#include <algorithm>
#include <math.h>

#include "support/cmpfit-1.3a/mpfit.hpp"
#include <string.h>

using namespace data_struct::xrf;


namespace fitting
{
namespace optimizers
{

int residuals_mpfit(int m, int params_size, real_t *params, real_t *dy, real_t **dvec, void *usr_data)
{
    //Get user passed data
    User_Data* ud = static_cast<User_Data*>(usr_data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(params, params_size);

    update_background_user_data(ud);

    //Model spectra based on new fit parameters
    ud->spectra_model = ud->fit_model->model_spectrum(ud->fit_parameters, ud->elements, ud->energy_range);
    ud->spectra_model += ud->spectra_background;

    //Calculate residuals
    for (int i=0; i<m; i++)
    {
		dy[i] = (ud->spectra[i] - ud->spectra_model[i]) * ud->weights[i];
    }
	
	return 0;
}

int gen_residuals_mpfit(int m, int params_size, real_t *params, real_t *dy, real_t **dvec, void *usr_data)
{
    //Get user passed data
    Gen_User_Data* ud = static_cast<Gen_User_Data*>(usr_data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(params, params_size);

    //Model spectra based on new fit parameters
    ud->func(ud->fit_parameters, &(ud->energy_range), &(ud->spectra_model));
    ud->spectra_model += ud->spectra_background;

    //Calculate residuals
    for (int i=0; i<m; i++)
    {
        dy[i] = ( ud->spectra[i] - ud->spectra_model[i] ) * ud->weights[i];
    }

    return 0;
}


// =====================================================================================================================


MPFit_Optimizer::MPFit_Optimizer() : Optimizer()
{

}

MPFit_Optimizer::~MPFit_Optimizer()
{


}

void MPFit_Optimizer::minimize(Fit_Parameters *fit_params,
                               const Spectra * const spectra,
                               const Fit_Element_Map_Dict * const elements_to_fit,
                               const Base_Model * const model,
                               const Range energy_range)
{
    User_Data ud;

    fill_user_data(ud, fit_params, spectra, elements_to_fit, model, energy_range);

    std::vector<real_t> fitp_arr = fit_params->to_array();
    std::vector<real_t> perror(fitp_arr.size());

    int info;

    /////// init config ////////////
    struct mp_config<real_t> config;
    config.ftol = (real_t)1.0e-10;       // Relative chi-square convergence criterium  Default: 1e-10
    config.xtol = (real_t)1.0e-10;       // Relative parameter convergence criterium   Default: 1e-10
    config.gtol = (real_t)1.0e-10;       // Orthogonality convergence criterium        Default: 1e-10
	//config.epsfcn = 1.0;  
    config.epsfcn = MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    config.stepfactor = (real_t)100.0;   // Initial step bound                         Default: 100.0
    config.covtol = (real_t)1.0e-14;     // Range tolerance for covariance calculation Default: 1e-14
    config.maxiter = MP_NO_ITER;     //    Maximum number of iterations.  If maxiter == MP_NO_ITER,
    //config.maxiter = 200;     //    Maximum number of iterations.  If maxiter == MP_NO_ITER,
                                    //    then basic error checking is done, and parameter
                                    //    errors/covariances are estimated based on input
                                    //    parameter values, but no fitting iterations are done.
                                    //    Default: 200

    config.maxfev = 0;
    //config.maxfev = 1000 *(fitp_arr.size()+1);        // Maximum number of function evaluations, or 0 for no limit
                                       // Default: 0 (no limit)
    config.nprint = 0;           // Default: 1
    config.douserscale = 0;      // Scale variables by user values?
                                    //    1 = yes, user scale values in diag;
                                    //    0 = no, variables scaled internally (Default)
    config.nofinitecheck = 0;    // Disable check for infinite quantities from user?
                                    //    0 = do not perform check (Default)
                                    //    1 = perform check

    config.iterproc = 0;         // Placeholder pointer - must set to 0


    /////////////// init params limits /////////////////////////
	std::vector<struct mp_par<real_t> > par;
	par.resize(fitp_arr.size());
    //struct mp_par<real_t> *par = new struct mp_par<real_t>[fitp_arr.size()];

    for(auto itr = fit_params->begin(); itr != fit_params->end(); itr++)
    {
        Fit_Param fit = (*fit_params)[itr->first];
        if (fit.opt_array_index > -1)
        {

            if(fit.value > fit.max_val)
            {
                fit.max_val = fit.value + (real_t)1.0;
                (*fit_params)[itr->first].max_val = fit.value + (real_t)1.0;
            }
            if(fit.value < fit.min_val)
            {
                fit.min_val = fit.value - (real_t)1.0;
                (*fit_params)[itr->first].min_val = fit.value - (real_t)1.0;
            }
            if(fit.bound_type == E_Bound_Type::LIMITED_HI
            || fit.bound_type == E_Bound_Type::LIMITED_LO
            || fit.bound_type == E_Bound_Type::LIMITED_LO_HI )
            {
                if(fit.max_val == fit.min_val)
                {
                    fit.max_val += (real_t)0.1;
                    fit.min_val -= (real_t)0.1;
                    (*fit_params)[itr->first].max_val += (real_t)1.0;
                    (*fit_params)[itr->first].min_val -= (real_t)1.0;
                }
            }


            par[fit.opt_array_index].fixed = 0;              // 1 = fixed; 0 = free
            switch (fit.bound_type)
            {
            case E_Bound_Type::LIMITED_HI:
                par[fit.opt_array_index].limited[0] = 0;   // 1 = low/upper limit; 0 = no limit
                par[fit.opt_array_index].limited[1] = 1;
                par[fit.opt_array_index].limits[0] = std::numeric_limits<real_t>::min();   // lower/upper limit boundary value
                par[fit.opt_array_index].limits[1] = fit.max_val;
                break;
            case E_Bound_Type::LIMITED_LO:
                par[fit.opt_array_index].limited[0] = 1;   // 1 = low/upper limit; 0 = no limit
                par[fit.opt_array_index].limited[1] = 0;
                par[fit.opt_array_index].limits[0] = fit.min_val;   // lower/upper limit boundary value
                par[fit.opt_array_index].limits[1] = std::numeric_limits<real_t>::max();
                break;
            case E_Bound_Type::LIMITED_LO_HI:
                par[fit.opt_array_index].limited[0] = 1;   // 1 = low/upper limit; 0 = no limit
                par[fit.opt_array_index].limited[1] = 1;
                par[fit.opt_array_index].limits[0] = fit.min_val;   // lower/upper limit boundary value
                par[fit.opt_array_index].limits[1] = fit.max_val;
                break;
            default:
                par[fit.opt_array_index].limited[0] = 0;   // 1 = low/upper limit; 0 = no limit
                par[fit.opt_array_index].limited[1] = 0;
                par[fit.opt_array_index].limits[0] = std::numeric_limits<real_t>::min();   // lower/upper limit boundary value
                par[fit.opt_array_index].limits[1] = std::numeric_limits<real_t>::max();
                break;
            }

			par[fit.opt_array_index].step = fit.step_size;      // Step size for finite difference
            par[fit.opt_array_index].parname = 0;
			par[fit.opt_array_index].relstep = fit.step_size * 0.1;   // Relative step size for finite difference
            par[fit.opt_array_index].side = 0;         // Sidedness of finite difference derivative
                     //     0 - one-sided derivative computed automatically
                     //     1 - one-sided derivative (f(x+h) - f(x)  )/h
                     //    -1 - one-sided derivative (f(x)   - f(x-h))/h
                     //     2 - two-sided derivative (f(x+h) - f(x-h))/(2*h)
                     // 3 - user-computed analytical derivatives

            par[fit.opt_array_index].deriv_debug = 0;  // Derivative debug mode: 1 = Yes; 0 = No;

                           //      If yes, compute both analytical and numerical
                           //      derivatives and print them to the console for
                           //      comparison.

                       //  NOTE: when debugging, do *not* set side = 3,
                       //  but rather to the kind of numerical derivative
                       //  you want to compare the user-analytical one to
                       //  (0, 1, -1, or 2).

            par[fit.opt_array_index].deriv_reltol = (real_t)0.00001; // Relative tolerance for derivative debug printout
            par[fit.opt_array_index].deriv_abstol = (real_t)0.00001; // Absolute tolerance for derivative debug printout
        }
    }


    mp_result<real_t> result;
    memset(&result,0,sizeof(result));
    result.xerror = &perror[0];

	//struct mp_config<real_t> *config2 = nullptr;
	//struct mp_par<real_t> * par2 = nullptr;
    //info = mpfit(residuals_mpfit, fitp_arr.size(), fitp_arr.size(), &fitp_arr[0], par2, config2, (void *) &ud, &result);

    info = mpfit(residuals_mpfit, energy_range.count(), fitp_arr.size(), &fitp_arr[0], &par[0], &config, (void *) &ud, &result);

    logit_s<<"*";


    switch(info)
    {
        case 0:
            logit<<"> Improper input parameters."<<"\n";
        break;
        case 1:
            logit<<"> Both actual and predicted relative reductions in the sum of squares are at most ftol. "<<"\n";
        break;
        case 2:
            logit<<"> Relative error between two consecutive iterates is at most xtol"<<"\n";
        break;
        case 3:
            logit<<"> Conditions for info = 1 and info = 2 both hold. "<<"\n";
        break;
        case 4:
            logit<<"> The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value."<<"\n";
        break;
        case 5:
            logit<<"> Number of calls to fcn has reached or exceeded maxfev."<<"\n";
        break;
        case 6:
            logit<<"> Ftol is too small. no further reduction in the sum of squares is possible."<<"\n";
        break;
        case 7:
            logit<<"> Xtol is too small. no further improvement in the approximate solution x is possible."<<"\n";
        break;
        case 8:
            logit<<"> Gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision."<<"\n";
        break;
        default:
            logit<<"!> Unknown info status"<<"\n";
        break;
    }

    fit_params->from_array(fitp_arr);
    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = result.nfev;
    }

    //delete [] mp_par;

}

void MPFit_Optimizer::minimize_func(Fit_Parameters *fit_params,
                                    const Spectra * const spectra,
                                    const Range energy_range,
									Gen_Func_Def gen_func)
{
    Gen_User_Data ud;

    fill_gen_user_data(ud, fit_params, spectra, energy_range, gen_func);
    

    std::vector<real_t> fitp_arr = fit_params->to_array();
    std::vector<real_t> perror(fitp_arr.size());

    int info;

    /////// init config ////////////
    struct mp_config<real_t> mp_config;
    mp_config.ftol = (real_t)1.0e-10;       // Relative chi-square convergence criterium  Default: 1e-10
    mp_config.xtol = (real_t)1.0e-10;       // Relative parameter convergence criterium   Default: 1e-10
    mp_config.gtol = (real_t)1.0e-10;       // Orthogonality convergence criterium        Default: 1e-10
    mp_config.epsfcn = MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    mp_config.stepfactor = (real_t)0.1;   // Initial step bound                         Default: 100.0
    mp_config.covtol = (real_t)1.0e-14;     // Range tolerance for covariance calculation Default: 1e-14
    mp_config.maxiter = 1000;     //    Maximum number of iterations.  If maxiter == MP_NO_ITER,
                                    //    then basic error checking is done, and parameter
                                    //    errors/covariances are estimated based on input
                                    //    parameter values, but no fitting iterations are done.
                                    //    Default: 200

    //mp_config.maxfev = 0;
    mp_config.maxfev = 1000 *(fitp_arr.size()+1);        // Maximum number of function evaluations, or 0 for no limit
                                       // Default: 0 (no limit)
    mp_config.nprint = 0;           // Default: 1
    mp_config.douserscale = 0;      // Scale variables by user values?
                                    //    1 = yes, user scale values in diag;
                                    //    0 = no, variables scaled internally (Default)
    mp_config.nofinitecheck = 1;    // Disable check for infinite quantities from user?
                                    //    0 = do not perform check (Default)
                                    //    1 = perform check

    mp_config.iterproc = 0;         // Placeholder pointer - must set to 0

/*
    /////////////// init params limits /////////////////////////
    struct mp_par_struct *mp_par = new struct mp_par_struct[fitp_arr.size()];
    for(const auto& itr : fit_params->_params)
    {
        Fit_Param fit = itr.second;
        if (fit.opt_array_index > -1)
        {

            if(fit.value > fit.max_val)
            {
                fit.max_val = fit.value + 1.0;
                fit_params->_params[itr.first].max_val = fit.value + 1.0;
            }
            if(fit.value < fit.min_val)
            {
                fit.min_val = fit.value - 1.0;
                fit_params->_params[itr.first].min_val = fit.value - 1.0;
            }
            if(fit.bound_type == E_Bound_Type::LIMITED_HI
            || fit.bound_type == E_Bound_Type::LIMITED_LO
            || fit.bound_type == E_Bound_Type::LIMITED_LO_HI )
            {
                if(fit.max_val == fit.min_val)
                {
                    fit.max_val += 0.1;
                    fit.min_val -= 0.1;
                    fit_params->_params[itr.first].max_val += 1.0;
                    fit_params->_params[itr.first].min_val -= 1.0;
                }
            }


            mp_par[fit.opt_array_index].fixed = 0;              // 1 = fixed; 0 = free
            switch (fit.bound_type)
            {
            case E_Bound_Type::LIMITED_HI:
                mp_par[fit.opt_array_index].limited[0] = 0;   // 1 = low/upper limit; 0 = no limit
                mp_par[fit.opt_array_index].limited[1] = 1;
                mp_par[fit.opt_array_index].limits[0] = std::numeric_limits<real_t>::min();   // lower/upper limit boundary value
                mp_par[fit.opt_array_index].limits[1] = fit.max_val;
                break;
            case E_Bound_Type::LIMITED_LO:
                mp_par[fit.opt_array_index].limited[0] = 1;   // 1 = low/upper limit; 0 = no limit
                mp_par[fit.opt_array_index].limited[1] = 0;
                mp_par[fit.opt_array_index].limits[0] = fit.min_val;   // lower/upper limit boundary value
                mp_par[fit.opt_array_index].limits[1] = std::numeric_limits<real_t>::max();
                break;
            case E_Bound_Type::LIMITED_LO_HI:
                mp_par[fit.opt_array_index].limited[0] = 1;   // 1 = low/upper limit; 0 = no limit
                mp_par[fit.opt_array_index].limited[1] = 1;
                mp_par[fit.opt_array_index].limits[0] = fit.min_val;   // lower/upper limit boundary value
                mp_par[fit.opt_array_index].limits[1] = fit.max_val;
                break;
            default:
                mp_par[fit.opt_array_index].limited[0] = 0;   // 1 = low/upper limit; 0 = no limit
                mp_par[fit.opt_array_index].limited[1] = 0;
                mp_par[fit.opt_array_index].limits[0] = std::numeric_limits<real_t>::min();   // lower/upper limit boundary value
                mp_par[fit.opt_array_index].limits[1] = std::numeric_limits<real_t>::max();
                break;
            }

            mp_par[fit.opt_array_index].step = fit.step_size;      // Step size for finite difference
            mp_par[fit.opt_array_index].parname = 0;
            mp_par[fit.opt_array_index].relstep = fit.step_size;   // Relative step size for finite difference
            mp_par[fit.opt_array_index].side = 0;         // Sidedness of finite difference derivative
                     //     0 - one-sided derivative computed automatically
                     //     1 - one-sided derivative (f(x+h) - f(x)  )/h
                     //    -1 - one-sided derivative (f(x)   - f(x-h))/h
                     //     2 - two-sided derivative (f(x+h) - f(x-h))/(2*h)
                     // 3 - user-computed analytical derivatives

            mp_par[fit.opt_array_index].deriv_debug = 0;  // Derivative debug mode: 1 = Yes; 0 = No;

                           //      If yes, compute both analytical and numerical
                           //      derivatives and print them to the console for
                           //      comparison.

                       //  NOTE: when debugging, do *not* set side = 3,
                       //  but rather to the kind of numerical derivative
                       //  you want to compare the user-analytical one to
                       //  (0, 1, -1, or 2).

            mp_par[fit.opt_array_index].deriv_reltol = 0.00001; // Relative tolerance for derivative debug printout
            mp_par[fit.opt_array_index].deriv_abstol = 0.00001; // Absolute tolerance for derivative debug printout
        }
    }
*/

    mp_result<real_t> result;
    memset(&result,0,sizeof(result));
    result.xerror = &perror[0];

    struct mp_par<real_t> *mp_par = nullptr;

    info = mpfit(gen_residuals_mpfit, energy_range.count(), fitp_arr.size(), &fitp_arr[0], mp_par, &mp_config, (void *) &ud, &result);
    logit_s<<"*";


  //  delete [] mp_par;
/*
    switch(info)
    {
        case 0:
            logit<<"> Improper input parameters."<<"\n";
        break;
        case 1:
            logit<<"> Both actual and predicted relative reductions in the sum of squares are at most ftol. "<<"\n";
        break;
        case 2:
            logit<<"> Relative error between two consecutive iterates is at most xtol"<<"\n";
        break;
        case 3:
            logit<<"> Conditions for info = 1 and info = 2 both hold. "<<"\n";
        break;
        case 4:
            logit<<"> The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value."<<"\n";
        break;
        case 5:
            logit<<"> Number of calls to fcn has reached or exceeded maxfev."<<"\n";
        break;
        case 6:
            logit<<"> Ftol is too small. no further reduction in the sum of squares is possible."<<"\n";
        break;
        case 7:
            logit<<"> Xtol is too small. no further improvement in the approximate solution x is possible."<<"\n";
        break;
        case 8:
            logit<<"> Gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision."<<"\n";
        break;
        default:
            logit<<"!> Unknown info status"<<"\n";
        break;
    }
*/
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = static_cast<real_t>(result.nfev);
    }

}

} //namespace optimizers
} //namespace fitting
