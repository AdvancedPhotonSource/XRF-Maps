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

#include "mpfit.h"
#include <string.h>

using namespace data_struct::xrf;


namespace fitting
{
namespace optimizers
{

int residuals_mpfit(int m, int params_size, double *params, double *dy, double **dvec, void *usr_data)
{
    //Get user passed data
    User_Data* ud = static_cast<User_Data*>(usr_data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(params, params_size);

    //Model spectra based on new fit parameters
    Spectra spectra_model = ud->fit_model->model_spectrum(ud->fit_parameters, ud->elements, *(ud->energy_range));

    //Calculate residuals
    std::valarray<real_t> residuals = ( (*ud->spectra) - spectra_model ) * (*ud->weights);

    for (int i=0; i<m; i++)
    {
        dy[i] = residuals[i];
    }

	return 0;
}

int gen_residuals_mpfit(int m, int params_size, double *params, double *dy, double **dvec, void *usr_data)
{
    //Get user passed data
    Gen_User_Data* ud = static_cast<Gen_User_Data*>(usr_data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(params, params_size);

    //Model spectra based on new fit parameters
    Spectra spectra_model = ud->func(ud->fit_parameters, ud->energy_range);

    //Calculate residuals
    std::valarray<real_t> residuals = ( (*ud->spectra) - spectra_model ) * (*ud->weights);

    for (int i=0; i<m; i++)
    {
        dy[i] = residuals[i];
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
                               const Base_Model * const model)
{
    User_Data ud;

    ud.fit_model = (Base_Model*)model;
    // set spectra to fit
    ud.spectra = (Spectra*)spectra;
    ud.fit_parameters = fit_params;
    ud.elements = (Fit_Element_Map_Dict *)elements_to_fit;

    //fitting::models::Range energy_range = fitting::models::get_energy_range(1.0, 11.0, spectra->size(), detector);
    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = spectra->size()-1;
    ud.energy_range = &energy_range;
/*
    Fit_Counts_Array counts_arr;
    counts_arr.resize(energy_range.count());
    ud.counts_arr = &counts_arr;
*/
    std::vector<real_t> fitp_arr = fit_params->to_array();
    std::vector<real_t> perror(fitp_arr.size());

    int info;


    std::valarray<real_t> weights = (real_t)1.0 / ( (real_t)1.0 + (*spectra) );
    weights = convolve1d(weights, 5);
    weights = std::abs(weights);
    weights /= weights.max();
    ud.weights = &weights;

    //std::valarray<real_t> weights = std::sqrt( *(spectra->buffer()) );
    //ud.weights = &weights;

    /////// init config ////////////
    struct mp_config_struct mp_config;
    mp_config.ftol = 1.0e-10;       // Relative chi-square convergence criterium  Default: 1e-10
    mp_config.xtol = 1.0e-10;       // Relative parameter convergence criterium   Default: 1e-10
    mp_config.gtol = 1.0e-10;       // Orthogonality convergence criterium        Default: 1e-10
    mp_config.epsfcn = MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    mp_config.stepfactor = 0.1;   // Initial step bound                         Default: 100.0
    mp_config.covtol = 1.0e-14;     // Range tolerance for covariance calculation Default: 1e-14
    mp_config.maxiter = 2000;     //    Maximum number of iterations.  If maxiter == MP_NO_ITER,
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


    /////////////// init params limits /////////////////////////
    struct mp_par_struct *mp_par = new struct mp_par_struct[fitp_arr.size()];

    for(auto itr = fit_params->begin(); itr != fit_params->end(); itr++)
    {
        Fit_Param fit = (*fit_params)[itr->first];
        if (fit.opt_array_index > -1)
        {

            if(fit.value > fit.max_val)
            {
                fit.max_val = fit.value + 1.0;
                (*fit_params)[itr->first].max_val = fit.value + 1.0;
            }
            if(fit.value < fit.min_val)
            {
                fit.min_val = fit.value - 1.0;
                (*fit_params)[itr->first].min_val = fit.value - 1.0;
            }
            if(fit.bound_type == E_Bound_Type::LIMITED_HI
            || fit.bound_type == E_Bound_Type::LIMITED_LO
            || fit.bound_type == E_Bound_Type::LIMITED_LO_HI )
            {
                if(fit.max_val == fit.min_val)
                {
                    fit.max_val += 0.1;
                    fit.min_val -= 0.1;
                    (*fit_params)[itr->first].max_val += 1.0;
                    (*fit_params)[itr->first].min_val -= 1.0;
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


    mp_result result;
    memset(&result,0,sizeof(result));
    result.xerror = &perror[0];

    //info = mpfit(residuals_mpfit, fitp_arr.size(), fitp_arr.size(), &fitp_arr[0], 0, 0, (void *) &ud, &result);

    //info = mpfit(residuals_mpfit, fitp_arr.size(), fitp_arr.size(), &fitp_arr[0], 0, &mp_config, (void *) &ud, &result);

    //info = mpfit(residuals_mpfit, fitp_arr.size(), fitp_arr.size(), &fitp_arr[0], &mp_par[0], 0, (void *) &ud, &result);

    info = mpfit(residuals_mpfit, spectra->size(), fitp_arr.size(), &fitp_arr[0], &mp_par[0], &mp_config, (void *) &ud, &result);

    std::cout<<"*** testlinfit status = "<<info<<std::endl;


  //  delete [] mp_par;

    switch(info)
    {
        case 0:
            std::cout<<"> Improper input parameters."<<std::endl;
        break;
        case 1:
            std::cout<<"> Both actual and predicted relative reductions in the sum of squares are at most ftol. "<<std::endl;
        break;
        case 2:
            std::cout<<"> Relative error between two consecutive iterates is at most xtol"<<std::endl;
        break;
        case 3:
            std::cout<<"> Conditions for info = 1 and info = 2 both hold. "<<std::endl;
        break;
        case 4:
            std::cout<<"> The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value."<<std::endl;
        break;
        case 5:
            std::cout<<"> Number of calls to fcn has reached or exceeded maxfev."<<std::endl;
        break;
        case 6:
            std::cout<<"> Ftol is too small. no further reduction in the sum of squares is possible."<<std::endl;
        break;
        case 7:
            std::cout<<"> Xtol is too small. no further improvement in the approximate solution x is possible."<<std::endl;
        break;
        case 8:
            std::cout<<"> Gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision."<<std::endl;
        break;
        default:
            std::cout<<"!> Unknown info status"<<std::endl;
        break;
    }

    fit_params->from_array(fitp_arr);
    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = result.nfev;
    }

}

void MPFit_Optimizer::minimize_func(Fit_Parameters *fit_params,
                                    const Spectra * const spectra,
                                    std::function<const Spectra(const Fit_Parameters* const, const struct Range* const)> gen_func)
{
    Gen_User_Data ud;


    ud.func = gen_func;
    // set spectra to fit
    ud.spectra = (Spectra*)spectra;
    ud.fit_parameters = fit_params;

    //fitting::models::Range energy_range = fitting::models::get_energy_range(1.0, 11.0, spectra->size(), detector);
    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = spectra->size()-1;
    ud.energy_range = &energy_range;

    std::vector<real_t> fitp_arr = fit_params->to_array();
    std::vector<real_t> perror(fitp_arr.size());

    int info;

    std::valarray<real_t> weights = (real_t)1.0 / ( (real_t)1.0 + (*spectra) );
    weights = convolve1d(weights, 5);
    weights = std::abs(weights);
    weights /= weights.max();
    ud.weights = &weights;

    //std::valarray<real_t> weights = std::sqrt( *(spectra->buffer()) );
    //ud.weights = &weights;

    /////// init config ////////////
    struct mp_config_struct mp_config;
    mp_config.ftol = 1.0e-10;       // Relative chi-square convergence criterium  Default: 1e-10
    mp_config.xtol = 1.0e-10;       // Relative parameter convergence criterium   Default: 1e-10
    mp_config.gtol = 1.0e-10;       // Orthogonality convergence criterium        Default: 1e-10
    mp_config.epsfcn = MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    mp_config.stepfactor = 0.1;   // Initial step bound                         Default: 100.0
    mp_config.covtol = 1.0e-14;     // Range tolerance for covariance calculation Default: 1e-14
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

    mp_result result;
    memset(&result,0,sizeof(result));
    result.xerror = &perror[0];
    //info = mpfit(residuals_mpfit, fitp_arr.size(), fitp_arr.size(), &fitp_arr[0], 0, 0, (void *) &ud, &result);
    //info = mpfit(residuals_mpfit, fitp_arr.size(), fitp_arr.size(), &fitp_arr[0], 0, &mp_config, (void *) &ud, &result);

    //info = mpfit(residuals_mpfit, fitp_arr.size(), fitp_arr.size(), &fitp_arr[0], &mp_par[0], 0, (void *) &ud, &result);

    info = mpfit(gen_residuals_mpfit, spectra->size(), fitp_arr.size(), &fitp_arr[0], 0, &mp_config, (void *) &ud, &result);
    std::cout<<"*** testlinfit status = "<<info<<std::endl;


  //  delete [] mp_par;

    switch(info)
    {
        case 0:
            std::cout<<"> Improper input parameters."<<std::endl;
        break;
        case 1:
            std::cout<<"> Both actual and predicted relative reductions in the sum of squares are at most ftol. "<<std::endl;
        break;
        case 2:
            std::cout<<"> Relative error between two consecutive iterates is at most xtol"<<std::endl;
        break;
        case 3:
            std::cout<<"> Conditions for info = 1 and info = 2 both hold. "<<std::endl;
        break;
        case 4:
            std::cout<<"> The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value."<<std::endl;
        break;
        case 5:
            std::cout<<"> Number of calls to fcn has reached or exceeded maxfev."<<std::endl;
        break;
        case 6:
            std::cout<<"> Ftol is too small. no further reduction in the sum of squares is possible."<<std::endl;
        break;
        case 7:
            std::cout<<"> Xtol is too small. no further improvement in the approximate solution x is possible."<<std::endl;
        break;
        case 8:
            std::cout<<"> Gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision."<<std::endl;
        break;
        default:
            std::cout<<"!> Unknown info status"<<std::endl;
        break;
    }

    fit_params->from_array(fitp_arr);

    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = result.nfev;
    }

}

} //namespace optimizers
} //namespace fitting
