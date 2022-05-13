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
#include <string.h>

using namespace data_struct;


namespace fitting
{
namespace optimizers
{

//-----------------------------------------------------------------------------
template<typename T_real>
int residuals_mpfit(int m, int params_size, T_real *params, T_real *dy, T_real **dvec, void *usr_data)
{
    // Get user passed data
    User_Data<T_real>* ud = static_cast<User_Data<T_real>*>(usr_data);

    // Update fit parameters from optimizer
    ud->fit_parameters->from_array(params, params_size);
    // Update background if fit_snip_width is set to fit
    update_background_user_data(ud);
    // Model spectra based on new fit parameters
    ud->spectra_model = ud->fit_model->model_spectrum_mp(ud->fit_parameters, ud->elements, ud->energy_range);
    // Add background
    ud->spectra_model += ud->spectra_background;
    // Remove nan's and inf's
    ud->spectra_model = (ArrayTr<T_real>)ud->spectra_model.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });

    //Calculate residuals
    for (int i=0; i<m; i++)
    {
		dy[i] = (ud->spectra[i] - ud->spectra_model[i]) * ud->weights[i];
		if (std::isfinite(dy[i]) == false)
		{
			logE << "\n\n\n";
			logE << ud->spectra[i] << "   :    " << ud->spectra_model[i] << "       ::      " << ud->weights[i];
			logE << "\n\n\n";
			//dy[i] = ud->spectra[i];
		}
    }
	
    ud->cur_itr++;
    if (ud->status_callback != nullptr)
    {
        (*ud->status_callback)(ud->cur_itr, ud->total_itr);
    }

	return 0;
}

//-----------------------------------------------------------------------------

template<typename T_real>
int gen_residuals_mpfit(int m, int params_size, T_real *params, T_real *dy, T_real **dvec, void *usr_data)
{
    // Get user passed data
    Gen_User_Data<T_real>* ud = static_cast<Gen_User_Data<T_real>*>(usr_data);

    // Update fit parameters from optimizer
    ud->fit_parameters->from_array(params, params_size);

    // Model spectra based on new fit parameters
    ud->func(ud->fit_parameters, &(ud->energy_range), &(ud->spectra_model));
    // Add background
    ud->spectra_model += ud->spectra_background;
    // Remove nan's and inf's
    ud->spectra_model = (ArrayTr<T_real>)ud->spectra_model.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });
    // Calculate residuals
    for (int i=0; i<m; i++)
    {
        dy[i] = ( ud->spectra[i] - ud->spectra_model[i] ) * ud->weights[i];
		if (std::isfinite(dy[i]) == false)
		{
			dy[i] = ud->spectra[i];
		}
    }

    return 0;
}


//-----------------------------------------------------------------------------

template<typename T_real>
int quantification_residuals_mpfit(int m, int params_size, T_real *params, T_real *dy, T_real **dvec, void *usr_data)
{
    ///(std::valarray<T_real> p, std::valarray<T_real> y, std::valarray<T_real> x)

    //y is array of elements standards
    //x is indexes of elements in standard
    //p is array size 2 but seems only first index is used
    ///return (y - this->fit_calibrationcurve(x, p));

    Quant_User_Data<T_real>* ud = (Quant_User_Data<T_real>*)(usr_data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(params, params_size);
    //Model spectra based on new fit parameters

    //Calculate residuals
    std::unordered_map<std::string, T_real> result_map = ud->quantification_model->model_calibrationcurve(ud->quant_map, params[0]);

    int idx = 0;
    for(auto& itr : ud->quant_map)
    {
		if (std::isfinite(result_map[itr.first]) == false)
		{
			dy[idx] = itr.second.e_cal_ratio;
		}
		else
		{
			dy[idx] = itr.second.e_cal_ratio - result_map[itr.first];
		}
        idx++;
    }

    return 0;
}

// =====================================================================================================================

template<typename T_real>
MPFit_Optimizer<T_real>::MPFit_Optimizer() : Optimizer<T_real>()
{
    //_options { 1e-10, 1e-10, 1e-10, MP_MACHEP0, 100.0, 1.0e-14, 2000, 0, 0, 0, 0, 0 };
    if (std::is_same<T_real, float>::value)
    {
        _options.ftol = 1.192e-10;       // Relative chi-square convergence criterium  Default: 1e-10
        _options.xtol = 1.192e-10;       // Relative parameter convergence criterium   Default: 1e-10
        _options.gtol = 1.192e-10;       // Orthogonality convergence criterium        Default: 1e-10
        _options.epsfcn = FP_MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    }
    else if (std::is_same<T_real, double>::value)
    {
        _options.ftol = 1.192e-10;       // Relative chi-square convergence criterium  Default: 1e-10
        _options.xtol = 1.192e-10;       // Relative parameter convergence criterium   Default: 1e-10
        _options.gtol = 1.192e-10;       // Orthogonality convergence criterium        Default: 1e-10
        _options.epsfcn = DP_MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    }
    _options.stepfactor = (T_real)100.0;   // Initial step bound                         Default: 100.0
    _options.covtol = (T_real)1.0e-14;     // Range tolerance for covariance calculation Default: 1e-14
    _options.maxiter = 1000;          //    Maximum number of iterations.  If maxiter == MP_NO_ITER,
                                    //    then basic error checking is done, and parameter
                                    //    errors/covariances are estimated based on input
                                    //    parameter values, but no fitting iterations are done.
                                    //    Default: 200

    _options.maxfev = 0;
    //_options.maxfev = 1000 *(fitp_arr.size()+1);        // Maximum number of function evaluations, or 0 for no limit
                                       // Default: 0 (no limit)
    _options.nprint = 0;           // Default: 1
    _options.douserscale = 0;      // Scale variables by user values?
                                    //    1 = yes, user scale values in diag;
                                    //    0 = no, variables scaled internally (Default)
    _options.nofinitecheck = 0;    // Disable check for infinite quantities from user?
                                    //    0 = do not perform check (Default)
                                    //    1 = perform check

    _options.iterproc = 0;         // Placeholder pointer - must set to 0


    this->_outcome_map[0] = OPTIMIZER_OUTCOME::FAILED;
    this->_outcome_map[1] = OPTIMIZER_OUTCOME::CONVERGED;
    this->_outcome_map[2] = OPTIMIZER_OUTCOME::CONVERGED;
    this->_outcome_map[3] = OPTIMIZER_OUTCOME::CONVERGED;
    this->_outcome_map[4] = OPTIMIZER_OUTCOME::TRAPPED;
    this->_outcome_map[5] = OPTIMIZER_OUTCOME::TRAPPED;
    this->_outcome_map[6] = OPTIMIZER_OUTCOME::F_TOL_LT_TOL;
    this->_outcome_map[7] = OPTIMIZER_OUTCOME::X_TOL_LT_TOL;
    this->_outcome_map[8] = OPTIMIZER_OUTCOME::G_TOL_LT_TOL;

}
// ----------------------------------------------------------------------------

template<typename T_real>
unordered_map<string, T_real> MPFit_Optimizer<T_real>::get_options()
{
    unordered_map<string, T_real> opts{
    {STR_OPT_FTOL, _options.ftol},
    {STR_OPT_XTOL, _options.xtol},
    {STR_OPT_GTOL, _options.gtol},
    {STR_OPT_EPSILON, _options.epsfcn},
    {STR_OPT_STEP, _options.stepfactor},
    {STR_OPT_COVTOL, _options.covtol},
    {STR_OPT_MAXITER, _options.maxiter}
    };

    return opts;
}

// ----------------------------------------------------------------------------

template<typename T_real>
void MPFit_Optimizer<T_real>::set_options(unordered_map<string, T_real> opt)
{
    if (opt.count(STR_OPT_FTOL) > 0)
    {
        _options.ftol = opt.at(STR_OPT_FTOL);
    }
    if (opt.count(STR_OPT_XTOL) > 0)
    {
        _options.xtol = opt.at(STR_OPT_XTOL);
    }
    if (opt.count(STR_OPT_GTOL) > 0)
    {
        _options.gtol = opt.at(STR_OPT_GTOL);
    }
    if (opt.count(STR_OPT_EPSILON) > 0)
    {
        _options.epsfcn = opt.at(STR_OPT_EPSILON);
    }
    if (opt.count(STR_OPT_STEP) > 0)
    {
        _options.stepfactor = opt.at(STR_OPT_STEP);
    }
    if (opt.count(STR_OPT_COVTOL) > 0)
    {
        _options.covtol = opt.at(STR_OPT_COVTOL);
    }
    if (opt.count(STR_OPT_MAXITER) > 0)
    {
        _options.maxiter = opt.at(STR_OPT_MAXITER);
    }
}

//-----------------------------------------------------------------------------


template<typename T_real>
void MPFit_Optimizer<T_real>::_fill_limits(Fit_Parameters<T_real> *fit_params , vector<struct mp_par<T_real> > &par)
{
	for (auto itr = fit_params->begin(); itr != fit_params->end(); itr++)
	{
		Fit_Param<T_real> fit = (*fit_params)[itr->first];
		if (fit.opt_array_index > -1)
		{

			if (fit.value > fit.max_val)
			{
				fit.max_val = fit.value + (T_real)1.0;
				(*fit_params)[itr->first].max_val = fit.value + (T_real)1.0;
			}
			if (fit.value < fit.min_val)
			{
				fit.min_val = fit.value - (T_real)1.0;
				(*fit_params)[itr->first].min_val = fit.value - (T_real)1.0;
			}
			if (fit.bound_type == E_Bound_Type::LIMITED_HI
				|| fit.bound_type == E_Bound_Type::LIMITED_LO
				|| fit.bound_type == E_Bound_Type::LIMITED_LO_HI)
			{
				if (fit.max_val == fit.min_val)
				{
					fit.max_val += (T_real)1.0;
					fit.min_val -= (T_real)1.0;
					(*fit_params)[itr->first].max_val += (T_real)1.0;
					(*fit_params)[itr->first].min_val -= (T_real)1.0;
				}
			}


			par[fit.opt_array_index].fixed = 0;              // 1 = fixed; 0 = free
			switch (fit.bound_type)
			{
			case E_Bound_Type::LIMITED_HI:
				par[fit.opt_array_index].limited[0] = 0;   // 1 = low/upper limit; 0 = no limit
				par[fit.opt_array_index].limited[1] = 1;
				par[fit.opt_array_index].limits[0] = std::numeric_limits<T_real>::min();   // lower/upper limit boundary value
				par[fit.opt_array_index].limits[1] = fit.max_val;
				break;
			case E_Bound_Type::LIMITED_LO:
				par[fit.opt_array_index].limited[0] = 1;   // 1 = low/upper limit; 0 = no limit
				par[fit.opt_array_index].limited[1] = 0;
				par[fit.opt_array_index].limits[0] = fit.min_val;   // lower/upper limit boundary value
				par[fit.opt_array_index].limits[1] = std::numeric_limits<T_real>::max();
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
				par[fit.opt_array_index].limits[0] = std::numeric_limits<T_real>::min();   // lower/upper limit boundary value
				par[fit.opt_array_index].limits[1] = std::numeric_limits<T_real>::max();
				break;
			}

			par[fit.opt_array_index].step = 0;      // 0 = auto ,Step size for finite difference
			par[fit.opt_array_index].parname = 0;
			par[fit.opt_array_index].relstep = 0;   // Relative step size for finite difference
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

			par[fit.opt_array_index].deriv_reltol = (T_real)0.00001; // Relative tolerance for derivative debug printout
			par[fit.opt_array_index].deriv_abstol = (T_real)0.00001; // Absolute tolerance for derivative debug printout
		}
	}
}

//-----------------------------------------------------------------------------

template<typename T_real>
void MPFit_Optimizer<T_real>::_print_info(int info)
{
	switch (info)
	{
	case 0:
		logI << "> Improper input parameters." << "\n";
		break;
	case 1:
		logI << "> Both actual and predicted relative reductions in the sum of squares are at most ftol. " << "\n";
		break;
	case 2:
		logI << "> Relative error between two consecutive iterates is at most xtol" << "\n";
		break;
	case 3:
		logI << "> Conditions for info = 1 and info = 2 both hold. " << "\n";
		break;
	case 4:
		logI << "> The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value." << "\n";
		break;
	case 5:
		logI << "> Number of calls to fcn has reached or exceeded maxfev." << "\n";
		break;
	case 6:
		logI << "> Ftol is too small. no further reduction in the sum of squares is possible." << "\n";
		break;
	case 7:
		logI << "> Xtol is too small. no further improvement in the approximate solution x is possible." << "\n";
		break;
	case 8:
		logI << "> Gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision." << "\n";
		break;
	default:
		logI << "!> Unknown info status" << "\n";
		break;
	}
}

//-----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME MPFit_Optimizer<T_real>::minimize(Fit_Parameters<T_real>*fit_params,
                                            const Spectra<T_real>* const spectra,
                                            const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                            const Base_Model<T_real>* const model,
                                            const Range energy_range,
                                            Callback_Func_Status_Def* status_callback)
{
    User_Data<T_real> ud;
    size_t num_itr = _options.maxiter;

    std::vector<T_real> fitp_arr = fit_params->to_array();
    std::vector<T_real> perror(fitp_arr.size());
    std::vector<T_real> resid(energy_range.count());

    size_t total_itr = num_itr * (fitp_arr.size() + 1);
    fill_user_data(ud, fit_params, spectra, elements_to_fit, model, energy_range, status_callback, total_itr);

    int info;
    /*
    /////// init config ////////////
    struct mp_config<T_real> config;
    config.ftol = 1e-10;       // Relative chi-square convergence criterium  Default: 1e-10
    config.xtol = 1e-10;       // Relative parameter convergence criterium   Default: 1e-10
    config.gtol = 1e-10;       // Orthogonality convergence criterium        Default: 1e-10
    config.epsfcn = MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    config.stepfactor = (T_real)100.0;   // Initial step bound                         Default: 100.0
    config.covtol = (T_real)1.0e-14;     // Range tolerance for covariance calculation Default: 1e-14
    config.maxiter = num_itr;          //    Maximum number of iterations.  If maxiter == MP_NO_ITER,
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

    */
    /////////////// init params limits /////////////////////////
	vector<struct mp_par<T_real> > par;
	par.resize(fitp_arr.size());

    _options.maxfev = _options.maxiter * (fitp_arr.size() + 1);

	_fill_limits(fit_params, par);

    mp_result<T_real> result;
    memset(&result,0,sizeof(result));
    result.xerror = &perror[0];
    result.resid = &resid[0];

    info = mpfit(residuals_mpfit<T_real>, energy_range.count(), fitp_arr.size(), &fitp_arr[0], &par[0], &_options, (void *) &ud, &result);

	_print_info(info);

    fit_params->from_array(fitp_arr);
    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = result.nfev;
    }
    if (fit_params->contains(STR_RESIDUAL))
    {
        T_real sum_resid = 0.0;
        for (int i = 0; i < energy_range.count(); i++)
        {
            sum_resid += resid[i];
        }
        (*fit_params)[STR_RESIDUAL].value = sum_resid;
    }

    if (this->_outcome_map.count(info) > 0)
        return this->_outcome_map[info];

    return OPTIMIZER_OUTCOME::FAILED;
}

//-----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME MPFit_Optimizer<T_real>::minimize_func(Fit_Parameters<T_real> *fit_params,
                                                const Spectra<T_real>* const spectra,
                                                const Range energy_range,
                                                const ArrayTr<T_real>* background,
									            Gen_Func_Def<T_real> gen_func)
{
    Gen_User_Data<T_real> ud;
    fill_gen_user_data(ud, fit_params, spectra, energy_range, background, gen_func);

    std::vector<T_real> fitp_arr = fit_params->to_array();
    std::vector<T_real> perror(fitp_arr.size());
    std::vector<T_real> resid(energy_range.count());

    int info;
    /*
    /////// init config ////////////
    struct mp_config<T_real> mp_config;
    mp_config.ftol = 1e-10;       // Relative chi-square convergence criterium  Default: 1e-10
    mp_config.xtol = 1e-10;       // Relative parameter convergence criterium   Default: 1e-10
    mp_config.gtol = 1e-10;       // Orthogonality convergence criterium        Default: 1e-10
    mp_config.epsfcn = MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    mp_config.stepfactor = (T_real)100.0;   // Initial step bound                         Default: 100.0
    mp_config.covtol = (T_real)1.0e-14;     // Range tolerance for covariance calculation Default: 1e-14
    mp_config.maxiter = 200;        //    Maximum number of iterations.  If maxiter == MP_NO_ITER,
                                    //    then basic error checking is done, and parameter
                                    //    errors/covariances are estimated based on input
                                    //    parameter values, but no fitting iterations are done.
                                    //    Default: 200


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
    */

    _options.maxfev = _options.maxiter * (fitp_arr.size() + 1);

	vector<struct mp_par<T_real> > par;
	par.resize(fitp_arr.size());
	_fill_limits(fit_params, par);

    mp_result<T_real> result;
    memset(&result,0,sizeof(result));
    result.xerror = &perror[0];
    result.resid = &resid[0];

    info = mpfit(gen_residuals_mpfit<T_real>, energy_range.count(), fitp_arr.size(), &fitp_arr[0], &par[0], &_options, (void*)&ud, &result);
/*
    
*/
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = static_cast<T_real>(result.nfev);
    }
    if (fit_params->contains(STR_RESIDUAL))
    {
        T_real sum_resid = 0.0;
        for (int i = 0; i< energy_range.count(); i++)
        {
             sum_resid += std::abs(resid[i]);
        }
        (*fit_params)[STR_RESIDUAL].value = sum_resid;
    }

    if (this->_outcome_map.count(info) > 0)
        return this->_outcome_map[info];

    return OPTIMIZER_OUTCOME::FAILED;
}

//-----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME MPFit_Optimizer<T_real>::minimize_quantification(Fit_Parameters<T_real> *fit_params,
                                                          std::unordered_map<std::string, Element_Quant<T_real>*> * quant_map,
                                                          quantification::models::Quantification_Model<T_real>* quantification_model)
{
    Quant_User_Data<T_real> ud;

    if (quant_map != nullptr)
    {
        for (const auto& itr : *quant_map)
        {
            ud.quant_map[itr.first] = *(itr.second);
        }
    }
    ud.quantification_model = quantification_model;
    ud.fit_parameters = fit_params;

    std::vector<T_real> fitp_arr = fit_params->to_array();
    std::vector<T_real> perror(fitp_arr.size());
    std::vector<T_real> resid(quant_map->size());

    int info;
    /*
    /////// init config ////////////
    struct mp_config<T_real> mp_config;
    mp_config.ftol = 1e-10;       // Relative chi-square convergence criterium  Default: 1e-10
    mp_config.xtol = 1e-10;       // Relative parameter convergence criterium   Default: 1e-10
    mp_config.gtol = 1e-10;       // Orthogonality convergence criterium        Default: 1e-10
    mp_config.epsfcn = MP_MACHEP0;  // Finite derivative step size                Default: MP_MACHEP0
    mp_config.stepfactor = (T_real)100.0;   // Initial step bound                         Default: 100.0
    mp_config.covtol = (T_real)1.0e-14;     // Range tolerance for covariance calculation Default: 1e-14
    mp_config.maxiter = 1000;        //    Maximum number of iterations.  If maxiter == MP_NO_ITER,
                                    //    then basic error checking is done, and parameter
                                    //    errors/covariances are estimated based on input
                                    //    parameter values, but no fitting iterations are done.
                                    //    Default: 200


    //mp_config.maxfev = 1000 *(fitp_arr.size()+1);        // Maximum number of function evaluations, or 0 for no limit
    mp_config.maxfev = 0;
                                       // Default: 0 (no limit)
    mp_config.nprint = 0;           // Default: 1
    mp_config.douserscale = 0;      // Scale variables by user values?
                                    //    1 = yes, user scale values in diag;
                                    //    0 = no, variables scaled internally (Default)
    mp_config.nofinitecheck = 1;    // Disable check for infinite quantities from user?
                                    //    0 = do not perform check (Default)
                                    //    1 = perform check

    mp_config.iterproc = 0;         // Placeholder pointer - must set to 0
    */

    _options.maxfev = _options.maxiter * (fitp_arr.size() + 1);

    mp_result<T_real> result;
    memset(&result,0,sizeof(result));
    result.xerror = &perror[0];
    result.resid = &resid[0];
//    struct mp_par<T_real> *mp_par = nullptr;
//	info = mpfit(quantification_residuals_mpfit, quant_map->size(), fitp_arr.size(), &fitp_arr[0], mp_par, &mp_config, (void *)&ud, &result);

	
	vector<struct mp_par<T_real> > par;
	par.resize(fitp_arr.size());
	_fill_limits(fit_params, par);

    info = mpfit(quantification_residuals_mpfit<T_real>, quant_map->size(), fitp_arr.size(), &fitp_arr[0], &par[0], &_options, (void *) &ud, &result);
    logI << "\nOutcome: " << optimizer_outcome_to_str(this->_outcome_map[info]) << "\nNum iter: " << result.niter << "\n Norm of the residue vector: " << *result.resid << "\n";

	_print_info(info);

    fit_params->from_array(fitp_arr);

    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = static_cast<T_real>(result.nfev);
    }
    if (fit_params->contains(STR_RESIDUAL))
    {
        T_real sum_resid = 0.0;
        for (int i = 0; i < quant_map->size(); i++)
        {
            sum_resid += resid[i];
        }
        (*fit_params)[STR_RESIDUAL].value = sum_resid;
    }
    if (this->_outcome_map.count(info) > 0)
        return this->_outcome_map[info];

    return OPTIMIZER_OUTCOME::FAILED;

}

//-----------------------------------------------------------------------------

} //namespace optimizers
} //namespace fitting
