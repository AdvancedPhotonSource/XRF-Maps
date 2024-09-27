/***
Copyright (c) 2024, UChicago Argonne, LLC. All rights reserved.

Copyright 2024. UChicago Argonne, LLC. This software was produced
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



#include "nlopt_optimizer.h"

#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>
#include <string.h>

using namespace data_struct;


namespace fitting
{
namespace optimizers
{

//-----------------------------------------------------------------------------

//int residuals_nlopt(int m, int params_size, T_real *params, T_real *dy, T_real **dvec, void *usr_data)
double residuals_nlopt(const std::vector<double> &x, std::vector<double> &grad, void *usr_data)
{
    bool first = true;
    // Get user passed data
    User_Data<double>* ud = static_cast<User_Data<double>*>(usr_data);

    // Debug to find which param changed last
    Fit_Parameters<double> prev_fit_p;
    prev_fit_p.update_and_add_values(ud->fit_parameters);

    // Update fit parameters from optimizer
    ud->fit_parameters->from_array_d(x);
    // Update background if fit_snip_width is set to fit
    update_background_user_data(ud);
    // Model spectra based on new fit parameters
    ud->spectra_model = ud->fit_model->model_spectrum_mp(ud->fit_parameters, ud->elements, ud->energy_range);
    // Add background
    ud->spectra_model += ud->spectra_background;
    // Remove nan's and inf's
    //ud->spectra_model = (ArrayTr<T_real>)ud->spectra_model.unaryExpr([ud](T_real v) { return std::isfinite(v) ? v : ud->normalizer; });

    double sum = 0.0;
    double dy = 0.;
    //Calculate residuals
    for (int i=0; i<ud->spectra.size(); i++)
    {
        dy += pow((ud->spectra[i] - ud->spectra_model[i]), 2.0) * ud->weights[i];

    	if (std::isfinite(dy) == false)
		{
            if(first)
            {
                first = false;
			    logE << "Spectra["<<i<<"] = "<< ud->spectra[i] << " ::spectra_model["<<i<<"] = " << ud->spectra_model[i] << "  ::weights["<<i<<"] = " << ud->weights[i]<<"\n";
                logE<<" \n Diff Param \n";
                for(auto &itr : prev_fit_p)
                {
                    if(itr.second.value != ud->fit_parameters->at(itr.first).value)
                    {
                        logE<<itr.first<<" : Old = "<<itr.second.value<<" ; New = "<< ud->fit_parameters->at(itr.first).value <<"\n";
                    }
                }
                logE<<" \n \n";
                ud->fit_parameters->print_non_fixed();
            }
            sum += ud->normalizer;
		}
        else
        {
            sum += dy;
        }
    }
    //logI << "f = " << sum << "\n";
    ud->cur_itr++;
    if (ud->status_callback != nullptr)
    {
        try
        {
            (*ud->status_callback)(ud->cur_itr, ud->total_itr);
        }
        catch (int e)
        {
            logI << "Cancel fitting" << std::endl;
            return -1;
        }
    }

	return sum;
}

//-----------------------------------------------------------------------------

template<typename T_real>
int gen_residuals_nlopt(int m, int params_size, T_real *params, T_real *dy, T_real **dvec, void *usr_data)
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
    //ud->spectra_model = (ArrayTr<T_real>)ud->spectra_model.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });
    // 
    
    // Calculate residuals
    for (int i=0; i<m; i++)
    {
        //dy[i] = std::pow(( ud->spectra[i] - ud->spectra_model[i] ), 2.0) * ud->weights[i];
		T_real n_raw = ud->spectra[i] / ud->normalizer;
        T_real n_model = ud->spectra_model[i] / ud->normalizer;
        dy[i] = pow((n_raw - n_model), (T_real)2.0) * ud->weights[i];
        /*
        if (std::isfinite(dy[i]) == false)
		{
			//dy[i] = ud->spectra[i] + ud->spectra_model[i];
            dy[i] = std::numeric_limits<T_real>::quiet_NaN();
		}
        */
    }

    return 0;
}


//-----------------------------------------------------------------------------

template<typename T_real>
int quantification_residuals_nlopt(int m, int params_size, T_real *params, T_real *dy, T_real **dvec, void *usr_data)
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
			dy[idx] = itr.second.e_cal_ratio * 10.0;
		}
		else
		{
			dy[idx] = pow((itr.second.e_cal_ratio - result_map[itr.first]), (T_real)2.0);
		}
        idx++;
    }

    return 0;
}

// =====================================================================================================================

template<typename T_real>
NLOPT_Optimizer<T_real>::NLOPT_Optimizer() : Optimizer<T_real>()
{

    this->_last_outcome = -1;
/*
    //_options { 1e-10, 1e-10, 1e-10, MP_MACHEP0, 100.0, 1.0e-14, 2000, 0, 0, 0, 0, 0 };
    if (std::is_same<T_real, float>::value)
    {
        _options.ftol = (T_real)1.192e-10;       // Relative chi-square convergence criterium  Default: 1e-10
        _options.xtol = (T_real)1.192e-10;       // Relative parameter convergence criterium   Default: 1e-10
        _options.gtol = (T_real)1.192e-10;       // Orthogonality convergence criterium        Default: 1e-10
        _options.epsfcn = (T_real)1.192e-14;  // Finite derivative step size                Default: MP_MACHEP0
    }
    else if (std::is_same<T_real, double>::value)
    {
        _options.ftol = (T_real)1.192e-10;       // Relative chi-square convergence criterium  Default: 1e-10
        _options.xtol = (T_real)1.192e-10;       // Relative parameter convergence criterium   Default: 1e-10
        _options.gtol = (T_real)1.192e-10;       // Orthogonality convergence criterium        Default: 1e-10
        _options.epsfcn = (T_real)1.192e-14;  // Finite derivative step size                Default: MP_MACHEP0
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

*/
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
std::unordered_map<std::string, T_real> NLOPT_Optimizer<T_real>::get_options()
{
    /*
    std::unordered_map<std::string, T_real> opts{
    {STR_OPT_FTOL, _options.ftol},
    {STR_OPT_XTOL, _options.xtol},
    {STR_OPT_GTOL, _options.gtol},
    {STR_OPT_EPSILON, _options.epsfcn},
    {STR_OPT_STEP, _options.stepfactor},
    {STR_OPT_COVTOL, _options.covtol},
    {STR_OPT_MAXITER, _options.maxiter}
    };
*/
std::unordered_map<std::string, T_real> opts;
    return opts;
}

// ----------------------------------------------------------------------------

template<typename T_real>
void NLOPT_Optimizer<T_real>::set_options(std::unordered_map<std::string, T_real> opt)
{
    /*
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
    */
}

//-----------------------------------------------------------------------------

/*
template<typename T_real>
bool NLOPT_Optimizer<T_real>::_fill_limits(Fit_Parameters<T_real> *fit_params , std::vector<struct mp_par<T_real> > &par)
{
	for (auto itr = fit_params->begin(); itr != fit_params->end(); itr++)
	{
		Fit_Param<T_real> fit = (*fit_params)[itr->first];
		if (fit.opt_array_index > -1)
		{

			if (fit.value > fit.max_val)
			{
                //logW << itr->first << " value (" << fit.value << ") > max_val(" << fit.max_val << ") : setting value = max_val\n";
                fit.value = fit.max_val - fit.step_size;
				(*fit_params)[itr->first].value = fit.value;

			}
			if (fit.value < fit.min_val)
			{
                //logW << itr->first << " value (" << fit.value << ") < min_val(" << fit.min_val << ") : setting value = min_val\n";
                fit.value = fit.min_val + fit.step_size;
				(*fit_params)[itr->first].value = fit.min_val;
			}
			if (fit.bound_type == E_Bound_Type::LIMITED_HI
				|| fit.bound_type == E_Bound_Type::LIMITED_LO
				|| fit.bound_type == E_Bound_Type::LIMITED_LO_HI)
			{
                if(fit.min_val == fit.max_val)
                {
                    logE<<"Fit parameter "<<fit.name<<" is limited and has min: "<<fit.min_val<<" == to max: "<<fit.max_val<<"\n";
                    return false;
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
    return true;
}
*/
//-----------------------------------------------------------------------------

template<typename T_real>
std::string NLOPT_Optimizer<T_real>::detailed_outcome(int info)
{
	switch (info)
	{
	case 0:
		return "Improper input parameters.";
	case 1:
        return "Both actual and predicted relative reductions in the sum of squares are at most ftol. ";
	case 2:
        return "Relative error between two consecutive iterates is at most xtol.";
	case 3:
        return "Conditions for info = 1 and info = 2 both hold.";
	case 4:
        return "The cosine of the angle between fvec and any column of the jacobian is at most gtol in absolute value.";
	case 5:
        return "Number of calls to fcn has reached or exceeded maxfev.";
	case 6:
        return "Ftol is too small. no further reduction in the sum of squares is possible.";	
	case 7:
        return "Xtol is too small. no further improvement in the approximate solution x is possible.";
	case 8:
        return "Gtol is too small. fvec is orthogonal to the columns of the jacobian to machine precision.";
        /*
    FAILURE = -1,         // generic failure code 
    INVALID_ARGS = -2,
    OUT_OF_MEMORY = -3,
    ROUNDOFF_LIMITED = -4,
    FORCED_STOP = -5,
    NUM_FAILURES = -6,    // not a result, just the number of possible failures 
    SUCCESS = 1,          // generic success code 
    STOPVAL_REACHED = 2,
    FTOL_REACHED = 3,
    XTOL_REACHED = 4,
    MAXEVAL_REACHED = 5,
    MAXTIME_REACHED = 6,
    NUM_RESULTS  

    case MP_ERR_NAN:
        return "User function produced non-finite values.";
    case MP_ERR_FUNC:
        return "No user function was supplied.";
    case MP_ERR_NPOINTS:
        return "No user data points were supplied.";
    case MP_ERR_NFREE:
        return "No free parameters.";
    case MP_ERR_MEMORY:
        return "Memory allocation error.";
    case MP_ERR_INITBOUNDS:
        return "Initial values inconsistent w constraints.";
    case MP_ERR_BOUNDS:
        return "Initial constraints inconsistent.";
    case MP_ERR_PARAM:
        return "General input parameter error.";
    case MP_ERR_DOF:
        return "Not enough degrees of freedom.";
        */
	default:
        return "Unknown info status";
	}
}

//-----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME NLOPT_Optimizer<T_real>::minimize(Fit_Parameters<T_real>*fit_params,
                                            const Spectra<T_real>* const spectra,
                                            const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                            const Base_Model<T_real>* const model,
                                            const Range energy_range,
                                            bool use_weights,
                                            Callback_Func_Status_Def* status_callback)
{
    User_Data<T_real> ud;

    std::vector<double> fitp_arr;
    std::vector<double> lb_arr;
    std::vector<double> ub_arr;
    
    fit_params->to_array_with_bounds(fitp_arr, lb_arr, ub_arr);
    if (fitp_arr.size() == 0)
    {
        return OPTIMIZER_OUTCOME::STOPPED;
    }
    std::vector<T_real> perror(fitp_arr.size());
    std::vector<T_real> resid(energy_range.count());
    std::vector<T_real> covar(fitp_arr.size() * fitp_arr.size());

    size_t total_itr = 1000; //num_itr
    fill_user_data(ud, fit_params, spectra, elements_to_fit, model, energy_range, status_callback, total_itr, use_weights);

    nlopt::opt opt(nlopt::LN_COBYLA, fitp_arr.size());
    opt.set_lower_bounds(lb_arr);
    opt.set_upper_bounds(ub_arr);
    opt.set_min_objective(residuals_nlopt, (void*)&ud);
    //my_constraint_data data[2] = { {2,0}, {-1,1} };
    //opt.add_inequality_constraint(myvconstraint, &data[0], 1e-8);
    //opt.add_inequality_constraint(myvconstraint, &data[1], 1e-8);
    opt.set_xtol_rel(1e-14);
    opt.set_maxeval(20000);

    double minf;

    try
    {
        nlopt::result result = opt.optimize(fitp_arr, minf);
        logI<< "resid = "<<minf<<"\n";
    }
    catch(std::exception &e) 
    {
        logW << "nlopt failed: " << e.what() << "\n";
    }


    /////////////// init params limits /////////////////////////
    /*
	std::vector<struct mp_par<T_real> > par;
	par.resize(fitp_arr.size());

    _options.maxfev = _options.maxiter * ((int)fitp_arr.size() + 1);

	if(false == _fill_limits(fit_params, par))
    {
        return OPTIMIZER_OUTCOME::FAILED;
    }

    this->_last_outcome = mpfit(residuals_mpfit<T_real>, (int)energy_range.count(), (int)fitp_arr.size(), &fitp_arr[0], &par[0], &_options, (void *) &ud, &result);

	logI<< detailed_outcome(this->_last_outcome) << "\n";
*/
    fit_params->from_array_d(fitp_arr);

    T_real sum_resid = 0.0;
    for (size_t i = 0; i < energy_range.count(); i++)
    {
        sum_resid += resid[i];
    }
/*
    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = result.nfev;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, result.nfev));
    }

    if (fit_params->contains(STR_RESIDUAL))
    {
        (*fit_params)[STR_RESIDUAL].value = sum_resid;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_RESIDUAL, sum_resid));
    }

    if (fit_params->contains(STR_CHISQUARE))
    {
        (*fit_params)[STR_CHISQUARE].value = result.bestnorm;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_CHISQUARE, result.bestnorm));
    }

    if (fit_params->contains(STR_CHISQRED))
    {
        (*fit_params)[STR_CHISQRED].value = result.bestnorm / fitp_arr.size();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_CHISQRED, result.bestnorm / fitp_arr.size()));
    }
*/
    if (fit_params->contains(STR_FREE_PARS))
    {
        (*fit_params)[STR_FREE_PARS].value = fitp_arr.size();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_FREE_PARS, fitp_arr.size()));
    }
    // add perror_ fit params
    Fit_Parameters<T_real> error_params;
    for (typename std::unordered_map<std::string, Fit_Param<T_real>>::const_iterator itr = fit_params->begin(); itr != fit_params->end(); itr++)
    {
        if (itr->second.opt_array_index > -1)
        {
            data_struct::Fit_Param<T_real> fp("perror_" + itr->first, 0.0);
            fp.opt_array_index = itr->second.opt_array_index;
            error_params.add_parameter(fp);
        }
    }
    error_params.from_array(perror);
    
    logI << "Covar Matrix : \n"; 
    for (int i = 0; i < fitp_arr.size(); ++i)
    {
        for (int j = 0; j < fitp_arr.size(); ++j)
        {
            logit_s << covar[i * fitp_arr.size() + j]<< "   ";
        }
        logit_s << "\n";
    }
    fit_params->append_and_update(error_params);

    if (this->_outcome_map.count(this->_last_outcome) > 0)
        return this->_outcome_map[this->_last_outcome];

    return OPTIMIZER_OUTCOME::FAILED;
}

//-----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME NLOPT_Optimizer<T_real>::minimize_func(Fit_Parameters<T_real> *fit_params,
                                                const Spectra<T_real>* const spectra,
                                                const Range energy_range,
                                                const ArrayTr<T_real>* background,
									            Gen_Func_Def<T_real> gen_func,
                                                bool use_weights)
{
    Gen_User_Data<T_real> ud;
    fill_gen_user_data(ud, fit_params, spectra, energy_range, background, gen_func, use_weights);

    std::vector<T_real> fitp_arr = fit_params->to_array();
    if (fitp_arr.size() == 0)
    {
        return OPTIMIZER_OUTCOME::STOPPED;
    }
    std::vector<T_real> perror(fitp_arr.size());
    std::vector<T_real> resid(energy_range.count());
/*
    int info;
    
    _options.maxfev = _options.maxiter * (fitp_arr.size() + 1);

	std::vector<struct mp_par<T_real> > par;
	par.resize(fitp_arr.size());
	if(false == _fill_limits(fit_params, par))
    {
        return OPTIMIZER_OUTCOME::FAILED;
    }

  
    info = mpfit(gen_residuals_mpfit<T_real>, (int)energy_range.count(), (int)fitp_arr.size(), &fitp_arr[0], &par[0], &_options, (void*)&ud, &result);
    
    this->_last_outcome = info;
*/
    fit_params->from_array(fitp_arr);

    T_real sum_resid = 0.0;
    for (size_t i = 0; i < energy_range.count(); i++)
    {
        sum_resid += resid[i];
    }
/*
    if (fit_params->contains(STR_NUM_ITR))
    {
        (*fit_params)[STR_NUM_ITR].value = result.nfev;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, result.nfev));
    }

    if (fit_params->contains(STR_RESIDUAL))
    {
        (*fit_params)[STR_RESIDUAL].value = sum_resid;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_RESIDUAL, sum_resid));
    }

    if (fit_params->contains(STR_CHISQUARE))
    {
        (*fit_params)[STR_CHISQUARE].value = result.bestnorm;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_CHISQUARE, result.bestnorm));
    }

    if (fit_params->contains(STR_CHISQRED))
    {
        (*fit_params)[STR_CHISQRED].value = result.bestnorm / fitp_arr.size();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_CHISQRED, result.bestnorm / fitp_arr.size()));
    }
*/
    if (fit_params->contains(STR_FREE_PARS))
    {
        (*fit_params)[STR_FREE_PARS].value = fitp_arr.size();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_FREE_PARS, fitp_arr.size()));
    }
    // add perror_ fit params
    Fit_Parameters<T_real> error_params;
    for (typename std::unordered_map<std::string, Fit_Param<T_real>>::const_iterator itr = fit_params->begin(); itr != fit_params->end(); itr++)
    {
        if (itr->second.opt_array_index > -1)
        {
            data_struct::Fit_Param<T_real> fp("perror_" + itr->first, 0.0);
            fp.opt_array_index = itr->second.opt_array_index;
            error_params.add_parameter(fp);
        }
    }
    error_params.from_array(perror);

    fit_params->append_and_update(error_params);

/*
    if (this->_outcome_map.count(info) > 0)
        return this->_outcome_map[info];
*/
    return OPTIMIZER_OUTCOME::FAILED;
}

//-----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME NLOPT_Optimizer<T_real>::minimize_quantification(Fit_Parameters<T_real> *fit_params,
                                                          std::unordered_map<std::string, Element_Quant<T_real>*> * quant_map,
                                                          quantification::models::Quantification_Model<T_real>* quantification_model)
{
    Quant_User_Data<T_real> ud;

    assert(quant_map != nullptr);
    
    for (const auto& itr : *quant_map)
    {
        ud.quant_map[itr.first] = *(itr.second);
    }
    
    ud.quantification_model = quantification_model;
    ud.fit_parameters = fit_params;

    std::vector<T_real> fitp_arr = fit_params->to_array();
    if (fitp_arr.size() == 0 || ud.quant_map.size() == 0)
    {
        return OPTIMIZER_OUTCOME::STOPPED;
    }
    std::vector<T_real> perror(fitp_arr.size());
    std::vector<T_real> resid(ud.quant_map.size());
/*
    _options.maxfev = _options.maxiter * (fitp_arr.size() + 1);
	
	std::vector<struct mp_par<T_real> > par;
	par.resize(fitp_arr.size());
	if(false == _fill_limits(fit_params, par))
    {
        return OPTIMIZER_OUTCOME::FAILED;
    }

    this->_last_outcome = mpfit(quantification_residuals_mpfit<T_real>, (int)ud.quant_map.size(), (int)fitp_arr.size(), &fitp_arr[0], &par[0], &_options, (void *) &ud, &result);
    logI << "\nOutcome: " << optimizer_outcome_to_str(this->_outcome_map[this->_last_outcome]) << "\nNum iter: " << result.niter << "\n Norm of the residue vector: " << *result.resid << "\n";

    logI << detailed_outcome(this->_last_outcome) << "\n";
*/
    fit_params->from_array(fitp_arr);

    T_real sum_resid = 0.0;
    for (size_t i = 0; i < quant_map->size(); i++)
    {
        sum_resid += resid[i];
    }
/*
    if (fit_params->contains(STR_NUM_ITR))
    {
        (*fit_params)[STR_NUM_ITR].value = result.nfev;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, result.nfev));
    }

    if (fit_params->contains(STR_RESIDUAL))
    {
        (*fit_params)[STR_RESIDUAL].value = sum_resid;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_RESIDUAL, sum_resid));
    }

    if (fit_params->contains(STR_CHISQUARE))
    {
        (*fit_params)[STR_CHISQUARE].value = result.bestnorm;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_CHISQUARE, result.bestnorm));
    }

    if (fit_params->contains(STR_CHISQRED))
    {
        (*fit_params)[STR_CHISQRED].value = result.bestnorm / fitp_arr.size();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_CHISQRED, result.bestnorm / fitp_arr.size()));
    }
*/
    if (fit_params->contains(STR_FREE_PARS))
    {
        (*fit_params)[STR_FREE_PARS].value = fitp_arr.size();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_FREE_PARS, fitp_arr.size()));
    }

    if (this->_outcome_map.count(this->_last_outcome) > 0)
        return this->_outcome_map[this->_last_outcome];

    return OPTIMIZER_OUTCOME::FAILED;

}

//-----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT NLOPT_Optimizer<float>;
TEMPLATE_CLASS_DLL_EXPORT NLOPT_Optimizer<double>;

} //namespace optimizers
} //namespace fitting
