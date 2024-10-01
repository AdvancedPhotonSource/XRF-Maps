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

double gen_residuals_nlopt(const std::vector<double> &x, std::vector<double> &grad, void *usr_data)
{
    // Get user passed data
    Gen_User_Data<double>* ud = static_cast<Gen_User_Data<double>*>(usr_data);

    // Update fit parameters from optimizer
    ud->fit_parameters->from_array_d(x);

    // Model spectra based on new fit parameters
    ud->func(ud->fit_parameters, &(ud->energy_range), &(ud->spectra_model));
    // Add background
    ud->spectra_model += ud->spectra_background;
    
    double sum = 0.0;
    // Calculate residuals
    for (int i=0; i<ud->spectra.size(); i++)
    {
        sum += pow((ud->spectra[i] - ud->spectra_model[i]), 2.0) * ud->weights[i];
    }

    return sum;
}


//-----------------------------------------------------------------------------

double quantification_residuals_nlopt(const std::vector<double> &x, std::vector<double> &grad, void *usr_data)
{
    ///(std::valarray<T_real> p, std::valarray<T_real> y, std::valarray<T_real> x)

    //y is array of elements standards
    //x is indexes of elements in standard
    //p is array size 2 but seems only first index is used
    ///return (y - this->fit_calibrationcurve(x, p));

    Quant_User_Data<double>* ud = (Quant_User_Data<double>*)(usr_data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array_d(x);
    
    //Calculate residuals
    std::unordered_map<std::string, double> result_map = ud->quantification_model->model_calibrationcurve(ud->quant_map, x[0]);

    double sum = 0.0;
    int idx = 0;
    for(auto& itr : ud->quant_map)
    {
		if (std::isfinite(result_map[itr.first]) == false)
		{
            logE<<"Quantification reuslted in NaN or Inf! "<< itr.first<<" : "<<result_map[itr.first]<<"\n";
			sum += itr.second.e_cal_ratio * 100.0;
		}
		else
		{
			sum += pow((itr.second.e_cal_ratio - result_map[itr.first]), 2.0);
		}
        idx++;
    }

    return sum;
}

// =====================================================================================================================

template<typename T_real>
NLOPT_Optimizer<T_real>::NLOPT_Optimizer() : Optimizer<T_real>()
{

    this->_last_outcome = -1;
/*
    //_options { 1e-10, 1e-10, 1e-10, MP_MACHEP0, 100.0, 1.0e-14, 2000, 0, 0, 0, 0, 0 };
    _options.maxiter = 1000;        
*/
    this->_outcome_map[nlopt::FAILURE] = OPTIMIZER_OUTCOME::FAILED;
    this->_outcome_map[nlopt::INVALID_ARGS] = OPTIMIZER_OUTCOME::FAILED;
    this->_outcome_map[nlopt::OUT_OF_MEMORY] = OPTIMIZER_OUTCOME::FAILED;
    this->_outcome_map[nlopt::ROUNDOFF_LIMITED] = OPTIMIZER_OUTCOME::FAILED;
    this->_outcome_map[nlopt::SUCCESS] = OPTIMIZER_OUTCOME::CONVERGED;
    this->_outcome_map[nlopt::STOPVAL_REACHED] = OPTIMIZER_OUTCOME::CONVERGED;
    this->_outcome_map[nlopt::FTOL_REACHED] = OPTIMIZER_OUTCOME::F_TOL_LT_TOL;
    this->_outcome_map[nlopt::XTOL_REACHED] = OPTIMIZER_OUTCOME::X_TOL_LT_TOL;
    this->_outcome_map[nlopt::MAXEVAL_REACHED] = OPTIMIZER_OUTCOME::EXHAUSTED;
    this->_outcome_map[nlopt::MAXTIME_REACHED] = OPTIMIZER_OUTCOME::EXHAUSTED;
    this->_outcome_map[nlopt::FORCED_STOP] = OPTIMIZER_OUTCOME::STOPPED;

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

template<typename T_real>
std::string NLOPT_Optimizer<T_real>::detailed_outcome(int info)
{
	switch (info)
	{
	case nlopt::FAILURE:
		return "Generic failure code .";
	case nlopt::INVALID_ARGS:
        return "Invalid Args. ";
	case nlopt::OUT_OF_MEMORY:
        return "Out of memory.";
	case nlopt::ROUNDOFF_LIMITED:
        return "Roundoff limited.";
	case nlopt::FORCED_STOP:
        return "Force Stop.";
	case nlopt::SUCCESS:
        return "Success.";	
	case nlopt::STOPVAL_REACHED:
        return "Stop val reached.";
	case nlopt::FTOL_REACHED:
        return "FTol reached.";
    case nlopt::XTOL_REACHED:
        return "XTol reached.";
    case nlopt::MAXEVAL_REACHED:
        return "Max function evaluations reached.";
    case nlopt::MAXTIME_REACHED:
        return "Max time reached.";
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
    std::vector<double> step_arr;
    
    fit_params->to_array_with_bounds(fitp_arr, lb_arr, ub_arr, step_arr);
    if (fitp_arr.size() == 0)
    {
        return OPTIMIZER_OUTCOME::STOPPED;
    }

    size_t total_itr = 20000; //num_itr
    fill_user_data(ud, fit_params, spectra, elements_to_fit, model, energy_range, status_callback, total_itr, use_weights);

//  LN_NEWUOA_BOUND
//  LN_COBYLA
// LN_NELDERMEAD
//  LN_SBPLX

    nlopt::opt opt(nlopt::LN_SBPLX, fitp_arr.size());
    opt.set_lower_bounds(lb_arr);
    opt.set_upper_bounds(ub_arr);
    opt.set_default_initial_step(step_arr);
    opt.set_min_objective(residuals_nlopt, (void*)&ud);
    opt.set_xtol_rel(1e-10);
    opt.set_maxeval(20000);

    double minf;

    try
    {
        nlopt::result result = opt.optimize(fitp_arr, minf);
        logI<<detailed_outcome(result)<< " : resid = "<<minf<<"\n\n";;
        this->_last_outcome = result;
    }
    catch(std::exception &e) 
    {
        logW << "nlopt failed: " << e.what() << "\n";
    }
    fit_params->from_array_d(fitp_arr);

    /////////////// init params limits /////////////////////////
    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = opt.get_numevals();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, opt.get_numevals()));
    }

    if (fit_params->contains(STR_RESIDUAL))
    {
        (*fit_params)[STR_RESIDUAL].value = minf;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_RESIDUAL, minf));
    }
  
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

    std::vector<double> fitp_arr;
    std::vector<double> lb_arr;
    std::vector<double> ub_arr;
    std::vector<double> step_arr;
    
    fit_params->to_array_with_bounds(fitp_arr, lb_arr, ub_arr, step_arr);
    if (fitp_arr.size() == 0)
    {
        return OPTIMIZER_OUTCOME::STOPPED;
    }

    size_t total_itr = 20000; 

    nlopt::opt opt(nlopt::LN_SBPLX, fitp_arr.size());
    opt.set_lower_bounds(lb_arr);
    opt.set_upper_bounds(ub_arr);
    opt.set_default_initial_step(step_arr);
    opt.set_min_objective(gen_residuals_nlopt, (void*)&ud);
    opt.set_xtol_rel(1e-10);
    opt.set_maxeval(20000);

    double minf;
    nlopt::result result;
    try
    {
        result = opt.optimize(fitp_arr, minf);
        logI<<detailed_outcome(result)<< " : resid = "<<minf<<"\n\n";;
        this->_last_outcome = result;
    }
    catch(std::exception &e) 
    {
        logW << "nlopt failed: " << e.what() << "\n";
    }
    fit_params->from_array_d(fitp_arr);

    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = opt.get_numevals();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, opt.get_numevals()));
    }

    if (fit_params->contains(STR_RESIDUAL))
    {
        (*fit_params)[STR_RESIDUAL].value = minf;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_RESIDUAL, minf));
    }
    
    if (fit_params->contains(STR_FREE_PARS))
    {
        (*fit_params)[STR_FREE_PARS].value = fitp_arr.size();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_FREE_PARS, fitp_arr.size()));
    }
    
    if (this->_outcome_map.count(result) > 0)
        return this->_outcome_map[result];

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

    std::vector<double> fitp_arr;
    std::vector<double> lb_arr;
    std::vector<double> ub_arr;
    std::vector<double> step_arr;
    
    fit_params->to_array_with_bounds(fitp_arr, lb_arr, ub_arr, step_arr);
    
    if (fitp_arr.size() == 0 || ud.quant_map.size() == 0)
    {
        return OPTIMIZER_OUTCOME::STOPPED;
    }
    
    nlopt::opt opt(nlopt::LN_SBPLX, fitp_arr.size());
    opt.set_lower_bounds(lb_arr);
    opt.set_upper_bounds(ub_arr);
    opt.set_default_initial_step(step_arr);
    opt.set_min_objective(quantification_residuals_nlopt, (void*)&ud);
    opt.set_xtol_rel(1e-10);
    opt.set_maxeval(20000);

    double minf;
    nlopt::result result;

    try
    {
        result = opt.optimize(fitp_arr, minf);
        logI<<detailed_outcome(result)<< " : resid = "<<minf<<"\n\n";;
        this->_last_outcome = result;
    }
    catch(std::exception &e) 
    {
        logW << "nlopt failed: " << e.what() << "\n";
    }
    fit_params->from_array_d(fitp_arr);

    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = opt.get_numevals();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, opt.get_numevals()));
    }

    if (fit_params->contains(STR_RESIDUAL))
    {
        (*fit_params)[STR_RESIDUAL].value = minf;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_RESIDUAL, minf));
    }
    
    if (fit_params->contains(STR_FREE_PARS))
    {
        (*fit_params)[STR_FREE_PARS].value = fitp_arr.size();
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_FREE_PARS, fitp_arr.size()));
    }

    if (this->_outcome_map.count(result) > 0)
        return this->_outcome_map[result];

    return OPTIMIZER_OUTCOME::FAILED;

}

//-----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT NLOPT_Optimizer<float>;
TEMPLATE_CLASS_DLL_EXPORT NLOPT_Optimizer<double>;

} //namespace optimizers
} //namespace fitting
