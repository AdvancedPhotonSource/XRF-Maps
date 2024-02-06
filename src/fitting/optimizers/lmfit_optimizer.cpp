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

using namespace data_struct;


namespace fitting
{
namespace optimizers
{

// ----------------------------------------------------------------------------

template<typename T_real>
void residuals_lmfit( const T_real *par, int m_dat, const void *data, T_real *fvec, int *userbreak )
{
    User_Data<T_real>* ud = (User_Data<T_real>*)(data);

    // Update fit parameters from optimizer
    ud->fit_parameters->from_array(par, m_dat);
    // Model spectra based on new fit parameters
    update_background_user_data(ud);
    ud->spectra_model = ud->fit_model->model_spectrum_mp(ud->fit_parameters, ud->elements, ud->energy_range);
    // Add background
    ud->spectra_model += ud->spectra_background;
    // Remove nan's and inf's
    // ud->spectra_model = (ArrayTr<T_real>)ud->spectra_model.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });

    // Calculate residuals
    for (int i = 0; i < m_dat; i++ )
    {
		fvec[i] = abs(ud->spectra[i] - ud->spectra_model[i]) * ud->weights[i];
		if (std::isfinite(fvec[i]) == false)
		{
			//logE << "\n\n\n";
			logE << "Spectra[i] = " << ud->spectra[i] << " :: spectra_model[i] = " << ud->spectra_model[i] << "  ::  weights[i] = " << ud->weights[i];
			//logE << "\n\n\n";
			//fvec[i] = ud->spectra[i];
            //fvec[i] = std::numeric_limits<T_real>::quiet_NaN();
		}
    }
    
    ud->cur_itr++;
    if (ud->status_callback != nullptr)
    {
        try
        {
            (*ud->status_callback)(ud->cur_itr, ud->total_itr);
        }
        catch (...)
        {
            logI << "Cancel fitting" << std::endl;
            *userbreak = 1;
        }
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
void general_residuals_lmfit( const T_real *par, int m_dat, const void *data, T_real *fvec, int *userbreak )
{

    Gen_User_Data<T_real>* ud = (Gen_User_Data<T_real>*)(data);

    // Update fit parameters from optimizer
    ud->fit_parameters->from_array(par, m_dat);
    // Model spectra based on new fit parameters
    ud->func(ud->fit_parameters, &(ud->energy_range), &(ud->spectra_model));
    // Add background
    ud->spectra_model += ud->spectra_background;
    // Remove nan's and inf's
	// Used to check for nan's here but there were some cases where the optimizer would return nan found. So moved to after subract of model
    // ud->spectra_model = (ArrayTr<T_real>)ud->spectra_model.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });
    
    // Calculate residuals
    for (int i = 0; i < m_dat; i++ )
    {
        fvec[i] = abs( ud->spectra[i] - ud->spectra_model[i] ) * ud->weights[i];
        /*
		if (std::isfinite(fvec[i]) == false)
		{
			fvec[i] = ud->spectra[i];
		}
        */
    }
    

}

//-----------------------------------------------------------------------------

template<typename T_real>
void quantification_residuals_lmfit( const T_real *par, int m_dat, const void *data, T_real *fvec, int *userbreak )
{
    ///(std::valarray<T_real> p, std::valarray<T_real> y, std::valarray<T_real> x)

    //y is array of elements standards
    //x is indexes of elements in standard
    //p is array size 2 but seems only first index is used
    ///return (y - this->fit_calibrationcurve(x, p));

    Quant_User_Data<T_real>* ud = (Quant_User_Data<T_real>*)(data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(par, m_dat);
    //Model spectra based on new fit parameters

    //Calculate residuals
    std::unordered_map<std::string, T_real> result_map = ud->quantification_model->model_calibrationcurve(ud->quant_map, par[0]);

    int idx = 0;
    for(auto& itr : ud->quant_map)
    {
		if (std::isfinite(result_map[itr.first]) == false)
		{
			fvec[idx] = itr.second.e_cal_ratio;
		}
		else
		{
			fvec[idx] = abs(itr.second.e_cal_ratio - result_map[itr.first]);
		}
        idx++;
    }
}


// =====================================================================================================================

template<typename T_real>
LMFit_Optimizer<T_real>::LMFit_Optimizer() : Optimizer<T_real>()
{
    this->_last_outcome = -1;

    if (std::is_same<T_real, float>::value)
    {
        _options.ftol = FP_LM_USERTOL; //LM_USERTOL; // Relative error desired in the sum of squares. Termination occurs when both the actualand predicted relative reductions in the sum of squares are at most ftol.
        _options.xtol = FP_LM_USERTOL; // Relative error between last two approximations. Termination occurs when the relative error between two consecutive iterates is at most xtol.
        _options.gtol = FP_LM_USERTOL;  //LM_USERTOL; // Orthogonality desired between fvec and its derivs. Termination occurs when the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value.
        _options.epsilon = FP_LM_EPSILON; // Step used to calculate the Jacobian, should be slightly larger than the relative error in the user-supplied functions.

    }
    else if (std::is_same<T_real, double>::value)
    {
        _options.ftol = DP_LM_USERTOL; //LM_USERTOL; // Relative error desired in the sum of squares. Termination occurs when both the actualand predicted relative reductions in the sum of squares are at most ftol.
        _options.xtol = DP_LM_USERTOL; // Relative error between last two approximations. Termination occurs when the relative error between two consecutive iterates is at most xtol.
        _options.gtol = DP_LM_USERTOL;  //LM_USERTOL; // Orthogonality desired between fvec and its derivs. Termination occurs when the cosine of the angle between fvec and any column of the Jacobian is at most gtol in absolute value.
        _options.epsilon = DP_LM_EPSILON; // Step used to calculate the Jacobian, should be slightly larger than the relative error in the user-supplied functions.
    }
	
    _options.stepbound = (T_real)100.; // Used in determining the initial step bound. This bound is set to the product of stepbound and the Euclidean norm of diag*x if nonzero, or else to stepbound itself. In most cases stepbound should lie in the interval (0.1,100.0). Generally, the value 100.0 is recommended.
    _options.patience = 1000; // Used to set the maximum number of function evaluations to patience*(number_of_parameters+1).
    _options.scale_diag = 1; // If 1, the variables will be rescaled internally. Recommended value is 1.
    _options.msgfile = NULL; //  Progress messages will be written to this file.
    _options.verbosity = 0; //  OR'ed: 1: print some messages; 2: print Jacobian. 
    _options.n_maxpri = -1; // -1, or max number of parameters to print.
    _options.m_maxpri = -1; // -1, or max number of residuals to print. 


    this->_outcome_map[0] = OPTIMIZER_OUTCOME::FOUND_ZERO;
    this->_outcome_map[1] = OPTIMIZER_OUTCOME::CONVERGED;
    this->_outcome_map[2] = OPTIMIZER_OUTCOME::CONVERGED;
    this->_outcome_map[3] = OPTIMIZER_OUTCOME::CONVERGED;
    this->_outcome_map[4] = OPTIMIZER_OUTCOME::TRAPPED;
    this->_outcome_map[5] = OPTIMIZER_OUTCOME::F_TOL_LT_TOL;
    this->_outcome_map[6] = OPTIMIZER_OUTCOME::X_TOL_LT_TOL;
    this->_outcome_map[7] = OPTIMIZER_OUTCOME::G_TOL_LT_TOL;
    this->_outcome_map[8] = OPTIMIZER_OUTCOME::CRASHED;
    this->_outcome_map[9] = OPTIMIZER_OUTCOME::EXPLODED;
    this->_outcome_map[10] = OPTIMIZER_OUTCOME::STOPPED;
    this->_outcome_map[11] = OPTIMIZER_OUTCOME::FOUND_NAN;

}

// ----------------------------------------------------------------------------

template<typename T_real>
std::unordered_map<std::string, T_real> LMFit_Optimizer<T_real>::get_options()
{
    std::unordered_map<std::string, T_real> opts{
        {STR_OPT_FTOL, _options.ftol},
        {STR_OPT_XTOL, _options.xtol},
        {STR_OPT_GTOL, _options.gtol},
        {STR_OPT_EPSILON, _options.epsilon},
        {STR_OPT_STEP, _options.stepbound},
        {STR_OPT_SCALE_DIAG, _options.scale_diag},
        {STR_OPT_MAXITER, _options.patience},
        {STR_OPT_VERBOSE_LEVEL, _options.verbosity}
    };
    return opts;
}

// ----------------------------------------------------------------------------

template<typename T_real>
void LMFit_Optimizer<T_real>::set_options(std::unordered_map<std::string, T_real> opt)
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
        _options.epsilon = opt.at(STR_OPT_EPSILON);
    }
    if (opt.count(STR_OPT_STEP) > 0)
    {
        _options.stepbound= opt.at(STR_OPT_STEP);
    }
    if (opt.count(STR_OPT_SCALE_DIAG) > 0)
    {
        _options.scale_diag = (int)opt.at(STR_OPT_SCALE_DIAG);
    }
    if (opt.count(STR_OPT_MAXITER) > 0)
    {
        _options.patience = (int)opt.at(STR_OPT_MAXITER);
    }
    if (opt.count(STR_OPT_VERBOSE_LEVEL) > 0)
    {
        _options.verbosity = (int)opt.at(STR_OPT_VERBOSE_LEVEL);
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME LMFit_Optimizer<T_real>::minimize(Fit_Parameters<T_real> *fit_params,
                                           const Spectra<T_real>* const spectra,
                                           const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                           const Base_Model<T_real>* const model,
                                           const Range energy_range,
                                           bool use_weights,
                                           Callback_Func_Status_Def* status_callback)
{

    User_Data<T_real> ud;
    std::vector<T_real> fitp_arr = fit_params->to_array();
    if (fitp_arr.size() == 0)
    {
        return OPTIMIZER_OUTCOME::STOPPED;
    }
    std::vector<T_real> perror(fitp_arr.size());

    //_options.verbosity = 2;
    size_t total_itr = _options.patience * (fitp_arr.size() + 1);
    fill_user_data(ud, fit_params, spectra, elements_to_fit, model, energy_range, status_callback, total_itr, use_weights);

    lm_status_struct<T_real> status;

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
    lmmin( (int)fitp_arr.size(), &fitp_arr[0], (int)energy_range.count(), (const void*) &ud, residuals_lmfit, &_options, &status );
    logI<< "Outcome: "<<lm_infmsg[status.outcome]<<"\nNum iter: "<<status.nfev<<"\n Norm of the residue vector: "<<status.fnorm<<"\n";
    this->_last_outcome = status.outcome;
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = static_cast<T_real>(status.nfev);
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, static_cast<T_real>(status.nfev)));
    }

    if (fit_params->contains(STR_RESIDUAL))
    {
        (*fit_params)[STR_RESIDUAL].value = status.fnorm;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_RESIDUAL, status.fnorm));
    }

    if (fit_params->contains(STR_OUTCOME))
    {
        if (this->_outcome_map.count(status.outcome) > 0)
            (*fit_params)[STR_OUTCOME].value = (T_real)(this->_outcome_map[status.outcome]);
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_OUTCOME, (T_real)(this->_outcome_map[status.outcome])));
    }

    if(this->_outcome_map.count(status.outcome)>0)
        return this->_outcome_map[status.outcome];

    return OPTIMIZER_OUTCOME::FAILED;

}

// ----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME LMFit_Optimizer<T_real>::minimize_func(Fit_Parameters<T_real>*fit_params,
                                                const Spectra<T_real>* const spectra,
                                                const Range energy_range,
                                                const ArrayTr<T_real>*background,
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

    lm_status_struct<T_real> status;
    _options.verbosity = 0;
    lmmin((int)fitp_arr.size(), &fitp_arr[0], (int)energy_range.count(), (const void*) &ud, general_residuals_lmfit, &_options, &status );
    this->_last_outcome = status.outcome;
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(STR_NUM_ITR) )
    {
        (*fit_params)[STR_NUM_ITR].value = status.nfev;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, static_cast<T_real>(status.nfev)));
    }
    if (fit_params->contains(STR_RESIDUAL))
    {
        //(*fit_params)[STR_RESIDUAL].value = status.fnorm;
        data_struct::ArrayTr<T_real> diff_arr = ud.spectra - ud.spectra_model;
        diff_arr = diff_arr.unaryExpr([](T_real v) { return std::abs(v); });
        (*fit_params)[STR_RESIDUAL].value = diff_arr.sum();
    }

    if (this->_outcome_map.count(status.outcome) > 0)
        return this->_outcome_map[status.outcome];

    return OPTIMIZER_OUTCOME::FAILED;

}

// ----------------------------------------------------------------------------

template<typename T_real>
OPTIMIZER_OUTCOME LMFit_Optimizer<T_real>::minimize_quantification(Fit_Parameters<T_real>*fit_params,
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
    //_options.verbosity = 2;
    std::vector<T_real> fitp_arr = fit_params->to_array();
    if (fitp_arr.size() == 0)
    {
        return OPTIMIZER_OUTCOME::STOPPED;
    }
    std::vector<T_real> perror(fitp_arr.size());

    lm_status_struct<T_real> status;
    lmmin((int)fitp_arr.size(), &fitp_arr[0], (int)ud.quant_map.size(), (const void*) &ud, quantification_residuals_lmfit, &_options, &status );
    logI << "\nOutcome: " << lm_infmsg[status.outcome] << "\nNum iter: " << status.nfev << "\nNorm of the residue vector: " << status.fnorm << "\n";
    this->_last_outcome = status.outcome;
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(STR_NUM_ITR))
    {
        (*fit_params)[STR_NUM_ITR].value = static_cast<T_real>(status.nfev);
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_NUM_ITR, static_cast<T_real>(status.nfev)));
    }

    if (fit_params->contains(STR_RESIDUAL))
    {
        (*fit_params)[STR_RESIDUAL].value = status.fnorm;
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_RESIDUAL, status.fnorm));
    }

    if (fit_params->contains(STR_OUTCOME))
    {
        if (this->_outcome_map.count(status.outcome) > 0)
            (*fit_params)[STR_OUTCOME].value = (T_real)(this->_outcome_map[status.outcome]);
    }
    else
    {
        fit_params->add_parameter(data_struct::Fit_Param<T_real>(STR_OUTCOME, (T_real)(this->_outcome_map[status.outcome])));
    }

    if (this->_outcome_map.count(status.outcome) > 0)
        return this->_outcome_map[status.outcome];

    return OPTIMIZER_OUTCOME::FAILED;
}

// ----------------------------------------------------------------------------

template<typename T_real>
std::string LMFit_Optimizer<T_real>::detailed_outcome(int outcome)
{
    if (outcome > -1 && outcome < 13)
    {
        return  lm_infmsg[outcome];
    }
    return "Unknown";
}

// ----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT LMFit_Optimizer<float>;
TEMPLATE_CLASS_DLL_EXPORT LMFit_Optimizer<double>;


} //namespace optimizers
} //namespace fitting
