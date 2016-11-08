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

#include "lmmin.h"


using namespace data_struct::xrf;


namespace fitting
{
namespace optimizers
{


static real_t saved_sum = 0.0;

void residuals_lmfit( const double *par, int m_dat, const void *data, double *fvec, int *userbreak )
{

    User_Data* ud = (User_Data*)(data);

    //Update fit parameters from optimizer
    ud->fit_parameters->from_array(par, m_dat);
    //Model spectra based on new fit parameters
    Spectra spectra_model = ud->fit_model->model_spectrum(ud->fit_parameters, ud->spectra, ud->detector, ud->elements, *(ud->energy_range));
    //Calculate residuals
    std::valarray<real_t> residuals = ( (*ud->spectra) - spectra_model ) * (*ud->weights);

    for (int i = 0; i < m_dat; i++ )
    {
        fvec[i] = residuals[i];
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
                               const Detector * const detector,
                               const Fit_Element_Map_Dict * const elements_to_fit,
                               Base_Model* model)
{

    User_Data ud;

    ud.fit_model = model;
    // set spectra to fit
    ud.spectra = (Spectra*)spectra;
    ud.fit_parameters = fit_params;
    ud.detector = (Detector*)detector;
    ud.elements = (Fit_Element_Map_Dict *)elements_to_fit;

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


    lm_status_struct status;
    lm_control_struct control = lm_control_double;

    //control.ftol = 1.0e-10;
    /* Relative error desired in the sum of squares.
                         Termination occurs when both the actual and
                         predicted relative reductions in the sum of squares
                         are at most ftol. */
    //control.xtol = 1.0e-10;
    /* Relative error between last two approximations.
                         Termination occurs when the relative error between
                         two consecutive iterates is at most xtol. */
    //control.gtol = 1.0e-10;
    /* Orthogonality desired between fvec and its derivs.
                         Termination occurs when the cosine of the angle
                         between fvec and any column of the Jacobian is at
                         most gtol in absolute value. */
    //control.epsilon = 2.2204460e-4;
    //control.epsilon = 1.0e-5;
    /* Step used to calculate the Jacobian, should be
                         slightly larger than the relative error in the
                         user-supplied functions. */
    //control.stepbound = 100.0;
    /* Used in determining the initial step bound. This
                         bound is set to the product of stepbound and the
                         Euclidean norm of diag*x if nonzero, or else to
                         stepbound itself. In most cases stepbound should lie
                         in the interval (0.1,100.0). Generally, the value
                         100.0 is recommended. */
    //control.patience = 1000;
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
    lmmin( fitp_arr.size(), &fitp_arr[0], spectra->size(), (const void*) &ud, residuals_lmfit, &control, &status );
	printf(".");
    /* print results */
    //printf( "\nResults:\n" );
    //printf( "status after %d function evaluations:\n  %s\n",  status.nfev, lm_infmsg[status.outcome] );
    //if(residuals_count_lmfit < 100)
     //   std::cout<<"=-=-=-=-=-=-=-=-=-=--==--=-=-==--==---==--==--== bad fit =-=-=-=-=-=-=-=-=-=-=-=-=-=-= "<<residuals_count_lmfit<<std::endl;

    //std::cout<<"residuals count = "<<residuals_count_lmfit<<std::endl;
    fit_params->from_array(fitp_arr);

    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = status.nfev;
    }

}


} //namespace optimizers
} //namespace fitting
