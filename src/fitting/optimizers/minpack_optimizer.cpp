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



#include "minpack_optimizer.h"

#include <iostream>
#include <algorithm>
#include <math.h>

#include "minpack.hpp"

#include <string.h>

using namespace data_struct::xrf;


namespace fitting
{
namespace optimizers
{


static int residuals_count_minpack = 0;

void residuals_minpack(void *usr_data, int params_size, real_t *params, real_t *fvec, int *iflag )
{
    // get user passed data
    User_Data* ud = static_cast<User_Data*>(usr_data);

    ud->fit_parameters->from_array(params, params_size);
    Spectra spectra_model = ud->fit_model->model_spectrum(ud->fit_parameters, ud->elements, *(ud->energy_range));

    std::valarray<real_t> err = ( (*ud->spectra) - spectra_model ) * (*ud->weights);

    for(size_t i=0; i<ud->spectra->size(); i++)
    {
        fvec[i] = err[i];
    }
    residuals_count_minpack ++;
}


void gen_residuals_minpack(void *usr_data, int params_size, real_t *params, real_t *fvec, int *iflag )
{
    // get user passed data
    Gen_User_Data* ud = static_cast<Gen_User_Data*>(usr_data);

    ud->fit_parameters->from_array(params, params_size);
    Spectra spectra_model = ud->func(ud->fit_parameters, ud->energy_range);

    std::valarray<real_t> err = ( (*ud->spectra) - spectra_model ) * (*ud->weights);

    for(size_t i=0; i<ud->spectra->size(); i++)
    {
        fvec[i] = err[i];
    }
    //gen_residuals_count_minpack ++;
}

// =====================================================================================================================


MinPack_Optimizer::MinPack_Optimizer() : Optimizer()
{

}

MinPack_Optimizer::~MinPack_Optimizer()
{


}

void MinPack_Optimizer::minimize(Fit_Parameters *fit_params,
                                const Spectra * const spectra,
                                const Fit_Element_Map_Dict * const elements_to_fit,
                                const Base_Model * const model)
{
    //const int params_size = 12;
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

    std::vector<real_t> fitp_arr = fit_params->to_array();
    std::vector<real_t> fvec;
    fvec.resize(spectra->size());
    //std::vector<real_t> perror(fitp_arr.size());


    real_t tol = 1.0e-10;
     //n * ( 3 * n + 13 ) ) / 2
    int lwa = fitp_arr.size() * ( 3 * fitp_arr.size() + 13) / 2 ;
    real_t *wa = new real_t[lwa];
    int info;

    std::valarray<real_t> weights = (real_t)1.0 / ( (real_t)1.0 + (*spectra) );
    weights = convolve1d(weights, 5);
    weights = std::abs(weights);
    weights /= weights.max();
    ud.weights = &weights;

    //std::valarray<real_t> weights = std::sqrt( *(spectra->buffer()) );
    //ud.weights = &weights;
    //      function, user data, fit parameters, result vec, tollerance
    info = hybrd1(&residuals_minpack, &ud, fitp_arr.size(), &fitp_arr[0], &fvec[0], tol, wa, lwa);
    switch(info)
    {
        case OPTIMIZER_INFO::IMPROPER_INPUT:
            std::cout<<"!> Improper input parameters."<<std::endl;
        break;
        case OPTIMIZER_INFO::MOST_TOL:
            std::cout<<"> Algorithm estimates that the relative error in the sum of squares is at most tol. "<<std::endl;
        break;
        case OPTIMIZER_INFO::EXCEED_CALL:
            std::cout<<"> Number of calls to fcn has reached or exceeded 200*(n+1)."<<std::endl;
        break;
        case OPTIMIZER_INFO::TOL_TOO_SMALL:
            std::cout<<">> Tol is too small. no further improvement in the approximate solution x is possible. "<<std::endl;
        break;
        case OPTIMIZER_INFO::NO_PROGRESS:
            std::cout<<"> Fiteration is not making good progress."<<std::endl;
        break;
        default:
            std::cout<<"!> Unknown info status"<<std::endl;
        break;
    }

    delete [] wa;
    std::cout<<"residuals count = "<<residuals_count_minpack<<std::endl;
    fit_params->from_array(fitp_arr);

}

void MinPack_Optimizer::minimize_func(Fit_Parameters *fit_params,
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
    std::vector<real_t> fvec;
    fvec.resize(spectra->size());
    //std::vector<real_t> perror(fitp_arr.size());


    real_t tol = 1.0e-10;
     //n * ( 3 * n + 13 ) ) / 2
    int lwa = fitp_arr.size() * ( 3 * fitp_arr.size() + 13) / 2 ;
    real_t *wa = new real_t[lwa];
    int info;

    std::valarray<real_t> weights = (real_t)1.0 / ( (real_t)1.0 + (*spectra) );
    weights = convolve1d(weights, 5);
    weights = std::abs(weights);
    weights /= weights.max();
    ud.weights = &weights;

    //std::valarray<real_t> weights = std::sqrt( *(spectra->buffer()) );
    //ud.weights = &weights;
    //      function, user data, fit parameters, result vec, tollerance
    //gen_residuals_count_minpack = 0;
    info = hybrd1(&gen_residuals_minpack, &ud, fitp_arr.size(), &fitp_arr[0], &fvec[0], tol, wa, lwa);
    switch(info)
    {
        case OPTIMIZER_INFO::IMPROPER_INPUT:
            std::cout<<"!> Improper input parameters."<<std::endl;
        break;
        case OPTIMIZER_INFO::MOST_TOL:
            std::cout<<"> Algorithm estimates that the relative error in the sum of squares is at most tol. "<<std::endl;
        break;
        case OPTIMIZER_INFO::EXCEED_CALL:
            std::cout<<"> Number of calls to fcn has reached or exceeded 200*(n+1)."<<std::endl;
        break;
        case OPTIMIZER_INFO::TOL_TOO_SMALL:
            std::cout<<">> Tol is too small. no further improvement in the approximate solution x is possible. "<<std::endl;
        break;
        case OPTIMIZER_INFO::NO_PROGRESS:
            std::cout<<"> Fiteration is not making good progress."<<std::endl;
        break;
        default:
            std::cout<<"!> Unknown info status"<<std::endl;
        break;
    }

    delete [] wa;
    //std::cout<<"residuals count = "<<gen_residuals_count_minpack<<std::endl;
    fit_params->from_array(fitp_arr);
/*
    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = gen_residuals_count_minpack;
    }
*/
}


} //namespace optimizers
} //namespace fitting
