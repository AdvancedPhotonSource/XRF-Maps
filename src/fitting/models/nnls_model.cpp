/***

Copyright (c) 2016 Arthur Glowacki

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

***/

#include "nnls_model.h"

//debug
#include <iostream>

namespace fitting
{
namespace models
{

NNLS_Model::NNLS_Model() : Gauss_Matrix_Model()
{

    _update_element_guess_value = false;
    _counts_log_10 = false;
    _snip_background = false;
    _fitmatrix = nullptr;
    _max_iter = 100;

}

// ----------------------------------------------------------------------------

NNLS_Model::NNLS_Model(size_t max_iter) : Gauss_Matrix_Model()
{

    _update_element_guess_value = false;
    _counts_log_10 = false;
    _snip_background = false;
    _fitmatrix = nullptr;
    _max_iter = max_iter;

}

// ----------------------------------------------------------------------------

NNLS_Model::~NNLS_Model()
{
    if(_fitmatrix != nullptr)
    {
        delete _fitmatrix;
        _fitmatrix = nullptr;
    }
}

// ----------------------------------------------------------------------------

void NNLS_Model::initialize(Fit_Parameters *fit_params,
                               const Calibration_Standard * const calibration,
                               const Fit_Element_Map_Dict * const elements_to_fit,
                               const struct Range energy_range)
{

    Base_Model::initialize(fit_params, calibration, elements_to_fit, energy_range);

    unordered_map<string, Spectra> element_models = _generate_element_models(fit_params, calibration, elements_to_fit, energy_range);

    _generate_fitmatrix(fit_params, &element_models, energy_range);

}

// ----------------------------------------------------------------------------

void NNLS_Model::_generate_fitmatrix(Fit_Parameters *fit_params,
                                        const unordered_map<string, Spectra> * const element_models,
                                        struct Range energy_range)
{

    if(_fitmatrix != nullptr)
    {
        delete _fitmatrix;
        _fitmatrix = nullptr;
    }
    _fitmatrix = new nsNNLS::denseMatrix(energy_range.count(), element_models->size());

    int i = 0;
    for(const auto& itr : *element_models)
    {
        //Spectra element_model = itr.second;
        for (int j=0; j<itr.second.size(); j++)
        {
            _fitmatrix->set(j,i,itr.second[j]);
        }
        //save element index for later
        (*fit_params)[itr.first].opt_array_index = i;
        i++;
    }

}

// ----------------------------------------------------------------------------

Spectra NNLS_Model::model_spectrum(const Fit_Parameters * const fit_params,
                                      const Spectra * const spectra,
                                      const Calibration_Standard * const calibration,
                                      const Fit_Element_Map_Dict * const elements_to_fit,
                                      const struct Range energy_range)
{
    //dummy function
    return *spectra;
}

// ----------------------------------------------------------------------------

void NNLS_Model::_fit_spectra(Fit_Parameters *fit_params,
                                 const Spectra * const spectra,
                                 const Calibration_Standard * const calibration,
                                 const Fit_Element_Map_Dict * const elements_to_fit)
{

    nsNNLS::nnls*   solver;
    nsNNLS::vector *result;
    int num_iter;

    nsNNLS::vector rhs(spectra->size(), (real_t*)&(*spectra)[0]);

    solver = new nsNNLS::nnls(_fitmatrix, &rhs, _max_iter);

    num_iter = solver->optimize();
    if (num_iter < 0)
    {
        std::cout<<"Error  NNLS_Model::_fit_spectra: in optimization routine"<<std::endl;
    }

    result = solver->getSolution();

    real_t *result_p = result->getData();

    for(const auto& itr : *elements_to_fit)
    {
        Fit_Param param = (*fit_params)[itr.first];
        (*fit_params)[itr.first].value = result_p[param.opt_array_index];
    }

    if (fit_params->contains(data_struct::xrf::STR_NUM_ITR) )
    {
        (*fit_params)[data_struct::xrf::STR_NUM_ITR].value = num_iter;
    }

    delete solver;

}

} //namespace models
} //namespace fitting
