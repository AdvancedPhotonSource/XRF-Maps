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



#ifndef Optimizer_H
#define Optimizer_H

#include <functional>
#include "data_struct/fit_parameters.h"
#include "fitting/models/base_model.h"
#include "quantification/models/quantification_model.h"


typedef std::function<void(size_t, size_t)> Callback_Func_Status_Def;


namespace fitting
{
namespace optimizers
{

using namespace data_struct;
using namespace fitting::models;

#define STR_OPT_FTOL "ftol"
#define STR_OPT_XTOL "xtol"
#define STR_OPT_GTOL "gtol"
#define STR_OPT_EPSILON "epsilon"
#define STR_OPT_STEP "stepbound"
//LM
#define STR_OPT_SCALE_DIAG "scale_diag"
#define STR_OPT_MAXITER "maxiter"
#define STR_OPT_VERBOSE_LEVEL "verbose_level"
//MP
#define STR_OPT_COVTOL "covtol"


//typedef std::function<void(const Fit_Parameters * const, const Range * const, Spectra*)> Gen_Func_Def;
template<typename T_real>
using Gen_Func_Def = std::function<void(const Fit_Parameters<T_real>* const, const Range* const, Spectra<T_real>*)>;

enum class OPTIMIZER_OUTCOME{ FOUND_ZERO, CONVERGED, TRAPPED,  EXHAUSTED, FAILED, CRASHED, EXPLODED, STOPPED, FOUND_NAN, F_TOL_LT_TOL, X_TOL_LT_TOL, G_TOL_LT_TOL};

DLL_EXPORT std::string optimizer_outcome_to_str(OPTIMIZER_OUTCOME outcome);

/**
 * @brief The User_Data struct : Structure used by minimize function for optimizers
 */
template<typename T_real>
struct User_Data
{
    User_Data()
    {
        fit_model = nullptr;
        fit_parameters = nullptr;
        elements = nullptr;
        orig_spectra = nullptr;
        status_callback = nullptr;
        cur_itr = 0;
        total_itr = 0;
    }

    Base_Model<T_real>* fit_model;
    Spectra<T_real> spectra;
	ArrayTr<T_real> weights;
    Fit_Parameters<T_real>*fit_parameters;
	ArrayTr<T_real> spectra_background;
    Fit_Element_Map_Dict<T_real> *elements;
    Range energy_range;
    Spectra<T_real>  spectra_model;
    const Spectra<T_real> *orig_spectra;
    T_real normalizer;
    Callback_Func_Status_Def* status_callback;
    size_t cur_itr;
    size_t total_itr;
};

template<typename T_real>
struct Gen_User_Data
{
    Gen_User_Data()
    {
        fit_parameters = nullptr;
    }

    Spectra<T_real> spectra;
	ArrayTr<T_real> weights;
    Fit_Parameters<T_real>*fit_parameters;
	ArrayTr<T_real> spectra_background;
    Range energy_range;
	Gen_Func_Def<T_real> func;
	Spectra<T_real>  spectra_model;
};

template<typename T_real>
struct Quant_User_Data
{
    Quant_User_Data()
    {
        quantification_model = nullptr;
        fit_parameters = nullptr;
    }

    quantification::models::Quantification_Model<T_real>* quantification_model;
    Fit_Parameters<T_real>* fit_parameters;
    std::unordered_map<std::string, Element_Quant<T_real>> quant_map;
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template<typename T_real>
void fill_user_data(User_Data<T_real> &ud,
                    Fit_Parameters<T_real>*fit_params,
                    const Spectra<T_real>* const spectra,
                    const Fit_Element_Map_Dict<T_real> * const elements_to_fit,
                    const Base_Model<T_real>* const model,
                    const Range energy_range,
                    Callback_Func_Status_Def* status_callback,
                    size_t total_itr,
                    bool use_weights = true)
{
    ud.fit_model = (Base_Model<T_real>*)model;
    // set spectra to fit
    ud.spectra = spectra->sub_spectra(energy_range.min, energy_range.count());
    data_struct::Spectra<T_real> sqr_spec = ud.spectra * ud.spectra; // square the spectra and sum it
    ud.normalizer = sqr_spec.sum();
    //not allocating memory. see https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
    //new (&(ud.spectra)) Eigen::Map<const ArrayTr<T_real>>(spectra->data() + energy_range.min, energy_range.count());
    ud.orig_spectra = spectra;
    ud.fit_parameters = fit_params;
    ud.elements = (Fit_Element_Map_Dict<T_real> *)elements_to_fit;
    ud.energy_range.min = energy_range.min;
    ud.energy_range.max = energy_range.max;

    ud.status_callback = status_callback;
    ud.cur_itr = 0;
    ud.total_itr = total_itr * fit_params->size_non_fixed();

    if (use_weights)
    {
        
        ArrayTr<T_real> weights = (T_real)1.0 / ((T_real)1.0 + (*spectra));
        weights = convolve1d(weights, 5);
        weights = Eigen::abs(weights);
        weights /= weights.maxCoeff();
        weights = weights.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)1.0; });
        ud.weights = weights.segment(energy_range.min, energy_range.count());
        /*
        ArrayTr<T_real> weights = (*spectra);
        weights = weights.log10();
        weights = weights.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });
        weights = convolve1d(weights, 5);
        weights = Eigen::abs(weights);
        weights /= weights.maxCoeff();
        ud.weights = weights.segment(energy_range.min, energy_range.count());
        */
    }
    else
    {
        ud.weights.resize(energy_range.count());
        ud.weights.fill(1.0);
    }

    ArrayTr<T_real> background(spectra->size());
    background.setZero(spectra->size());
    if (fit_params->contains(STR_SNIP_WIDTH))
    {
        background = snip_background<T_real>(spectra,
            fit_params->value(STR_ENERGY_OFFSET),
            fit_params->value(STR_ENERGY_SLOPE),
            fit_params->value(STR_ENERGY_QUADRATIC),
            fit_params->value(STR_SNIP_WIDTH),
            energy_range.min,
            energy_range.max);
    }
    ud.spectra_background = background.segment(energy_range.min, energy_range.count());
    ud.spectra_background = ud.spectra_background.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });
    ud.spectra_model.resize(energy_range.count());
}

//----------------------------------------------------------------------------

template<typename T_real>
void fill_gen_user_data(Gen_User_Data<T_real>& ud,
                        Fit_Parameters<T_real>* fit_params,
                        const Spectra<T_real>* const spectra,
                        const Range energy_range,
                        const ArrayTr<T_real>* background,
                        Gen_Func_Def<T_real> gen_func,
                        bool use_weights = true)
{
    ud.func = gen_func;
    // set spectra to fit
    ud.spectra = spectra->sub_spectra(energy_range.min, energy_range.count());;
    //not allocating memory. see https://eigen.tuxfamily.org/dox/group__TutorialMapClass.html
    //new (&ud.spectra) Eigen::Map<const ArrayTr<T_real>>(spectra->data() + energy_range.min, energy_range.count());
    ud.fit_parameters = fit_params;
    ud.energy_range.min = energy_range.min;
    ud.energy_range.max = energy_range.max;

    if (use_weights)
    {
        ArrayTr<T_real> weights = (T_real)1.0 / ((T_real)1.0 + (*spectra));
        weights = convolve1d(weights, 5);
        weights = Eigen::abs(weights);
        weights /= weights.maxCoeff();
        weights = weights.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });
        ud.weights = weights.segment(energy_range.min, energy_range.count());
    }
    else
    {
        ud.weights.resize(energy_range.count());
        ud.weights.fill(1.0);
    }
    ud.spectra_background = *background;
    ud.spectra_background = ud.spectra_background.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });

    ud.spectra_model.resize(energy_range.count());
}

//----------------------------------------------------------------------------

template<typename T_real>
void update_background_user_data(User_Data<T_real> *ud)
{
    if (ud->fit_parameters->contains(STR_SNIP_WIDTH))
    {
        Fit_Param<T_real> fit_snip_width = ud->fit_parameters->at(STR_SNIP_WIDTH);
        if (fit_snip_width.bound_type != E_Bound_Type::FIXED && ud->orig_spectra != nullptr)
        {
            //ud->spectra_background = snip_background(ud->orig_spectra,
            ArrayTr<T_real> background = snip_background<T_real>(ud->orig_spectra,
                ud->fit_parameters->value(STR_ENERGY_OFFSET),
                ud->fit_parameters->value(STR_ENERGY_SLOPE),
                ud->fit_parameters->value(STR_ENERGY_QUADRATIC),
                fit_snip_width.value,
                ud->energy_range.min,
                ud->energy_range.max);

            ud->spectra_background = background.segment(ud->energy_range.min, ud->energy_range.count());
            ud->spectra_background = ud->spectra_background.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });

        }
    }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

/**
 * @brief The Optimizer class : Base class for error minimization to find optimal specta model
 */
template<typename T_real>
class DLL_EXPORT Optimizer
{
public:
    Optimizer()
    {
        _last_outcome = -1;
    }

    virtual ~Optimizer(){}

    virtual OPTIMIZER_OUTCOME minimize(Fit_Parameters<T_real> *fit_params,
                          const Spectra<T_real>* const spectra,
                          const Fit_Element_Map_Dict<T_real> * const elements_to_fit,
                          const Base_Model<T_real>* const model,
                          const Range energy_range,
                          bool use_weights,
                          Callback_Func_Status_Def* status_callback = nullptr) = 0;

    virtual OPTIMIZER_OUTCOME minimize_func(Fit_Parameters<T_real>*fit_params,
                               const Spectra<T_real>* const spectra,
                               const Range energy_range,
                               const ArrayTr<T_real>* background,
                               Gen_Func_Def<T_real> gen_func,
                               bool use_weights) = 0;


    virtual OPTIMIZER_OUTCOME minimize_quantification(Fit_Parameters<T_real>*fit_params,
                                         std::unordered_map<std::string, Element_Quant<T_real>*> * quant_map,
                                         quantification::models::Quantification_Model<T_real>* quantification_model) = 0;

    virtual std::unordered_map<std::string, T_real> get_options() = 0;

    virtual void set_options(std::unordered_map<std::string, T_real> opt) = 0;

    virtual std::string detailed_outcome(int outcome) = 0;

    std::string get_last_detailed_outcome() {return detailed_outcome(_last_outcome);}

protected:

    std::map<int, OPTIMIZER_OUTCOME> _outcome_map;

    int _last_outcome;

};

} //namespace optimizers

} //namespace fitting

#endif // Optimizer
