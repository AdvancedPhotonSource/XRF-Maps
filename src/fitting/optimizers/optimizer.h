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

#include "fit_parameters.h"
#include "detector.h"
#include "base_model.h"
#include "quantification_model.h"


namespace fitting
{
namespace optimizers
{

using namespace std;
using namespace data_struct::xrf;
using namespace fitting::models;


enum OPTIMIZER_INFO { IMPROPER_INPUT, MOST_TOL, EXCEED_CALL, TOL_TOO_SMALL, NO_PROGRESS };

typedef std::function<void(const Fit_Parameters * const, const Range * const, Spectra*)> Gen_Func_Def;


/**
 * @brief The User_Data struct : Structure used by minimize function for optimizers
 */
struct User_Data
{
    Base_Model* fit_model;
    Spectra *spectra;
    std::valarray<real_t> weights;
    Fit_Parameters *fit_parameters;
    std::valarray<real_t> *spectra_background;
    Fit_Element_Map_Dict *elements;
    Range energy_range;
	std::valarray<real_t> residuals;
	Spectra  spectra_model;
};

struct Gen_User_Data
{
    Spectra *spectra;
    std::valarray<real_t> weights;
    Fit_Parameters *fit_parameters;
    Range energy_range;
	Gen_Func_Def func;
	std::valarray<real_t> residuals;
	Spectra  spectra_model;
};

struct Quant_User_Data
{
    quantification::models::Quantification_Model * quantification_model;
    Fit_Parameters * fit_parameters;
    std::unordered_map<std::string, Element_Quant> * quant_map;
};

void fill_user_data(User_Data &ud,
	Fit_Parameters *fit_params,
	const Spectra * const spectra,
	const Fit_Element_Map_Dict * const elements_to_fit,
	const Base_Model * const model);

	void fill_gen_user_data(Gen_User_Data &ud,
		Fit_Parameters *fit_params,
		const Spectra * const spectra,
		Gen_Func_Def gen_func);

/**
 * @brief The Optimizer class : Base class for error minimization to find optimal specta model
 */
class DLL_EXPORT Optimizer
{
public:
    Optimizer(){}

    ~Optimizer(){}

    virtual void minimize(Fit_Parameters *fit_params,
                          const Spectra * const spectra,
                          const Fit_Element_Map_Dict * const elements_to_fit,
                          const Base_Model * const model) = 0;

    virtual void minimize_func(Fit_Parameters *fit_params,
                               const Spectra * const spectra,
								Gen_Func_Def gen_func) = 0;


    virtual void minimize_quantification(Fit_Parameters *fit_params,
                                         std::unordered_map<std::string, Element_Quant> * quant_map,
                                         quantification::models::Quantification_Model * quantification_model) = 0;

private:


};

} //namespace optimizers

} //namespace fitting

#endif // Optimizer
