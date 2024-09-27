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



#ifndef NLOPT_Optimizer_H
#define NLOPT_Optimizer_H

#include "fitting/optimizers/optimizer.h"
#include <math.h>
#include <nlopt.hpp>

namespace fitting
{
namespace optimizers
{

using namespace data_struct;

template<typename T_real>
class DLL_EXPORT NLOPT_Optimizer: public Optimizer<T_real>
{
public:
    NLOPT_Optimizer();

    virtual ~NLOPT_Optimizer() {}

    virtual OPTIMIZER_OUTCOME minimize(Fit_Parameters<T_real>*fit_params,
                                        const Spectra<T_real>* const spectra,
                                        const Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                                        const Base_Model<T_real>* const model,
                                        const Range energy_range,
                                        bool use_weights,
                                        Callback_Func_Status_Def* status_callback = nullptr);

    virtual OPTIMIZER_OUTCOME minimize_func(Fit_Parameters<T_real>*fit_params,
                                            const Spectra<T_real>* const spectra,
                                            const Range energy_range,
                                            const ArrayTr<T_real>* background,
                                            Gen_Func_Def<T_real> gen_func,
                                            bool use_weights);

    virtual OPTIMIZER_OUTCOME minimize_quantification(Fit_Parameters<T_real>*fit_params,
                                                     std::unordered_map<std::string, Element_Quant<T_real>*> * quant_map,
                                                     quantification::models::Quantification_Model<T_real>* quantification_model);

    virtual std::unordered_map<std::string, T_real> get_options();

    virtual void set_options(std::unordered_map<std::string, T_real> opt);

    virtual std::string detailed_outcome(int outcome);

private:

	//bool _fill_limits(Fit_Parameters<T_real> *fit_params, std::vector<struct mp_par<T_real> > &par);

  //struct mp_config<T_real> _options;
  nlopt_opt _options;
};

} //namespace optimizers

} //namespace fitting

#endif // NLOPT_Optimizer_H
