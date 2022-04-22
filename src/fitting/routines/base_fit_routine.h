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


#ifndef Base_Fit_Routine_H
#define Base_Fit_Routine_H

#include <unordered_map>

#include "fitting/optimizers/optimizer.h"
#include "data_struct/spectra.h"
#include "fitting/models/base_model.h"
#include "data_struct/fit_element_map.h"

namespace fitting
{
namespace routines
{

using namespace data_struct;
using namespace std;


/**
 * @brief The Base_Fit_Routine class: base class for modeling spectra and fitting elements
 */
template<typename T_real>
class DLL_EXPORT Base_Fit_Routine
{
public:
    /**
     * @brief Base_Fit_Routine : Constructor
     */
    Base_Fit_Routine() {}

    /**
     * @brief ~Base_Fit_Routine : Destructor
     */
    virtual ~Base_Fit_Routine() {}
    
    /**
     * @brief fit_spectra : Fit a single specra ( typically 2048 in size )
     * @param fit_params : Fitting parameters required by the routine
     * @param spectra : Pointer to the spectra we are fitting to
     * @param calibration : Energy calibration
     * @param elements_to_fit : List of elemetns to fit to the spectra. This is an out variable also. Must be allocated to saved fitted value to using row_idx and col_idx
     * @param row_idx : row index used to save the fitted value back into elements_to_fit class
     * @param col_idx : column index used to save the fitted value back into elements_to_fit class
     */
    virtual optimizers::OPTIMIZER_OUTCOME fit_spectra(const models::Base_Model<T_real> * const model,
                                                      const Spectra<T_real>* const spectra,
                                                      const Fit_Element_Map_Dict<T_real> * const elements_to_fit,
                                                      std::unordered_map<std::string, T_real>& out_counts) = 0;

    /**
     * @brief get_name : Returns fit routine name
     * @return
     */
    virtual std::string get_name() = 0;

    /**
     * @brief initialize : Initialize the model
     * @param fit_params
     * @param calibration
     * @param elements_to_fit
     * @param energy_range
     */
    virtual void initialize(models::Base_Model<T_real>* const model,
                            const Fit_Element_Map_Dict<T_real> * const elements_to_fit,
                            const struct Range energy_range) = 0;


protected:


private:


};

TEMPLATE_CLASS_DLL_EXPORT Base_Fit_Routine<float>;
TEMPLATE_CLASS_DLL_EXPORT Base_Fit_Routine<double>;

} //namespace routines

} //namespace fitting

#endif // Base_Fit_Routine_H
