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

#ifndef Fit_Parameters_H
#define Fit_Parameters_H

#include <algorithm>
#include <unordered_map>
#include <string>
#include <vector>
#include <valarray>
#include <limits>
#include <math.h>

#include "defines.h"

namespace data_struct
{
namespace xrf
{

    /**
     * @brief String defines for fit parameters string value pair.
     */
    const std::string STR_FWHM_OFFSET = "FWHM_OFFSET";
    const std::string STR_FWHM_FANOPRIME = "FWHM_FANOPRIME";

    const std::string STR_COHERENT_SCT_ENERGY = "COHERENT_SCT_ENERGY";
    const std::string STR_COHERENT_SCT_AMPLITUDE = "COHERENT_SCT_AMPLITUDE";

    const std::string STR_COMPTON_ANGLE = "COMPTON_ANGLE";
    const std::string STR_COMPTON_FWHM_CORR = "COMPTON_FWHM_CORR";
    const std::string STR_COMPTON_AMPLITUDE = "COMPTON_AMPLITUDE";
    const std::string STR_COMPTON_F_STEP = "COMPTON_F_STEP";
    const std::string STR_COMPTON_F_TAIL = "COMPTON_F_TAIL";
    const std::string STR_COMPTON_GAMMA = "COMPTON_GAMMA";
    const std::string STR_COMPTON_HI_F_TAIL = "COMPTON_HI_F_TAIL";
    const std::string STR_COMPTON_HI_GAMMA = "COMPTON_HI_GAMMA";

    const std::string STR_SNIP_WIDTH = "SNIP_WIDTH";
    const std::string STR_SI_ESCAPE = "SI_ESCAPE";
    const std::string STR_GE_ESCAPE = "GE_ESCAPE";
    const std::string STR_ESCAPE_LINEAR = "ESCAPE_LINEAR";
    const std::string STR_PILEUP0 = "PILEUP0";
    const std::string STR_PILEUP1 = "PILEUP1";
    const std::string STR_PILEUP2 = "PILEUP2";
    const std::string STR_PILEUP3 = "PILEUP3";
    const std::string STR_PILEUP4 = "PILEUP4";
    const std::string STR_PILEUP5 = "PILEUP5";
    const std::string STR_PILEUP6 = "PILEUP6";
    const std::string STR_PILEUP7 = "PILEUP7";
    const std::string STR_PILEUP8 = "PILEUP8";

    const std::string STR_F_STEP_OFFSET = "F_STEP_OFFSET";
    const std::string STR_F_STEP_LINEAR = "F_STEP_LINEAR";
    const std::string STR_F_STEP_QUADRATIC = "F_STEP_QUADRATIC";
    const std::string STR_F_TAIL_OFFSET = "F_TAIL_OFFSET";
    const std::string STR_F_TAIL_LINEAR = "F_TAIL_LINEAR";
    const std::string STR_F_TAIL_QUADRATIC = "F_TAIL_QUADRATIC";
    const std::string STR_GAMMA_OFFSET = "GAMMA_OFFSET";
    const std::string STR_GAMMA_LINEAR = "GAMMA_LINEAR";
    const std::string STR_GAMMA_QUADRATIC = "GAMMA_QUADRATIC";
    const std::string STR_KB_F_TAIL_OFFSET = "KB_F_TAIL_OFFSET";
    const std::string STR_KB_F_TAIL_LINEAR = "KB_F_TAIL_LINEAR";
    const std::string STR_KB_F_TAIL_QUADRATIC = "KB_F_TAIL_QUADRATIC";

const std::string STR_NUM_ITR = "Num_Iter";
const std::string STR_DETECTOR_ELEMENT = "DETECTOR_ELEMENT";

const std::vector<std::string> k_elements {"Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
                                             "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
                                             "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
                                             "In", "Sn", "Sb", "Te", "I"};

const std::vector<std::string> l_elements {"Mo_L", "Tc_L", "Ru_L", "Rh_L", "Pd_L", "Ag_L", "Cd_L", "In_L", "Sn_L", "Sb_L", "Te_L", "I_L", "Xe_L,", "Cs_L", "Ba_L", "La_L", "Ce_L", "Pr_L", "Nd_L", "Pm_L", "Sm_L",
                                             "Eu_L", "Gd_L", "Tb_L", "Dy_L", "Ho_L", "Er_L", "Tm_L", "Yb_L", "Lu_L", "Hf_L", "Ta_L", "W_L", "Re_L", "Os_L", "Ir_L", "Pt_L", "Au_L", "Hg_L", "Tl_L",
                                             "Pb_L", "Bi_L", "Po_L", "At_L", "Rn_L", "Fr_L", "Ac_L", "Th_L", "Pa_L", "U_L", "Np_L", "Pu_L", "Am_L"};

const std::vector<std::string> m_elements {"Au_M", "Pb_M", "U_M", "noise", "Pt_M"};

//-----------------------------------------------------------------------------

enum E_Bound_Type {NOT_INIT=0, FIXED=1, LIMITED_LO_HI=2, LIMITED_LO=3, LIMITED_HI=4, FIT=5};

//-----------------------------------------------------------------------------
/**
 * @brief The Fit_Counts_Array struct: stucture used to seperate fit counts to element lines for debugging or more details. Not yet implemeneted.
 */
struct DLL_EXPORT Fit_Counts_Array
{
    Fit_Counts_Array()
    {

    }

    Fit_Counts_Array(size_t rows, size_t cols)
    {
        resize(rows, cols);
    }

    std::valarray<real_t>& operator [](size_t idx) { return counts[idx]; }

    const std::valarray<real_t>& operator [](size_t idx) const { return counts[idx]; }

    void resize(size_t rows, size_t cols)
    {
        counts.resize(rows);
        for(size_t i=0; i<counts.size(); i++)
        {
            counts[i].resize(cols);
            counts[i] = 0.0;
        }
    }

    const size_t rows() const {return counts.size(); }

    const size_t cols() const {return counts[0].size(); }

    std::valarray< std::valarray<real_t> > counts;
};

typedef std::unordered_map<std::string, Fit_Counts_Array > Fit_Count_Dict;

//-----------------------------------------------------------------------------
/**
 * @brief The Fit_Param struct : Structure that holds a parameter which consists of a value, min, max, and if it should be used in the fit routine.
 *                                Many fit routines use arrays so there are convert to and from array functions.
 */
struct DLL_EXPORT Fit_Param
{
    Fit_Param()
    {
        name = "N/A";
        min_val = std::numeric_limits<real_t>::quiet_NaN();
        max_val = std::numeric_limits<real_t>::quiet_NaN();
        value = std::numeric_limits<real_t>::quiet_NaN();
        step_size = std::numeric_limits<real_t>::quiet_NaN();
        bound_type = E_Bound_Type::NOT_INIT;
        opt_array_index = -1;
    }

    Fit_Param(const Fit_Param& param)
    {
        name = param.name;
        min_val = param.min_val;
        max_val = param.max_val;
        value = param.value;
        step_size = param.step_size;
        bound_type = param.bound_type;
        opt_array_index = param.opt_array_index;
    }

    Fit_Param(std::string name_)
    {
        name = name_;
        min_val = std::numeric_limits<real_t>::quiet_NaN();
        max_val = std::numeric_limits<real_t>::quiet_NaN();
        value = std::numeric_limits<real_t>::quiet_NaN();
        step_size = std::numeric_limits<real_t>::quiet_NaN();
        bound_type = E_Bound_Type::NOT_INIT;
        opt_array_index = -1;
    }

    Fit_Param(std::string name_, real_t val_)
    {
        name = name_;
        min_val = std::numeric_limits<real_t>::min();
        max_val = std::numeric_limits<real_t>::max();
        step_size = 0.000001;
        value = val_;
        bound_type = E_Bound_Type::FIXED;
        opt_array_index = -1;
    }

    Fit_Param(std::string name_, real_t min_, real_t max_, real_t val_, real_t step_size_, E_Bound_Type b_type)
    {
        name = name_;
        min_val = min_;
        max_val = max_;
        value = val_;
        bound_type = b_type;
        step_size = step_size_;
        opt_array_index = -1;
    }

    std::string name;
    real_t min_val;
    real_t max_val;
    real_t value;
    real_t step_size;
    E_Bound_Type bound_type;
    int opt_array_index;
};


//-----------------------------------------------------------------------------
/**
 * @brief The Fit_Parameters class: Dictionary of fit parameters. Many fit routines use arrays so there are convert to and from array functions. TODO: Maybe derive from unordered map?
 */
class DLL_EXPORT Fit_Parameters
{
public:

    Fit_Parameters(){}

    Fit_Parameters(const Fit_Parameters& fit_pars);

    ~Fit_Parameters(){}

    Fit_Param& operator [](std::string name) { return _params[name]; }

    //const Fit_Param& operator [](std::string name) const { return _params[name]; }

    void add_parameter(std::string name, Fit_Param param);

    bool contains(std::string name) const { return ( _params.find(name) != _params.end()); }

    std::vector<real_t> to_array();

    void from_array(std::vector<real_t> arr);

    void from_array(const real_t* arr, size_t arr_size);

    void set_all_value(real_t value, data_struct::xrf::E_Bound_Type btype);

    void set_all(data_struct::xrf::E_Bound_Type btype);

    void update_values(Fit_Parameters override_fit_params);

    void print();

    void print_non_fixed();

    void update_value_by_idx(double * val, int idx);

    const Fit_Param at(std::string name) const {return _params.at(name); }

    void pow10values();

private:

    std::unordered_map<std::string, Fit_Param> _params;

};


} //namespace xrf

} //namespace data_struct

#endif // Fit_Parameters_H
