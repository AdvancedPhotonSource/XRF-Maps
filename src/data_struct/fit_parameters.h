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
#include <limits>
#include <math.h>
#include <vector>
#include <Eigen/Core>

#include "core/defines.h"

namespace data_struct
{
	
//-----------------------------------------------------------------------------

enum class Fitting_Routines { ROI=1 , GAUSS_TAILS=2, GAUSS_MATRIX=4, SVD=8, NNLS=16 };

enum class E_Bound_Type {NOT_INIT=0, FIXED=1, LIMITED_LO_HI=2, LIMITED_LO=3, LIMITED_HI=4, FIT=5};

 const static std::unordered_map<Fitting_Routines, std::string> Fitting_Routine_To_Str = { {Fitting_Routines::ROI, STR_FIT_ROI},
    {Fitting_Routines::GAUSS_TAILS,STR_FIT_GAUSS_TAILS},
    {Fitting_Routines::GAUSS_MATRIX,STR_FIT_GAUSS_MATRIX},
    {Fitting_Routines::SVD,STR_FIT_SVD},
    {Fitting_Routines::NNLS,STR_FIT_NNLS} };

//-----------------------------------------------------------------------------

typedef Eigen::Array<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ArrayXXr;

typedef std::unordered_map<std::string, ArrayXXr > Fit_Count_Dict;

/**
* @brief The Range struct to determine size of spectra we want to fit or model
*/
struct Range
{
	Range() { min = 0; max = 0; }
	Range(const Range& r) { min = r.min; max = r.max; }
    Range(size_t rmin, size_t rmax) { min = rmin; max = rmax; }
	size_t count() const { return (max - min) + 1; }
    size_t min;
    size_t max;
};

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
        step_size = (real_t)0.000001;
        value = val_;
        bound_type = E_Bound_Type::FIXED;
        opt_array_index = -1;
    }

	Fit_Param(std::string name_, real_t val_, E_Bound_Type b_type)
	{
		name = name_;
		min_val = std::numeric_limits<real_t>::min();
		max_val = std::numeric_limits<real_t>::max();
		step_size = (real_t)0.000001;
		value = val_;
		bound_type = b_type;
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

    const std::string bound_type_str() const; 

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
 * @brief The Fit_Parameters class: Dictionary of fit parameters. Many fit routines use arrays so there are convert to and from array functions.
 */
class DLL_EXPORT Fit_Parameters
{
public:

    Fit_Parameters(){}

    Fit_Parameters(const Fit_Parameters& fit_pars);

    ~Fit_Parameters(){_params.clear();}

    Fit_Param& operator [](std::string name) { return _params[name]; }

    //const Fit_Param& operator [](std::string name) const { return _params[name]; }

    void add_parameter(Fit_Param param);

	void append_and_update(Fit_Parameters* fit_params);

    inline auto begin() { return _params.begin(); }

    inline auto end() { return _params.end(); }

    void sum_values(Fit_Parameters fit_params);

    void divide_fit_values_by(real_t divisor);

    bool contains(std::string name) const { return ( _params.find(name) != _params.end()); }

    std::vector<real_t> to_array();

    std::vector<std::string> names_to_array();

    void from_array(std::vector<real_t> &arr);

    void from_array(const real_t* arr, size_t arr_size);

    void set_all_value(real_t value, E_Bound_Type btype);

    void set_all(E_Bound_Type btype);

    void update_values(Fit_Parameters *override_fit_params);

    void update_and_add_values(Fit_Parameters *override_fit_params);

    void update_and_add_values_gt_zero(Fit_Parameters *override_fit_params);

    void remove(Fit_Parameters* override_fit_params);

    void remove(std::string key);

    inline const real_t& value(std::string key) const { return _params.at(key).value; }

    void print();

    void print_non_fixed();

    const Fit_Param& at(std::string name) const {return _params.at(name); }

    size_t size() const { return _params.size(); }
	
	std::unordered_map<std::string, Fit_Param>* Params() { return &_params; }

private:

    std::unordered_map<std::string, Fit_Param> _params;

};

//-----------------------------------------------------------------------------

DLL_EXPORT Range get_energy_range(size_t spectra_size, Fit_Parameters* params);

/**
* @brief get_energy_range: genereates a range which consists of min and max. This represents the min energy and max enegry of the spectra to fit.
* @param min_energy
* @param max_energy
* @param spectra_size
* @param calibration: energy calibration
* @return Range structure with the min energy and max enegry of the spectra to fit.
*/
DLL_EXPORT Range get_energy_range(real_t min_energy, real_t max_energy, size_t spectra_size, real_t energy_offset, real_t energy_slope);

} //namespace data_struct

#endif // Fit_Parameters_H
