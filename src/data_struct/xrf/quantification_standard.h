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



#ifndef Quantification_Standard_H
#define Quantification_Standard_H

#include "defines.h"

#include <string>
#include <unordered_map>

#include "element_info.h"
#include "element_quant.h"
#include "optimizer.h"


namespace data_struct
{
namespace xrf
{

using namespace std;

//-----------------------------------------------------------------------------

struct Calibration_Curve
{
    //enum C_KEYS {K_SHELL, L_SHELL, M_SHELL};

    Calibration_Curve()
    {
        shell_curves.resize(3);
        shell_curves_labels.resize(3);
    }

    Calibration_Curve(size_t i)
    {
        shell_curves.resize(3);
        shell_curves_labels.resize(3);
        resize(i);
    }

    void resize(size_t i)
    {
        for(auto& curve : shell_curves)
        {
            curve.resize(i);
        }
    }

    string quantifier_name;
    vector<vector<real_t> > shell_curves;
    vector<vector<std::string> > shell_curves_labels;

};

//-----------------------------------------------------------------------------

struct Quantifiers
{
    enum Q_KEYS {CURRENT, US_IC, DS_IC};

    Quantifiers()
    {
        calib_curves.resize(3);
        calib_curves[CURRENT].quantifier_name = "Current";
        calib_curves[US_IC].quantifier_name = "US_IC";
        calib_curves[DS_IC].quantifier_name = "DS_IC";
    }

    Quantifiers(size_t i)
    {
        calib_curves.resize(3);
        calib_curves[CURRENT].quantifier_name = "Current";
        calib_curves[US_IC].quantifier_name = "US_IC";
        calib_curves[DS_IC].quantifier_name = "DS_IC";
        resize(i);
    }

    void resize(size_t i)
    {
        for(auto& curve : calib_curves)
        {
            curve.resize(i);
        }
    }

    vector<Calibration_Curve> calib_curves;
    //Calibration_Curve current;
    //Calibration_Curve us_ic;
    //Calibration_Curve ds_ic;
};

//-----------------------------------------------------------------------------

///
/// \brief The Quantification_Standard class:
///
class DLL_EXPORT Quantification_Standard
{

public:
    Quantification_Standard();

    ~Quantification_Standard();

    void append_element(string name, real_t weight);

    const real_t& element_weight(string element_symb) const { return _element_quants.at(element_symb).weight; }

    const unordered_map<string, Element_Quant>& element_weights() const { return _element_quants; }

    void standard_filename(string standard_filename) { _standard_filename = standard_filename; }

    const string& standard_filename() const { return _standard_filename; }

    void element_counts(string, unordered_map<string, real_t> map);

    const unordered_map<string, unordered_map<string, real_t> >& element_counts() const { return _element_counts; }

    void integrated_spectra(const Spectra &spec) {_integrated_spectra = spec; }

    const Spectra& integrated_spectra() const { return _integrated_spectra; }

    const real_t& sr_current() {return _sr_current;}

    void sr_current(real_t val) {_sr_current = val;}

    const real_t& US_IC() {return _US_IC;}

    void US_IC(real_t val) {_US_IC = val;}

    const real_t& DS_IC() {return _DS_IC;}

    void DS_IC(real_t val) {_DS_IC = val;}

    bool processed() {return _processed;}

    bool quantifiy(fitting::optimizers::Optimizer * optimizer,
                   string proc_type_str,
                   unordered_map<string, real_t>  *element_counts,
                   real_t incident_energy,
                   Element_Info* detector_element,
                   bool airpath,
                   real_t detector_chip_thickness,
                   real_t beryllium_window_thickness,
                   real_t germanium_dead_layer);

    //           proc_type  quantifier
    unordered_map<string, Quantifiers> calibration_curves;

protected:

    bool _processed;

    //          element    quant
    unordered_map<string, Element_Quant> _element_quants;

    //          proc_type               Element   Counts
    unordered_map<string, unordered_map<string, real_t> >  _element_counts;

    Spectra _integrated_spectra;

    std::string _standard_filename;

    real_t _sr_current;
    real_t _US_IC;
    real_t _DS_IC;



};

//-----------------------------------------------------------------------------


} //namespace xrf

} //namespace data_struct

#endif // Quantification_Standard_H
