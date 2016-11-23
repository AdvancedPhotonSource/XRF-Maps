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
#include <valarray>
#include <vector>
#include <string>
#include <unordered_map>
#include "spectra.h"
#include "element_info.h"

#include "optimizer.h"

namespace data_struct
{
namespace xrf
{


enum Electron_Shell {K_SHELL, L_SHELL, M_SHELL, N_SHELL, O_SHELL, P_SHELL, Q_SHELL};

struct Element_Quant
{
    Element_Quant()
    {
        zero();
    }
    Element_Quant(real_t weight_)
    {
        zero();
        weight = weight_;
    }
    void zero()
    {
        weight = 0.0;
        absorption = 0.0;
        transmission_Be = 0.0;
        transmission_Ge = 0.0; // or Si dead layer
        yield = 0.0;
        transmission_through_Si_detector = 0.0;
        transmission_through_air = 0.0;// (N2)
        index = -1;
    }

    real_t weight;  // in ug/cm2
    real_t absorption;
    real_t transmission_Be;
    real_t transmission_Ge; // or Si dead layer
    real_t yield;
    real_t transmission_through_Si_detector;
    real_t transmission_through_air;// (N2)
    int index; //used to map to and from fitting function
};

//-----------------------------------------------------------------------------

///
/// \brief The Quantification_Standard class: Class of resulting fit per element
///
class DLL_EXPORT Quantification_Standard
{

public:
    Quantification_Standard();

    ~Quantification_Standard();

    void append_element(std::string name, real_t weight);

    const real_t& element_weight(std::string element_symb) const { return _element_quants.at(element_symb).weight; }

    void standard_filename(std::string standard_filename) { _standard_filename = standard_filename; }

    const std::string& standard_filename() { return _standard_filename; }

    bool quantifiy(fitting::optimizers::Optimizer * optimizer,
                   real_t incident_energy,
                   Element_Info* detector_element,
                   bool airpath,
                   real_t detector_chip_thickness,
                   real_t beryllium_window_thickness,
                   real_t germanium_dead_layer);

    Element_Quant generate_element_quant(real_t incident_energy,
                                        Element_Info* detector_element,
                                        Electron_Shell shell,
                                        bool airpath,
                                        real_t detector_chip_thickness,
                                        real_t beryllium_window_thickness,
                                        real_t germanium_dead_layer,
                                        size_t z_number);

    std::unordered_map<std::string, Element_Quant> generate_quant_map(real_t incident_energy,
                                                                      Element_Info* detector_element,
                                                                      Electron_Shell shell,
                                                                      bool airpath = false,
                                                                      real_t detector_chip_thickness = 0.0,
                                                                      real_t beryllium_window_thickness = 0.0,
                                                                      real_t germanium_dead_layer = 0.0,
                                                                      size_t start_z = 0,
                                                                      size_t end_z = 95);

    real_t transmission(real_t thickness, real_t beta, real_t llambda) const;

    real_t absorption(real_t thickness, real_t beta, real_t llambda, real_t shell_factor=1) const;

    std::unordered_map<std::string, real_t> model_calibrationcurve(std::unordered_map<std::string, Element_Quant> quant_map, real_t p);

protected:

    std::string _standard_filename;

    real_t _sr_current;
    real_t _IC_US;
    real_t _IC_DS;

    //std::unordered_map<std::string, real_t> _element_weights; // in ug/cm2

    std::unordered_map<std::string, Element_Quant> _element_quants;

    std::unordered_map<std::string, Element_Quant> _calibration_curve;

    //std::valarray<real_t> _us_amp;
    //std::valarray<real_t> _ds_amp;

};

//-----------------------------------------------------------------------------


} //namespace xrf

} //namespace data_struct

#endif // Quantification_Standard_H
