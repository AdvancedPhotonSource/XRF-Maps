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

/// Initial Author <2017>: Arthur Glowacki



#ifndef Stream_Block_H
#define Stream_Block_H

#include "defines.h"
#include "element_info.h"
#include "base_fit_routine.h"
#include <vector>
#include <string>

namespace data_struct
{
namespace xrf
{

///
/// \brief The Stream_Fitting_Block struct
///
struct Stream_Fitting_Block
{
    fitting::routines::Base_Fit_Routine * fit_routine;
    std::unordered_map<std::string, real_t> fit_counts;
};

//-----------------------------------------------------------------------------

///
/// \brief The Stream_Block class
///
class DLL_EXPORT Stream_Block
{

public:

    Stream_Block(size_t row, size_t col);

    ~Stream_Block();

    void init_fitting_blocks(std::unordered_map<int, fitting::routines::Base_Fit_Routine *> *fit_routines, Fit_Element_Map_Dict * elements_to_fit_);

    const size_t& row() { return _row; }

    const size_t& col() { return _col; }

    //by Processing_Type
    std::unordered_map<int, Stream_Fitting_Block> fitting_blocks;

    size_t detector_number;

    Spectra * spectra;

    Fit_Element_Map_Dict * elements_to_fit;

    fitting::models::Base_Model * model;

protected:

    size_t _row;

    size_t _col;

};


} //namespace xrf

} //namespace data_struct

#endif // Stream_Block_H
