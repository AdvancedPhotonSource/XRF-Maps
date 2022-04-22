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

#include "core/defines.h"
#include "data_struct/element_info.h"
#include "fitting/routines/base_fit_routine.h"

namespace data_struct
{

//-----------------------------------------------------------------------------

///
/// \brief The Stream_Fitting_Block struct
///
template<typename T_real>
struct Stream_Fitting_Block
{
    fitting::routines::Base_Fit_Routine<T_real>* fit_routine;
    std::unordered_map<std::string, T_real> fit_counts;
};

//-----------------------------------------------------------------------------

///
/// \brief The Stream_Block class
///
template<typename T_real>
class DLL_EXPORT Stream_Block
{

public:

	Stream_Block();

    Stream_Block(int detector, size_t row, size_t col, size_t height, size_t width);

	Stream_Block(const Stream_Block& stream_block);

	Stream_Block& operator=(const Stream_Block&);

    ~Stream_Block();

    void init_fitting_blocks(std::unordered_map<Fitting_Routines, fitting::routines::Base_Fit_Routine<T_real>*> *fit_routines, Fit_Element_Map_Dict<T_real>* elements_to_fit_);

    const size_t& row() { return _row; }

    const size_t& col() { return _col; }

    const size_t& height() { return _height; }

    const size_t& width() { return _width; }

    const int& detector_number() { return _detector; }

    inline bool is_end_of_row() { return (_col == _width-1); }

	inline bool is_end_block() { return (_detector == -1 && _row == -1 && _height == -1 && _col == -1 && _width == -1); }


    //by Fitting_Routines
    std::unordered_map<Fitting_Routines, Stream_Fitting_Block<T_real>> fitting_blocks;

    size_t dataset_hash();

    std::string *dataset_directory;

    std::string *dataset_name;

    Spectra<T_real>* spectra;

    Fit_Element_Map_Dict<T_real> * elements_to_fit;
    //data_struct::Params_Override *fit_params_override_dict;

    fitting::models::Fit_Params_Preset optimize_fit_params_preset;

    fitting::models::Base_Model<T_real> * model;

    float theta;

    bool del_str_ptr;

protected:

    size_t _row;

    size_t _col;

    size_t _height;

    size_t _width;

    int _detector;

};

TEMPLATE_CLASS_DLL_EXPORT Stream_Block<float>;
TEMPLATE_CLASS_DLL_EXPORT Stream_Block<double>;

} //namespace data_struct

#endif // Stream_Block_H
