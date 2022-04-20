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



#include "spectra_line.h"

namespace data_struct
{

// ----------------------------------------------------------------------------

template<typename T_real>
Spectra_Line<T_real>::Spectra_Line()
{

}

// ----------------------------------------------------------------------------

template<typename T_real>
Spectra_Line<T_real>::~Spectra_Line()
{
    /*
    for(auto spec : _data_line)
    {
        delete spec;
    }
    _data_line.clear();
    */
}

// ----------------------------------------------------------------------------

/*
void Spectra_Line<T_real>::append_spectra(Spectra* spectra)
{
   _data_line.push_back(spectra);
}
*/

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Line<T_real>::resize_and_zero(size_t cols, size_t samples)
{
    alloc_row_size(cols);
    _alloc_spectra_size(samples);
}

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Line<T_real>::alloc_row_size(size_t n)
{
    _data_line.resize(n);
}

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Line<T_real>::_alloc_spectra_size(size_t n)
{
    for(size_t i=0; i<_data_line.size(); i++)
    {
        _data_line[i].resize(n);
        _data_line[i].setZero(n);
    }
    /*
    for(Spectra s : _data_line)
    {
        s.alloc_buffer(n);
    }*/
    //_data_line.resize(n);
}

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Line<T_real>::recalc_elapsed_livetime()
{
    for(size_t i=0; i<_data_line.size(); i++)
    {
        _data_line[i].recalc_elapsed_livetime();
    }
}

} //namespace data_struct
