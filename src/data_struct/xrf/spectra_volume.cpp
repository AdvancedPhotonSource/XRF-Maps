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



#include "spectra_volume.h"

namespace data_struct
{
namespace xrf
{

Spectra_Volume::Spectra_Volume() //: Base_Dataset()
{
    //_data = nullptr;
}

Spectra_Volume::~Spectra_Volume()
{
    /*
    for(auto spec_line : _data_vol)
    {
        delete spec_line;
    }
    _data_vol.clear();

    if(_data != nullptr)
        delete [] _data;
    _data = nullptr;
*/
}
/*
void Spectra_Volume::append_spectra_line(Spectra_Line* spec_line)
{
  // _data_vol.push_back(spec_line);

}

void Spectra_Volume::alloc(unsigned short rank, int *dims)
{
    long total = 0;
    for (unsigned short i=0; i<rank; i++)
    {
        total *= (long)dims[i];
    }

    if (_data != nullptr)
        delete [] _data;
    _data = new real_t[total];
}
*/
void Spectra_Volume::resize(size_t rows, size_t cols, size_t samples)
{

    _data_vol.resize(rows);
    for(size_t i=0; i<_data_vol.size(); i++)
    {
        _data_vol[i].alloc(cols, samples);
    }
/*
    for(Spectra_Line line : _data_vol)
    {
        line.alloc(rows, samples);
    }
    */
/*
    _rows = rows;
    _cols = cols;
    _samples = samples;

    long total = _rows * _cols * _samples;
    if (_data != nullptr)
        delete [] _data;
    _data = new real_t[total];
*/
/*
    array3D.resize(_samples);
      for (int i = 0; i < _samples; ++i) {
        array3D[i].resize(_rows);

        for (int j = 0; j < _rows; ++j)
          array3D[i][j].resize(_cols);
      }

      array3D[115][0][0] = 10;
      */
/*
    _data = new real_t**[_cols];
    for (unsigned int i=0; i<_cols; i++)
    {
        _data[i] = new real_t*[_rows];
        for (unsigned int j=0; j<_rows; j++)
        {
            _data[i][j] = new real_t[_samples];
        }
    }
*/
/*
    _data = new real_t**[_samples];
    for (unsigned int s=0; s<_samples; s++)
    {
        _data[s] = new real_t*[_rows];
        for (unsigned int r=0; r<_rows; r++)
        {
            _data[s][r] = new real_t[_cols];
        }
    }
    _data[115][0][0] = 10;
*/
}
/*
void Spectra_Volume::del_data()
{

    if (_data != nullptr)
    {
        for (unsigned int i=0; i<_cols; i++)
        {
            for (unsigned int j=0; j<_rows; j++)
            {
                delete [] _data[i][j];
            }
            delete [] _data[i];
        }
        delete [] _data;
    }
    _data = nullptr;

}
*/

const Spectra Spectra_Volume::integrate()
{

    Spectra i_spectra(_data_vol[0][0].size());
    real_t elt = 0.0;
    real_t ert = 0.0;
    real_t in_cnt = 0.0;
    real_t out_cnt = 0.0;
    for(size_t i = 0; i < _data_vol.size(); i++)
    {
        for(size_t j = 0; j < _data_vol[0].size(); j++)
        {
            for(size_t k = 0; k < _data_vol[0][0].size(); k++)
            {
                i_spectra[k] += _data_vol[i][j][k];
                elt += _data_vol[i][j].elapsed_lifetime();
                ert += _data_vol[i][j].elapsed_realtime();
                in_cnt += _data_vol[i][j].input_counts();
                out_cnt += _data_vol[i][j].output_counts();
            }
        }
    }

    i_spectra.elapsed_lifetime(elt);
    i_spectra.elapsed_realtime(ert);
    i_spectra.input_counts(in_cnt);
    i_spectra.output_counts(out_cnt);

    i_spectra.recalc_elapsed_lifetime();

    return i_spectra;
}
/*
real_t* Spectra_Volume::get_spectra(unsigned int row, unsigned int col)
{

    return 0;//&_data[ (_samples * col * _rows)+(_samples * row) ];

}
*/
void Spectra_Volume::recalc_elapsed_lifetime()
{

    for(size_t i=0; i<_data_vol.size(); i++)
    {
        _data_vol[i].recalc_elapsed_lifetime();
    }

}
/*
void Spectra_Volume::transpose(real_t *buf)
{

    int col_row = _rows * _cols;
    for(unsigned int c = 0; c< _cols; c++)
    {
        int cTrow = c*_rows;
        for(unsigned int r = 0; r< _rows; r++)
        {
            int rTcol = r*_cols;
            for(unsigned int s = 0; s< _samples; s++)
            {
                //real_t tmp = _data[ s*col_row + rTcol + c ];
                //_data[ s*col_row + rTcol + c ] = _data[s + (_samples * r)+(_samples * cTrow)];
                //_data[s + (_samples * r)+(_samples * cTrow)] = tmp;
                _data[s + (_samples * r)+(_samples * cTrow)] = buf[ s*col_row + rTcol + c ];
            }
        }
    }
/////////////
    int col_row = _rows * _cols;

    for(unsigned int s = 0; s< _samples; s++)
    {
        for(unsigned int r = 0; r< _rows; r++)
        {
            for(unsigned int c = 0; c< _cols; c++)
            {

                real_t tmp = _data[ (s*col_row) + (r*_cols) + c ];
                _data[ (s*col_row) + (r*_cols) + c ] = _data[s + (_samples * r)+(_samples * c*_rows)];
                _data[s + (_samples * r)+(_samples * c*_rows)] = tmp;
            }
        }
    }

}
*/
} //namespace xrf
} //namespace data_struct
