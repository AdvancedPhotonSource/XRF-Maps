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

// ----------------------------------------------------------------------------

template<typename T_real>
Spectra_Volume<T_real>::Spectra_Volume()
{

}

// ----------------------------------------------------------------------------

template<typename T_real>
Spectra_Volume<T_real>::~Spectra_Volume()
{

}

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Volume<T_real>::resize_and_zero(size_t rows, size_t cols, size_t samples)
{

    _data_vol.resize(rows);
    for(size_t i=0; i<_data_vol.size(); i++)
    {
        _data_vol[i].resize_and_zero(cols, samples);
    }

}

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Volume<T_real>::resize_samples(size_t samples)
{
    for(size_t i=0; i<_data_vol.size(); i++)
    {
        _data_vol[i].resize_samples(samples);
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
Spectra<T_real> Spectra_Volume<T_real>::integrate()
{

    if(_data_vol.size() > 0)
    {
        Spectra<T_real> i_spectra(_data_vol[0][0].size());
        
        T_real elt = 0.0;
        T_real ert = 0.0;
        T_real in_cnt = 0.0;
        T_real out_cnt = 0.0;
        for(size_t i = 0; i < _data_vol.size(); i++)
        {
            for(size_t j = 0; j < _data_vol[0].size(); j++)
            {
                i_spectra += _data_vol[i][j];
                elt += _data_vol[i][j].elapsed_livetime();
                ert += _data_vol[i][j].elapsed_realtime();
                in_cnt += _data_vol[i][j].input_counts();
                out_cnt += _data_vol[i][j].output_counts();
            }
        }

        i_spectra.elapsed_livetime(elt);
        i_spectra.elapsed_realtime(ert);
        i_spectra.input_counts(in_cnt);
        i_spectra.output_counts(out_cnt);

        i_spectra.recalc_elapsed_livetime();

        return i_spectra;
    }
    else
    {
        Spectra<T_real> i_spectra(1);
        return i_spectra;
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
bool Spectra_Volume<T_real>::integrate_polar(Spectra<T_real> &lhs, Spectra<T_real> &rhs)
{

    if(_data_vol.size() > 0)
    {
        lhs.resize(_data_vol[0][0].size());
        rhs.resize(_data_vol[0][0].size());
        
        T_real lhs_elt = 0.0;
        T_real lhs_ert = 0.0;
        T_real lhs_in_cnt = 0.0;
        T_real lhs_out_cnt = 0.0;

        T_real rhs_elt = 0.0;
        T_real rhs_ert = 0.0;
        T_real rhs_in_cnt = 0.0;
        T_real rhs_out_cnt = 0.0;
        bool lhs_now = true;
        for(size_t i = 0; i < _data_vol.size(); i++)
        {
            for(size_t j = 0; j < _data_vol[0].size(); j++)
            {
                if(lhs_now)
                {
                    lhs += _data_vol[i][j];
                    lhs_elt += _data_vol[i][j].elapsed_livetime();
                    lhs_ert += _data_vol[i][j].elapsed_realtime();
                    lhs_in_cnt += _data_vol[i][j].input_counts();
                    lhs_out_cnt += _data_vol[i][j].output_counts();
                }
                else
                {
                    rhs += _data_vol[i][j];
                    rhs_elt += _data_vol[i][j].elapsed_livetime();
                    rhs_ert += _data_vol[i][j].elapsed_realtime();
                    rhs_in_cnt += _data_vol[i][j].input_counts();
                    rhs_out_cnt += _data_vol[i][j].output_counts();
                }
            }
            lhs_now = !lhs_now;
        }

        lhs.elapsed_livetime(lhs_elt);
        lhs.elapsed_realtime(lhs_ert);
        lhs.input_counts(lhs_in_cnt);
        lhs.output_counts(lhs_out_cnt);

        rhs.elapsed_livetime(rhs_elt);
        rhs.elapsed_realtime(rhs_ert);
        rhs.input_counts(rhs_in_cnt);
        rhs.output_counts(rhs_out_cnt);

        lhs.recalc_elapsed_livetime();
        rhs.recalc_elapsed_livetime();
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Volume<T_real>::recalc_elapsed_livetime()
{

    for(size_t i=0; i<_data_vol.size(); i++)
    {
        _data_vol[i].recalc_elapsed_livetime();
    }

}

// ----------------------------------------------------------------------------

template<typename T_real>
void Spectra_Volume<T_real>::generate_scaler_maps(std::vector<Scaler_Map<T_real>> *scaler_maps)
{
    if (scaler_maps != nullptr)
    {
        data_struct::Scaler_Map<T_real> elt_map, ert_map, in_cnt_map, out_cnt_map, dead_time_map;

        elt_map.name = STR_ELT;
        ert_map.name = STR_ERT;
        in_cnt_map.name = "INCNT";
        out_cnt_map.name = "OUTCNT";
        dead_time_map.name = STR_DEAD_TIME;
        elt_map.unit = "seconds";
        ert_map.unit = "seconds";
        in_cnt_map.unit = "cts/s";
        out_cnt_map.unit = "cts/s";
        dead_time_map.unit = "%";

        elt_map.values.resize(_data_vol.size(), _data_vol[0].size());
        ert_map.values.resize(_data_vol.size(), _data_vol[0].size());
        in_cnt_map.values.resize(_data_vol.size(), _data_vol[0].size());
        out_cnt_map.values.resize(_data_vol.size(), _data_vol[0].size());
        dead_time_map.values.resize(_data_vol.size(), _data_vol[0].size());

        for (size_t i = 0; i < _data_vol.size(); i++)
        {
            for (size_t j = 0; j < _data_vol[0].size(); j++)
            {
                elt_map.values(i, j) = _data_vol[i][j].elapsed_livetime();
                ert_map.values(i, j) = _data_vol[i][j].elapsed_realtime();
                in_cnt_map.values(i, j) = _data_vol[i][j].input_counts();
                out_cnt_map.values(i, j) = _data_vol[i][j].output_counts();
                dead_time_map.values(i, j) = (1.0 - (out_cnt_map.values(i, j) / in_cnt_map.values(i, j))) * 100.0;
            }
        }
        
        scaler_maps->push_back(elt_map);
        scaler_maps->push_back(ert_map);
        scaler_maps->push_back(in_cnt_map);
        scaler_maps->push_back(out_cnt_map);
        scaler_maps->push_back(dead_time_map);
    }
}

// ----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT Spectra_Volume<float>;
TEMPLATE_CLASS_DLL_EXPORT Spectra_Volume<double>;

} //namespace data_struct
