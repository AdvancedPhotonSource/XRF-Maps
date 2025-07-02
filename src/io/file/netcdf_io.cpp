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


#include "netcdf_io.h"

#include <iostream>
#include <string>

#include <chrono>
#include <ctime>
#include <thread>

namespace io
{
namespace file
{

template<typename T_real>
std::mutex NetCDF_IO<T_real>::_mutex;

template<typename T_real>
NetCDF_IO<T_real>* NetCDF_IO<T_real>::_this_inst(nullptr);


#define ELAPSED_REALTIME_OFFSET 32
#define ELAPSED_LIVETIME_OFFSET 34
#define INPUT_COUNTS_OFFSET 36
#define OUTPUT_COUNTS_OFFSET 38

#define MAX_NUM_SUPPORTED_DETECOTRS_PER_COL 4

//-----------------------------------------------------------------------------

template<typename T_real>
NetCDF_IO<T_real>::NetCDF_IO()
{

}

//-----------------------------------------------------------------------------

template<typename T_real>
NetCDF_IO<T_real>* NetCDF_IO<T_real>::inst()
{
    std::lock_guard<std::mutex> lock(_mutex);

    if (_this_inst == nullptr)
    {
        _this_inst = new NetCDF_IO();
    }
    return _this_inst;
}

//-----------------------------------------------------------------------------

template<typename T_real>
size_t NetCDF_IO<T_real>::_load_spectra(E_load_type ltype,
                                std::string path,
                                size_t detector,
                                data_struct::Spectra_Line<T_real>* spec_line,
                                size_t line_size,
                                data_struct::Spectra<T_real>* spectra,
                                size_t cur_row,
                                size_t max_rows,
                                data_struct::IO_Callback_Func_Def<T_real> *callback_fun,
                                void* user_data)
{
    std::lock_guard<std::mutex> lock(_mutex);

    size_t header_size = 256;
    int ncid = 0, varid = 0, retval = 0;
    size_t start[] = {0, 0, 0};
    size_t count[] = {1, 1, 1047808};
    ptrdiff_t stride[] = {1, 1, 1};
    T_real data_in[1][1][1047808];
    size_t spectra_size = 0;
    nc_type rh_type;
    int rh_ndims = 0;
    int  rh_dimids[NC_MAX_VAR_DIMS] = {0};
    int rh_natts = 0;
    T_real elapsed_livetime = 0.;
    T_real elapsed_realtime = 0.;
    T_real input_counts = 0.;
    T_real output_counts = 0.;
    size_t cols_before_inc = 124; // default is 124 but read in also
    data_struct::Spectra<T_real>* callback_spectra = nullptr;

    size_t dim2size[NC_MAX_VAR_DIMS] = {0};

    if( (retval = nc_open(path.c_str(), NC_NOWRITE, &ncid)) != 0)
    {
        logE<<path<<" :: "<< nc_strerror(retval)<<"\n";
        return 0;
    }


    if( (retval = nc_inq_varid(ncid, "array_data", &varid)) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return 0;
    }

    if( (retval = nc_inq_var (ncid, varid, nullptr, &rh_type, &rh_ndims, rh_dimids, &rh_natts) ) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return 0;
    }

    for (int i=0; i <  rh_ndims; i++)
    {
        if( (retval = nc_inq_dimlen(ncid, rh_dimids[i], &dim2size[i]) ) != 0)
        {
            logE<< path << " :: " << nc_strerror(retval)<<"\n";
            nc_close(ncid);
            return 0;
        }
    }

    if (detector > 3)
    {
        if (dim2size[1] != 2)
        {
            logE << "NetCDF dims: [" << dim2size[0] <<"]["<< dim2size[1] <<"]["<< dim2size[2] <<"] needs to be [x][2][x] for detector "<<detector<<" " << path << "\n";
            nc_close(ncid);
            return -1;
        }
        start[1] = 1;
    }

    if( (retval = _nc_get_vars_real(ncid, varid, start, count, stride, &data_in[0][0][0]) ) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return 0;
    }

    if (data_in[0][0][0] != 21930 || data_in[0][0][1] != -21931)
    {
        logE<<"NetCDF header [0][0][0]  not found! Stopping load : "<<path<<"\n";
        nc_close(ncid);
        return 0;
    }

    
    int d_idx = 12;
    if (detector > 3)
    {
        d_idx += 2 * (detector-4);
    }
    else
    {
        d_idx += 2 * detector;
    }
    
    size_t dset_det = size_t(data_in[0][0][d_idx]);
    if (dset_det != detector)
    {
        logE << "detector not found! "<< dset_det <<" != "<<detector<<" Stopping load : " << path << "\n";
        nc_close(ncid);
        return -1;
    }
    

    header_size = data_in[0][0][2];
    cols_before_inc = data_in[0][0][8];  //sum all across the first dim looking at value 8
    spectra_size = data_in[0][0][20];


    if (detector > 3)
    {
        detector -= MAX_NUM_SUPPORTED_DETECOTRS_PER_COL; // 4,5,6,7 = 0,1,2,3
    }

    size_t spec_cntr = 0;
    if (ltype == E_load_type::LINE)
    {
        spec_cntr = spec_line->size();
    }
    else if (ltype == E_load_type::INTEGRATED || ltype == E_load_type::CALLBACKF)
    {
        spec_cntr = line_size;
    }

    size_t l = header_size;
    //size_t max_spec_per_netcdf_idx = (1047808 - header_size) / (header_size + (spectra_size * MAX_NUM_SUPPORTED_DETECOTRS_PER_COL))

    for(size_t j=0; j<spec_cntr; j++)
    {
		if (ltype == E_load_type::LINE)
		{
			(*spec_line)[j].resize(spectra_size); // should be renames to resize
		}
        else if (ltype == E_load_type::CALLBACKF)
        {
            callback_spectra = new data_struct::Spectra<T_real>(spectra_size);
        }

        //if ( j>0 &&  (j % max_spec_per_netcdf_idx) == 0)
        if (j > cols_before_inc || l > count[2])
        {
            l = header_size;
            start[0]++;
            //read header
            if( (retval = _nc_get_vars_real(ncid, varid, start, count, stride, &data_in[0][0][0]) ) != 0)
            {
                logE<< nc_strerror(retval)<<"\n";
                nc_close(ncid);
                return j;
            }
            cols_before_inc = data_in[0][0][8];
        }
        if (data_in[0][0][l] != 13260 || data_in[0][0][l+1] != -13261)
        {
            if(j < spec_cntr -2)
            {
                logE<<"NetCDF sub header not found! Stopping load at Col: "<<j<<" path :"<<path<<"\n";
                nc_close(ncid);
                return j;
            }
            //last two may not be filled with data
            //TODO: send end of row stream_block down pipeline
            nc_close(ncid);
            return j;
        }

        
        unsigned short i1 = data_in[0][0][l+ELAPSED_LIVETIME_OFFSET+(detector*8)];
        unsigned short i2 = data_in[0][0][l+ELAPSED_LIVETIME_OFFSET+(detector*8)+1];
        unsigned int ii = i1 | i2<<16;
        if (ltype == E_load_type::LINE)
        {
            elapsed_livetime = ((float)ii) * 320e-9f; // need to multiply by this value becuase of the way it is saved
            if (elapsed_livetime == 0)
            {
                if (j > 0 && j < spec_cntr - 2) // copy the previous value
                {
                    logW << "Reading in elapsed lifetime for Col:" << j << " is 0. Setting it to " << j - 1 << ". path :" << path << "\n";
                    elapsed_livetime = (*spec_line)[j - 1].elapsed_livetime();
                }
                else if (j < spec_cntr - 2) // usually the last two are missing which spams the log ouput.
                {
                    logW << "Reading in elapsed lifetime for Col:" << j << " is 0. Setting it to 1.0. path :" << path << "\n";
                    elapsed_livetime = 1.0;
                }
            }
            else
            {
                (*spec_line)[j].elapsed_livetime(elapsed_livetime);
            }
        }
        else if (ltype == E_load_type::INTEGRATED)
        {
            elapsed_livetime += ((float)ii) * 320e-9f; 
        }
        else if(ltype == E_load_type::CALLBACKF && callback_spectra != nullptr)
        {
            callback_spectra->elapsed_livetime( ((float)ii) * 320e-9f);
        }

        i1 = data_in[0][0][l+ELAPSED_REALTIME_OFFSET+(detector*8)];
        i2 = data_in[0][0][l+ELAPSED_REALTIME_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        if (ltype == E_load_type::LINE)
        {
            elapsed_realtime = ((float)ii) * 320e-9f; // need to multiply by this value becuase of the way it is saved
            if (elapsed_realtime == 0)
            {
                if (j > 0 && j < spec_cntr - 2) // copy the previous value
                {
                    logW << "Reading in elapsed realtime for Col:" << j << " is 0. Setting it to " << j - 1 << ". path :" << path << "\n";
                    elapsed_realtime = (*spec_line)[j - 1].elapsed_realtime();
                }
                else if (j < spec_cntr - 2) // usually the last two are missing which spams the log ouput.
                {
                    logW << "Reading in elapsed realtime for Col:" << j << " is 0. Setting it to 1.0. path :" << path << "\n";
                    elapsed_realtime = 1.0;
                }
            }
            else
            {
                (*spec_line)[j].elapsed_realtime(elapsed_realtime);
            }
        }
        else if (ltype == E_load_type::INTEGRATED)
        {
            elapsed_realtime += ((float)ii) * 320e-9f;
        }
        else if(ltype == E_load_type::CALLBACKF && callback_spectra != nullptr)
        {
            callback_spectra->elapsed_realtime(((float)ii) * 320e-9f);
        }

        i1 = data_in[0][0][l+INPUT_COUNTS_OFFSET+(detector*8)];
        i2 = data_in[0][0][l+INPUT_COUNTS_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        if (ltype == E_load_type::LINE || ltype == E_load_type::CALLBACKF)
        {
            input_counts = ((float)ii) / elapsed_livetime;
            if (input_counts == 0)
            {
                if (j > 0 && j < spec_cntr - 2) // copy the previous value
                {
                    logW << "Reading in elapsed input_counts for Col:" << j << " is 0. Setting it to " << j - 1 << ". path :" << path << "\n";
                    input_counts = (*spec_line)[j - 1].input_counts();
                }
                else if (j < spec_cntr - 2) // usually the last two are missing which spams the log ouput.
                {
                    logW << "Reading in elapsed input_counts for Col:" << j << " is 0. Setting it to 1.0. path :" << path << "\n";
                    input_counts = 1.0;
                }
            }
            else
            {
                (*spec_line)[j].input_counts(input_counts);
            }
        }
        else if (ltype == E_load_type::INTEGRATED)
        {
            input_counts += ((float)ii) / elapsed_livetime;
        }
        else if(ltype == E_load_type::CALLBACKF && callback_fun != nullptr)
        {
            callback_spectra->input_counts(((float)ii) / elapsed_livetime);
        }


        i1 = data_in[0][0][l+OUTPUT_COUNTS_OFFSET+(detector*8)];
        i2 = data_in[0][0][l+OUTPUT_COUNTS_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        if (ltype == E_load_type::LINE || ltype == E_load_type::CALLBACKF)
        {
            output_counts = ((float)ii) / elapsed_realtime;
            if (output_counts == 0)
            {
                if (j > 0 && j < spec_cntr - 2) // copy the previous value
                {
                    logW << "Reading in elapsed output_counts for Col:" << j << " is 0. Setting it to " << j - 1 << ". path :" << path << "\n";
                    output_counts = (*spec_line)[j - 1].output_counts();
                }
                else if (j < spec_cntr - 2) // usually the last two are missing which spams the log ouput.
                {
                    logW << "Reading in elapsed output_counts for Col:" << j << " is 0. Setting it to 1.0. path :" << path << "\n";
                    output_counts = 1.0;
                }
            }
            else
            {
                (*spec_line)[j].output_counts(output_counts);
            }
        }
        else if(ltype == E_load_type::INTEGRATED)
        {
            output_counts += ((float)ii) / elapsed_realtime;
        }
        else if(ltype == E_load_type::CALLBACKF && callback_fun != nullptr)
        {
            callback_spectra->output_counts(((float)ii) / elapsed_realtime);
        }

        if (ltype == E_load_type::LINE)
        {
            (*spec_line)[j].elapsed_livetime(elapsed_livetime);
            (*spec_line)[j].elapsed_realtime(elapsed_realtime);
            (*spec_line)[j].input_counts(input_counts);
            (*spec_line)[j].output_counts(output_counts);
            // recalculate elapsed lifetime
            (*spec_line)[j].recalc_elapsed_livetime();
        }

        l += header_size + (spectra_size * detector);
        if (ltype == E_load_type::LINE)
        {
            for (size_t k = 0; k < spectra_size; k++)
            {
                (*spec_line)[j][k] = data_in[0][0][l+k];
            }
        }
        else if (ltype == E_load_type::INTEGRATED)
        {
            for (size_t k = 0; k < spectra_size; k++)
            {
                (*spectra)(k) += data_in[0][0][l+k];
            }
        }
        else if(ltype == E_load_type::CALLBACKF && callback_fun != nullptr)
        {
            callback_spectra->recalc_elapsed_livetime();
            (*callback_fun)(cur_row, j, max_rows, line_size, detector, callback_spectra, user_data);
        }

        l+=spectra_size * (MAX_NUM_SUPPORTED_DETECOTRS_PER_COL - detector);
    }

    if (ltype == E_load_type::INTEGRATED)
    {
        if(spectra->elapsed_livetime() == default_time_and_io_counts) //first spectra being loaded
        {
            spectra->elapsed_livetime(elapsed_livetime);
            spectra->elapsed_realtime(elapsed_realtime);
            spectra->input_counts(input_counts);
            spectra->output_counts(output_counts);
        }
        else
        {
            spectra->elapsed_livetime(spectra->elapsed_livetime() + elapsed_livetime);
            spectra->elapsed_realtime(spectra->elapsed_realtime() + elapsed_realtime);
            spectra->input_counts(spectra->input_counts() + input_counts);
            spectra->output_counts(spectra->output_counts() + output_counts);
        }
    }

    if ((retval = nc_close(ncid)))
    {
        logE<<" path :"<<path<<" : "<< nc_strerror(retval)<<"\n";
        return spec_cntr;
    }
    return spec_cntr;
}

//-----------------------------------------------------------------------------

template<typename T_real>
size_t NetCDF_IO<T_real>::load_spectra_line(std::string path, size_t detector, data_struct::Spectra_Line<T_real>* spec_line)
{
    return _load_spectra(E_load_type::LINE, path, detector, spec_line, -1, nullptr, 0, 0, nullptr, nullptr);
}

//-----------------------------------------------------------------------------

template<typename T_real>
size_t NetCDF_IO<T_real>::load_spectra_line_integrated(std::string path, size_t detector, size_t line_size, data_struct::Spectra<T_real>* spectra)
{
    return _load_spectra(E_load_type::INTEGRATED, path, detector, nullptr, line_size, spectra, 0, 0, nullptr, nullptr);
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool NetCDF_IO<T_real>::load_spectra_line_with_callback(std::string path,
												const std::vector<size_t>& detector_num_arr,
                                                int row,
                                                size_t max_rows,
                                                size_t max_cols,
                                                data_struct::IO_Callback_Func_Def<T_real> callback_fun,
                                                void* user_data)
{
    if (detector_num_arr.size() > 0)
    {
        size_t detector = detector_num_arr[0];
        return _load_spectra(E_load_type::CALLBACKF, path, detector, nullptr, max_cols, nullptr, row, max_rows, &callback_fun, nullptr);
    }
    return false;
}

//-----------------------------------------------------------------------------
 
template<typename T_real>
size_t NetCDF_IO<T_real>::load_scalers_line(std::string path, std::string tag, size_t row, data_struct::Scan_Info<T_real>* scan_info)
{
    std::lock_guard<std::mutex> lock(_mutex);

    size_t header_size = 256;
    int ncid = 0, varid = 0, retval = 0;
    size_t start[] = {0, 0, 0};
    size_t count[] = {1, 1, 1};
    ptrdiff_t stride[] = {1, 1, 1};
    //T_real data_in[10000][1][11];
    T_real *data_in = nullptr;
    size_t dim2size[NC_MAX_VAR_DIMS] = {0};
    size_t col_size = 0;
    size_t scalers_size = 0;
    nc_type rh_type;
    int rh_ndims = 0;
    int  rh_dimids[NC_MAX_VAR_DIMS] = {0};
    int rh_natts = 0;

    if(scan_info == nullptr)
    {
        logW<<"Scan info is null. Stopping load : "<<path<<"\n";
        return 0;
    }

    if( (retval = nc_open(path.c_str(), NC_NOWRITE, &ncid)) != 0)
    {
        logE<<path<<" :: "<< nc_strerror(retval)<<"\n";
        return 0;
    }

    if( (retval = nc_inq_varid(ncid, "array_data", &varid)) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return 0;
    }

    if( (retval = nc_inq_var (ncid, varid, nullptr, &rh_type, &rh_ndims, rh_dimids, &rh_natts) ) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return 0;
    }

    for (int i=0; i <  rh_ndims; i++)
    {
        if( (retval = nc_inq_dimlen(ncid, rh_dimids[i], &dim2size[i]) ) != 0)
        {
            logE<< path << " :: " << nc_strerror(retval)<<"\n";
            nc_close(ncid);
            return 0;
        }
    }


    count[0] = dim2size[0]; // num cols
    count[2] = 1;
    data_in = new T_real[dim2size[0]];

    //count[2] = dim2size[2]; // num scalers
    //data_in = new T_real[dim2size[0] * dim2size[1] * dim2size[2]];

    for(size_t i=0; i < dim2size[2]; i++)
    {
        start[2] = i;
        if( (retval = _nc_get_vars_real(ncid, varid, start, count, stride, data_in) ) != 0)
        {
            delete []data_in;
            logE<< path << " :: " << nc_strerror(retval)<<"\n";
            nc_close(ncid);
            return 0;
        }
        std::string search_name = tag+std::to_string(i);
        for(auto& scaler_map : scan_info->scaler_maps)
        {
            if(scaler_map.unit == search_name)
            {
                for(size_t j=0; j < dim2size[0]; j++)
                {
                    scaler_map.values(row,j) = data_in[j];
                }
                break;
            }
        }
    }

    delete []data_in;
    return count[0];
}

//-----------------------------------------------------------------------------

TEMPLATE_CLASS_DLL_EXPORT NetCDF_IO<float>;
TEMPLATE_CLASS_DLL_EXPORT NetCDF_IO<double>;

} //end namespace file
}// end namespace io
