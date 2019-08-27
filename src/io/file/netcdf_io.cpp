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

std::mutex NetCDF_IO::_mutex;

NetCDF_IO* NetCDF_IO::_this_inst(nullptr);


#define ELAPSED_REALTIME_OFFSET 32
#define ELAPSED_LIVETIME_OFFSET 34
#define INPUT_COUNTS_OFFSET 36
#define OUTPUT_COUNTS_OFFSET 38

#define MAX_NUM_SUPPORTED_DETECOTRS 4

//-----------------------------------------------------------------------------

NetCDF_IO::NetCDF_IO()
{

}

//-----------------------------------------------------------------------------

NetCDF_IO* NetCDF_IO::inst()
{
    std::lock_guard<std::mutex> lock(_mutex);

    if (_this_inst == nullptr)
    {
        _this_inst = new NetCDF_IO();
    }
    return _this_inst;
}

//-----------------------------------------------------------------------------

size_t NetCDF_IO::load_spectra_line(std::string path, size_t detector, data_struct::Spectra_Line* spec_line)
{
    std::lock_guard<std::mutex> lock(_mutex);

    size_t header_size = 256;
    int ncid, varid, retval;
    size_t start[] = {0, 0, 0};
    size_t count[] = {1, 1, header_size};
    ptrdiff_t stride[] = {1, 1, 1};
    real_t data_in[1][1][10000];
    size_t spectra_size;

    nc_type rh_type;
    int rh_ndims;
    int  rh_dimids[NC_MAX_VAR_DIMS] = {0};
    int rh_natts;
    real_t elapsed_livetime;
    real_t elapsed_realtime;
    real_t input_counts;
    real_t output_counts;

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

    if( (retval = nc_get_vars_real(ncid, varid, start, count, stride, &data_in[0][0][0]) ) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return 0;
    }

    if (data_in[0][0][0] != 21930 || data_in[0][0][1] != -21931)
    {
        logE<<"NetCDF header not found! Stopping load : "<<path<<"\n";
        nc_close(ncid);
        return 0;
    }

    header_size = data_in[0][0][2];
    //num_cols = data_in[][0][8];  //sum all across the first dim looking at value 8
    spectra_size = data_in[0][0][20];

    /*
    if( num_cols != spec_line->size() )
    {
        logW<<"Number of columns in NetCDF are "<<num_cols<<". Number of columns in spectra line are "<<spec_line->size()<< "\n";
    }
    */
    start[2] += header_size;
    count[2] = header_size;
    size_t j=0;

    for(; j<spec_line->size(); j++)
    {
        (*spec_line)[j].resize(spectra_size); // should be renames to resize

        //read header
        if( (retval = nc_get_vars_real(ncid, varid, start, count, stride, &data_in[0][0][0]) ) != 0)
        {
            logE<< nc_strerror(retval)<<"\n";
            nc_close(ncid);
            return j;
        }

        if (data_in[0][0][0] != 13260 || data_in[0][0][1] != -13261)
        {
            if(j < spec_line->size() -2)
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

        unsigned short i1 = data_in[0][0][ELAPSED_LIVETIME_OFFSET+(detector*8)];
        unsigned short i2 = data_in[0][0][ELAPSED_LIVETIME_OFFSET+(detector*8)+1];
        unsigned int ii = i1 | i2<<16;
        elapsed_livetime = ((float)ii) * 320e-9f; // need to multiply by this value becuase of the way it is saved
        if(elapsed_livetime == 0)
        {
            if(j < spec_line->size()-2) // usually the last two are missing which spams the log ouput.
            {
                logW<<"Reading in elapsed lifetime for Col:"<<j<<" is 0. Setting it to 1.0. path :"<<path<<"\n";
                elapsed_livetime = 1.0;
            }
        }
        (*spec_line)[j].elapsed_livetime(elapsed_livetime);

        i1 = data_in[0][0][ELAPSED_REALTIME_OFFSET+(detector*8)];
        i2 = data_in[0][0][ELAPSED_REALTIME_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        elapsed_realtime = ((float)ii) * 320e-9f; // need to multiply by this value becuase of the way it is saved
        if(elapsed_realtime == 0)
        {
            if(j < spec_line->size()-2) // usually the last two are missing which spams the log ouput.
            {
                logW<<"Reading in elapsed realtime for Col:"<<j<<" is 0. Setting it to 1.0. path :"<<path<<"\n";
                elapsed_realtime = 1.0;
            }
        }
        (*spec_line)[j].elapsed_realtime(elapsed_realtime);


        i1 = data_in[0][0][INPUT_COUNTS_OFFSET+(detector*8)];
        i2 = data_in[0][0][INPUT_COUNTS_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        input_counts = ((float)ii) / elapsed_livetime;
        (*spec_line)[j].input_counts(input_counts);

        i1 = data_in[0][0][OUTPUT_COUNTS_OFFSET+(detector*8)];
        i2 = data_in[0][0][OUTPUT_COUNTS_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        output_counts = ((float)ii) / elapsed_realtime;
        (*spec_line)[j].output_counts(output_counts);

        // recalculate elapsed lifetime
        (*spec_line)[j].recalc_elapsed_livetime();

        start[2] += header_size + (spectra_size * detector);
        count[2] = spectra_size;

        if( (retval = nc_get_vars_real(ncid, varid, start, count, stride, &data_in[0][0][0]) ) != 0)
        {
            logE<<" path :"<<path<<" : "<< nc_strerror(retval)<<"\n";
            nc_close(ncid);
            return j;
        }

        for(size_t k=0; k<spectra_size; k++)
        {
            (*spec_line)[j][k] = data_in[0][0][k];
        }

        start[2] += spectra_size * (4 - detector);
        count[2] = header_size;

        if(start[2] >= dim2size[2])
        {
            start[0]++;
            start[2] = header_size;
        }

    }

    if ((retval = nc_close(ncid)))
    {
        logE<<" path :"<<path<<" : "<< nc_strerror(retval)<<"\n";
        return j;
    }
    return j;
}

//-----------------------------------------------------------------------------

bool NetCDF_IO::load_spectra_line_with_callback(std::string path,
												const std::vector<size_t>& detector_num_arr,
                                                int row,
                                                size_t max_rows,
                                                size_t max_cols,
                                                data_struct::IO_Callback_Func_Def callback_fun,
                                                void* user_data)
{

    std::lock_guard<std::mutex> lock(_mutex);

    size_t header_size = 256;
    int ncid, varid, retval;
    size_t start[] = {0, 0, 0};
    size_t count[] = {1, 1, header_size};
    ptrdiff_t stride[] = {1, 1, 1};
    real_t data_in[1][1][10000];
    size_t spectra_size;
	size_t num_detectors = detector_num_arr.size();

    nc_type rh_type;
    int rh_ndims;
    int  rh_dimids[NC_MAX_VAR_DIMS] = {0};
    int rh_natts;

    real_t elapsed_livetime[MAX_NUM_SUPPORTED_DETECOTRS];
    real_t elapsed_realtime[MAX_NUM_SUPPORTED_DETECOTRS];
    real_t input_counts[MAX_NUM_SUPPORTED_DETECOTRS];
    real_t output_counts[MAX_NUM_SUPPORTED_DETECOTRS];

    size_t dim2size[NC_MAX_VAR_DIMS] = {0};

    if( num_detectors > 4)
    {
        logE<<"Max detectors supported is 4, requesting :"<< num_detectors<<"\n";
        return false;
    }

    if( (retval = nc_open(path.c_str(), NC_NOWRITE, &ncid)) != 0 )
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        return false;
    }

    if( (retval = nc_inq_varid(ncid, "array_data", &varid)) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return false;
    }

    if( (retval = nc_inq_var (ncid, varid, 0, &rh_type, &rh_ndims, rh_dimids, &rh_natts) ) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return false;
    }

    for (int i=0; i <  rh_ndims; i++)
    {
        if( (retval = nc_inq_dimlen(ncid, rh_dimids[i], &dim2size[i]) ) != 0)
        {
            logE<< path << " :: " << nc_strerror(retval)<<"\n";
            nc_close(ncid);
            return false;
        }
    }

    //read in last col sector to get total number of cols
    //start[0] = dim2size[0] - 1;
    if( (retval = nc_get_vars_real(ncid, varid, start, count, stride, &data_in[0][0][0]) ) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        nc_close(ncid);
        return false;
    }

    if (data_in[0][0][0] != 21930 || data_in[0][0][1] != -21931)
    {
        logE<<"NetCDF header not found! Stopping load : "<<path<<"\n";
        nc_close(ncid);
        return false;
    }

    //can't read from file because it can change inbetween rows ...
    //max_cols = (124 * (dim2size[0] - 1) ) + data_in[0][0][8];
    header_size = data_in[0][0][2];
    //num_cols = data_in[][0][8];  //sum all across the first dim looking at value 8
    spectra_size = data_in[0][0][20];

    start[2] += count[2];
    count[2] = header_size + (spectra_size * MAX_NUM_SUPPORTED_DETECOTRS); //only 4 element detector supported
    //loop through col sectors
    for(size_t j = 0; j < max_cols; j++)
    {
        //read header
        if( (retval = nc_get_vars_real(ncid, varid, start, count, stride, &data_in[0][0][0]) ) != 0)
        {
            logE<< path << " :: " << nc_strerror(retval)<<"\n";
            nc_close(ncid);
            return false;
        }

        if (data_in[0][0][0] != 13260 || data_in[0][0][1] != -13261)
        {
            if(j < max_cols -2)
            {
                logE<<"NetCDF sub header not found! Stopping load at Col: "<<j<<" path :"<<path<<"\n";
                nc_close(ncid);
                return false;
            }
            //last two may not be filled with data
            //TODO: send end of row stream_block down pipeline
            nc_close(ncid);
            return true;
        }

        for(size_t detector_num : detector_num_arr)
        {
            unsigned short i1 = data_in[0][0][ELAPSED_LIVETIME_OFFSET+(detector_num*8)];
            unsigned short i2 = data_in[0][0][ELAPSED_LIVETIME_OFFSET+(detector_num*8)+1];
            unsigned int ii = i1 | i2<<16;
            elapsed_livetime[detector_num] = ((float)ii) * 320e-9f; // need to multiply by this value becuase of the way it is saved
            if(elapsed_livetime[detector_num] == 0)
            {
                if(j < max_cols-2) // usually the last two are missing which spams the log ouput.
                {
                    logW<<"Reading in elapsed lifetime for Col:"<<j<<" is 0. Setting it to 1.0 " << path <<"\n";
                    elapsed_livetime[detector_num] = 1.0;
                }
            }

            i1 = data_in[0][0][ELAPSED_REALTIME_OFFSET+(detector_num*8)];
            i2 = data_in[0][0][ELAPSED_REALTIME_OFFSET+(detector_num*8)+1];
            ii = i1 | i2<<16;
            elapsed_realtime[detector_num] = ((float)ii) * 320e-9f; // need to multiply by this value becuase of the way it is saved
            if(elapsed_realtime[detector_num] == 0)
            {
                if(j < max_cols-2) // usually the last two are missing which spams the log ouput.
                {
                    logW<<"Reading in elapsed realtime for Col:"<<j<<" is 0. Setting it to 1.0 "<<path<<"\n";
                    elapsed_realtime[detector_num] = 1.0;
                }
            }

            i1 = data_in[0][0][INPUT_COUNTS_OFFSET+(detector_num*8)];
            i2 = data_in[0][0][INPUT_COUNTS_OFFSET+(detector_num*8)+1];
            ii = i1 | i2<<16;
            input_counts[detector_num] = ((float)ii) / elapsed_livetime[detector_num];

            i1 = data_in[0][0][OUTPUT_COUNTS_OFFSET+(detector_num*8)];
            i2 = data_in[0][0][OUTPUT_COUNTS_OFFSET+(detector_num*8)+1];
            ii = i1 | i2<<16;
            output_counts[detector_num] = ((float)ii) / elapsed_realtime[detector_num];

            data_struct::Spectra * spectra = new data_struct::Spectra(spectra_size);

            spectra->elapsed_livetime(elapsed_livetime[detector_num]);
            spectra->elapsed_realtime(elapsed_realtime[detector_num]);
            spectra->input_counts(input_counts[detector_num]);
            spectra->output_counts(output_counts[detector_num]);
            spectra->recalc_elapsed_livetime();

            int idx = header_size + (detector_num*spectra_size);

            for(size_t k=0; k<spectra_size; k++)
            {
                (*spectra)[k] = data_in[0][0][idx+k];
            }

            callback_fun(row, j, max_rows, max_cols, detector_num, spectra, user_data);

        }

        start[2] += count[2];

        if(start[2] >= dim2size[2])
        {
            start[0]++;
            start[2] = header_size;
        }
    }

    if((retval = nc_close(ncid)) != 0)
    {
        logE<< path << " :: " << nc_strerror(retval)<<"\n";
        return false;
    }

    return true;

}

//-----------------------------------------------------------------------------

} //end namespace file
}// end namespace io
