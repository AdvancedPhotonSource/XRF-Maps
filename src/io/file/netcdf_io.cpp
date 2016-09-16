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

#define ERR(e) {std::cout<<"Error: "<< nc_strerror(e)<<std::endl; return false;}

#define ELAPSED_REALTIME_OFFSET 32
#define ELAPSED_LIFETIME_OFFSET 34
#define INPUT_COUNTS_OFFSET 36
#define OUTPUT_COUNTS_OFFSET 38

NetCDF_IO::NetCDF_IO() : Base_File_IO()
{

}

void NetCDF_IO::lazy_load()
{

}

bool NetCDF_IO::load_dataset(std::string path, Base_Dataset *dset)
{
    return false;

}

bool NetCDF_IO::load_spectra_line(std::string path, size_t detector, data_struct::xrf::Spectra_Line* spec_line)
{

    size_t header_size = 256;
    int ncid, varid, retval;
    size_t start[] = {0, 0, 0};
    size_t count[] = {1, 1, header_size};
    ptrdiff_t stride[] = {1, 1, 1};
    double data_in[1][1][5000];
    size_t num_cols;
    size_t spectra_size;

    nc_type rh_type;
    int rh_ndims;
    int  rh_dimids[NC_MAX_VAR_DIMS] = {0};
    int rh_natts;
    real_t elapsed_lifetime;
    real_t elapsed_realtime;
    real_t input_counts;
    real_t output_counts;

    size_t dim2size[NC_MAX_VAR_DIMS] = {0};

    if( (retval = nc_open(path.c_str(), NC_NOWRITE, &ncid)) )
       ERR(retval);

    if( (retval = nc_inq_varid(ncid, "array_data", &varid)) )
        ERR(retval);

    if( retval = nc_inq_var (ncid, varid, 0, &rh_type, &rh_ndims, rh_dimids, &rh_natts) )
       ERR(retval);

    for (int i=0; i <  rh_ndims; i++)
    {
        if( retval = nc_inq_dimlen(ncid, rh_dimids[i], &dim2size[i]) )
            ERR(retval);
    }

    if( retval = nc_get_vars_double(ncid, varid, start, count, stride, &data_in[0][0][0]) )
       ERR(retval);

    if (data_in[0][0][0] != 21930 || data_in[0][0][1] != -21931)
    {
        std::cout<<"Error: NetCDF header not found! Stopping load : "<<path<<std::endl;
        return false;
    }

    header_size = data_in[0][0][2];
    num_cols = data_in[0][0][8];
    spectra_size = data_in[0][0][20];

    /*
    if( num_cols != spec_line->size() )
    {
        std::cout<<"Warning: Number of columns in NetCDF are "<<num_cols<<". Number of columns in spectra line are "<<spec_line->size()<< std::endl;
    }
    */
    start[2] += header_size;

    for(size_t j=0; j<spec_line->size(); j++)
    {
        (*spec_line)[j].resize(spectra_size); // should be renames to resize

        count[2] = header_size;
        //read header
        if( retval = nc_get_vars_double(ncid, varid, start, count, stride, &data_in[0][0][0]) )
           ERR(retval);

        header_size = data_in[0][0][2];

        unsigned short i1 = data_in[0][0][ELAPSED_LIFETIME_OFFSET+(detector*8)];
        unsigned short i2 = data_in[0][0][ELAPSED_LIFETIME_OFFSET+(detector*8)+1];
        unsigned int ii = i1 | i2<<16;
        elapsed_lifetime = ((float)ii) * 320e-9; // need to multiply by this value becuase of the way it is saved
        (*spec_line)[j].elapsed_lifetime(elapsed_lifetime);

        i1 = data_in[0][0][ELAPSED_REALTIME_OFFSET+(detector*8)];
        i2 = data_in[0][0][ELAPSED_REALTIME_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        elapsed_realtime = ((float)ii) * 320e-9; // need to multiply by this value becuase of the way it is saved
        (*spec_line)[j].elapsed_realtime(elapsed_realtime);

        i1 = data_in[0][0][INPUT_COUNTS_OFFSET+(detector*8)];
        i2 = data_in[0][0][INPUT_COUNTS_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        input_counts = ((float)ii) / elapsed_lifetime;
        (*spec_line)[j].input_counts(input_counts);

        i1 = data_in[0][0][OUTPUT_COUNTS_OFFSET+(detector*8)];
        i2 = data_in[0][0][OUTPUT_COUNTS_OFFSET+(detector*8)+1];
        ii = i1 | i2<<16;
        output_counts = ((float)ii) / elapsed_realtime;
        (*spec_line)[j].output_counts(output_counts);

        // recalculate elapsed lifetime
        (*spec_line)[j].recalc_elapsed_lifetime();

        start[2] += header_size + (spectra_size * detector);
        count[2] = spectra_size;

        if( retval = nc_get_vars_double(ncid, varid, start, count, stride, &data_in[0][0][0]) )
           ERR(retval);

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
       ERR(retval);

}

} //end namespace file
}// end namespace io
