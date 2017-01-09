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


#include "hdf5_io.h"

#include <iostream>
#include <string>

#include <chrono>
#include <ctime>
#include <thread>


#include "element_info.h"

#define HDF5_SAVE_VERSION 10.0

namespace io
{
namespace file
{

//-----------------------------------------------------------------------------

HDF5_IO::HDF5_IO() : Base_File_IO()
{
	//disable hdf print to std err
	hid_t status;
	status = H5Eset_auto(H5E_DEFAULT, NULL, NULL);
}

//-----------------------------------------------------------------------------

void HDF5_IO::lazy_load()
{
   _is_loaded = ERROR_LOADING;

   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   //std::cout<<"lazy_load "<< _filename <<std::endl;

   hid_t    file_id, dset_id, dataspace_id, memoryspace, datatype;
   herr_t   error;

    H5T_class_t dtype_class;
   //H5T_order_t order;
   //size_t      size;

    //std::cout<<"pre open file"<< std::endl;
    file_id = H5Fopen(_filename.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    //std::cout<<"pre open dset "<<std::endl;
    dset_id = H5Dopen2(file_id, "/MAPS_RAW/data_a", H5P_DEFAULT);


    datatype = H5Dget_type(dset_id);
    dtype_class = H5Tget_class(datatype);
    //if (dtype_class == H5T_INTEGER)
    //   printf("Data set has INTEGER type \n");
    //order = H5Tget_order(datatype);
    //if (order == H5T_ORDER_LE)
    //   printf("Little endian order \n");
    dataspace_id = H5Dget_space(dset_id);
    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    hsize_t* dims_out = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    hsize_t* sel_dims = new hsize_t[rank];
    //std::cout<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_out[0], NULL);

    unsigned long total = 1;
    for (int i=0; i < rank; i++)
    {
       std::cout<<"dims ["<<i<<"] ="<<dims_out[i]<< std::endl;
       total *= dims_out[i];
       offset[i] = 0;
       sel_dims[i] = count[i] = dims_out[i];
    }
    //std::cout<<"total = "<<total<<std::endl;

    //unsigned short *buffer = new unsigned short[total];
    float *buffer = new float[total];
    //std::cout<<"allocated "<<std::endl;
    memoryspace = H5Screate_simple(rank, sel_dims, NULL);

    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    //error = H5Dread(dset_id, H5T_NATIVE_USHORT, memoryspace, dataspace_id, H5P_DEFAULT, buffer);
    error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, buffer);
    printf("read in: error = %d \n", error);
/*
    for(unsigned long i=0; i<total; i++)
    {
       std::cout<<buffer[i]<<" ";
       if( i%2048 == 1)
          std::cout<<std::endl;
    }
*/
    delete [] dims_out;
    delete [] offset;
    delete [] count;
    delete [] sel_dims;

    H5Dclose(dset_id);
    H5Sclose(memoryspace);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
                 << "elapsed time: " << elapsed_seconds.count() << "s\n";

    _is_loaded = LAZY_LOAD;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_dataset(std::string path, Base_Dataset *dset)
{
    return false;

}

//-----------------------------------------------------------------------------

/*
void HDF5_IO::load_dataset2(std::string path, HDF5_Spectra_Layout layout, HDF5_Range range, data_struct::xrf::Spectra_Volume* spec_vol)
{
     _is_loaded = ERROR_LOADING;

     std::chrono::time_point<std::chrono::system_clock> start, end;
     start = std::chrono::system_clock::now();

     hid_t    file_id, dset_id, dataspace_id, memoryspace, datatype;
     herr_t   error;

     H5T_class_t dtype_class;

     file_id = H5Fopen(_filename.c_str(), H5P_DEFAULT, H5P_DEFAULT);
     dset_id = H5Dopen2(file_id, "/MAPS_RAW/data_a", H5P_DEFAULT);

     datatype = H5Dget_type(dset_id);
     dtype_class = H5Tget_class(datatype);

     dataspace_id = H5Dget_space(dset_id);
     int rank = H5Sget_simple_extent_ndims(dataspace_id);
     hsize_t* dims_out = new hsize_t[rank];
     hsize_t* offset = new hsize_t[rank];
     hsize_t* count = new hsize_t[rank];
     hsize_t* sel_dims = new hsize_t[rank];
     //std::cout<<"rank = "<<rank<< std::endl;
     unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_out[0], NULL);

     spec_vol->resize(dims_out[layout.row_dim],dims_out[layout.col_dim],dims_out[layout.spectrum_dim]);

//float bufs[2048][11][9];
//float* pbufs = new float[2048 * 11 * 9];

float *** ppbufs = new float**[dims_out[0]];
for (int s=0; s<dims_out[0]; s++)
{
    ppbufs[s] = new float*[dims_out[1]];
    for (int r =0; r<dims_out[1]; r++)
    {
        ppbufs[s][r] = new float[dims_out[2]];
    }
}


//std::vector<std::vector<std::vector<float> > > *vec = spec_vol->get_vector();
     for (int i=0; i < rank; i++)
     {
        std::cout<<"dims ["<<i<<"] ="<<dims_out[i]<< std::endl;
        offset[i] = 0;
        sel_dims[i] = count[i] = dims_out[i];
     }
     //int total = (dims_out[layout.row_dim]*dims_out[layout.col_dim]*dims_out[layout.spectrum_dim]);
     sel_dims[0] = dims_out[0];
     sel_dims[1] = dims_out[1];
     sel_dims[2] = dims_out[2];

     count[0] = dims_out[0];
     count[1] = dims_out[1];
     count[2] = dims_out[2];

     memoryspace = H5Screate_simple(rank, sel_dims, NULL);

     H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, sel_dims, NULL);

  //   H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, sel_dims, NULL);
     H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

     real_t *buf = spec_vol->get_buffer_ptr();
     //float *buf = new float[dims_out[0] * dims_out[1] * dims_out[2]];
     //  error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, buf);
#ifdef _REAL_FLOAT
       error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, ppbufs);
#else
       error = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memoryspace, dataspace_id, H5P_DEFAULT, ppbufs);
#endif
     //     spec_vol->transpose(buf);


     //delete [] buf;
     //error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, buf[0][0][0]);
//          error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, &(bufs[0][0][0]));
 //         error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, pbufs);


     //    error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, &(ppbufs[0][0][0]));

 //         error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, &((*vec)[0][0][0]));
    // printf("read in: error = %d \n", error);


    // std::cout<<"A[115] "<<buf[115][0][0]<<std::endl;
     //std::cout<<"A[115] "<<(*vec)[115][0][0]<<std::endl;
     //i*ny*nz + j*nz + k

//     int i = 115;
//     int j = 0;
//     int k = 0;
//     int ny = 11;
//     int nz = 9;

 //    std::cout<<"A[115] "<<pbufs[ i*ny*nz + j*nz + k  ]<<std::endl;
     //std::cout<<"B[115] "<<ppbufs[0][0][115]<<std::endl;
     //std::cout<<"B[115] "<<ppbufs[0][1][1]<<std::endl;
     //std::cout<<"B[115] "<<ppbufs[0][0]<<std::endl;

//     for(int i=0; i<2048; i++)
//     {
//        std::cout<<ppbufs[0][0][i]<<std::endl;
//     }

 //    std::cout<<"C[115] "<<bufs[115][0][0]<<std::endl;

    for(unsigned int c=0; c<dims_out[layout.col_dim]; c++)
    {
        for(unsigned int r=0; r<dims_out[layout.row_dim]; r++)
        {
            //float *spec_buf = spec_vol->get_spectra(k, j);
            std::cout<<"|||||"<<std::endl;
            for(unsigned int s=0; s<dims_out[layout.spectrum_dim]; s++)
            {
     //           std::cout<<spec_buf[i]<<" ";
                  std::cout<<buf[s + (dims_out[layout.spectrum_dim] * dims_out[layout.row_dim] * c)+(dims_out[layout.spectrum_dim] * r) ]<<" ";
                //std::cout<<bufs[i][k][j]<<" ";

            }
            std::cout<<std::endl<<std::endl;
        }
    }




     float * spec_buf = spec_vol->get_spectra(0,0);
     for (int i=1004; i<1083; i++)
        std::cout<<"buff["<<i<<"] = "<<spec_buf[i]<<std::endl;

       for (int i=1004; i<1083; i++)
          std::cout<<"buff["<<i<<"] = "<<ppbufs[i][0][0]<<std::endl;

     delete [] dims_out;
     delete [] offset;
     delete [] count;
     delete [] sel_dims;

     H5Dclose(dset_id);
     H5Sclose(memoryspace);
     H5Sclose(dataspace_id);
     H5Fclose(file_id);

     end = std::chrono::system_clock::now();
     std::chrono::duration<double> elapsed_seconds = end-start;
     std::time_t end_time = std::chrono::system_clock::to_time_t(end);

     std::cout << "finished computation at " << std::ctime(&end_time)
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";

     _is_loaded = LAZY_LOAD;
}
*/

//-----------------------------------------------------------------------------

//void HDF5_IO::parse_group_info(hid_t h5file, hid_t id)
//{

//}

////-----------------------------------------------------------------------------

//void HDF5_IO::parse_dataset_info(hid_t h5file, hid_t id)
//{
//    /*
//    DataSet dataset = h5file.openDataSet(DATASET_NAME);

//    H5T_class_t type_class = dataset.getTypeClass();

//    DataSpace dataspace = dataset.getSpace();

//    int rank = dataspace.getSimpleExtentNdims();

//    hsize_t *dims_out = new hsize_t[rank];
//    int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
//*/
//}

////-----------------------------------------------------------------------------

//void HDF5_IO::parse_attr_info(hid_t h5file, hid_t id)
//{

//}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume(std::string path, size_t detector_num, data_struct::xrf::Spectra_Volume* spec_vol)
{

   //_is_loaded = ERROR_LOADING;

   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   std::cout<<"load_spectra_volume "<< path <<" detector : "<<detector_num<<std::endl;

   hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
   hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
   herr_t   error;
   std::string detector_path;
   real_t * buffer;
   hsize_t offset_row[2] = {0,0};
   hsize_t count_row[2] = {0,0};
   hsize_t offset_meta[3] = {0,0,0};
   hsize_t count_meta[3] = {1,1,1};


   switch(detector_num)
   {
   case 0:
       detector_path = "data_a";
       break;
   case 1:
       detector_path = "data_b";
       break;
   case 2:
       detector_path = "data_c";
       break;
   case 3:
       detector_path = "data_d";
       break;
   default:
       detector_path = "";
       break;
   }

    file_id = H5Fopen(path.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    if(file_id < 0)
    {

        std::cout<<"Error opening file "<<path<<std::endl;
        return false;
    }

    maps_grp_id = H5Gopen(file_id, "MAPS_RAW", H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        std::cout<<"Error opening group MAPS_RAW"<<std::endl;
        return false;
    }


    std::cout<<"pre open dset "<<std::endl;
    //std::cout<<"pre open dset "<<std::endl;
    dset_id = H5Dopen2(maps_grp_id, detector_path.c_str(), H5P_DEFAULT);
    if(dset_id < 0)
    {
        std::cout<<"Error opening dataset /MAPS_RAW/"<<detector_path<<std::endl;
        return false;
    }
    dataspace_id = H5Dget_space(dset_id);

    dset_lt_id = H5Dopen2(maps_grp_id, "livetime", H5P_DEFAULT);
    if(dset_lt_id < 0)
    {
        std::cout<<"Error opening dataset /MAPS_RAW/livetime"<<std::endl;
        return false;
    }
    dataspace_lt_id = H5Dget_space(dset_lt_id);

    dset_rt_id = H5Dopen2(maps_grp_id, "realtime", H5P_DEFAULT);
    if(dset_rt_id < 0)
    {
        std::cout<<"Error opening dataset /MAPS_RAW/realtime"<<std::endl;
        return false;
    }
    dataspace_rt_id = H5Dget_space(dset_rt_id);

    dset_incnt_id = H5Dopen2(maps_grp_id, "inputcounts", H5P_DEFAULT);
    if(dset_incnt_id < 0)
    {
        std::cout<<"Error opening dataset /MAPS_RAW/inputcounts"<<std::endl;
        return false;
    }
    dataspace_inct_id = H5Dget_space(dset_incnt_id);

    dset_outcnt_id = H5Dopen2(maps_grp_id, "ouputcounts", H5P_DEFAULT);
    if(dset_outcnt_id < 0)
    {
        std::cout<<"Error opening dataset /MAPS_RAW/ouputcounts"<<std::endl;
        return false;
    }
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);


    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        std::cout<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    std::cout<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
    if(status_n < 0)
    {
         std::cout<<"Error getting dataset rank for MAPS_RAW/"<< detector_path<<std::endl;
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       std::cout<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
    count_row[0] = dims_in[0];
    count_row[1] = dims_in[2];

    if(spec_vol->rows() < dims_in[1] || spec_vol->cols() < dims_in[2] || spec_vol->samples_size() < dims_in[0])
    {
        spec_vol->resize(dims_in[1], dims_in[2], dims_in[0]);
    }

    count[1] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, NULL);
    memoryspace_meta_id = H5Screate_simple(1, count_meta, NULL);
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    offset_meta[0] = detector_num;
    for (size_t row=0; row < dims_in[1]; row++)
    {
         offset[1] = row;
         offset_meta[1] = row;

         H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
         error = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

         if (error > -1 )
         {
             for(size_t col=0; col<count_row[1]; col++)
             {
                 offset_meta[2] = col;
                 data_struct::xrf::Spectra *spectra = &((*spec_vol)[row][col]);

                 H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_lt_id, H5T_NATIVE_DOUBLE, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                 spectra->elapsed_lifetime(live_time);

                 H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_rt_id, H5T_NATIVE_DOUBLE, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                 spectra->elapsed_realtime(real_time);

                 H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_incnt_id, H5T_NATIVE_DOUBLE, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                 spectra->input_counts(in_cnt);

                 H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                 error = H5Dread(dset_outcnt_id, H5T_NATIVE_DOUBLE, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                 spectra->output_counts(out_cnt);

                 for(size_t s=0; s<count_row[0]; s++)
                 {
                     (*spectra)[s] = buffer[(count_row[1] * s) + col];
                 }
                 //std::cout<<"saved col "<<col<<std::endl;
             }

            std::cout<<"read row "<<row<<std::endl;
         }
         else
         {
            std::cout<<"Error: reading row "<<row<<std::endl;
         }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;
    delete [] buffer;

    H5Dclose(dset_id);
    H5Dclose(dset_incnt_id);
    H5Dclose(dset_outcnt_id);
    H5Dclose(dset_rt_id);
    H5Dclose(dset_lt_id);
    H5Sclose(memoryspace_meta_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_lt_id);
    H5Sclose(dataspace_rt_id);
    H5Sclose(dataspace_inct_id);
    H5Sclose(dataspace_outct_id);
    H5Sclose(dataspace_id);
    H5Fclose(file_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
                 << "elapsed time: " << elapsed_seconds.count() << "s\n";

}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_and_integrate_spectra_volume(std::string path, size_t detector_num, data_struct::xrf::Spectra* spectra)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::cout<<"load_spectra_volume "<< path <<" detector : "<<detector_num<<std::endl;

    hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
    hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
    herr_t   error;
    std::string detector_path;
    real_t * buffer;
    hsize_t offset_row[2] = {0,0};
    hsize_t count_row[2] = {0,0};
    hsize_t offset_meta[3] = {0,0,0};
    hsize_t count_meta[3] = {1,1,1};


    switch(detector_num)
    {
    case 0:
        detector_path = "data_a";
        break;
    case 1:
        detector_path = "data_b";
        break;
    case 2:
        detector_path = "data_c";
        break;
    case 3:
        detector_path = "data_d";
        break;
    default:
        detector_path = "";
        break;
    }

     file_id = H5Fopen(path.c_str(), H5P_DEFAULT, H5P_DEFAULT);
     if(file_id < 0)
     {

         std::cout<<"Error opening file "<<path<<std::endl;
         return false;
     }

     maps_grp_id = H5Gopen(file_id, "MAPS_RAW", H5P_DEFAULT);
     if(maps_grp_id < 0)
     {
         std::cout<<"Error opening group MAPS_RAW"<<std::endl;
         return false;
     }


     std::cout<<"pre open dset "<<std::endl;
     //std::cout<<"pre open dset "<<std::endl;
     dset_id = H5Dopen2(maps_grp_id, detector_path.c_str(), H5P_DEFAULT);
     if(dset_id < 0)
     {
         std::cout<<"Error opening dataset /MAPS_RAW/"<<detector_path<<std::endl;
         return false;
     }
     dataspace_id = H5Dget_space(dset_id);

     dset_lt_id = H5Dopen2(maps_grp_id, "livetime", H5P_DEFAULT);
     if(dset_lt_id < 0)
     {
         std::cout<<"Error opening dataset /MAPS_RAW/livetime"<<std::endl;
         return false;
     }
     dataspace_lt_id = H5Dget_space(dset_lt_id);

     dset_rt_id = H5Dopen2(maps_grp_id, "realtime", H5P_DEFAULT);
     if(dset_rt_id < 0)
     {
         std::cout<<"Error opening dataset /MAPS_RAW/realtime"<<std::endl;
         return false;
     }
     dataspace_rt_id = H5Dget_space(dset_rt_id);

     dset_incnt_id = H5Dopen2(maps_grp_id, "inputcounts", H5P_DEFAULT);
     if(dset_incnt_id < 0)
     {
         std::cout<<"Error opening dataset /MAPS_RAW/inputcounts"<<std::endl;
         return false;
     }
     dataspace_inct_id = H5Dget_space(dset_incnt_id);

     dset_outcnt_id = H5Dopen2(maps_grp_id, "ouputcounts", H5P_DEFAULT);
     if(dset_outcnt_id < 0)
     {
         std::cout<<"Error opening dataset /MAPS_RAW/ouputcounts"<<std::endl;
         return false;
     }
     dataspace_outct_id = H5Dget_space(dset_outcnt_id);


     int rank = H5Sget_simple_extent_ndims(dataspace_id);
     if (rank != 3)
     {
         std::cout<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<std::endl;
         return false;
        //throw exception ("Dataset is not a volume");
     }
     hsize_t* dims_in = new hsize_t[rank];
     hsize_t* offset = new hsize_t[rank];
     hsize_t* count = new hsize_t[rank];
     std::cout<<"rank = "<<rank<< std::endl;
     unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
     if(status_n < 0)
     {
          std::cout<<"Error getting dataset rank for MAPS_RAW/"<< detector_path<<std::endl;
          return false;
     }

     for (int i=0; i < rank; i++)
     {
        std::cout<<"dims ["<<i<<"] ="<<dims_in[i]<< std::endl;
        offset[i] = 0;
        count[i] = dims_in[i];
     }

     buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
     count_row[0] = dims_in[0];
     count_row[1] = dims_in[2];

     count[1] = 1; //1 row

     memoryspace_id = H5Screate_simple(2, count_row, NULL);
     memoryspace_meta_id = H5Screate_simple(1, count_meta, NULL);
     H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, NULL, count_row, NULL);
     H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);

     real_t live_time = 1.0;
     real_t real_time = 1.0;
     real_t in_cnt = 1.0;
     real_t out_cnt = 1.0;

     real_t live_time_total = 0.0;
     real_t real_time_total = 0.0;
     real_t in_cnt_total = 0.0;
     real_t out_cnt_total = 0.0;

     offset_meta[0] = detector_num;
     for (size_t row=0; row < dims_in[1]; row++)
     {
          offset[1] = row;
          offset_meta[1] = row;

          H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
          error = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

          if (error > -1 )
          {
              for(size_t col=0; col<count_row[1]; col++)
              {
                  offset_meta[2] = col;

                  H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                  error = H5Dread(dset_lt_id, H5T_NATIVE_DOUBLE, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                  live_time_total += live_time;

                  H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                  error = H5Dread(dset_rt_id, H5T_NATIVE_DOUBLE, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                  real_time_total += real_time;

                  H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                  error = H5Dread(dset_incnt_id, H5T_NATIVE_DOUBLE, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                  in_cnt_total += in_cnt;

                  H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, NULL, count_meta, NULL);
                  error = H5Dread(dset_outcnt_id, H5T_NATIVE_DOUBLE, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                  out_cnt_total += out_cnt;

                  for(size_t s=0; s<count_row[0]; s++)
                  {
                      (*spectra)[s] += buffer[(count_row[1] * s) + col];
                  }
                  //std::cout<<"saved col "<<col<<std::endl;
              }

             std::cout<<"read row "<<row<<std::endl;
          }
          else
          {
             std::cout<<"Error: reading row "<<row<<std::endl;
          }
     }

     spectra->elapsed_lifetime(live_time_total);
     spectra->elapsed_realtime(real_time_total);
     spectra->input_counts(in_cnt_total);
     spectra->output_counts(out_cnt_total);

     delete [] dims_in;
     delete [] offset;
     delete [] count;
     delete [] buffer;

     H5Dclose(dset_id);
     H5Dclose(dset_incnt_id);
     H5Dclose(dset_outcnt_id);
     H5Dclose(dset_rt_id);
     H5Dclose(dset_lt_id);
     H5Sclose(memoryspace_meta_id);
     H5Sclose(memoryspace_id);
     H5Sclose(dataspace_lt_id);
     H5Sclose(dataspace_rt_id);
     H5Sclose(dataspace_inct_id);
     H5Sclose(dataspace_outct_id);
     H5Sclose(dataspace_id);
     H5Fclose(file_id);

     end = std::chrono::system_clock::now();
     std::chrono::duration<double> elapsed_seconds = end-start;
     std::time_t end_time = std::chrono::system_clock::to_time_t(end);

     std::cout << "finished computation at " << std::ctime(&end_time)
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";

     return true;
}

//-----------------------------------------------------------------------------

hid_t HDF5_IO::start_save_seq(const std::string filename)
{
    hid_t    file_id;

    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (file_id < 1)
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(file_id < 0)
    {
        std::cout<<"Error opening file "<<filename<<std::endl;
        return -1;
    }

    return file_id;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::end_save_seq(const hid_t file_id)
{

    H5Fclose(file_id);
    H5close(); //test first, might break h5
    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_spectra_volume(const hid_t file_id,
                                  const std::string path,
                                  data_struct::xrf::Spectra_Volume * spectra_volume,
                                  real_t energy_offset,
                                  real_t energy_slope,
                                  real_t energy_quad,
                                  size_t row_idx_start,
                                  int row_idx_end,
                                  size_t col_idx_start,
                                  int col_idx_end)
{
    std::unique_lock<std::mutex> lock(_mutex);

    std::cout<<"HDF5_IO::save_spectra_volume()"<<std::endl;

//herr_t (*old_func)(void*);
//void *old_client_data;
//hid_t error_stack = H5Eget_current_stack();
//H5Eget_auto2(error_stack, &old_func, &old_client_data);

//H5Eset_auto2(error_stack, NULL, NULL);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    dset_id, spec_grp_id, int_spec_grp_id, dataspace_id, memoryspace_id, memoryspace_time_id, dataspace_time_id, file_time_id, filespace_id, status, maps_grp_id, dcpl_id;
    hid_t   dset_rt_id, dset_lt_id, incnt_dset_id, outcnt_dset_id;
    herr_t   error;

    hsize_t chunk_dims[3];
    hsize_t dims_out[3];
    hsize_t offset[3];
    hsize_t count[3];
    hsize_t dims_time_out[2];
    hsize_t offset_time[2];
    hsize_t count_time[2];

    if (row_idx_end < row_idx_start || row_idx_end > spectra_volume->rows() -1)
    {
        row_idx_end = spectra_volume->rows();
    }
    if (col_idx_end < col_idx_start || col_idx_end > spectra_volume->cols() -1)
    {
        col_idx_end = spectra_volume->cols();
    }

    //get one element
    data_struct::xrf::Fit_Element_Map* element;

    //H5T_FLOAT
    dims_out[0] = spectra_volume->samples_size();
    dims_out[1] = spectra_volume->rows();
    dims_out[2] = spectra_volume->cols();
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    count[0] = dims_out[0];
    count[1] = 1;
    count[2] = 1;
    chunk_dims[0] = dims_out[0];
    chunk_dims[1] = 1;
    chunk_dims[2] = 1;


    dims_time_out[0] = spectra_volume->rows();
    dims_time_out[1] = spectra_volume->cols();
    offset_time[0] = 0;
    offset_time[1] = 0;
    count_time[0] = 1;
    count_time[1] = 1;


    memoryspace_id = H5Screate_simple(3, count, NULL);
    filespace_id = H5Screate_simple(3, dims_out, NULL);

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    dataspace_id = H5Screate_simple (3, dims_out, NULL);

    memoryspace_time_id = H5Screate_simple(2, count_time, NULL);
    file_time_id = H5Screate_simple(2, dims_time_out, NULL);
    dataspace_time_id = H5Screate_simple (2, dims_time_out, NULL);

    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

    maps_grp_id = H5Gopen(file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        std::cout<<"Error creating group MAPS"<<std::endl;
        return false;
    }

    spec_grp_id = H5Gopen(maps_grp_id, "Spectra", H5P_DEFAULT);
    if(spec_grp_id < 0)
        spec_grp_id = H5Gcreate(maps_grp_id, "Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(spec_grp_id < 0)
    {
        std::cout<<"Error creating group MAPS/Spectra"<<std::endl;
        return false;
    }
/*
    scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
    if(scalers_grp_id < 0)
        scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(scalers_grp_id < 0)
    {
        std::cout<<"Error creating group MAPS/Scalers"<<std::endl;
        return false;
    }
*/
    dset_id = H5Dcreate (spec_grp_id, path.c_str(), H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

    dset_rt_id = H5Dcreate (spec_grp_id, "Elapsed_Realtime", H5T_INTEL_F64, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dset_lt_id = H5Dcreate (spec_grp_id, "Elapsed_Livetime", H5T_INTEL_F64, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    incnt_dset_id = H5Dcreate (spec_grp_id, "Output_Counts", H5T_INTEL_F64, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    outcnt_dset_id = H5Dcreate (spec_grp_id, "Input_Counts", H5T_INTEL_F64, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Sselect_hyperslab (memoryspace_time_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);

    real_t real_time;
    real_t life_time;
    real_t in_cnt;
    real_t out_cnt;
    for(size_t row=row_idx_start; row < row_idx_end; row++)
    {
        offset[1] = row;
        offset_time[1] = row;
        for(size_t col=col_idx_start; col < col_idx_end; col++)
        {
            const data_struct::xrf::Spectra *spectra = &((*spectra_volume)[row][col]);
            offset[2] = col;
            offset_time[0] = col;
            H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

            status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);

            H5Sselect_hyperslab (file_time_id, H5S_SELECT_SET, offset_time, NULL, count_time, NULL);

            real_time = spectra->elapsed_realtime();
            life_time = spectra->elapsed_lifetime();
            in_cnt = spectra->input_counts();
            out_cnt = spectra->output_counts();
            status = H5Dwrite (dset_rt_id, H5T_NATIVE_DOUBLE, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&real_time);
            status = H5Dwrite (dset_lt_id, H5T_NATIVE_DOUBLE, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&life_time);
            status = H5Dwrite (incnt_dset_id, H5T_NATIVE_DOUBLE, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&in_cnt);
            status = H5Dwrite (outcnt_dset_id, H5T_NATIVE_DOUBLE, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&out_cnt);
        }
    }


    H5Dclose(dset_id);
    H5Dclose(dset_rt_id);
    H5Dclose(dset_lt_id);
    H5Dclose(incnt_dset_id);
    H5Dclose(outcnt_dset_id);
    H5Sclose(memoryspace_time_id);
    H5Sclose(memoryspace_id);
    H5Sclose(file_time_id);
    H5Sclose(filespace_id);
    H5Sclose(dataspace_time_id);
    H5Sclose(dataspace_id);
    H5Pclose(dcpl_id);


    int_spec_grp_id = H5Gopen(spec_grp_id, "Integrated_Spectra", H5P_DEFAULT);
    if(int_spec_grp_id < 0)
        int_spec_grp_id = H5Gcreate(spec_grp_id, "Integrated_Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //save quantification_standard integrated spectra
    data_struct::xrf::Spectra spectra = spectra_volume->integrate();
    count[0] = spectra.size();
    memoryspace_id = H5Screate_simple(1, count, NULL);
    filespace_id = H5Screate_simple(1, count, NULL);
    dataspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (int_spec_grp_id, "Spectra", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&spectra[0]);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(filespace_id);
    H5Sclose(dataspace_id);


    //save real_time
    count[0] = 1;
    real_t save_val = spectra.elapsed_realtime();
    dataspace_id = H5Screate_simple (1, count, NULL);
    memoryspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (int_spec_grp_id, "Elapsed_Realtime", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save life_time
    save_val = spectra.elapsed_lifetime();
    dataspace_id = H5Screate_simple (1, count, NULL);
    memoryspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (int_spec_grp_id, "Elapsed_Livetime", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save energy vector
    std::vector<real_t> out_vec;
    data_struct::xrf::gen_energy_vector(spectra.size(), energy_offset, energy_slope, &out_vec);
    count[0] = out_vec.size();
    memoryspace_id = H5Screate_simple(1, count, NULL);
    filespace_id = H5Screate_simple(1, count, NULL);
    dataspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (spec_grp_id, "Energy", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&out_vec[0]);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(filespace_id);
    H5Sclose(dataspace_id);

    //save energy calibration
    count[0] = 3;
    dataspace_id = H5Screate_simple (1, count, NULL);
    filespace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (spec_grp_id, "Energy_Calibration", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    count[0] = 1;
    memoryspace_id = H5Screate_simple(1, count, NULL);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_offset);
    offset[0] = 1;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_slope);
    offset[0] = 2;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_quad);
    H5Dclose(dset_id);
    H5Sclose(filespace_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save file version
    save_val = HDF5_SAVE_VERSION;
    dataspace_id = H5Screate_simple (1, count, NULL);
    dset_id = H5Dcreate (maps_grp_id, "version", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(dataspace_id);

    H5Gclose(int_spec_grp_id);
    H5Gclose(spec_grp_id);
    //H5Gclose(scalers_grp_id);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "\n\n save mca_arr elapsed time: " << elapsed_seconds.count() << "s\n\n\n";

    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_element_fits(const hid_t file_id,
                                std::string path,
                                const data_struct::xrf::Fit_Count_Dict * const element_counts,
                                size_t row_idx_start,
                                int row_idx_end,
                                size_t col_idx_start,
                                int col_idx_end)
{
    std::unique_lock<std::mutex> lock(_mutex);

    std::cout<<"HDF5_IO::save_element_fits()"<<std::endl;

//hid_t error_stack = H5Eget_current_stack();
//H5Eset_auto2(error_stack, NULL, NULL);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::string xrf_grp_name = "XRF_Analyzed";
    hid_t   dset_id, dset_ch_id;
    hid_t   memoryspace, filespace, dataspace_id, dataspace_ch_id, dataspace_ch_off_id;
    hid_t   filetype, memtype, status;
    hid_t   dcpl_id;
    hid_t   xrf_grp_id, fit_grp_id, maps_grp_id;
    herr_t   error;

    //hid_t q_dataspace_id, q_dataspace_ch_id, q_memoryspace_id, q_filespace_id, q_dset_id, q_dset_ch_id, q_grp_id, q_fit_grp_id;

    dset_id = -1;
    dset_ch_id = -1;
    hsize_t dims_out[3];
    hsize_t offset[3];
    hsize_t count[3];
    hsize_t chunk_dims[3];

    //hsize_t q_dims_out[2];

    //get one element
    data_struct::xrf::Fit_Counts_Array element;
    //fix this
    for(const auto& iter : *element_counts)
    {
        element = iter.second;
        break;
    }
    //H5T_FLOAT
    dims_out[0] = element_counts->size();
    dims_out[1] = element.rows();
    dims_out[2] = element.cols();
    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    count[0] = 1;
    count[1] = 1;
    count[2] = dims_out[2];
    chunk_dims[0] = 1;
    chunk_dims[1] = 1;
    chunk_dims[2] = dims_out[2];

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    memoryspace = H5Screate_simple(3, count, NULL);
    filespace = H5Screate_simple(3, dims_out, NULL);

    dataspace_id = H5Screate_simple (3, dims_out, NULL);

    dataspace_ch_id = H5Screate_simple (1, dims_out, NULL);
    dataspace_ch_off_id = H5Screate_simple (1, dims_out, NULL);

    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    maps_grp_id = H5Gopen(file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        std::cout<<"Error creating group 'MAPS'"<<std::endl;
        return false;
    }

    xrf_grp_id = H5Gopen(maps_grp_id, xrf_grp_name.c_str(), H5P_DEFAULT);
    if(xrf_grp_id < 0)
        xrf_grp_id = H5Gcreate(maps_grp_id, xrf_grp_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(xrf_grp_id < 0)
    {
        std::cout<<"Error creating group MAPS/"<<xrf_grp_name<<std::endl;
        return false;
    }

    fit_grp_id = H5Gopen(xrf_grp_id, path.c_str(), H5P_DEFAULT);
    if(fit_grp_id < 0)
        fit_grp_id = H5Gcreate(xrf_grp_id, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(fit_grp_id < 0)
    {
        std::cout<<"Error creating group MAPS/"<<xrf_grp_name<<"/"<<path<<std::endl;
        return false;
    }

    //dset_id = H5Dopen (maps_grp_id, path.c_str(), H5P_DEFAULT);
    if(dset_id < 0)
        dset_id = H5Dcreate (fit_grp_id, "Counts_Per_Sec", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    if(dset_id < 0)
    {
        std::cout<<"Error creating dataset MAPS/"<<xrf_grp_name<<"/"<<path<<"/Counts_Per_Sec"<<std::endl;
        return false;
    }

    filetype = H5Tcopy (H5T_FORTRAN_S1);
    H5Tset_size (filetype, 255);
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, 256);

    dset_ch_id = H5Dopen (fit_grp_id, "Channel_Names", H5P_DEFAULT);
    if(dset_ch_id < 0)
        dset_ch_id = H5Dcreate (fit_grp_id, "Channel_Names", filetype, dataspace_ch_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dset_ch_id < 0)
    {
        std::cout<<"Error creating dataset MAPS/"<<xrf_grp_name<<"/"<<path<<"/Channel_Names"<<std::endl;
        return false;
    }

/*
   if (row_idx_end < row_idx_start || row_idx_end > spectra_volume->rows() -1)
    {
        row_idx_end = spectra_volume->rows();
    }
    if (col_idx_end < col_idx_start || col_idx_end > spectra_volume->cols() -1)
    {
        col_idx_end = spectra_volume->cols();
    }
*/

    //create save ordered vector by element Z number with K , L, M lines
    std::vector<std::string> element_lines;
    for (std::string el_name : data_struct::xrf::Element_Symbols)
    {
        element_lines.push_back(el_name);
    }
    for (std::string el_name : data_struct::xrf::Element_Symbols)
    {
        element_lines.push_back(el_name+"_L");
    }
    for (std::string el_name : data_struct::xrf::Element_Symbols)
    {
        element_lines.push_back(el_name+"_M");
    }

    element_lines.push_back(data_struct::xrf::STR_COHERENT_SCT_AMPLITUDE);
    element_lines.push_back(data_struct::xrf::STR_COMPTON_AMPLITUDE);
    element_lines.push_back(data_struct::xrf::STR_NUM_ITR);


    int i=0;
    //save by element Z order
    //for(const auto& iter : *element_counts)
    for (std::string el_name : element_lines)
    {
        if(element_counts->count(el_name) < 1 )
        {
            continue;
        }
        offset[0] = i;
        element = element_counts->at(el_name);

        H5Sselect_hyperslab (dataspace_ch_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Sselect_hyperslab (dataspace_ch_off_id, H5S_SELECT_SET, &offset[2], NULL, count, NULL);
        status = H5Dwrite (dset_ch_id, memtype, dataspace_ch_off_id, dataspace_ch_id, H5P_DEFAULT, (void*)(el_name.c_str()));

        for(size_t row = 0; row < element.rows(); row++)
        {
            offset[1] = row;
            H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

            status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace, filespace, H5P_DEFAULT, (void*)&(element[row][0]));
        }
        i++;
    }

    H5Dclose(dset_id);
    H5Dclose(dset_ch_id);
    H5Sclose(memoryspace);
    H5Sclose(filespace);
    H5Sclose(dataspace_ch_off_id);
    H5Sclose(dataspace_ch_id);
    //h5Tclose(filetype);
    //h5Tclose(memtype);
    H5Pclose(dcpl_id);    
    H5Sclose(dataspace_id);
    H5Gclose(fit_grp_id);
    H5Gclose(xrf_grp_id);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "\n\n save channels elapsed time: " << elapsed_seconds.count() << "s\n\n\n";



    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_quantification(const hid_t file_id,
                                  data_struct::xrf::Quantification_Standard * quantification_standard,
                                  size_t row_idx_start,
                                  int row_idx_end,
                                  size_t col_idx_start,
                                  int col_idx_end)
{
    std::unique_lock<std::mutex> lock(_mutex);

    std::cout<<"HDF5_IO::save_quantification()"<<std::endl;

//hid_t error_stack = H5Eget_current_stack();
//H5Eset_auto2(error_stack, NULL, NULL);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    dset_id, memoryspace_id, filespace_id, dataspace_id, filetype, dataspace_ch_id, memtype,  dset_ch_id, q_int_spec_grp_id;
//    hid_t   xrf_grp_id, fit_grp_id;
    herr_t   error;

    hid_t q_dataspace_id, q_memoryspace_id, q_filespace_id, q_dset_id, q_grp_id, q_fit_grp_id, maps_grp_id, status, scalers_grp_id, xrf_fits_grp_id;

//    dset_id = -1;
//    dset_ch_id = -1;
//    hsize_t dims_out[3];
    hsize_t offset[3];
    hsize_t count[3];
//    hsize_t chunk_dims[3];

    hsize_t q_dims_out[2];


    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    count[0] = 1;
    count[1] = 0;
    count[2] = 0;
//    chunk_dims[0] = 1;
//    chunk_dims[1] = 1;
//    chunk_dims[2] = dims_out[2];

//    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
//    H5Pset_chunk(dcpl_id, 3, chunk_dims);
//    H5Pset_deflate (dcpl_id, 7);

//    memoryspace = H5Screate_simple(3, count, NULL);
//    filespace = H5Screate_simple(3, dims_out, NULL);

//    dataspace_id = H5Screate_simple (3, dims_out, NULL);

//    dataspace_ch_id = H5Screate_simple (1, dims_out, NULL);


//    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    maps_grp_id = H5Gopen(file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        std::cout<<"Error creating group 'MAPS'"<<std::endl;
        return false;
    }

    filetype = H5Tcopy (H5T_FORTRAN_S1);
    H5Tset_size (filetype, 255);
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, 256);



    //--                        save calibration curve                  --

    if( quantification_standard != nullptr && quantification_standard->processed() == true)
    {

        q_grp_id = H5Gopen(maps_grp_id, "Quantification", H5P_DEFAULT);
        if(q_grp_id < 0)
            q_grp_id = H5Gcreate(maps_grp_id, "Quantification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(q_grp_id < 0)
        {
            std::cout<<"Error creating group MAPS/Quantification"<<std::endl;
            return false;
        }

        scalers_grp_id = H5Gopen(q_grp_id, "Scalers", H5P_DEFAULT);
        if(scalers_grp_id < 0)
            scalers_grp_id = H5Gcreate(q_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(scalers_grp_id < 0)
        {
            std::cout<<"Error creating group MAPS/Quantification/Scalers"<<std::endl;
            return false;
        }

        xrf_fits_grp_id = H5Gopen(q_grp_id, "XRF_Analyzed", H5P_DEFAULT);
        if(xrf_fits_grp_id < 0)
            xrf_fits_grp_id = H5Gcreate(q_grp_id, "XRF_Analyzed", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(xrf_fits_grp_id < 0)
        {
            std::cout<<"Error creating group MAPS/Quantification/XRF_Analyzed"<<std::endl;
            return false;
        }

        //save quantification_standard element weights
        const std::unordered_map<std::string, data_struct::xrf::Element_Quant> e_weights = quantification_standard->element_weights();
        count[0] = e_weights.size();
        memoryspace_id = H5Screate_simple(1, count, NULL);
        filespace_id = H5Screate_simple(1, count, NULL);
        dataspace_id = H5Screate_simple (1, count, NULL);
        dataspace_ch_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_grp_id, "Element_Weights", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        dset_ch_id = H5Dcreate (q_grp_id, "Element_Weights_Names", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        offset[0] = 0;
        count[0] = 1;
        H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        int offset_idx = 0;
        for(auto itr: e_weights)
        {
            offset[0] = offset_idx;
            H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

            status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(itr.second.weight));
            status = H5Dwrite (dset_ch_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)(itr.first.c_str()));

            offset_idx++;
        }
        H5Dclose(dset_id);
        H5Dclose(dset_ch_id);
        H5Sclose(memoryspace_id);
        H5Sclose(filespace_id);
        H5Sclose(dataspace_id);
        H5Sclose(dataspace_ch_id);


        q_int_spec_grp_id = H5Gopen(q_grp_id, "Integrated_Spectra", H5P_DEFAULT);
        if(q_int_spec_grp_id < 0)
           q_int_spec_grp_id  = H5Gcreate(q_grp_id, "Integrated_Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(q_int_spec_grp_id < 0)
        {
            std::cout<<"Error creating group MAPS/Quantification/Integrated_Spectra"<<std::endl;
            return false;
        }

        //save quantification_standard integrated spectra
        data_struct::xrf::Spectra spectra = quantification_standard->integrated_spectra();
        count[0] = spectra.size();
        memoryspace_id = H5Screate_simple(1, count, NULL);
        filespace_id = H5Screate_simple(1, count, NULL);
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Spectra", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        offset[0] = 0;
        H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&spectra[0]);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(filespace_id);
        H5Sclose(dataspace_id);


        //save standard name
        count[0] = 1;
        dataspace_id = H5Screate_simple (1, count, NULL);
        memoryspace_id = H5Screate_simple (1, count, NULL);
        dset_ch_id = H5Dcreate (q_grp_id, "Standard_Name", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_ch_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)(quantification_standard->standard_filename().c_str()));
        H5Dclose(dset_ch_id);
        H5Sclose(dataspace_id);

        //save sr_current
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (scalers_grp_id, "SR_Current", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->sr_current()));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save us_ic
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (scalers_grp_id, "US_IC", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->US_IC()));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save ds_ic
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (scalers_grp_id, "DS_IC", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->DS_IC()));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);


        //save real_time
        real_t save_val = spectra.elapsed_realtime();
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Elapsed_Realtime", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save life_time
        save_val = spectra.elapsed_lifetime();
        dataspace_id = H5Screate_simple (1, count, NULL);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Elapsed_Livetime", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(dataspace_id);


        //save calibration curves
        //data_struct::xrf::Quantifiers quantifiers = quantification_standard->calibration_curves.at(path);
        //auto shell_itr = quantification_standard->_calibration_curves.begin();

        q_dims_out[0] = 3;// shells K, L, and M
        q_dims_out[1] = 93; // elements 1 - 92

        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;

        count[0] = 1;
        count[1] = q_dims_out[1];
        count[2] = 0;

        q_memoryspace_id = H5Screate_simple(2, count, NULL);
        q_filespace_id = H5Screate_simple(2, q_dims_out, NULL);

        q_dataspace_id = H5Screate_simple (2, q_dims_out, NULL);
        //q_dataspace_ch_id = H5Screate_simple (1, q_dims_out, NULL);


        for (auto& qitr: quantification_standard->calibration_curves)
        {

            q_fit_grp_id = H5Gopen(xrf_fits_grp_id, qitr.first.c_str(), H5P_DEFAULT);
            if(q_fit_grp_id < 0)
                q_fit_grp_id = H5Gcreate(xrf_fits_grp_id, qitr.first.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if(q_fit_grp_id < 0)
            {
                std::cout<<"Error creating group MAPS/Quantification/"<<qitr.first<<std::endl;
                return false;
            }

            //save quantification_standard counts
            const std::unordered_map<std::string, std::unordered_map<std::string, real_t> > all_element_counts = quantification_standard->element_counts();
            std::unordered_map<std::string, real_t> element_counts = all_element_counts.at(qitr.first);
            count[0] = element_counts.size();
            dataspace_id = H5Screate_simple (1, count, NULL);
            memoryspace_id = H5Screate_simple(1, count, NULL);
            filespace_id = H5Screate_simple(1, count, NULL);
            dset_id = H5Dcreate (q_fit_grp_id, "Counts_Per_Sec", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_ch_id = H5Dcreate (q_fit_grp_id, "Channel_Names", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            //create save ordered vector by element Z number with K , L, M lines
            std::vector<std::string> element_lines;
            for (std::string el_name : data_struct::xrf::Element_Symbols)
            {
                element_lines.push_back(el_name);
            }
            for (std::string el_name : data_struct::xrf::Element_Symbols)
            {
                element_lines.push_back(el_name+"_L");
            }
            for (std::string el_name : data_struct::xrf::Element_Symbols)
            {
                element_lines.push_back(el_name+"_M");
            }

            element_lines.push_back(data_struct::xrf::STR_NUM_ITR);

            offset_idx = 0;
            offset[0] = 0;
            count[0] = 1;
            H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            //save by element Z order
            for (std::string el_name : element_lines)
            {
                if(element_counts.count(el_name) < 1 )
                {
                    continue;
                }
                offset[0] = offset_idx;
                real_t val = element_counts.at(el_name);

                H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                status = H5Dwrite (dset_ch_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)(el_name.c_str()));
                status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&val);

                offset_idx++;
            }
            H5Dclose(dset_id);
            H5Dclose(dset_ch_id);
            H5Sclose(memoryspace_id);
            H5Sclose(filespace_id);
            H5Sclose(dataspace_id);


            for (size_t i =0; i< qitr.second.calib_curves.size(); i++)
            {

                std::string q_dset_name = "Calibration_Curve_" + qitr.second.calib_curves[i].quantifier_name;

                //dset_id = H5Dopen (fit_grp_id, "Calibration_Curve", H5P_DEFAULT);
                //if(q_dset_id < 0)
                    q_dset_id = H5Dcreate (q_fit_grp_id, q_dset_name.c_str(), H5T_INTEL_F64, q_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if(q_dset_id < 0)
                {
                    std::cout<<"Error creating dataset MAPS/Quantification/"<<q_dset_name<<"/"<<q_dset_name<<std::endl;
                    continue;
                }


                //for(auto& shell_itr : quant_itr.second)
                for(size_t j=0; j<3; j++)
                {
                    //int element_offset = 0;
                    //create dataset for different shell curves
                    std::vector<real_t> calibration_curve = qitr.second.calib_curves[i].shell_curves[j];

                    offset[0] = j;
                    offset[1] = 0;

                    //H5Sselect_hyperslab (q_dataspace_ch_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    //H5Sselect_hyperslab (q_dataspace_ch_off_id, H5S_SELECT_SET, &offset[2], NULL, count, NULL);
                    //status = H5Dwrite (q_dset_ch_id, memtype, q_dataspace_ch_off_id, q_dataspace_ch_id, H5P_DEFAULT, (void*)(el_name.c_str()));

                    H5Sselect_hyperslab (q_filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);

                    status = H5Dwrite (q_dset_id, H5T_NATIVE_DOUBLE, q_memoryspace_id, q_filespace_id, H5P_DEFAULT, (void*)&calibration_curve[0]);
                }
                H5Dclose(q_dset_id);
            }
        }

        //H5Dclose(q_dset_ch_id);
        H5Sclose(q_dataspace_id);
        //H5Sclose(q_dataspace_ch_id);
        H5Gclose(xrf_fits_grp_id);
        H5Gclose(scalers_grp_id);
        H5Gclose(q_int_spec_grp_id);
        H5Gclose(q_fit_grp_id);
        H5Gclose(q_grp_id);

    }

//    H5Dclose(dset_id);
//    H5Dclose(dset_ch_id);
//    H5Sclose(memoryspace);
//    H5Sclose(filespace);
//    H5Sclose(dataspace_ch_off_id);
//    H5Sclose(dataspace_ch_id);
//    H5Pclose(dcpl_id);
//    H5Sclose(dataspace_id);
//    H5Gclose(xrf_grp_id);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "\n\n save channels elapsed time: " << elapsed_seconds.count() << "s\n\n\n";



    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::_save_scan_meta_data(hid_t scan_grp_id, struct mda_file *mda_scalers)
{
    hid_t dataspace_id = -1, memoryspace_id = -1, filespace_id = -1;
	hid_t status;
	hid_t filetype, memtype;
    hid_t dset_id = -1;
	
	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };

	bool single_row_scan = false;

    std::cout<<"HDF5_IO::_save_scan_meta_data()"<<std::endl;

	try
	{

		filetype = H5Tcopy(H5T_FORTRAN_S1);
		H5Tset_size(filetype, 255);
		memtype = H5Tcopy(H5T_C_S1);
		status = H5Tset_size(memtype, 256);


		memoryspace_id = H5Screate_simple(1, count, NULL);

		//save scan positions
		if (mda_scalers->scan->scan_rank > 1)
		{
			if (mda_scalers->header->data_rank == 2)
			{
				if (mda_scalers->header->dimensions[1] == 2000)
				{
					single_row_scan = true;
				}
			}

			if (single_row_scan)
			{
				count[0] = 1;
				dataspace_id = H5Screate_simple(1, count, NULL);
				dset_id = H5Dcreate(scan_grp_id, "y_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
					std::cout << "Error creating dataset 'y_axis'" << std::endl;
					return false;
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);

				count[0] = mda_scalers->scan->last_point;
				dataspace_id = H5Screate_simple(1, count, NULL);
				filespace_id = H5Screate_simple(1, count, NULL);
				dset_id = H5Dcreate(scan_grp_id, "x_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
					std::cout << "Error creating dataset 'x_axis'" << std::endl;
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->positioners_data[0][i];
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
					status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);
			}
			else
			{
				count[0] = mda_scalers->scan->last_point;

				dataspace_id = H5Screate_simple(1, count, NULL);
				filespace_id = H5Screate_simple(1, count, NULL);
				dset_id = H5Dcreate(scan_grp_id, "y_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
					std::cout << "Error creating dataset 'y_axis'" << std::endl;
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->positioners_data[0][i];
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
					status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);

				count[0] = mda_scalers->scan->sub_scans[0]->last_point;
				dataspace_id = H5Screate_simple(1, count, NULL);
				filespace_id = H5Screate_simple(1, count, NULL);
				dset_id = H5Dcreate(scan_grp_id, "x_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
					std::cout << "Error creating dataset 'x_axis'" << std::endl;
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->sub_scans[0]->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->sub_scans[0]->positioners_data[0][i];
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
					status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);
			}

		}

		//save write date
		count[0] = 1;

		dset_id = H5Dcreate(scan_grp_id, "scan_time_stamp", filetype, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
			std::cout << "Error creating dataset 'scan_time_stamp'" << std::endl;
			return false;
		}
		if (mda_scalers->scan->time != nullptr)
		{
			std::string str_time = std::string(mda_scalers->scan->time);
			status = H5Dwrite(dset_id, memtype, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)str_time.c_str());
		}
        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
		dset_id = H5Dcreate(scan_grp_id, "name", filetype, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
			std::cout << "Error creating dataset 'name'" << std::endl;
			return false;
		}
		if(mda_scalers->scan->name != nullptr)
		{
			std::string str_name = std::string(mda_scalers->scan->name);
			status = H5Dwrite(dset_id, memtype, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)str_name.c_str());
		}

        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
	}
	catch (...)
	{
		std::cout << "Error saving MAPS/Scan meta data" << std::endl;
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_save_extras(hid_t scan_grp_id, struct mda_file *mda_scalers)
{

    hid_t filespace_id = -1, filespace_id2 = -1, filespace_id3 = -1, filespace_id4 = -1;
    hid_t memoryspace_id = -1, memoryspace_id2 = -1, memoryspace_id3 = -1, memoryspace_id4 = -1;
    hid_t status = -1;
	hid_t filetype, memtype;
    hid_t dset_desc_id = -1, dset_unit_id = -1, dset_id = -1, dset_val_id = -1;
    hid_t extra_grp_id = -1;
	
	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };

    std::cout<<"HDF5_IO::_save_extras()"<<std::endl;

	if (mda_scalers->extra == nullptr)
	{
		return false;
	}

	try
	{

        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 512);
		memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 512);

        //std::cout<<"HDF5_IO::_save_extras() - extra_pvs group"<<std::endl;
		extra_grp_id = H5Gopen(scan_grp_id, "Extra_PVs", H5P_DEFAULT);
		if (extra_grp_id < 0)
			extra_grp_id = H5Gcreate(scan_grp_id, "Extra_PVs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (extra_grp_id < 0)
		{
			std::cout << "Error creating group MAPS/Scan/Extra_PVs" << std::endl;
			return false;
		}

		//save extra pv's
		count[0] = (size_t)mda_scalers->extra->number_pvs;
		filespace_id = H5Screate_simple(1, count, NULL);
        //std::cout<<"HDF5_IO::_save_extras() - names dataset"<<std::endl;
        dset_id = H5Dcreate(extra_grp_id, "Names", filetype, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
            std::cout << "Error creating dataset MAPS/Scan/Extra_PVs/Names" << std::endl;
            if(filespace_id > -1)
            {
                H5Sclose(filespace_id);
                filespace_id = -1;
            }
            if(extra_grp_id > -1)
            {
                H5Gclose(extra_grp_id);
                extra_grp_id = -1;
            }
            return false;
		}

        filespace_id2 = H5Screate_simple(1, count, NULL);
        //std::cout<<"HDF5_IO::_save_extras() - values dataset"<<std::endl;
        dset_val_id = H5Dcreate(extra_grp_id, "Values", filetype, filespace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_val_id < 0)
		{
            std::cout << "Error creating dataset MAPS/Scan/Extra_PVs/Values" << std::endl;
            if(filespace_id > -1)
            {
                H5Sclose(filespace_id);
                filespace_id = -1;
            }
            if(filespace_id2 > -1)
            {
                H5Sclose(filespace_id2);
                filespace_id2 = -1;
            }
            if(extra_grp_id > -1)
            {
                H5Gclose(extra_grp_id);
                extra_grp_id = -1;
            }
            return false;
		}

        filespace_id3 = H5Screate_simple(1, count, NULL);
        //std::cout<<"HDF5_IO::_save_extras() - desc dataset"<<std::endl;
        dset_desc_id = H5Dcreate(extra_grp_id, "Description", filetype, filespace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_desc_id < 0)
		{
            std::cout << "Error creating dataset MAPS/Scan/Extra_PVs/Description" << std::endl;
            if(filespace_id > -1)
            {
                H5Sclose(filespace_id);
                filespace_id = -1;
            }
            if(filespace_id2 > -1)
            {
                H5Sclose(filespace_id2);
                filespace_id2 = -1;
            }
            if(filespace_id3 > -1)
            {
                H5Sclose(filespace_id3);
                filespace_id3 = -1;
            }
            if(extra_grp_id > -1)
            {
                H5Gclose(extra_grp_id);
                extra_grp_id = -1;
            }
            return false;
		}

        filespace_id4 = H5Screate_simple(1, count, NULL);
        //std::cout<<"HDF5_IO::_save_extras() - unit dataset"<<std::endl;
        dset_unit_id = H5Dcreate(extra_grp_id, "Unit", filetype, filespace_id4, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_unit_id < 0)
		{
            std::cout << "Error creating dataset MAPS/Scan/Extra_PVs/Unit" << std::endl;
            if(filespace_id > -1)
            {
                H5Sclose(filespace_id);
                filespace_id = -1;
            }
            if(filespace_id2 > -1)
            {
                H5Sclose(filespace_id2);
                filespace_id2 = -1;
            }
            if(filespace_id3 > -1)
            {
                H5Sclose(filespace_id3);
                filespace_id3 = -1;
            }
            if(filespace_id4 > -1)
            {
                H5Sclose(filespace_id4);
                filespace_id4 = -1;
            }
            if(extra_grp_id > -1)
            {
                H5Gclose(extra_grp_id);
                extra_grp_id = -1;
            }
            return false;
		}

        count[0] = 1;

		std::string str_val;
		short* s_val;
		int* i_val;
		float* f_val;
		double* d_val;

        //std::cout<<"HDF5_IO::_save_extras() - memspaces"<<std::endl;
		count[0] = 1;
		memoryspace_id = H5Screate_simple(1, count, NULL);
		H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
		memoryspace_id2 = H5Screate_simple(1, count, NULL);
		H5Sselect_hyperslab(memoryspace_id2, H5S_SELECT_SET, offset, NULL, count, NULL);
		memoryspace_id3 = H5Screate_simple(1, count, NULL);
		H5Sselect_hyperslab(memoryspace_id3, H5S_SELECT_SET, offset, NULL, count, NULL);
		memoryspace_id4 = H5Screate_simple(1, count, NULL);
		H5Sselect_hyperslab(memoryspace_id4, H5S_SELECT_SET, offset, NULL, count, NULL);

		for (int16_t i = 0; i < mda_scalers->extra->number_pvs; i++)
		{
			offset[0] = i;
			struct mda_pv * pv = mda_scalers->extra->pvs[i];
            if(pv == nullptr)
            {
                continue;
            }
			switch (pv->type)
			{

			case EXTRA_PV_STRING:
				str_val = std::string(pv->values);
				break;
				//case EXTRA_PV_INT8:

				//    break;
			case EXTRA_PV_INT16:
				s_val = (short*)pv->values;
				str_val = std::to_string(*s_val);
				break;
			case EXTRA_PV_INT32:
				i_val = (int*)pv->values;
				str_val = std::to_string(*i_val);
				break;
			case EXTRA_PV_FLOAT:
				f_val = (float*)pv->values;
				str_val = std::to_string(*f_val);
				break;
			case EXTRA_PV_DOUBLE:
				d_val = (double*)pv->values;
				str_val = std::to_string(*d_val);
				break;

			}

			H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id2, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id3, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id4, H5S_SELECT_SET, offset, NULL, count, NULL);

			if (pv->name != nullptr)
			{
				//need this becuase it crashes on windows if we just pass the char *
				std::string name = std::string(pv->name);
                //std::cout<<"HDF5_IO::_save_extras() - save name : "<<name<<" "<<name.length()<<std::endl;
                //BUG: causes crash randomly when writing name to hdf. not sure why yet
                name = std::string("a");
                status = H5Dwrite(dset_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)name.c_str());
			}

            //std::cout<<"HDF5_IO::_save_extras() - save val : "<<str_val<<" "<<str_val.length()<<std::endl;
            status = H5Dwrite(dset_val_id, memtype, memoryspace_id2, filespace_id2, H5P_DEFAULT, (void*)str_val.c_str());
			if (pv->description != nullptr)
			{
                std::string description = std::string(pv->description);
                //std::cout<<"HDF5_IO::_save_extras() - save desc : "<<description<<" "<<description.length()<<std::endl;
                status = H5Dwrite(dset_desc_id, memtype, memoryspace_id3, filespace_id3, H5P_DEFAULT, (void*)description.c_str());
			}

			if (pv->unit != nullptr)
			{
                std::string unit = std::string(pv->unit);
                //std::cout<<"HDF5_IO::_save_extras() - save unit : "<<unit<<" "<<unit.length()<<std::endl;
                status = H5Dwrite(dset_unit_id, memtype, memoryspace_id4, filespace_id4, H5P_DEFAULT, (void*)unit.c_str());
			}
		}

        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
        if(dset_val_id > -1)
        {
            H5Dclose(dset_val_id);
            dset_val_id = -1;
        }
        if(dset_desc_id > -1)
        {
            H5Dclose(dset_desc_id);
            dset_desc_id = -1;
        }
        if(dset_unit_id > -1)
        {
            H5Dclose(dset_unit_id);
            dset_unit_id = -1;
        }
        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
        if(memoryspace_id2 > -1)
        {
            H5Sclose(memoryspace_id2);
            memoryspace_id2 = -1;
        }
        if(memoryspace_id3 > -1)
        {
            H5Sclose(memoryspace_id3);
            memoryspace_id3 = -1;
        }
        if(memoryspace_id4 > -1)
        {
            H5Sclose(memoryspace_id4);
            memoryspace_id4 = -1;
        }
        if(filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if(filespace_id2 > -1)
        {
            H5Sclose(filespace_id2);
            filespace_id2 = -1;
        }
        if(filespace_id3 > -1)
        {
            H5Sclose(filespace_id3);
            filespace_id3 = -1;
        }
        if(filespace_id4 > -1)
        {
            H5Sclose(filespace_id4);
            filespace_id4 = -1;
        }
        if(extra_grp_id > -1)
        {
            H5Gclose(extra_grp_id);
            extra_grp_id = -1;
        }

	}
	catch (...)
	{
        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
        if(dset_val_id > -1)
        {
            H5Dclose(dset_val_id);
            dset_val_id = -1;
        }
        if(dset_desc_id > -1)
        {
            H5Dclose(dset_desc_id);
            dset_desc_id = -1;
        }
        if(dset_unit_id > -1)
        {
            H5Dclose(dset_unit_id);
            dset_unit_id = -1;
        }
        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
        if(memoryspace_id2 > -1)
        {
            H5Sclose(memoryspace_id2);
            memoryspace_id2 = -1;
        }
        if(memoryspace_id3 > -1)
        {
            H5Sclose(memoryspace_id3);
            memoryspace_id3 = -1;
        }
        if(memoryspace_id4 > -1)
        {
            H5Sclose(memoryspace_id4);
            memoryspace_id4 = -1;
        }
        if(filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if(filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if(filespace_id2 > -1)
        {
            H5Sclose(filespace_id2);
            filespace_id2 = -1;
        }
        if(filespace_id3 > -1)
        {
            H5Sclose(filespace_id3);
            filespace_id3 = -1;
        }
        if(filespace_id4 > -1)
        {
            H5Sclose(filespace_id4);
            filespace_id4 = -1;
        }
        if(extra_grp_id > -1)
        {
            H5Gclose(extra_grp_id);
            extra_grp_id = -1;
        }
		std::cout << "Error creating group MAPS/Scan/Extra_PVs" << std::endl;
		return false;
	}
    //std::cout<<"HDF5_IO::_save_extras() - done"<<std::endl;
	return true;
}

bool HDF5_IO::_save_scalers(hid_t maps_grp_id, struct mda_file *mda_scalers, size_t detector_num, std::unordered_map< std::string, std::string > *extra_override_values)
{
    hid_t dataspace_id = -1, memoryspace_id = -1, filespace_id = -1, filespace_name_id = -1;
    hid_t filetype, memtype;
    hid_t dset_names_id = -1;
    hid_t scalers_grp_id = -1;
    hid_t dcpl_id = -1, status;

    hsize_t offset[1] = { 0 };
    hsize_t count[1] = { 1 };

    hsize_t offset_3d[3] = { 0, 0, 0 };
    hsize_t count_3d[3] = { 1, 1, 1 };

    MDA_IO mda_io;

    bool single_row_scan = false;

    std::cout<<"HDF5_IO::_save_scalers()"<<std::endl;

    //don't save these scalers
    std::list<std::string> ignore_scaler_strings = { "ELT1", "ERT1", "ICR1", "OCR1" };
    std::list<struct scaler_struct> scalers;

    hid_t dset_cps_id;


    filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 255);
    memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(memtype, 256);


    try
    {
        count[0] = 1;
        memoryspace_id = H5Screate_simple(1, count, NULL);

        real_t val;
        bool save_cfg_abs = false;
        if (extra_override_values != nullptr && mda_scalers->scan->scan_rank > 1)
        {
            int us_ic_idx = -1;
            int ds_ic_idx = -1;
            int cfg_2_idx = -1;
            int cfg_3_idx = -1;
            int cfg_4_idx = -1;
            int cfg_5_idx = -1;

            int hdf_idx = 0;
            for (auto itr : *extra_override_values)
            {
                //don't save ELT1, ERT1, ICR1, OCR1. these are saved elsewhere
                std::list<std::string>::iterator s_itr = std::find(ignore_scaler_strings.begin(), ignore_scaler_strings.end(), itr.first);
                if (s_itr != ignore_scaler_strings.end())
                    continue;

                int mda_idx = mda_io.find_2d_detector_index(mda_scalers, itr.second, detector_num, val);
                if (mda_idx > -1)
                {
                    scalers.push_back(scaler_struct(itr.first, mda_idx, hdf_idx));
                    hdf_idx++;
                    if (itr.first == "US_IC")
                        us_ic_idx = mda_idx;
                    else if (itr.first == "DS_IC")
                        ds_ic_idx = mda_idx;
                    else if (itr.first == "CFG_2")
                        cfg_2_idx = mda_idx;
                    else if (itr.first == "CFG_3")
                        cfg_3_idx = mda_idx;
                    else if (itr.first == "CFG_4")
                        cfg_4_idx = mda_idx;
                    else if (itr.first == "CFG_5")
                        cfg_5_idx = mda_idx;
                }
            }

            if (scalers.size() > 0)
            {

                scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
                if (scalers_grp_id < 0)
                    scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (scalers_grp_id < 0)
                {
                    std::cout << "Error creating group MAPS/Scalers" << std::endl;
                    return false;
                }

                if (mda_scalers->header->data_rank == 2)
                {
                    if (mda_scalers->header->dimensions[1] == 2000)
                    {
                        count_3d[0] = 1;
                        count_3d[1] = 1;
                        count_3d[2] = mda_scalers->scan->last_point;
                        single_row_scan = true;
                    }
                }
                else if (mda_scalers->header->data_rank == 3)
                {
                    count_3d[0] = 1;
                    count_3d[1] = mda_scalers->scan->last_point;
                    count_3d[2] = mda_scalers->scan->sub_scans[0]->last_point;
                }
                else
                {
                    std::cout << "Unsupported rank " << mda_scalers->header->data_rank << " . Skipping scalers" << std::endl;
                    //TODO: return / throw exception
                }
                dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(dcpl_id, 3, count_3d);
                H5Pset_deflate(dcpl_id, 7);

                if (us_ic_idx > -1 && ds_ic_idx > -1 && cfg_2_idx > -1 && cfg_2_idx > -1 && cfg_2_idx > -1 && cfg_2_idx > -1)
                {
                    count_3d[0] = scalers.size() + 6; //abs_ic, abs_cfg, H_dpc_cfg, V_dpc_cfg, dia1_dpc_cfg, dia2_dpc_cfg
                    save_cfg_abs = true;
                }
                else
                {
                    count_3d[0] = scalers.size();
                }
                dataspace_id = H5Screate_simple(3, count_3d, NULL);
                filespace_id = H5Screate_simple(3, count_3d, NULL);

                count[0] = count_3d[0];
                filespace_name_id = H5Screate_simple(1, count, NULL);

                dset_cps_id = H5Dcreate(scalers_grp_id, "Values", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
                dset_names_id = H5Dcreate(scalers_grp_id, "Names", filetype, filespace_name_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                count_3d[0] = 1;
                count_3d[1] = 1;
                count_3d[2] = 1;
                count[0] = 1;
                //TODO: create Names, and Units to go with the Values
                //save scalers



                for (auto &itr : scalers)
                {
                    if (single_row_scan)
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            double val = mda_scalers->scan->detectors_data[itr.mda_idx][i];
                            offset_3d[0] = itr.hdf_idx;
                            offset_3d[1] = 1;
                            offset_3d[2] = i;
                            H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                            status = H5Dwrite(dset_cps_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&val);
                        }
                    }
                    else
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            for (int32_t j = 0; j < mda_scalers->scan->sub_scans[0]->last_point; j++)
                            {
                                double val = mda_scalers->scan->sub_scans[i]->detectors_data[itr.mda_idx][j];
                                offset_3d[0] = itr.hdf_idx;
                                offset_3d[1] = i;
                                offset_3d[2] = j;
                                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                                status = H5Dwrite(dset_cps_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&val);
                            }
                        }
                    }
                    offset[0] = itr.hdf_idx;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_id, filespace_name_id, H5P_DEFAULT, (void*)itr.hdf_name.c_str());
                }

                if (save_cfg_abs)
                {
                    std::string tmp_name = "abs_ic";
                    offset[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_name.c_str());

                    tmp_name = "abs_cfg";
                    offset[0] = hdf_idx + 1;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_name.c_str());

                    tmp_name = "H_dpc_cfg";
                    offset[0] = hdf_idx + 2;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_name.c_str());

                    tmp_name = "V_dpc_cfg";
                    offset[0] = hdf_idx + 3;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_name.c_str());

                    tmp_name = "dia1_dpc_cfg";
                    offset[0] = hdf_idx + 4;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_name.c_str());

                    tmp_name = "dia2_dpc_cfg";
                    offset[0] = hdf_idx + 5;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_name.c_str());

                    if (single_row_scan)
                    {
                        //TODO: finish single row scan scalers
                    }
                    else
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            for (int32_t j = 0; j < mda_scalers->scan->sub_scans[0]->last_point; j++)
                            {
                                double us_ic = mda_scalers->scan->sub_scans[i]->detectors_data[us_ic_idx][j];
                                double ds_ic = mda_scalers->scan->sub_scans[i]->detectors_data[ds_ic_idx][j];
                                double t_2 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_2_idx][j];
                                double t_3 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_3_idx][j];
                                double t_4 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_4_idx][j];
                                double t_5 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_5_idx][j];


                                offset_3d[1] = i;
                                offset_3d[2] = j;

                                double t_abs = t_2 + t_3 + t_4 + t_5;
                                double abs_ic = ds_ic / us_ic;

                                offset_3d[0] = hdf_idx;
                                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                                status = H5Dwrite(dset_cps_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&abs_ic);

                                double abs_cfg = t_abs / us_ic;

                                offset_3d[0] = hdf_idx + 1;
                                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                                status = H5Dwrite(dset_cps_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&abs_cfg);

                                double H_dpc_cfg = (t_2 - t_3 - t_4 + t_5) / t_abs;

                                offset_3d[0] = hdf_idx + 2;
                                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                                status = H5Dwrite(dset_cps_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&H_dpc_cfg);

                                double V_dpc_cfg = (t_2 + t_3 - t_4 - t_5) / t_abs;

                                offset_3d[0] = hdf_idx + 3;
                                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                                status = H5Dwrite(dset_cps_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&V_dpc_cfg);

                                double dia1_dpc_cfg = (t_2 - t_4) / t_abs;

                                offset_3d[0] = hdf_idx + 4;
                                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                                status = H5Dwrite(dset_cps_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&dia1_dpc_cfg);

                                double dia2_dpc_cfg = (t_3 - t_5) / t_abs;

                                offset_3d[0] = hdf_idx + 5;
                                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                                status = H5Dwrite(dset_cps_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&dia2_dpc_cfg);
                            }
                        }
                    }
                }
            }
        }

        if(dset_names_id > -1)
        {
            H5Dclose(dset_names_id);
            dset_names_id = -1;
        }
        if(dset_cps_id > -1)
        {
            H5Dclose(dset_cps_id);
            dset_cps_id = -1;
        }
        if(dataspace_id > -1)
        {
            H5Sclose(dataspace_id);
            dataspace_id = -1;
        }
        if(filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if(filespace_name_id > -1)
        {
            H5Sclose(filespace_name_id);
            filespace_name_id = -1;
        }
        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
        if(dcpl_id > -1)
        {
            H5Pclose(dcpl_id);
            dcpl_id = -1;
        }
        if(scalers_grp_id > -1)
        {
            H5Gclose(scalers_grp_id);
            scalers_grp_id = -1;
        }
    }
    catch (...)
    {
        if(memoryspace_id > -1)
            H5Sclose(memoryspace_id);
        if(dset_names_id > -1)
            H5Dclose(dset_names_id);
        if(dset_cps_id > -1)
            H5Dclose(dset_cps_id);
        if(dataspace_id > -1)
            H5Sclose(dataspace_id);
        if(filespace_id > -1)
            H5Sclose(filespace_id);
        if(filespace_name_id > -1)
            H5Sclose(filespace_name_id);
        if(dcpl_id > -1)
            H5Pclose(dcpl_id);
        if(scalers_grp_id > -1)
            H5Gclose(scalers_grp_id);

        std::cout << "Error creating group MAPS/Scalers" << std::endl;
        return false;
    }

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_scan_scalers(const std::string filename,
	size_t detector_num,
	struct mda_file *mda_scalers,
	std::unordered_map< std::string, std::string > *extra_override_values,
	size_t row_idx_start,
	int row_idx_end,
	size_t col_idx_start,
	int col_idx_end)
{

    hid_t scan_grp_id, maps_grp_id;

	if (mda_scalers == nullptr)
    {
        std::cout << "Warning: mda_scalers == nullptr. Not returning from save_scan_scalers" << std::endl;
		return false;
    }

	std::cout << "Saving scalers to hdf5: " << filename << std::endl;

	hid_t file_id = start_save_seq(filename);

	maps_grp_id = H5Gopen(file_id, "MAPS", H5P_DEFAULT);
	if (maps_grp_id < 0)
		maps_grp_id = H5Gcreate(file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (maps_grp_id < 0)
	{
		std::cout << "Error creating group 'MAPS'" << std::endl;
        end_save_seq(file_id);
		return false;
	}

	scan_grp_id = H5Gopen(maps_grp_id, "Scan", H5P_DEFAULT);
	if (scan_grp_id < 0)
		scan_grp_id = H5Gcreate(maps_grp_id, "Scan", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (scan_grp_id < 0)
	{
		std::cout << "Error creating group MAPS/Scan" << std::endl;
        end_save_seq(file_id);
		return false;
	}

    _save_scan_meta_data(scan_grp_id, mda_scalers);
	
    _save_extras(scan_grp_id, mda_scalers);
	
    _save_scalers(maps_grp_id, mda_scalers, detector_num, extra_override_values);

    H5Gclose(scan_grp_id);
    H5Gclose(maps_grp_id);

    end_save_seq(file_id);

    return true;
}

//-----------------------------------------------------------------------------

//This should probably be in a different section instead of hdf_io
// -----------------------------------------------------------------------------
//  Applies 2-D Fourier integration to reconstruct dpc images,
//  presents fluorescence maps for select elements if present in data directory
// VARIABLE DECLARATIONS / MEANINGS:
// info				  structure of info about illumination etc
// delta, beta		  optical parameters
// nrml				  right-minus-left normalised to transmission signal
// ntmb				  top-minus-bottom normalised to transmission signal
//
// gxdt				   x component of the gradient of the delta.thickness
// gydt				   y component of the gradient of the delta.thickness
//
// ngxdt				   x component of the gradient of the delta.thickness, normalised
//						  so as to have zero mean
// ngydt				   y component of the gradient of the delta.thickness, normalised
//						  so as to have zero mean
// =============================================
bool simple_dpc_integration(std::valarray<real_t> nrml, std::valarray<real_t> ntmb, bool no_int)//=true)
{

    real_t hc = 0.001239842;  // wavelength-energy relationship, microns / keV
/*
    if nrml.ndim < 1:
        nrml = 0
        ntmb = 0
        rdt = 0
        return false;
        */
/*
    sz = nrml.shape
    nx = sz[0]
    ny = sz[1]

    // "what goes up must come down"
    //	  - can be used to remove beam intensity variations AND
    //	  - removes first order effect of detector misalignment
    //			 (i.e., removes 'constant gradient')

    ylo = 0
    yhi = ny - 1

    // find the vertical lines with the smalles
    // spread, which hopefully are the background

    for i in range(ylo, yhi)
    {
        nrml[:, i] = nrml[:, i] - nrml[:, i].mean(axis=0)
        // added this for the other direction, too.
        ntmb[:, i] = ntmb[:, i] - ntmb[:, i].mean(axis=0)
    }
    // remove first order effect of detector misalignment in vertical
    //			 (i.e., remove 'constant gradient')
    ntmb = ntmb - ntmb.mean(axis=0)
    // added this for the other direction, too.
    nrml = nrml - nrml.mean(axis=0)

    rdt = 0

    if (no_int == False)
    {
        cs_d = 40.0
        zp_d = 160.0
        zp_dr = 50.0 / 1000.
        zp_f = 18.26 * 1000.
        energy = 10.1

        zz = 82.17

        hx = 0.1
        hy = 0.1

        xlin = np.arange(float(nx)) * hx
        ylin = np.arange(float(ny)) * hy
        ylin = ylin[::-1]
        xax = xlin  // (ylin*0+1)
        yax = (xlin * 0 + 1)  // ylin

        // =============================================
        // calculate as gradient of t
        // gxdt, gydt refers to the gradient of the delta.thickness

        // extra factor of 2 comes from use of diameters, not radii...
        ngxdt = (np.math.pi * (zp_d + cs_d)) / (8. * zp_f) * nrml
        ngydt = (np.math.pi * (zp_d + cs_d)) / (8. * zp_f) * ntmb

        // =============================================
        // implement FFT reconstruction

        dpc = ngxdt + 1j * ngydt

        fx = (np.arange(float(nx)) - nx / 2) / ((nx - 1) * hx)
        fy = (np.arange(float(ny)) - ny / 2) / ((ny - 1) * hy)
        fy = fy[::-1]
        fxt = np.outer(fx, (fy * 0. + 1.))
        fy = np.outer((fx * 0. + 1.), fy)
        fx = fxt
        fxt = 0

        xy = 2j * np.math.pi * (fx + 1j * fy)
        xy[(nx - 1) / 2, (ny - 1) / 2] = 1  // to avoid 0/0 error
        xy = np.fft.fftshift(xy)

        Fdpc = np.fft.fft2(dpc)

        Fdpc[0, 0] = 0  // dc level information not available, leads to numerical error

        dt = np.fft.ifft2(Fdpc / xy)

        // is dt the magnitude of the complex value or the real part?
        // note that the real part dominates, and the magnitude
        // loses dynamic range due to abs(-1) = 1, i.e. real part is positive only

        idt = dt.imag
        rdt = dt.real

        temp = np.concatenate((rdt[0, 1:ny - 1].flatten(), rdt[nx - 1, 1:ny - 1].flatten(), \ rdt[0:nx, 0].flatten(), rdt[0:nx, ny - 1].flatten()))

        // set the average of the perimetric values to be zero
        rdt = rdt - np.mean(temp)
    }
    return nrml, ntmb, rdt
                 */
    return true;
}

//-----------------------------------------------------------------------------


} //end namespace file
}// end namespace io
