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

namespace io
{
namespace file
{

HDF5_IO::HDF5_IO() : Base_File_IO()
{

}

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

bool HDF5_IO::load_dataset(std::string path, Base_Dataset *dset)
{
    return false;

}
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
void HDF5_IO::parse_group_info(hid_t h5file, hid_t id)
{

}

void HDF5_IO::parse_dataset_info(hid_t h5file, hid_t id)
{
    /*
    DataSet dataset = h5file.openDataSet(DATASET_NAME);

    H5T_class_t type_class = dataset.getTypeClass();

    DataSpace dataspace = dataset.getSpace();

    int rank = dataspace.getSimpleExtentNdims();

    hsize_t *dims_out = new hsize_t[rank];
    int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
*/
}

void HDF5_IO::parse_attr_info(hid_t h5file, hid_t id)
{

}

void HDF5_IO::load_spectra_volume(std::string path, HDF5_Spectra_Layout layout, data_struct::xrf::Spectra_Volume* spec_vol)
{

   _is_loaded = ERROR_LOADING;
/*
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   //std::cout<<"load_spectra_volume "<< _filename <<" : "<<path<<std::endl;

   hid_t    file_id, dset_id, dataspace_id, memoryspace, datatype;
   herr_t   error;

    H5T_class_t dtype_class;
   //H5T_order_t order;
   //size_t      size;

    file_id = H5Fopen(_filename.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    std::cout<<"pre open dset "<<std::endl;
    //std::cout<<"pre open dset "<<std::endl;
    dset_id = H5Dopen2(file_id, path.c_str(), H5P_DEFAULT);


    datatype = H5Dget_type(dset_id);
    dtype_class = H5Tget_class(datatype);
    //if (dtype_class == H5T_INTEGER)
    //   printf("Data set has INTEGER type \n");
    if (dtype_class == H5T_INTEGER)
       printf("Data set has INTEGER type \n");

    dataspace_id = H5Dget_space(dset_id);
    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_out = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    hsize_t* sel_dims = new hsize_t[rank];
    std::cout<<"rank = "<<rank<< std::endl;
    unsigned int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_out[0], NULL);

    unsigned long total = 1;
    for (int i=0; i < rank; i++)
    {
       std::cout<<"dims ["<<i<<"] ="<<dims_out[i]<< std::endl;
       total *= dims_out[i];
       offset[i] = 0;
       sel_dims[i] = count[i] = dims_out[i];
    }
    std::cout<<"total = "<<total<<std::endl;

    //std::cout<<"allocated "<<std::endl;
    sel_dims[layout.spectrum_dim] = dims_out[layout.spectrum_dim];
    sel_dims[layout.row_dim] = 1;
    sel_dims[layout.col_dim] = 1;
    count[layout.spectrum_dim] = dims_out[layout.spectrum_dim];
    count[layout.row_dim] = 1;
    count[layout.col_dim] = 1;
    unsigned short *buffer = new unsigned short[total];
    std::cout<<"allocated "<<std::endl;
    memoryspace = H5Screate_simple(rank, sel_dims, NULL);
    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    //float *buffer = new float[total];
    offset[layout.col_dim] = 0;
    for (int col=0; col < dims_out[layout.col_dim]; col++)
    {
      offset[layout.row_dim] = 0;
      data_struct::xrf::Spectra_Line *spectra_line = new data_struct::xrf::Spectra_Line();
      for (int row=0; row < dims_out[layout.row_dim]; row++)
      {
         data_struct::xrf::Spectra *spectra = new data_struct::xrf::Spectra(dims_out[layout.spectrum_dim]);
         float* spectra_buffer = spectra->get_buffer();
         H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
         error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace, dataspace_id, H5P_DEFAULT, spectra_buffer);

//         for(unsigned long i=0; i<dims_out[layout.spectrum_dim]; i++)
//         {
//            std::cout<<spectra_buffer[i]<<" ";
//         }
//         std::cout<<std::endl;

         spectra_line->append_spectra(spectra);

         offset[layout.row_dim] ++;
      }
      std::cout<<"read col "<<col<<std::endl;
      spec_vol->append_spectra_line(spectra_line);
      offset[layout.col_dim] ++;
    }

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
    error = H5Dread(dset_id, H5T_NATIVE_USHORT, memoryspace, dataspace_id, H5P_DEFAULT, buffer);
    printf("read in \n");
*/
}

bool HDF5_IO::save_spectra_volume(const std::string filename,
                                  const std::string path,
                                  const data_struct::xrf::Spectra_Volume * const spectra_volume,
                                  size_t row_idx_start,
                                  int row_idx_end,
                                  size_t col_idx_start,
                                  int col_idx_end)
{
    std::unique_lock<std::mutex> lock(_mutex);

//herr_t (*old_func)(void*);
//void *old_client_data;
//hid_t error_stack = H5Eget_current_stack();
//H5Eget_auto2(error_stack, &old_func, &old_client_data);

//H5Eset_auto2(error_stack, NULL, NULL);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    file_id, dset_id, dataspace_id, memoryspace, filespace, status, maps_grp_id, dcpl_id;
    herr_t   error;

    hsize_t chunk_dims[3];
    hsize_t dims_out[3];
    hsize_t offset[3];
    hsize_t count[3];

    if (row_idx_end < row_idx_start || row_idx_end > spectra_volume->rows() -1)
    {
        row_idx_end = spectra_volume->rows();
    }
    if (col_idx_end < col_idx_start || col_idx_end > spectra_volume->cols() -1)
    {
        col_idx_end = spectra_volume->cols();
    }

    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    if (file_id < 1)
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

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
    memoryspace = H5Screate_simple(3, count, NULL);
    filespace = H5Screate_simple(3, dims_out, NULL);

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    dataspace_id = H5Screate_simple (3, dims_out, NULL);

    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, NULL, count, NULL);

    maps_grp_id = H5Gopen(file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dset_id = H5Dcreate (maps_grp_id, path.c_str(), H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

    for(size_t row=row_idx_start; row < row_idx_end; row++)
    {
        offset[1] = row;
        for(size_t col=col_idx_start; col < col_idx_end; col++)
        {
            const data_struct::xrf::Spectra *spectra = &((*spectra_volume)[row][col]);
            offset[2] = col;
            H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

            status = H5Dwrite (dset_id, H5T_NATIVE_DOUBLE, memoryspace, filespace, H5P_DEFAULT, (void*)&(*spectra)[0]);
        }
    }



    H5Dclose(dset_id);
    H5Sclose(memoryspace);
    H5Sclose(filespace);
    H5Sclose(dataspace_id);
    H5Pclose(dcpl_id);
    H5Gclose(maps_grp_id);
    H5Fclose(file_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "\n\n save mca_arr elapsed time: " << elapsed_seconds.count() << "s\n\n\n";

    return true;

}

bool HDF5_IO::save_element_fits(std::string filename,
                                std::string path,
                                const data_struct::xrf::Fit_Count_Dict * const element_counts,
                                data_struct::xrf::Quantification_Standard * quantification_standard,
                                size_t row_idx_start,
                                int row_idx_end,
                                size_t col_idx_start,
                                int col_idx_end)
{
    std::unique_lock<std::mutex> lock(_mutex);

//hid_t error_stack = H5Eget_current_stack();
//H5Eset_auto2(error_stack, NULL, NULL);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    file_id, dset_id, memoryspace, filespace, dataspace_id, filetype, dataspace_ch_id, dataspace_ch_off_id, memtype, status, maps_grp_id, dset_ch_id, dcpl_id;
    herr_t   error;

    dset_id = -1;
    dset_ch_id = -1;
    hsize_t dims_out[3];
    hsize_t offset[3];
    hsize_t count[3];
    hsize_t chunk_dims[3];

    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    if (file_id < 1)
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(file_id < 0)
    {
        std::cout<<"Error opening file "<<filename<<std::endl;
        return false;
    }

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

    //dset_id = H5Dopen (maps_grp_id, path.c_str(), H5P_DEFAULT);
    if(dset_id < 0)
        dset_id = H5Dcreate (maps_grp_id, path.c_str(), H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    if(dset_id < 0)
    {
        std::cout<<"Error creating dataset "<<path<<std::endl;
        return false;
    }

    filetype = H5Tcopy (H5T_FORTRAN_S1);
    H5Tset_size (filetype, 255);
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, 256);

    dset_ch_id = H5Dopen (maps_grp_id, "channel_names", H5P_DEFAULT);
    if(dset_ch_id < 0)
        dset_ch_id = H5Dcreate (maps_grp_id, "channel_names", filetype, dataspace_ch_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dset_ch_id < 0)
    {
        std::cout<<"Error creating 'channel_names' dataset"<<std::endl;
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

    //save quantification
    for(auto& quant_itr : quantification_standard->_calibration_curves)
    {
        for(auto& shell_itr : quant_itr.second)
        {
            //create dataset for different shell curves
            std::unordered_map<std::string, real_t> calibration_curve = shell_itr.second;
            for (std::string el_name : element_lines)
            {
                if(calibration_curve.count(el_name) < 1 )
                {
                    continue;
                }
                //save the values
            }
        }
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
    H5Gclose(maps_grp_id);
    H5Fclose(file_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "\n\n save channels elapsed time: " << elapsed_seconds.count() << "s\n\n\n";



    return true;

}


} //end namespace file
}// end namespace io
