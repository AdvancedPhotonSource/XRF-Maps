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
#include <mutex>
#include <algorithm>
#include <cctype>

#include "data_struct/element_info.h"

#include "csv_io.h"

#define HDF5_SAVE_VERSION 10.0

#define HDF5_EXCHANGE_VERSION 1.0

const std::vector<std::string> hdf5_copy_dset_names = {"Element_Weights",
                                                       "Element_Weights_Names",
                                                       STR_DS_IC,
                                                       STR_US_IC,
                                                       "Standard_Name",
                                                       "Channel_Names",
													   "Channel_Units",
                                                       "Names",
                                                       "version"
                                                      };

const std::vector<std::string> hdf5_copy_grp_names = {"Scan"
                                                      };

const std::vector<std::string> xrf_analysis_save_names = { STR_FIT_ROI,
                                                          STR_FIT_GAUSS_MATRIX,
                                                          STR_FIT_SVD,
                                                          STR_FIT_NNLS
                                                         };





namespace io
{
namespace file
{

std::mutex HDF5_IO::_mutex;

HDF5_IO* HDF5_IO::_this_inst(nullptr);


struct Detector_HDF5_Struct
{
    hid_t    dset_id;
    hid_t    dataspace_id;
    real_t * buffer;
};

//-----------------------------------------------------------------------------

HDF5_IO::HDF5_IO()
{
	//disable hdf print to std err
	hid_t status;
    status = H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
    _cur_file_id = -1;
}

//-----------------------------------------------------------------------------

HDF5_IO* HDF5_IO::inst()
{
    std::lock_guard<std::mutex> lock(_mutex);

    if (_this_inst == nullptr)
    {
        _this_inst = new HDF5_IO();
    }
    return _this_inst;
}

//-----------------------------------------------------------------------------

HDF5_IO::~HDF5_IO()
{
	_cur_file_id = -1;
	_cur_filename = "";
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_open_h5_object(hid_t &id, H5_OBJECTS obj, std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map, std::string s1, hid_t id2, bool log_error, bool close_on_fail)
{
    if (obj == H5O_FILE)
    {
        id = H5Fopen(s1.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if(id < 0)
        {
			if (close_on_fail)
			{
				_close_h5_objects(close_map);
			}
			if (log_error)
			{
                logW << "Failed to open file " << s1 << "\n";
			}
            return false;
        }
    }
    else if (obj == H5O_GROUP)
    {
        id = H5Gopen(id2, s1.c_str(), H5P_DEFAULT);
        if(id < 0)
        {
			if (close_on_fail)
			{
				_close_h5_objects(close_map);
			}
			if (log_error)
			{
                logW << "Failed to open group " << s1 << "\n";
			}
            return false;
        }
    }
    else if (obj ==  H5O_DATASET)
    {
        id = H5Dopen2(id2, s1.c_str(), H5P_DEFAULT);
        if(id < 0)
        {
			if (close_on_fail)
			{
				_close_h5_objects(close_map);
			}
			if (log_error)
			{
                logW << "Failed to open dataset " << s1 << "\n";
			}
            return false;
        }
    }
    else
    {
		if (close_on_fail)
		{
			_close_h5_objects(close_map);
		}
		if (log_error)
		{
            logW << "Unknown H5_OBJECT type " << obj << "\n";
		}
        return false;
    }

    close_map.push({id, obj});
    return true;
}

void HDF5_IO::_close_h5_objects(std::stack<std::pair<hid_t, H5_OBJECTS> >  &close_map)
{

    herr_t err;
    while (false == close_map.empty())
    {
        auto& obj = close_map.top();
        close_map.pop();

        if (obj.second == H5O_FILE)
        {
            err = H5Fclose(obj.first);
            if (err > 0)
            {
                logW<<"Could not close h5 file id "<<obj.first<<"\n";
            }
        }
        else if (obj.second == H5O_GROUP)
        {
            err = H5Gclose(obj.first);
            if (err > 0)
            {
                logW<<"Could not close h5 group id "<<obj.first<<"\n";
            }
        }
        else if (obj.second == H5O_DATASPACE)
        {
            err = H5Sclose(obj.first);
            if (err > 0)
            {
                logW<<"Could not close h5 dataspace id "<<obj.first<<"\n";
            }
        }
        else if (obj.second == H5O_DATASET)
        {
            err = H5Dclose(obj.first);
            if (err > 0)
            {
                logW<<"Could not close h5 dataset id "<<obj.first<<"\n";
            }
        }
        else if(obj.second == H5O_ATTRIBUTE)
        {
            err = H5Aclose(obj.first);
            if (err > 0)
            {
                logW<<"Could not close h5 dataset id "<<obj.first<<"\n";
            }
        }
    }
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume(std::string path, size_t detector_num, data_struct::Spectra_Volume* spec_vol)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   logI<< path <<" detector : "<<detector_num<<"\n";

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

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id) )
       return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id) )
       return false;
    dataspace_rt_id = H5Dget_space(dset_rt_id);
    close_map.push({dataspace_rt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id) )
       return false;
    dataspace_inct_id = H5Dget_space(dset_incnt_id);
    close_map.push({dataspace_inct_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id) )
       return false;
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);
    close_map.push({dataspace_outct_id, H5O_DATASPACE});

    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logW<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<"\n";
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logE<<"Could not get dataset rank for MAPS_RAW/"<< detector_path<<"\n";
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       //logI<<"dims ["<<i<<"] ="<<dims_in[i]<< "\n";
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
    count_row[0] = dims_in[0];
    count_row[1] = dims_in[2];

/* TODO: maybe use greatest size (like xpress) becaue x_axis and y_axis will be diff than images size
    size_t greater_rows = std::max(spec_vol->rows() , dims_in[1]);
    size_t greater_cols = std::max(spec_vol->cols() , dims_in[2]);
    size_t greater_channels = std::max(spec_vol->samples_size() , dims_in[0]);
*/
/* Bug fix: mda files has dims_in[2] as one value, hdf5 has a larger value but all zeros, use mda size.
    if(spec_vol->rows() < dims_in[1] || spec_vol->cols() < dims_in[2] || spec_vol->samples_size() < dims_in[0])
    {
        spec_vol->resize_and_zero(dims_in[1], dims_in[2], dims_in[0]);
    }
*/
    count[1] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, nullptr);
    close_map.push({memoryspace_id, H5O_DATASPACE});
    memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
    close_map.push({memoryspace_meta_id, H5O_DATASPACE});
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    offset_meta[0] = detector_num;
    for (size_t row=0; row < spec_vol->rows(); row++)
    {
         offset[1] = row;
         offset_meta[1] = row;

         H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
         error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

         if (error > -1 )
         {
             for(size_t col=0; col<spec_vol->cols(); col++)
             {
                 offset_meta[2] = col;
                 data_struct::Spectra *spectra = &((*spec_vol)[row][col]);

                 H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                 error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                 spectra->elapsed_livetime(live_time);

                 H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                 error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                 spectra->elapsed_realtime(real_time);

                 H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                 error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                 spectra->input_counts(in_cnt);

                 H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                 error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                 spectra->output_counts(out_cnt);

                 spectra->recalc_elapsed_livetime();

                 for(size_t s=0; s<count_row[0]; s++)
                 {
                     (*spectra)[s] = buffer[(count_row[1] * s) + col];
                 }
                 //logD<<"saved col "<<col<<"\n";
             }

            //logD<<"read row "<<row<<"\n";
         }
         else
         {
            logE<<"reading row "<<row<<"\n";
         }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;
    delete [] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";
    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_line_xspress3(std::string path, size_t detector_num, data_struct::Spectra_Line* spec_row)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   logI<< path <<" detector : "<<detector_num<<"\n";

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
   hid_t    file_id, dset_id, dataspace_id, maps_grp_id, scaler_grp_id, memoryspace_id, memoryspace_meta_id, dset_lt_id;
   //hid_t    dset_incnt_id, dset_outcnt_id, dset_rt_id;
   hid_t    dataspace_lt_id;
  //hid_t dataspace_inct_id, dataspace_outct_id;
   herr_t   error;
   std::string detector_path;
   real_t * buffer;
   hsize_t offset_row[3] = {0,0,0};
   hsize_t count_row[3] = {0,0,0};
   hsize_t offset_meta[1] = {0};
   hsize_t count_meta[1] = {1};

   std::string live_time_dataset_name = "CHAN" + std::to_string(detector_num+1) + "SCA0";
   //std::string real_time_dataset_name = "CHAN" + std::to_string(detector_num+1) + "EventWidth";


    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
       return false;

    if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/entry/data", file_id) )
      return false;

    if ( false == _open_h5_object(scaler_grp_id, H5O_GROUP, close_map, "/entry/instrument/NDAttributes", file_id) )
      return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, "data", maps_grp_id) )
      return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, live_time_dataset_name, scaler_grp_id) )
      return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logW<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<"\n";
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logE<<"Could not get dataset rank for MAPS_RAW/"<< detector_path<<"\n";
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    buffer = new real_t [dims_in[0] * dims_in[2]]; // cols x spectra_size
    count_row[0] = dims_in[0];
    count_row[1] = 1;
    count_row[2] = dims_in[2];

    size_t greater_cols = std::max(spec_row->size() , (size_t)dims_in[0]);
    size_t greater_channels = std::max(spec_row[0].size() , (size_t)dims_in[2]);

    if( spec_row->size() < dims_in[0] || spec_row[0].size() < dims_in[2])
    {
        spec_row->resize_and_zero(greater_cols, greater_channels);
    }

    memoryspace_id = H5Screate_simple(3, count_row, nullptr);
    close_map.push({memoryspace_id, H5O_DATASPACE});
    memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
    close_map.push({memoryspace_meta_id, H5O_DATASPACE});
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
    //H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

    real_t live_time = 1.0;
    //real_t real_time = 1.0;
    //real_t in_cnt = 1.0;
    //real_t out_cnt = 1.0;

    //offset[1] = row;

    offset_row[1] = detector_num;

    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
    error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);


    if (error > -1 )
    {

        for(size_t col=0; col < dims_in[0]; col++)
        {
            offset_meta[0] = col;
            data_struct::Spectra *spectra = &((*spec_row)[col]);

            H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
            error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
            spectra->elapsed_livetime(live_time * 0.000000125 );
            
            //H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
            //error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
            //spectra->elapsed_realtime(real_time);

            //H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
            //error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
            //spectra->input_counts(in_cnt);

            //H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
            //error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
            //spectra->output_counts(out_cnt);

            //spectra->recalc_elapsed_livetime();
            


            for(size_t s=0; s<count_row[2]; s++) // num samples
            {
                //(*spectra)[s] = buffer[(count_row[1] * s) + col];
                (*spectra)[s] = buffer[(count_row[2] * col) + s];
            }
            //logD<<"saved col "<<col<<"\n";
        }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;
    delete [] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";
    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume_confocal(std::string path, size_t detector_num, data_struct::Spectra_Volume* spec_vol, bool log_error)
{
    std::lock_guard<std::mutex> lock(_mutex);

    //_is_loaded = ERROR_LOADING;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    bool confocal_ver_2020 = false;

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
    if (log_error)
    {
        if (detector_num == -1)
        {
            logI << path <<  "\n";
        }
        else
        {
            logI << path << " detector : " << detector_num << "\n";
        }
    }
    hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_detectors_id, memoryspace_id2;
    //hid_t    dset_xpos_id, dset_ypos_id, dataspace_xpos_id, dataspace_ypos_id;
    hid_t    dataspace_detectors_id;
    hid_t    attr_detector_names_id, attr_timebase_id;
    hid_t   elt_id = -1;
    hid_t   incnt_id = -1;
    hid_t   outcnt_id = -1;

    herr_t   error;
    std::string detector_path;
    char* detector_names[256];
    real_t time_base = 1.0f;
    real_t el_time = 1.0f;
    real_t* buffer;
    hsize_t offset_row[2] = { 0,0 };
    hsize_t count_row[2] = { 0,0 };
    hsize_t offset2[2] = { 0,0 };
    hsize_t count2[2] = { 1,1 };
    hsize_t offset_meta[3] = { 0,0,0 };
    hsize_t count_meta[3] = { 1,1,1 };
    std::unordered_map<std::string, int> detector_lookup;
    std::string elt_str = "Timer";
    std::string incnt_str = "ICR Ch ";
    std::string outcnt_str = "OCR Ch ";

    switch (detector_num)
    {
    case 0:
        detector_path = "MCA 1";
        incnt_str += "1";
        outcnt_str += "1";
        break;
    case 1:
        detector_path = "MCA 2";
        incnt_str += "2";
        outcnt_str += "2";
        break;
    case 2:
        detector_path = "MCA 3";
        incnt_str += "3";
        outcnt_str += "3";
        break;
    case 3:
        detector_path = "MCA 4";
        incnt_str += "4";
        outcnt_str += "4";
        break;
    default:
        detector_path = "";
        break;
    }

    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1, log_error))
        return false;

    if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "2D Scan", file_id))
    {
        _close_h5_objects(close_map);
        return false;
    }

    if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({ dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(dset_detectors_id, H5O_DATASET, close_map, "Detectors", maps_grp_id, false, false))
    {
        // try Detectors as group since it changed in 2020

        if (false == _open_h5_object(dset_detectors_id, H5O_GROUP, close_map, "Detectors", maps_grp_id))
        {
            _close_h5_objects(close_map);
            return false;
        }
        else
        {
            confocal_ver_2020 = true;
        }
    }

    if (confocal_ver_2020)
    {
        hsize_t nobj;
        H5Gget_num_objs(dset_detectors_id, &nobj);
        for (hsize_t i = 0; i < nobj; i++)
        {
            char str_grp_name[2048] = { 0 };
            hsize_t len = H5Gget_objname_by_idx(dset_detectors_id, i, str_grp_name, 2048);
            std::string scaler_name = std::string(str_grp_name, len);

            if (scaler_name == elt_str)
            {
                if (false == _open_h5_object(elt_id, H5O_DATASET, close_map, elt_str, dset_detectors_id))
                {
                    elt_id = -1;
                }
            }
            if (scaler_name == incnt_str)
            {
                if (false == _open_h5_object(incnt_id, H5O_DATASET, close_map, incnt_str, dset_detectors_id))
                {
                    incnt_id = -1;
                }
            }
            if (scaler_name == outcnt_str)
            {
                if (false == _open_h5_object(outcnt_id, H5O_DATASET, close_map, outcnt_str, dset_detectors_id))
                {
                    outcnt_id = -1;
                }
            }
        }
    }
    else
    {
        dataspace_detectors_id = H5Dget_space(dset_detectors_id);
        close_map.push({ dataspace_detectors_id, H5O_DATASPACE });
        attr_detector_names_id = H5Aopen(dset_detectors_id, "Detector Names", H5P_DEFAULT);
        close_map.push({ attr_detector_names_id, H5O_ATTRIBUTE });
    }
    attr_timebase_id = H5Aopen(dset_detectors_id, "TIMEBASE", H5P_DEFAULT);
    close_map.push({ attr_timebase_id, H5O_ATTRIBUTE });
    
//    if ( false == _open_h5_object(dset_xpos_id, H5O_DATASET, close_map, "X Positions", maps_grp_id) )
//    {
//        _close_h5_objects(close_map);
        //       return false;
//    }
//    dataspace_xpos_id = H5Dget_space(dset_xpos_id);
//    close_map.push({dataspace_xpos_id, H5O_DATASPACE});

//    if ( false == _open_h5_object(dset_ypos_id, H5O_DATASET, close_map, "Y Positions", maps_grp_id) )
    //{
    //_close_h5_objects(close_map);
//       return false;
    //}
//    dataspace_ypos_id = H5Dget_space(dset_ypos_id);
//    close_map.push({dataspace_ypos_id, H5O_DATASPACE});

    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logW<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<"\n";
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
		delete[] dims_in;
		delete[] offset;
		delete[] count;
         logE<<"Could not get dataset rank for MAPS_RAW/"<< detector_path<<"\n";
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       //logI<<"dims ["<<i<<"] ="<<dims_in[i]<< "\n";
       offset[i] = 0;
       count[i] = dims_in[i];
    }


    //chunking is 1 x col x samples
    buffer = new real_t [dims_in[1] * dims_in[2]]; //  cols x spectra_size
    count_row[0] = dims_in[1];
    count_row[1] = dims_in[2];

// TODO: maybe use greatest size (like xpress) becaue x_axis and y_axis will be diff than images size
//    size_t greater_rows = std::max(spec_vol->rows() , dims_in[1]);
//    size_t greater_cols = std::max(spec_vol->cols() , dims_in[2]);
//    size_t greater_channels = std::max(spec_vol->samples_size() , dims_in[0]);
	if (dims_in[0] == 0 && dims_in[1] == 0 && dims_in[2] == 0)
	{
		_close_h5_objects(close_map);
		delete[] buffer;
		delete[] dims_in;
		delete[] offset;
		delete[] count;
		return false;
	}

    if(spec_vol->rows() < dims_in[0] || spec_vol->cols() < dims_in[1] || spec_vol->samples_size() < dims_in[2])
    {
        spec_vol->resize_and_zero(dims_in[0], dims_in[1], dims_in[2]);
    }

    if (false == confocal_ver_2020)
    {
        int det_rank = H5Sget_simple_extent_ndims(dataspace_detectors_id);
        hsize_t* det_dims_in = new hsize_t[det_rank];
        H5Sget_simple_extent_dims(dataspace_detectors_id, &det_dims_in[0], nullptr);

        hid_t ftype = H5Aget_type(attr_detector_names_id);
        hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
        error = H5Aread(attr_detector_names_id, type, detector_names);
        if (error == 0)
        {
            //generate lookup map
            for (size_t z = 0; z < det_dims_in[2]; z++)
            {
                detector_lookup[std::string(detector_names[z])] = z;
                //free(detector_names[z]);
            }
        }
        delete[] det_dims_in;
    }
    error = H5Aread(attr_timebase_id, H5T_NATIVE_REAL, &time_base);
    

    count[0] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, nullptr);
    close_map.push({memoryspace_id, H5O_DATASPACE});

    memoryspace_id2 = H5Screate_simple(2, count2, nullptr);
    close_map.push({ memoryspace_id2, H5O_DATASPACE });

    memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
    close_map.push({memoryspace_meta_id, H5O_DATASPACE});
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

    real_t live_time = 1.0;
    //real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;


    for (size_t row=0; row < dims_in[0]; row++)
    {
         offset[0] = row;
         offset_meta[0] = row;

         H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
         error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

         if (error > -1 )
         {
             for(size_t col=0; col<dims_in[1]; col++)
             {
                 offset_meta[1] = col;
                 
                 data_struct::Spectra *spectra = &((*spec_vol)[row][col]);

                 

                 if (confocal_ver_2020)
                 {
                     offset2[0] = row;
                     offset2[1] = col;
                     H5Sselect_hyperslab(memoryspace_id2, H5S_SELECT_SET, offset2, nullptr, count2, nullptr);
                     if (elt_id > -1)
                     {
                         error = H5Dread(elt_id, H5T_NATIVE_REAL, memoryspace_id2, memoryspace_id2, H5P_DEFAULT, &live_time);
                         el_time = live_time / time_base;
                         spectra->elapsed_livetime(el_time);
                     }
                     if (incnt_id > -1)
                     {
                         error = H5Dread(incnt_id, H5T_NATIVE_REAL, memoryspace_id2, memoryspace_id2, H5P_DEFAULT, &in_cnt);
                         in_cnt *= 1000.0;
                         spectra->input_counts(in_cnt);
                     }
                     if (outcnt_id > -1)
                     {
                         error = H5Dread(outcnt_id, H5T_NATIVE_REAL, memoryspace_id2, memoryspace_id2, H5P_DEFAULT, &out_cnt);
                         out_cnt *= 1000.0;
                         spectra->output_counts(out_cnt);
                     }
                 }
                 else
                 {
                     offset_meta[2] = detector_lookup[elt_str];
                     H5Sselect_hyperslab(dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                     error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &live_time);
                     el_time = live_time / time_base;
                     spectra->elapsed_livetime(el_time);

                     offset_meta[2] = detector_lookup[incnt_str];
                     H5Sselect_hyperslab(dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                     error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &in_cnt);
                     spectra->input_counts(in_cnt * 1000.0);

                     offset_meta[2] = detector_lookup[outcnt_str];
                     H5Sselect_hyperslab(dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                     error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &out_cnt);
                     spectra->output_counts(out_cnt * 1000.0);
                 }

                 for(size_t s=0; s<dims_in[2]; s++)
                 {
                     (*spectra)[s] = buffer[(col * dims_in[2] ) + s];
                 }
             }
         }
         else
         {
            logE<<"Could not read row "<<row<<"\n";
         }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;
    delete [] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";
    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume_gsecars(std::string path, size_t detector_num, data_struct::Spectra_Volume* spec_vol, bool log_error)
{
    std::lock_guard<std::mutex> lock(_mutex);

    //_is_loaded = ERROR_LOADING;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
    if (log_error)
    {
        logI << path << " detector : " << detector_num << "\n";
    }
    hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_detectors_id;
    hid_t    dset_xypos_id, dataspace_xypos_id;
    hid_t	 livetime_id, realtime_id, inpcounts_id, outcounts_id;
    hid_t    livetime_dataspace_id, realtime_dataspace_id, inpcounts_dataspace_id, outcounts_dataspace_id;
    herr_t   error;
    std::string detector_path;
    char* detector_names[256];
    real_t time_base = 1.0f;
    real_t el_time = 1.0f;
    real_t* buffer;
    hsize_t offset_row[2] = { 0,0 };
    hsize_t count_row[2] = { 0,0 };
    hsize_t offset_meta[2] = { 0,0 };
    hsize_t count_meta[2] = { 1,1 };

    std::string counts_path;
    std::string incnt_path;
    std::string outcnt_path;
    std::string livetime_path;
    std::string realtime_path;
    GSE_CARS_SAVE_VER version = GSE_CARS_SAVE_VER::UNKNOWN;


    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1, log_error))
        return false;

    // check if this is old save style /xrfmap or new save style /xrmmap
    if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "xrmmap", file_id, false, false))
    {
        if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "xrfmap", file_id))
        {
            _close_h5_objects(close_map);
            return false;
        }
        else
        {
            version = GSE_CARS_SAVE_VER::XRFMAP;
        }
    }
    else
    {
        version = GSE_CARS_SAVE_VER::XRMMAP;
    }

    if (version == GSE_CARS_SAVE_VER::XRMMAP)
    {
        switch (detector_num)
        {
        case 0:
            detector_path = "mca1";
            break;
        case 1:
            detector_path = "mca2";
            break;
        case 2:
            detector_path = "mca3";
            break;
        case 3:
            detector_path = "mca4";
            break;
        default:
            detector_path = "mcasum";
            break;
        }
    }
    else if (version == GSE_CARS_SAVE_VER::XRFMAP)
    {
        switch (detector_num)
        {
        case 0:
            detector_path = "det1";
            break;
        case 1:
            detector_path = "det2";
            break;
        case 2:
            detector_path = "det3";
            break;
        case 3:
            detector_path = "det4";
            break;
        default:
            detector_path = "detsum";
            break;
        }
    }

    counts_path = detector_path + "/counts";
    incnt_path = detector_path + "/inpcounts";
    outcnt_path = detector_path + "/outcounts";
    realtime_path = detector_path + "/realtime";
    livetime_path = detector_path + "/livetime";

    if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, counts_path, maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({ dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(inpcounts_id, H5O_DATASET, close_map, incnt_path, maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    inpcounts_dataspace_id = H5Dget_space(inpcounts_id);
    close_map.push({ inpcounts_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(outcounts_id, H5O_DATASET, close_map, outcnt_path, maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    outcounts_dataspace_id = H5Dget_space(outcounts_id);
    close_map.push({ outcounts_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(realtime_id, H5O_DATASET, close_map, realtime_path, maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    realtime_dataspace_id = H5Dget_space(realtime_id);
    close_map.push({ realtime_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(livetime_id, H5O_DATASET, close_map, livetime_path, maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
	livetime_dataspace_id = H5Dget_space(livetime_id);
	close_map.push({ livetime_dataspace_id, H5O_DATASPACE });

	int rank = H5Sget_simple_extent_ndims(dataspace_id);
	if (rank != 3)
	{
		_close_h5_objects(close_map);
		logW << "Dataset /MAPS_RAW/" << detector_path << " rank != 3. rank = " << rank << ". Can't load dataset. returning" << "\n";
		return false;
		//throw exception ("Dataset is not a volume");
	}
	hsize_t* dims_in = new hsize_t[rank];
	hsize_t* offset = new hsize_t[rank];
	hsize_t* count = new hsize_t[rank];
	int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
	if (status_n < 0)
	{
		_close_h5_objects(close_map);
		logE << "Could not get dataset rank for MAPS_RAW/" << detector_path << "\n";
		return false;
	}

	for (int i = 0; i < rank; i++)
	{
		//logI<<"dims ["<<i<<"] ="<<dims_in[i]<< "\n";
		offset[i] = 0;
		count[i] = dims_in[i];
	}


	//chunking is 1 x col x samples
	buffer = new real_t[dims_in[1] * dims_in[2]]; //  cols x spectra_size
	count_row[0] = dims_in[1];
	count_row[1] = dims_in[2];

	if (dims_in[0] == 0 && dims_in[1] == 0 && dims_in[2] == 0)
	{
		_close_h5_objects(close_map);
		delete[] buffer;
		delete[] dims_in;
		delete[] offset;
		delete[] count;
		return false;
	}

	if (spec_vol->rows() < dims_in[0] || spec_vol->cols() < dims_in[1] || spec_vol->samples_size() < dims_in[2])
	{
		spec_vol->resize_and_zero(dims_in[0], dims_in[1], dims_in[2]);
	}

	count[0] = 1; //1 row

	memoryspace_id = H5Screate_simple(2, count_row, nullptr);
	close_map.push({ memoryspace_id, H5O_DATASPACE });
	memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
	close_map.push({ memoryspace_meta_id, H5O_DATASPACE });
	H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
	H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

	real_t live_time = 1.0;
	real_t real_time = 1.0;
	real_t in_cnt = 1.0;
	real_t out_cnt = 1.0;


	for (size_t row = 0; row < dims_in[0]; row++)
	{
		offset[0] = row;
		offset_meta[0] = row;

		H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
		error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

		if (error > -1) //no error
		{
			for (size_t col = 0; col < dims_in[1]; col++)
			{
				offset_meta[1] = col;

				data_struct::Spectra* spectra = &((*spec_vol)[row][col]);
				
				H5Sselect_hyperslab(livetime_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
				H5Sselect_hyperslab(realtime_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
				H5Sselect_hyperslab(inpcounts_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
				H5Sselect_hyperslab(outcounts_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

				error = H5Dread(realtime_id, H5T_NATIVE_REAL, memoryspace_meta_id, realtime_dataspace_id, H5P_DEFAULT, &real_time);
				if (error > -1)
				{
					spectra->elapsed_realtime(real_time);
				}
				error = H5Dread(livetime_id, H5T_NATIVE_REAL, memoryspace_meta_id, livetime_dataspace_id, H5P_DEFAULT, &live_time);
				if (error > -1)
				{
					spectra->elapsed_livetime(live_time);
				}
				error = H5Dread(inpcounts_id, H5T_NATIVE_REAL, memoryspace_meta_id, inpcounts_dataspace_id, H5P_DEFAULT, &in_cnt);
				if (error > -1)
				{
					spectra->input_counts(in_cnt);
				}
				error = H5Dread(outcounts_id, H5T_NATIVE_REAL, memoryspace_meta_id, outcounts_dataspace_id, H5P_DEFAULT, &out_cnt);
				if (error > -1)
				{
					spectra->output_counts(out_cnt);
				}
				
				//spectra->recalc_elapsed_livetime();

				for (size_t s = 0; s < dims_in[2]; s++)
				{
					(*spectra)[s] = buffer[(col * dims_in[2]) + s];
				}
			}
		}
		else
		{
			logE << "Could not read row " << row << "\n";
		}
	}

	delete[] dims_in;
	delete[] offset;
	delete[] count;
	delete[] buffer;

	_close_h5_objects(close_map);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	//std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";
	return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume_bnl(std::string path, size_t detector_num, data_struct::Spectra_Volume* spec_vol, bool log_error)
{
    std::lock_guard<std::mutex> lock(_mutex);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
    if (log_error)
    {
        logI << path << " detector : " << detector_num << "\n";
    }
    hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_detectors_id;
    hid_t    dset_xypos_id, dataspace_xypos_id;
    hid_t	 livetime_id, realtime_id, inpcounts_id, outcounts_id;
    hid_t    livetime_dataspace_id, realtime_dataspace_id, inpcounts_dataspace_id, outcounts_dataspace_id;
    herr_t   error;
    std::string detector_path;
    char* detector_names[256];
    real_t time_base = 1.0f;
    real_t el_time = 1.0f;
    real_t* buffer;
    hsize_t offset_row[2] = { 0,0 };
    hsize_t count_row[2] = { 0,0 };
    hsize_t offset_meta[2] = { 0,0 };
    hsize_t count_meta[2] = { 1,1 };

    std::string counts_path;
    std::string incnt_path;
    std::string outcnt_path;
    std::string livetime_path;
    std::string realtime_path;


    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1, log_error))
        return false;

    if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "xrfmap", file_id, false, false))
    {
        _close_h5_objects(close_map);
        return false;
    }

    switch (detector_num)
    {
    case -1:
        detector_path = "detsum";
        break;
    case 0:
        detector_path = "det1";
        break;
    case 1:
        detector_path = "det2";
        break;
    case 2:
        detector_path = "det3";
        break;
    case 3:
        detector_path = "det4";
        break;
    default:
        detector_path = "detsum";
        break;
    }

    counts_path = detector_path + "/counts";
    //incnt_path = detector_path + "/inpcounts";
    //outcnt_path = detector_path + "/outcounts";
    //realtime_path = detector_path + "/realtime";
    //livetime_path = detector_path + "/livetime";

    if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, counts_path, maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({ dataspace_id, H5O_DATASPACE });
    /*
    if (false == _open_h5_object(inpcounts_id, H5O_DATASET, close_map, incnt_path, maps_grp_id))
        return false;
    inpcounts_dataspace_id = H5Dget_space(inpcounts_id);
    close_map.push({ inpcounts_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(outcounts_id, H5O_DATASET, close_map, outcnt_path, maps_grp_id))
        return false;
    outcounts_dataspace_id = H5Dget_space(outcounts_id);
    close_map.push({ outcounts_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(realtime_id, H5O_DATASET, close_map, realtime_path, maps_grp_id))
        return false;
    realtime_dataspace_id = H5Dget_space(realtime_id);
    close_map.push({ realtime_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(livetime_id, H5O_DATASET, close_map, livetime_path, maps_grp_id))
        return false;
    livetime_dataspace_id = H5Dget_space(livetime_id);
    close_map.push({ livetime_dataspace_id, H5O_DATASPACE });
    */
    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logW << "Dataset /MAPS_RAW/" << detector_path << " rank != 3. rank = " << rank << ". Can't load dataset. returning" << "\n";
        return false;
        //throw exception ("Dataset is not a volume");
    }

    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];

    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
    if (status_n < 0)
    {
        _close_h5_objects(close_map);
        logE << "Could not get dataset rank for MAPS_RAW/" << detector_path << "\n";
        return false;
    }

    // width x height x samples
    for (int i = 0; i < rank; i++)
    {
        //logI<<"dims ["<<i<<"] ="<<dims_in[i]<< "\n";
        offset[i] = 0;
        count[i] = dims_in[i];
    }


    //chunking is 1 x col x samples
    buffer = new real_t[dims_in[1] * dims_in[2]]; //  cols x spectra_size
    count_row[0] = dims_in[1];
    count_row[1] = dims_in[2];

    if (dims_in[0] == 0 && dims_in[1] == 0 && dims_in[2] == 0)
    {
        _close_h5_objects(close_map);
        delete[] buffer;
        delete[] dims_in;
        delete[] offset;
        delete[] count;
        return false;
    }

    if (spec_vol->rows() < dims_in[0] || spec_vol->cols() < dims_in[1] || spec_vol->samples_size() < dims_in[2])
    {
        spec_vol->resize_and_zero(dims_in[0], dims_in[1], dims_in[2]);
    }

    count[0] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, nullptr);
    close_map.push({ memoryspace_id, H5O_DATASPACE });
    memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
    close_map.push({ memoryspace_meta_id, H5O_DATASPACE });
    H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
    H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;


    for (size_t row = 0; row < dims_in[0]; row++)
    {
        offset[0] = row;
        offset_meta[0] = row;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

        if (error > -1) //no error
        {
            for (size_t col = 0; col < dims_in[1]; col++)
            {
                offset_meta[1] = col;

                data_struct::Spectra* spectra = &((*spec_vol)[row][col]);
                /*
                H5Sselect_hyperslab(livetime_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                H5Sselect_hyperslab(realtime_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                H5Sselect_hyperslab(inpcounts_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                H5Sselect_hyperslab(outcounts_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                
                error = H5Dread(realtime_id, H5T_NATIVE_REAL, memoryspace_meta_id, realtime_dataspace_id, H5P_DEFAULT, &real_time);
                if (error > -1)
                {
                    spectra->elapsed_realtime(real_time);
                }
                error = H5Dread(livetime_id, H5T_NATIVE_REAL, memoryspace_meta_id, livetime_dataspace_id, H5P_DEFAULT, &live_time);
                if (error > -1)
                {
                    spectra->elapsed_livetime(live_time);
                }
                error = H5Dread(inpcounts_id, H5T_NATIVE_REAL, memoryspace_meta_id, inpcounts_dataspace_id, H5P_DEFAULT, &in_cnt);
                if (error > -1)
                {
                    spectra->input_counts(in_cnt);
                }
                error = H5Dread(outcounts_id, H5T_NATIVE_REAL, memoryspace_meta_id, outcounts_dataspace_id, H5P_DEFAULT, &out_cnt);
                if (error > -1)
                {
                    spectra->output_counts(out_cnt);
                }

                //spectra->recalc_elapsed_livetime();
                */
                for (size_t s = 0; s < dims_in[2]; s++)
                {
                    (*spectra)[s] = buffer[(col * dims_in[2]) + s];
                }
            }
        }
        else
        {
            logE << "Could not read row " << row << "\n";
        }
    }

    delete[] dims_in;
    delete[] offset;
    delete[] count;
    delete[] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";
    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_integrated_spectra_bnl(std::string path, size_t detector_num, data_struct::Spectra* spec, bool log_error)
{
    std::lock_guard<std::mutex> lock(_mutex);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
    if (log_error)
    {
        logI << path << " detector : " << detector_num << "\n";
    }
    hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_detectors_id;
    hid_t    dset_xypos_id, dataspace_xypos_id;
    hid_t	 livetime_id, realtime_id, inpcounts_id, outcounts_id;
    hid_t    livetime_dataspace_id, realtime_dataspace_id, inpcounts_dataspace_id, outcounts_dataspace_id;
    herr_t   error;
    std::string detector_path;
    char* detector_names[256];
    real_t time_base = 1.0f;
    real_t el_time = 1.0f;
    real_t* buffer;
    hsize_t offset_row[2] = { 0,0 };
    hsize_t count_row[2] = { 0,0 };
    hsize_t offset_meta[2] = { 0,0 };
    hsize_t count_meta[2] = { 1,1 };

    std::string counts_path;
    std::string incnt_path;
    std::string outcnt_path;
    std::string livetime_path;
    std::string realtime_path;


    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1, log_error))
        return false;

    if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "xrfmap", file_id, false, false))
    {
        _close_h5_objects(close_map);
        return false;
    }

    switch (detector_num)
    {
    case -1:
        detector_path = "detsum";
        break;
    case 0:
        detector_path = "det1";
        break;
    case 1:
        detector_path = "det2";
        break;
    case 2:
        detector_path = "det3";
        break;
    case 3:
        detector_path = "det4";
        break;
    default:
        detector_path = "detsum";
        break;
    }

    counts_path = detector_path + "/counts";
    //incnt_path = detector_path + "/inpcounts";
    //outcnt_path = detector_path + "/outcounts";
    //realtime_path = detector_path + "/realtime";
    //livetime_path = detector_path + "/livetime";

    if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, counts_path, maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({ dataspace_id, H5O_DATASPACE });
    /*
    if (false == _open_h5_object(inpcounts_id, H5O_DATASET, close_map, incnt_path, maps_grp_id))
        return false;
    inpcounts_dataspace_id = H5Dget_space(inpcounts_id);
    close_map.push({ inpcounts_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(outcounts_id, H5O_DATASET, close_map, outcnt_path, maps_grp_id))
        return false;
    outcounts_dataspace_id = H5Dget_space(outcounts_id);
    close_map.push({ outcounts_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(realtime_id, H5O_DATASET, close_map, realtime_path, maps_grp_id))
        return false;
    realtime_dataspace_id = H5Dget_space(realtime_id);
    close_map.push({ realtime_dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(livetime_id, H5O_DATASET, close_map, livetime_path, maps_grp_id))
        return false;
    livetime_dataspace_id = H5Dget_space(livetime_id);
    close_map.push({ livetime_dataspace_id, H5O_DATASPACE });
    */
    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logW << "Dataset /MAPS_RAW/" << detector_path << " rank != 3. rank = " << rank << ". Can't load dataset. returning" << "\n";
        return false;
        //throw exception ("Dataset is not a volume");
    }

    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];

    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
    if (status_n < 0)
    {
        _close_h5_objects(close_map);
        logE << "Could not get dataset rank for MAPS_RAW/" << detector_path << "\n";
        return false;
    }

    // width x height x samples
    for (int i = 0; i < rank; i++)
    {
        //logI<<"dims ["<<i<<"] ="<<dims_in[i]<< "\n";
        offset[i] = 0;
        count[i] = dims_in[i];
    }


    //chunking is 1 x col x samples
    buffer = new real_t[dims_in[1] * dims_in[2]]; //  cols x spectra_size
    count_row[0] = dims_in[1];
    count_row[1] = dims_in[2];

    if (dims_in[0] == 0 && dims_in[1] == 0 && dims_in[2] == 0)
    {
        _close_h5_objects(close_map);
        delete[] buffer;
        delete[] dims_in;
        delete[] offset;
        delete[] count;
        return false;
    }

    if (spec->size() < dims_in[2] )
    {
        spec->resize(dims_in[2]);
        spec->setZero(dims_in[2]);
    }

    count[0] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, nullptr);
    close_map.push({ memoryspace_id, H5O_DATASPACE });
    memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
    close_map.push({ memoryspace_meta_id, H5O_DATASPACE });
    H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
    H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;


    for (size_t row = 0; row < dims_in[0]; row++)
    {
        offset[0] = row;
        offset_meta[0] = row;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

        if (error > -1) //no error
        {
            for (size_t col = 0; col < dims_in[1]; col++)
            {
                offset_meta[1] = col;
                /*
                H5Sselect_hyperslab(livetime_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                H5Sselect_hyperslab(realtime_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                H5Sselect_hyperslab(inpcounts_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                H5Sselect_hyperslab(outcounts_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

                error = H5Dread(realtime_id, H5T_NATIVE_REAL, memoryspace_meta_id, realtime_dataspace_id, H5P_DEFAULT, &real_time);
                if (error > -1)
                {
                    spectra->elapsed_realtime(real_time);
                }
                error = H5Dread(livetime_id, H5T_NATIVE_REAL, memoryspace_meta_id, livetime_dataspace_id, H5P_DEFAULT, &live_time);
                if (error > -1)
                {
                    spectra->elapsed_livetime(live_time);
                }
                error = H5Dread(inpcounts_id, H5T_NATIVE_REAL, memoryspace_meta_id, inpcounts_dataspace_id, H5P_DEFAULT, &in_cnt);
                if (error > -1)
                {
                    spectra->input_counts(in_cnt);
                }
                error = H5Dread(outcounts_id, H5T_NATIVE_REAL, memoryspace_meta_id, outcounts_dataspace_id, H5P_DEFAULT, &out_cnt);
                if (error > -1)
                {
                    spectra->output_counts(out_cnt);
                }

                //spectra->recalc_elapsed_livetime();
                */
                for (size_t s = 0; s < dims_in[2]; s++)
                {
                    (*spec)[s] += buffer[(col * dims_in[2]) + s];
                }
            }
        }
        else
        {
            logE << "Could not read row " << row << "\n";
        }
    }

    delete[] dims_in;
    delete[] offset;
    delete[] count;
    delete[] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";
    return true;
}

//-----------------------------------------------------------------------------

int parse_str_val_to_int(std::string start_delim, std::string end_delim, std::string lookup_str)
{
    int find_idx = lookup_str.find(start_delim);
    if(find_idx < 0)
    {
        logW<<"Could not find string property start delimeter "<<start_delim<<"\n";
        return -1;
    }
    std::string leftover = lookup_str.substr(find_idx+start_delim.length());
    find_idx = leftover.find(end_delim); // find end double quote
    if(find_idx < 0)
    {
        logW<<"Could not find spectra sample count end delimeter \"\n";
        return -1;
    }
    std::string str_val = leftover.substr(0, find_idx);
    return std::stoi(str_val.c_str());
}


//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume_emd(std::string path,
                                        size_t frame_num,
                                        data_struct::Spectra_Volume *spec_vol,
                                        bool logerr)
{
    std::lock_guard<std::mutex> lock(_mutex);

    const std::string STR_BINCOUNT = "\"bincount\": \"";
    const std::string STR_WIDTH = "\"Width\": \"";
    const std::string STR_HEIGHT = "\"Height\": \"";

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
    int height = -1;
    int width = -1;

    hsize_t start_offset = 0;

    logI << path << " frames : " << frame_num << "\n";

    hid_t    file_id, maps_grp_id, memoryspace_id, spectra_grp_id, dataset_id, acqui_id, frame_id, meta_id;
    hid_t    image_grp_id, image_hash_grp_id, image_id, image_ds_id;
    hid_t    dataspace_id, dataspace_acqui_id, dataspace_frame_id, dataspace_meta_id;
    herr_t   error;
    std::string detector_path;
    hsize_t image_dims[3] = {0,0,0};
    char str_grp_name[2048] = { 0 };

    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
        return false;


    if (false == _open_h5_object(image_grp_id, H5O_GROUP, close_map, "/Data/Image", file_id))
        return false;

    // get hash group name for first image
    error = H5Gget_objname_by_idx(image_grp_id, 0, str_grp_name, 2048);

    if (false == _open_h5_object(image_hash_grp_id, H5O_GROUP, close_map, std::string(str_grp_name), image_grp_id))
        return false;

    if (false == _open_h5_object(image_id, H5O_DATASET, close_map, "Data", image_hash_grp_id))
        return false;
    image_ds_id = H5Dget_space(image_id);
    close_map.push({ image_ds_id, H5O_DATASPACE });

    error = H5Sget_simple_extent_dims(image_ds_id, &image_dims[0], nullptr);
    if(error > -1)
    {
        width = image_dims[0];
        height = image_dims[1];
    }
    else
    {
        if (logerr)
        {
            logE << "Could not get image dimensions from /Data/Image/" << str_grp_name << "/Data\n ";
        }
        _close_h5_objects(close_map);
        return false;
    }

    if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/Data/SpectrumStream", file_id))
        return false;

    error = H5Gget_objname_by_idx(maps_grp_id, 0, str_grp_name, 2048);

    if (false == _open_h5_object(spectra_grp_id, H5O_GROUP, close_map, std::string(str_grp_name), maps_grp_id))
        return false;

    if (false == _open_h5_object(dataset_id, H5O_DATASET, close_map, "Data", spectra_grp_id))
        return false;
    dataspace_id = H5Dget_space(dataset_id);
    close_map.push({ dataspace_id, H5O_DATASPACE });

    if (false == _open_h5_object(acqui_id, H5O_DATASET, close_map, "AcquisitionSettings", spectra_grp_id))
        return false;
    dataspace_acqui_id = H5Dget_space(acqui_id);
    close_map.push({ dataspace_acqui_id, H5O_DATASPACE });

    if(frame_num > 0)
    {
        if (false == _open_h5_object(frame_id, H5O_DATASET, close_map, "FrameLocationTable", spectra_grp_id))
            return false;
        dataspace_frame_id = H5Dget_space(frame_id);
        close_map.push({ dataspace_frame_id, H5O_DATASPACE });

        hsize_t frame_dims[2] = {0,0};
        error = H5Sget_simple_extent_dims(dataspace_frame_id, &frame_dims[0], nullptr);

        if(frame_num < frame_dims[0])
        {
            int *frames = new int[frame_dims[0]];
            //read frames
            hid_t frame_memoryspace_id = H5Screate_simple(2, frame_dims, nullptr);
            error = H5Dread(frame_id, H5T_NATIVE_INT, frame_memoryspace_id, dataspace_frame_id, H5P_DEFAULT, &frames[0]);
            start_offset = frames[frame_num];
            delete [] frames;
            close_map.push({frame_memoryspace_id, H5O_DATASPACE});
        }
        else
        {
            logW<<"frame_num "<<frame_num<<" > Max Frames"<<frame_dims[0]<<". Setting frame_num = 0\n";
        }

    }
    if (false == _open_h5_object(meta_id, H5O_DATASET, close_map, "Metadata", spectra_grp_id))
        return false;
    dataspace_meta_id = H5Dget_space(meta_id);
    close_map.push({ dataspace_meta_id, H5O_DATASPACE });

    char *acquisition[1];
    hid_t ftype = H5Dget_type(acqui_id);
    hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
    error = H5Dread(acqui_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, acquisition);

    //prase bincount/samples width and height
    std::string str_acqui = std::string(acquisition[0]);
    int samples = parse_str_val_to_int(STR_BINCOUNT, "\"", str_acqui);

    if( samples < 0 || height < 0 || width < 0)
    {
        logE<<"Unknown spectra volume size Width: "<<width<<" Height: "<<height<<" Samples: "<<samples<<"\n";
        _close_h5_objects(close_map);
        return false;
    }

    spec_vol->resize_and_zero(height, width, samples);

    hsize_t rank;
    rank = H5Sget_simple_extent_ndims(dataspace_id);
    unsigned short *buffer;
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
    hsize_t* chunk_dims = new hsize_t[rank];
    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
    if (status_n < 0)
    {
        _close_h5_objects(close_map);
        logE << "Getting dataset rank for MAPS_RAW/" << detector_path << "\n";
        return false;
    }

    for (hsize_t i = 0; i < rank; i++)
    {
        offset[i] = 0;
        count[i] = dims_in[i];
    }

    hid_t dcpl_id = H5Dget_create_plist(dataset_id);
    error = H5Pget_chunk(dcpl_id, rank, chunk_dims);
    buffer = new unsigned short[chunk_dims[0]];

    if(start_offset > count[0])
    {
        logW<<"frame start offset: "<<start_offset<<" > dataset count: "<<count[0]<<". Setting start offset = 0\n.";
        start_offset = 0;
    }

    memoryspace_id = H5Screate_simple(2, chunk_dims, nullptr);
    close_map.push({memoryspace_id, H5O_DATASPACE});
    H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);

    size_t read_leftover = (count[0] - start_offset) % chunk_dims[0];
    size_t read_amt = (count[0] - start_offset) / chunk_dims[0];

    int row = 0;
    int col = 0;

    offset[0] = start_offset; // start of frame offset
    real_t incnt = 0.0;
    real_t outcnt = 0.0;
    real_t elt = 1.0; // TODO: read from metadata
    data_struct::Spectra *spectra = &((*spec_vol)[row][col]);
    for (size_t i = 0; i < read_amt; i++)
    {

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
        error = H5Dread(dataset_id, H5T_NATIVE_USHORT, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
        offset[0] += chunk_dims[0];

        if (error > -1)
        {
            for (hsize_t j = 0; j < chunk_dims[0]; j++)
            {
                //check for delim to move to new spectra
                if (buffer[j] == 65535)
                {
                    col++;
                    if(col >= width)
                    {
                        col = 0;
                        row++;
                        logI<<"Reading row "<<row<<" of "<<height<<"\n";
                    }
                    if(row >= height)
                    {
                        delete[] dims_in;
                        delete[] offset;
                        delete[] count;
                        delete[] chunk_dims;
                        delete[] buffer;

                        _close_h5_objects(close_map);

                        end = std::chrono::system_clock::now();
                        std::chrono::duration<double> elapsed_seconds = end - start;
                        //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

                        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

                        return true;
                    }
                    spectra->input_counts(incnt);
                    spectra->output_counts(outcnt);
                    spectra->elapsed_livetime(elt);
                    incnt = 0.0;
                    outcnt = 0.0;
                    spectra = &((*spec_vol)[row][col]);
                }
                else
                {
                    (*spectra)[buffer[j]] += 1.0;
                    incnt += 1.0;
                    outcnt += 1.0;
                }
            }
        }
        else
        {
            logE << "reading chunk " << i << " of "<< read_amt <<"\n";
        }

    }

    chunk_dims[0] = read_leftover;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
    //offset[0] = 0;
    //H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
    error = H5Dread(dataset_id, H5T_NATIVE_USHORT, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
    if (error > -1)
    {
        for (hsize_t j = 0; j < read_leftover; j++)
        {
            //check for delim to move to new spectra
            if (buffer[j] == 65535)
            {
                col++;
                if(col >= width)
                {
                    col = 0;
                    row++;
                    logI<<"Reading row "<<row<<" of "<<height<<"\n";
                }
                if(row >= height)
                {
                    col = 0;
                    row = 0;
                }
                spectra = &((*spec_vol)[row][col]);
            }
            else
            {
                (*spectra)[buffer[j]] += 1.0;
            }
        }
    }
    else
    {
        logE << "reading last chunk \n";
    }

    delete[] dims_in;
    delete[] offset;
    delete[] count;
    delete[] chunk_dims;
    delete[] buffer;

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume_emd_with_callback(std::string path,
	const std::vector<size_t>& detector_num_arr,
	data_struct::IO_Callback_Func_Def callback_func,
	void* user_data)
{
	std::lock_guard<std::mutex> lock(_mutex);
//TDOO: change to unique_lock so we can unlock and allow stream saver to be called


    const std::string STR_BINCOUNT = "\"bincount\": \"";
    const std::string STR_WIDTH = "\"Width\": \"";
    const std::string STR_HEIGHT = "\"Height\": \"";

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
    int height = -1;
    int width = -1;

    hsize_t start_offset = 0;

    hid_t    file_id, maps_grp_id, memoryspace_id, spectra_grp_id, dataset_id, acqui_id, frame_id, meta_id;
    hid_t    image_grp_id, image_hash_grp_id, image_id, image_ds_id;
	hid_t    dataspace_id, dataspace_acqui_id, dataspace_frame_id, dataspace_meta_id;
	herr_t   error;
	std::string detector_path;
    hsize_t image_dims[3] = {0,0,0};
	char str_grp_name[2048] = { 0 };
    size_t frame = detector_num_arr[0];
	size_t frame_idx = 0;

	if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
		return false;


    if (false == _open_h5_object(image_grp_id, H5O_GROUP, close_map, "/Data/Image", file_id))
        return false;

    // get hash group name for first image
    error = H5Gget_objname_by_idx(image_grp_id, 0, str_grp_name, 2048);

    if (false == _open_h5_object(image_hash_grp_id, H5O_GROUP, close_map, std::string(str_grp_name), image_grp_id))
        return false;

    if (false == _open_h5_object(image_id, H5O_DATASET, close_map, "Data", image_hash_grp_id))
        return false;
    image_ds_id = H5Dget_space(image_id);
    close_map.push({ image_ds_id, H5O_DATASPACE });

    error = H5Sget_simple_extent_dims(image_ds_id, &image_dims[0], nullptr);
    if(error > -1)
    {
        width = image_dims[0];
        height = image_dims[1];
    }
    else
    {
        logE<<"Could not get image dimensions from /Data/Image/"<<str_grp_name<<"/Data\n ";
        _close_h5_objects(close_map);
        return false;
    }


	if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/Data/SpectrumStream", file_id))
		return false;

	error = H5Gget_objname_by_idx(maps_grp_id, 0, str_grp_name, 2048);

	if (false == _open_h5_object(spectra_grp_id, H5O_GROUP, close_map, std::string(str_grp_name), maps_grp_id))
		return false;

	if (false == _open_h5_object(dataset_id, H5O_DATASET, close_map, "Data", spectra_grp_id))
		return false;
	dataspace_id = H5Dget_space(dataset_id);
	close_map.push({ dataspace_id, H5O_DATASPACE });

	if (false == _open_h5_object(acqui_id, H5O_DATASET, close_map, "AcquisitionSettings", spectra_grp_id))
		return false;
	dataspace_acqui_id = H5Dget_space(acqui_id);
	close_map.push({ dataspace_acqui_id, H5O_DATASPACE });

    if(frame > 0)
    {
        if (false == _open_h5_object(frame_id, H5O_DATASET, close_map, "FrameLocationTable", spectra_grp_id))
            return false;
        dataspace_frame_id = H5Dget_space(frame_id);
        close_map.push({ dataspace_frame_id, H5O_DATASPACE });

        hsize_t frame_dims[2] = {0,0};
        error = H5Sget_simple_extent_dims(dataspace_frame_id, &frame_dims[0], nullptr);

        if(frame < frame_dims[0])
        {
            int *frames = new int[frame_dims[0]];
            //read frames
            hid_t frame_memoryspace_id = H5Screate_simple(2, frame_dims, nullptr);
            error = H5Dread(frame_id, H5T_NATIVE_INT, frame_memoryspace_id, dataspace_frame_id, H5P_DEFAULT, &frames[0]);
            start_offset = frames[frame];
            delete [] frames;
            close_map.push({frame_memoryspace_id, H5O_DATASPACE});
        }
        else
        {
            logW<<" detector_num_arr[0] "<<frame<<" > Max Frames"<<frame_dims[0]<<". Setting frame_num_start = 0\n";
        }

    }
	if (false == _open_h5_object(meta_id, H5O_DATASET, close_map, "Metadata", spectra_grp_id))
		return false;
	dataspace_meta_id = H5Dget_space(meta_id);
	close_map.push({ dataspace_meta_id, H5O_DATASPACE });

	char *acquisition[1];
	hid_t ftype = H5Dget_type(acqui_id);
	hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
	error = H5Dread(acqui_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, acquisition);

    //prase bincount/samples width and height
	std::string str_acqui = std::string(acquisition[0]);
    int samples = parse_str_val_to_int(STR_BINCOUNT, "\"", str_acqui);
    //int height = parse_str_val_to_int(STR_HEIGHT, "\"", str_acqui);; //wrong place, this data is incorrect
    //int width = parse_str_val_to_int(STR_WIDTH, "\"", str_acqui);; // need to read image size in from /Data/Image/...

    if( samples < 0 || height < 0 || width < 0)
    {
        logE<<"Unknown spectra volume size Width: "<<width<<" Height: "<<height<<" Samples: "<<samples<<"\n";
        _close_h5_objects(close_map);
        return false;
    }


	hsize_t rank;
	rank = H5Sget_simple_extent_ndims(dataspace_id);
    unsigned short *buffer;
	hsize_t* dims_in = new hsize_t[rank];
	hsize_t* offset = new hsize_t[rank];
	hsize_t* count = new hsize_t[rank];
	hsize_t* chunk_dims = new hsize_t[rank];
	
    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
	if (status_n < 0)
	{
		_close_h5_objects(close_map);
		logE << "getting dataset rank for MAPS_RAW/" << detector_path << "\n";
		return false;
	}

    for (size_t i = 0; i < rank; i++)
	{
		
		offset[i] = 0;
		count[i] = dims_in[i];
	}

	hid_t dcpl_id = H5Dget_create_plist(dataset_id);
	H5Pget_chunk(dcpl_id, rank, chunk_dims);
    buffer = new unsigned short[chunk_dims[0]];

    if(start_offset > count[0])
    {
        logW<<" frame start offset: "<<start_offset<<" > dataset count: "<<count[0]<<". Setting start offset = 0\n.";
        start_offset = 0;
    }

    memoryspace_id = H5Screate_simple(2, chunk_dims, nullptr);
    close_map.push({memoryspace_id, H5O_DATASPACE});
    H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);

    size_t read_leftover = (count[0] - start_offset) % chunk_dims[0];
    size_t read_amt = (count[0] - start_offset) / chunk_dims[0];

	int row = 0;
	int col = 0;

    real_t incnt = 0.0;
    real_t outcnt = 0.0;
    real_t elt = 1.0;
    offset[0] = start_offset; // start of frame offset
	data_struct::Spectra *spectra = new data_struct::Spectra(samples);
	for (size_t i = 0; i < read_amt; i++)
	{
		
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
        error = H5Dread(dataset_id, H5T_NATIVE_USHORT, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
		offset[0] += chunk_dims[0];

		if (error > -1)
		{
            for (hsize_t j = 0; j < chunk_dims[0]; j++)
			{
				//check for delim to move to new spectra
				if (buffer[j] == 65535)
				{
                    spectra->input_counts(incnt);
                    spectra->output_counts(outcnt);
                    spectra->elapsed_livetime(elt);
                    incnt = 0.0;
                    outcnt = 0.0;

                    callback_func(row, col, height, width, frame, spectra, user_data);
                    col++;
                    if(col >= width)
                    {
                        col = 0;
                        row++;
                        logI<<"Reading row "<<row<<" of "<<height<<"\n";
                    }
                    if(row >= height)
                    {
                        col = 0;
                        row = 0;
						frame_idx++;
						frame = detector_num_arr[frame_idx];
                        if(frame_idx > detector_num_arr.size())
                        {
                            delete[] dims_in;
                            delete[] offset;
                            delete[] count;
                            delete[] chunk_dims;
                            delete[] buffer;

                            _close_h5_objects(close_map);

                            end = std::chrono::system_clock::now();
                            std::chrono::duration<double> elapsed_seconds = end - start;
                            //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

                            logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

                            return true;
                        }
                    }
                    spectra = new data_struct::Spectra(samples);
				}
				else
				{
					(*spectra)[buffer[j]] += 1.0;
                    incnt += 1.0;
                    outcnt += 1.0;
				}
			}
		}
		else
		{
			logE << "reading chunk " << i << " of "<< read_amt <<"\n";
		}
		
	}

    chunk_dims[0] = read_leftover;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
    //offset[0] = 0;
    //H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
    error = H5Dread(dataset_id, H5T_NATIVE_USHORT, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
	if (error > -1)
	{
        for (hsize_t j = 0; j < read_leftover; j++)
        {
            //check for delim to move to new spectra
            if (buffer[j] == 65535)
            {
                callback_func(row, col, height, width, frame, spectra, user_data);
                col++;
                if(col >= width)
                {
                    col = 0;
                    row++;
                    logI<<"Reading row "<<row<<" of "<<height<<"\n";
                }
                if(row >= height)
                {
                    col = 0;
                    row = 0;
                    frame++;
                }
                spectra = new data_struct::Spectra(samples);
            }
            else
            {
                (*spectra)[buffer[j]] += 1.0;
            }
        }
	}
	else
	{
		logE << "reading last chunk \n";
	}

	delete[] dims_in;
	delete[] offset;
	delete[] count;
	delete[] chunk_dims;
	delete[] buffer;

	_close_h5_objects(close_map);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	//std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_spectra_volume_with_callback(std::string path,
												const std::vector<size_t>& detector_num_arr,
												data_struct::IO_Callback_Func_Def callback_func,
                                                void* user_data)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();
   std::map<size_t, struct Detector_HDF5_Struct> detector_hid_map;

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   hid_t    file_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
   hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
   herr_t   error = -1;
   std::string detector_path;
   hsize_t offset_row[2] = {0,0};
   hsize_t count_row[2] = {0,0};
   hsize_t offset_meta[3] = {0,0,0};
   hsize_t count_meta[3] = {1,1,1};

   if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
       return false;

   if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id) )
      return false;
    for(size_t detector_num : detector_num_arr)
    {
        detector_hid_map.insert( {detector_num, Detector_HDF5_Struct()} );

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


        if ( false == _open_h5_object(detector_hid_map[detector_num].dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id) )
           return false;
        detector_hid_map[detector_num].dataspace_id = H5Dget_space(detector_hid_map[detector_num].dset_id);
        close_map.push({detector_hid_map[detector_num].dataspace_id, H5O_DATASPACE});
    }

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id) )
       return false;
    dataspace_rt_id = H5Dget_space(dset_rt_id);
    close_map.push({dataspace_rt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id) )
       return false;
    dataspace_inct_id = H5Dget_space(dset_incnt_id);
    close_map.push({dataspace_inct_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id) )
       return false;
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);
    close_map.push({dataspace_outct_id, H5O_DATASPACE});


    int rank = H5Sget_simple_extent_ndims(detector_hid_map[detector_num_arr[0]].dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logE<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<"\n";
        return false;
       //throw exception ("Dataset is not a volume");
    }
    hsize_t* dims_in = new hsize_t[rank];
    hsize_t* offset = new hsize_t[rank];
    hsize_t* count = new hsize_t[rank];
   
    int status_n = H5Sget_simple_extent_dims(detector_hid_map[detector_num_arr[0]].dataspace_id, &dims_in[0], nullptr);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logE<<"getting dataset rank for MAPS_RAW/"<< detector_path<<"\n";
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    for(size_t detector_num : detector_num_arr)
    {
        detector_hid_map[detector_num].buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
    }
    count_row[0] = dims_in[0];
    count_row[1] = dims_in[2];

    //if(spec_vol->rows() < dims_in[1] || spec_vol->cols() < dims_in[2] || spec_vol->samples_size() < dims_in[0])
    //{
    //    spec_vol->resize(dims_in[1], dims_in[2], dims_in[0]);
    //}

    count[1] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, nullptr);
    memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    for (size_t row=0; row < dims_in[1]; row++)
    {
         offset[1] = row;
         offset_meta[1] = row;

         for(size_t detector_num : detector_num_arr)
         {
            H5Sselect_hyperslab (detector_hid_map[detector_num].dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            error = H5Dread(detector_hid_map[detector_num].dset_id, H5T_NATIVE_REAL, memoryspace_id, detector_hid_map[detector_num].dataspace_id, H5P_DEFAULT, detector_hid_map[detector_num].buffer);
         }

         if (error > -1 )
         {
             for(size_t col=0; col<count_row[1]; col++)
             {
                 offset_meta[2] = col;

                 for(size_t detector_num : detector_num_arr)
                 {
                     offset_meta[0] = detector_num;
                     data_struct::Spectra * spectra = new data_struct::Spectra(dims_in[0]);

                     H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                     error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                     spectra->elapsed_livetime(live_time);

                     H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                     error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                     spectra->elapsed_realtime(real_time);

                     H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                     error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                     spectra->input_counts(in_cnt);

                     H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                     error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                     spectra->output_counts(out_cnt);

                     for(size_t s=0; s<count_row[0]; s++)
                     {
                         (*spectra)[s] = detector_hid_map[detector_num].buffer[(count_row[1] * s) + col];
                     }
                     callback_func(row, col, dims_in[1], count_row[1], detector_num, spectra, user_data);
                }
             }

         }
         else
         {
            logE<<"reading row "<<row<<"\n";
         }
    }

    delete [] dims_in;
    delete [] offset;
    delete [] count;

    for(size_t detector_num : detector_num_arr)
    {
        delete [] detector_hid_map[detector_num].buffer;
    }

    detector_hid_map.clear();

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

	return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_and_integrate_spectra_volume(std::string path, size_t detector_num, data_struct::Spectra* spectra)
{
    std::lock_guard<std::mutex> lock(_mutex);


    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    logI<< path <<" detector : "<<detector_num<<"\n";

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

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id) )
       return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id) )
       return false;
    dataspace_rt_id = H5Dget_space(dset_rt_id);
    close_map.push({dataspace_rt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id) )
       return false;
    dataspace_inct_id = H5Dget_space(dset_incnt_id);
    close_map.push({dataspace_inct_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id) )
       return false;
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);
    close_map.push({dataspace_outct_id, H5O_DATASPACE});


     int rank = H5Sget_simple_extent_ndims(dataspace_id);
     if (rank != 3)
     {
         _close_h5_objects(close_map);
         logE<<"Dataset /MAPS_RAW/"<<detector_path<<" rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<"\n";
         return false;
        //throw exception ("Dataset is not a volume");
     }
     hsize_t* dims_in = new hsize_t[rank];
     hsize_t* offset = new hsize_t[rank];
     hsize_t* count = new hsize_t[rank];
    
     int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
     if(status_n < 0)
     {
         _close_h5_objects(close_map);
          logE<<"getting dataset rank for MAPS_RAW/"<< detector_path<<"\n";
          return false;
     }

     for (int i=0; i < rank; i++)
     {
        
        offset[i] = 0;
        count[i] = dims_in[i];
     }

     buffer = new real_t [dims_in[0] * dims_in[2]]; // spectra_size x cols
     count_row[0] = dims_in[0];
     count_row[1] = dims_in[2];

     count[1] = 1; //1 row

     if((hsize_t)spectra->size() != dims_in[0])
     {
        spectra->resize(dims_in[0]);
        spectra->setZero(dims_in[0]);
     }

     memoryspace_id = H5Screate_simple(2, count_row, nullptr);
     memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
     H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
     H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

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

          H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
          error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

          if (error > -1 )
          {
              for(size_t col=0; col<count_row[1]; col++)
              {
                  offset_meta[2] = col;

                  H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                  error = H5Dread(dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                  live_time_total += live_time;

                  H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                  error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                  real_time_total += real_time;

                  H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                  error = H5Dread(dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                  in_cnt_total += in_cnt;

                  H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                  error = H5Dread(dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                  out_cnt_total += out_cnt;

                  for(size_t s=0; s<count_row[0]; s++)
                  {
                      (*spectra)[s] += buffer[(count_row[1] * s) + col];
                  }
              }

             //logI<<"read row "<<row<<"\n";
          }
          else
          {
             logE<<"reading row "<<row<<"\n";
          }
     }

     spectra->elapsed_livetime(live_time_total);
     spectra->elapsed_realtime(real_time_total);
     spectra->input_counts(in_cnt_total);
     spectra->output_counts(out_cnt_total);

     delete [] dims_in;
     delete [] offset;
     delete [] count;
     delete [] buffer;

     _close_h5_objects(close_map);

     end = std::chrono::system_clock::now();
     std::chrono::duration<double> elapsed_seconds = end-start;
     //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

     logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

     return true;
}

//-----------------------------------------------------------------------------


bool HDF5_IO::load_spectra_vol_analyzed_h5(std::string path,
                                           data_struct::Spectra_Volume* spectra_volume,
                                           int row_idx_start,
                                           int row_idx_end,
                                           int col_idx_start,
                                           int col_idx_end)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   logI<< path <<"\n";

   hid_t    file_id, dset_id, dataspace_id, spec_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
   hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
   herr_t   error;
   hsize_t dims_in[3] = {0,0,0};
   hsize_t offset[3] = {0,0,0};
   hsize_t count[3] = {1,1,1};
   hsize_t offset_time[2] = {0,0};
   hsize_t count_time[2] = {1,1};

   if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1, false) )
       return false;

   if ( false == _open_h5_object(spec_grp_id, H5O_GROUP, close_map, "/MAPS/Spectra", file_id) )
      return false;

   if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, "mca_arr", spec_grp_id) )
      return false;
   dataspace_id = H5Dget_space(dset_id);
   close_map.push({dataspace_id, H5O_DATASPACE});

   if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "Elapsed_Livetime", spec_grp_id) )
      return false;
   dataspace_lt_id = H5Dget_space(dset_lt_id);
   close_map.push({dataspace_lt_id, H5O_DATASPACE});

   if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "Elapsed_Realtime", spec_grp_id) )
      return false;
   dataspace_rt_id = H5Dget_space(dset_rt_id);
   close_map.push({dataspace_rt_id, H5O_DATASPACE});

   if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "Input_Counts", spec_grp_id) )
      return false;
   dataspace_inct_id = H5Dget_space(dset_incnt_id);
   close_map.push({dataspace_inct_id, H5O_DATASPACE});

   if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "Output_Counts", spec_grp_id) )
      return false;
   dataspace_outct_id = H5Dget_space(dset_outcnt_id);
   close_map.push({dataspace_outct_id, H5O_DATASPACE});


    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 3)
    {
        _close_h5_objects(close_map);
        logE<<"Dataset /MAPS/Spectra/mca_arr  rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<"\n";
        return false;
       //throw exception ("Dataset is not a volume");
    }
   
    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logE<<"getting dataset dims for /MAPS/Spectra/mca_arr"<<"\n";
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    if(row_idx_end < row_idx_start)
    {
        row_idx_end = dims_in[1];
    }

    if(col_idx_end < col_idx_start)
    {
        col_idx_end = dims_in[2];
    }

    spectra_volume->resize_and_zero(dims_in[1], dims_in[2], dims_in[0]);

    //buffer = new real_t [dims_in[0] * dims_in[2]]; // cols x spectra_size
    count[0] = dims_in[0];
    count[1] = 1;
    count[2] = 1;

    memoryspace_id = H5Screate_simple(3, count, nullptr);
    memoryspace_meta_id = H5Screate_simple(2, count_time, nullptr);
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    //offset[1] = row;

    for(size_t row=(size_t)row_idx_start; row < (size_t)row_idx_end; row++)
    {
        offset[1] = row;
        offset_time[0] = row;
        for(size_t col=(size_t)col_idx_start; col < (size_t)col_idx_end; col++)
        {
            data_struct::Spectra *spectra = &((*spectra_volume)[row][col]);
            offset[2] = col;
            offset_time[1] = col;
            H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

            //error = H5Dread (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);
			error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)spectra->data());
			if (error > 0)
			{
				logW << "Counld not read row " << row << " col " << col << "\n";
			}

            H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
            H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
            H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
            H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

            error = H5Dread (dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, (void*)&real_time);
            error = H5Dread (dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, (void*)&live_time);
            error = H5Dread (dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, (void*)&in_cnt);
            error = H5Dread (dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, (void*)&out_cnt);

            spectra->elapsed_livetime(live_time);
            spectra->elapsed_realtime(real_time);
            spectra->input_counts(in_cnt);
            spectra->output_counts(out_cnt);
        }
    }

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    logI << "elapsed time: " << elapsed_seconds.count() << "s\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_quantification_scalers_analyzed_h5(std::string path,
                                                      data_struct::Params_Override *override_values)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   std::chrono::time_point<std::chrono::system_clock> start, end;
   start = std::chrono::system_clock::now();

   std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

   logI<< path <<"\n";

   real_t srcurrent, us_ic, ds_ic;
   hid_t    file_id, ds_ic_id, us_ic_id, srcurrent_id;
   //herr_t   error;
   //hsize_t offset[1] = {0};
   hsize_t count[1] = {1};
   hid_t readwrite_space = H5Screate_simple(1, &count[0], &count[0]);

   if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
       return false;

   if ( false == _open_h5_object(ds_ic_id, H5O_DATASET, close_map, "/MAPS/Quantification/Scalers/DS_IC", file_id) )
      return false;

   if ( false == _open_h5_object(us_ic_id, H5O_DATASET, close_map, "/MAPS/Quantification/Scalers/US_IC", file_id) )
      return false;

   if ( false == _open_h5_object(srcurrent_id, H5O_DATASET, close_map, "/MAPS/Quantification/Scalers/SR_Current", file_id) )
      return false;


   //read in scaler
   hid_t d_space = H5Dget_space(srcurrent_id);
   hid_t d_type = H5Dget_type(srcurrent_id);
   hid_t status = H5Dread(srcurrent_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)&srcurrent);
   if(status > -1)
   {
        override_values->sr_current = (srcurrent);
   }
   d_space = H5Dget_space(us_ic_id);
   d_type = H5Dget_type(us_ic_id);
   status = H5Dread(us_ic_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)&us_ic);
   if(status > -1)
   {
        override_values->US_IC = (us_ic);
   }
   d_space = H5Dget_space(ds_ic_id);
   d_type = H5Dget_type(ds_ic_id);
   status = H5Dread(ds_ic_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)&ds_ic);
   if(status > -1)
   {
        override_values->DS_IC = (ds_ic);
   }

   _close_h5_objects(close_map);

   end = std::chrono::system_clock::now();
   std::chrono::duration<double> elapsed_seconds = end-start;
   //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

   logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

   return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_quantification_scalers_gsecars(std::string path, data_struct::Params_Override *override_values)
{
    std::lock_guard<std::mutex> lock(_mutex);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    int status = 0;
    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    logI<< path <<"\n";

    real_t us_ic, ds_ic;
    hid_t file_id, ds_ic_id, us_ic_id;
    hsize_t offset3[3] = { 0,0,0 };
    hsize_t count3[3] = { 1,1,1 };
    hsize_t dims3[3] = { 1,1,1 };
    hsize_t offset[1] = {1};
    hsize_t count[1] = {1};
    hid_t readwrite_space = H5Screate_simple(1, &count[0], &count[0]);

    GSE_CARS_SAVE_VER version = GSE_CARS_SAVE_VER::UNKNOWN;

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(ds_ic_id, H5O_DATASET, close_map, "/xrmmap/scalars/I1_raw", file_id, false, false) )
    {
        if ( false == _open_h5_object(ds_ic_id, H5O_DATASET, close_map, "/xrfmap/roimap/det_name", file_id) )
        {
            return false;
        }
        version = GSE_CARS_SAVE_VER::XRFMAP;
    }
    else
    {
        version = GSE_CARS_SAVE_VER::XRMMAP;
    }

    if ( false == _open_h5_object(us_ic_id, H5O_DATASET, close_map, "/xrmmap/scalars/I0_raw", file_id, false, false) )
    {
        if ( false == _open_h5_object(us_ic_id, H5O_DATASET, close_map, "/xrfmap/roimap/det_raw", file_id) )
        {
            return false;
        }
    }

    //read in scaler
    hid_t d_space = H5Dget_space(us_ic_id);
    //hid_t d_type = H5Dget_type(us_ic_id);
    status = H5Sget_simple_extent_dims(d_space, &dims3[0], nullptr);
    if(status < 0)
    {
        _close_h5_objects(close_map);
         logE<<"getting dataset dims"<<"\n";
         return false;
    }

    if(version == GSE_CARS_SAVE_VER::XRMMAP)
    {
        for(offset3[0] = 0; offset3[0] < dims3[0]; offset3[0]++)
        {
            for(offset3[1] = 0; offset3[1] < dims3[1]; offset3[1]++)
            {
                H5Sselect_hyperslab (d_space, H5S_SELECT_SET, offset3, nullptr, count3, nullptr);

                status = H5Dread(us_ic_id, H5T_NATIVE_REAL, readwrite_space, d_space, H5P_DEFAULT, (void*)&us_ic);
                if(status > -1)
                {
                    override_values->US_IC += (us_ic);
                }

                status = H5Dread(ds_ic_id, H5T_NATIVE_REAL, readwrite_space, d_space, H5P_DEFAULT, (void*)&ds_ic);
                if(status > -1)
                {
                    override_values->DS_IC += (ds_ic);
                }
            }
        }
    }
    else if(version == GSE_CARS_SAVE_VER::XRFMAP)
    {
        int usIDX = -1;
        int dsIDX = -1;

        hid_t dtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(dtype, 255);

        hid_t name_space = H5Dget_space(ds_ic_id);
        status = H5Sget_simple_extent_dims(name_space, &count[0], nullptr);
        int name_amt = count[0];
        count[0] = 1;

        for(offset[0] = 0; offset[0] < name_amt; offset[0]++)
        {
            char tmp_char[256] = {0};
            //read the detector names and find I0 and I1 indicies
            H5Sselect_hyperslab (name_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            status = H5Dread(ds_ic_id, dtype, readwrite_space, name_space, H5P_DEFAULT, (void*)tmp_char);
            if (status > -1)
            {
                std::string sname(tmp_char);
                if (sname == "I0")
                {
                    usIDX = offset[0];
                }
                if (sname == "I1")
                {
                    dsIDX = offset[0];
                }
                if(usIDX > -1 && dsIDX > -1)
                {
                    break;
                }
            }
        }

        if(usIDX == -1 || dsIDX == -1)
        {
            logW<<"Could not find up stream ion or down stream ion chamber indicies\n";
        }

        for(offset3[0] = 0; offset3[0] < dims3[0]; offset3[0]++)
        {
            for(offset3[1] = 0; offset3[1] < dims3[1]; offset3[1]++)
            {
                offset3[2] = usIDX;
                H5Sselect_hyperslab (d_space, H5S_SELECT_SET, offset3, nullptr, count3, nullptr);

                status = H5Dread(us_ic_id, H5T_NATIVE_REAL, readwrite_space, d_space, H5P_DEFAULT, (void*)&us_ic);
                if(status > -1)
                {
                    override_values->US_IC += (us_ic);
                }

                offset3[2] = dsIDX;
                H5Sselect_hyperslab (d_space, H5S_SELECT_SET, offset3, nullptr, count3, nullptr);
                status = H5Dread(us_ic_id, H5T_NATIVE_REAL, readwrite_space, d_space, H5P_DEFAULT, (void*)&ds_ic);
                if(status > -1)
                {
                    override_values->DS_IC += (ds_ic);
                }
            }
        }
    }

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_quantification_scalers_BNL(std::string path, data_struct::Params_Override* override_values)
{
    std::lock_guard<std::mutex> lock(_mutex);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    logI << path << "\n";

    hid_t scan_grp_id, maps_grp_id, scalers_grp_id, status, error;
    hid_t    file_id, src_maps_grp_id;
    hid_t    dataspace_detectors_id, dset_detectors_id;
    hid_t   xypos_dataspace_id, xypos_id;
    hid_t x_dataspace_id, y_dataspace_id, x_dataset_id, y_dataset_id;
    char* detector_names[256];
    int det_rank;
    hsize_t* det_dims_in = nullptr;
    hsize_t* val_dims_in = nullptr;;
    hsize_t scaler_offset[3] = { 0,0,0 };
    hsize_t single_offset[1] = { 0 };
    hsize_t single_count[1] = { 1 };
    hsize_t mem_offset[2] = { 0,0 };
    hsize_t mem_count[2] = { 1,1 };
    hid_t ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
    double* buffer = nullptr;

    std::string ds_ic_search = "";

    if (override_values == nullptr)
    {
        return false;
    }

    if (override_values->scaler_pvs.count(STR_DS_IC) == 0)
    {
        logW << "Need to set " << STR_DS_IC << ":scaler_name in maps_fit_parameter_override.txt\n";
        return false;
    }
    else
    {
        ds_ic_search = override_values->scaler_pvs.at(STR_DS_IC);
    }

    
    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
        return false;

    if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "xrfmap", file_id))
    {
        return false;
    }

    //Save scalers
    // name, val
    hid_t scaler_name_id, scaler_val_id, scaler_grp_id;
    if (false == _open_h5_object(scaler_grp_id, H5O_GROUP, close_map, "scalers", src_maps_grp_id))
    {
        return false;
    }
    if (false == _open_h5_object(scaler_name_id, H5O_DATASET, close_map, "name", scaler_grp_id))
    {
        return false;
    }
    if (false == _open_h5_object(scaler_val_id, H5O_DATASET, close_map, "val", scaler_grp_id))
    {
        return false;
    }

    hid_t name_type = H5Dget_type(scaler_name_id);
    hid_t scaler_type = H5Dget_type(scaler_val_id);
    hid_t scaler_val_space = H5Dget_space(scaler_val_id);
    hid_t val_rank = H5Sget_simple_extent_ndims(scaler_val_space);
    val_dims_in = new hsize_t[val_rank];
    H5Sget_simple_extent_dims(scaler_val_space, &val_dims_in[0], NULL);

    mem_count[0] = val_dims_in[0];
    mem_count[1] = val_dims_in[1];
    //names
    hid_t mem_single_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
    single_count[0] = val_dims_in[2];
    hid_t name_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
    single_count[0] = 1;

    char tmp_char[255] = { 0 };
    buffer = new double[val_dims_in[0] * val_dims_in[1]];
    hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);
    close_map.push({ mem_space, H5O_DATASPACE });
    size_t scaler_cnt = val_dims_in[2];
    val_dims_in[2] = 1;
    for (hsize_t i = 0; i < scaler_cnt; i++)
    {
        single_offset[0] = i;
        scaler_offset[2] = i;

        H5Sselect_hyperslab(scaler_val_space, H5S_SELECT_SET, scaler_offset, NULL, val_dims_in, NULL);
        H5Sselect_hyperslab(name_space, H5S_SELECT_SET, single_offset, nullptr, single_count, nullptr);

        status = H5Dread(scaler_name_id, name_type, mem_single_space, name_space, H5P_DEFAULT, (void*)tmp_char);
        if (status > -1)
        {
            std::string read_name = std::string(tmp_char, 255);
            read_name.erase(std::remove_if(read_name.begin(), read_name.end(), ::isspace), read_name.end());
            read_name.erase(std::find(read_name.begin(), read_name.end(), '\0'), read_name.end());
            if (read_name == ds_ic_search)
            {
                status = H5Dread(scaler_val_id, scaler_type, mem_space, scaler_val_space, H5P_DEFAULT, (void*)buffer);
                override_values->DS_IC = 0.0;
                if (status > -1)
                {
                    for (hsize_t x = 0; x < val_dims_in[0] * val_dims_in[1]; x++)
                    {
                        override_values->DS_IC += (real_t)buffer[x];
                    }
                    
                }
                break;
            }
        }
    }

    _close_h5_objects(close_map);

    if (det_dims_in != nullptr)
        delete[] det_dims_in;
    if (val_dims_in != nullptr)
        delete[] val_dims_in;
    if (buffer != nullptr)
        delete[] buffer;

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::load_integrated_spectra_analyzed_h5(std::string path,
                                           data_struct::Spectra* spectra, bool log_error)
{
    std::lock_guard<std::mutex> lock(_mutex);

   //_is_loaded = ERROR_LOADING;
   //std::chrono::time_point<std::chrono::system_clock> start, end;
   //start = std::chrono::system_clock::now();
/*
   
   if (detector_num > -1)
   {
	   path += ".h5" + std::to_string(detector_num);
   }
 */
   hid_t    file_id;

    file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if(file_id < 0)
    {
		if (log_error)
		{
			logE << "opening file " << path << "\n";
		}
        return false;
    }

    return _load_integrated_spectra_analyzed_h5(file_id, spectra);
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_load_integrated_spectra_analyzed_h5(hid_t file_id, data_struct::Spectra* spectra)
{

    hid_t    dset_id, dataspace_id, spec_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
    hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
    herr_t   error;
    hsize_t dims_in[1] = {0};
    hsize_t offset[1] = {0};
    hsize_t count[1] = {1};
    hsize_t offset_time[1] = {0};
    hsize_t count_time[1] = {1};

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    close_map.push({file_id, H5O_FILE});

    if ( false == _open_h5_object(spec_grp_id, H5O_GROUP, close_map, "/MAPS/Spectra/Integrated_Spectra", file_id) )
       return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, "Spectra", spec_grp_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "Elapsed_Livetime", spec_grp_id) )
       return false;
    dataspace_lt_id = H5Dget_space(dset_lt_id);
    close_map.push({dataspace_lt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "Elapsed_Realtime", spec_grp_id) )
       return false;
    dataspace_rt_id = H5Dget_space(dset_rt_id);
    close_map.push({dataspace_rt_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "Input_Counts", spec_grp_id) )
       return false;
    dataspace_inct_id = H5Dget_space(dset_incnt_id);
    close_map.push({dataspace_inct_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "Output_Counts", spec_grp_id) )
       return false;
    dataspace_outct_id = H5Dget_space(dset_outcnt_id);
    close_map.push({dataspace_outct_id, H5O_DATASPACE});


    int rank = H5Sget_simple_extent_ndims(dataspace_id);
    if (rank != 1)
    {
        _close_h5_objects(close_map);
        logE<<"Dataset /MAPS/Spectra/mca_arr  rank != 3. rank = "<<rank<<". Can't load dataset. returning"<<"\n";
        return false;
       //throw exception ("Dataset is not a volume");
    }
   
    int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
    if(status_n < 0)
    {
        _close_h5_objects(close_map);
         logE<<"getting dataset dims for /MAPS/Spectra/mca_arr"<<"\n";
         return false;
    }

    for (int i=0; i < rank; i++)
    {
       
       offset[i] = 0;
       count[i] = dims_in[i];
    }

    spectra->resize(dims_in[0]);
    spectra->setZero(dims_in[0]);

    count[0] = dims_in[0];

    memoryspace_id = H5Screate_simple(1, count, nullptr);
    memoryspace_meta_id = H5Screate_simple(1, count_time, nullptr);
    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

    real_t live_time = 1.0;
    real_t real_time = 1.0;
    real_t in_cnt = 1.0;
    real_t out_cnt = 1.0;

    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    error = H5Dread (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);

    H5Sselect_hyperslab (dataspace_lt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
    H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
    H5Sselect_hyperslab (dataspace_inct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
    H5Sselect_hyperslab (dataspace_outct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

    error = H5Dread (dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, (void*)&real_time);
    error = H5Dread (dset_lt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, (void*)&live_time);
    error = H5Dread (dset_incnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, (void*)&in_cnt);
    error = H5Dread (dset_outcnt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, (void*)&out_cnt);

    spectra->elapsed_livetime(live_time);
    spectra->elapsed_realtime(real_time);
    spectra->input_counts(in_cnt);
    spectra->output_counts(out_cnt);

    _close_h5_objects(close_map);

    //end = std::chrono::system_clock::now();
    //std::chrono::duration<double> elapsed_seconds = end-start;
    ////std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    //logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::get_scalers_and_metadata_emd(std::string path, data_struct::Scan_Info* scan_info)
{
    std::lock_guard<std::mutex> lock(_mutex);
    hid_t    file_id, src_maps_grp_id, detectors_grp_id, hash_grp_id, data_id, spectrumstream_grp_id, hash2_grp_id, data2_id;
    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    if (scan_info == nullptr)
    {
        return false;
    }

    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
    {
        return false;
    }

    if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "Data", file_id))
    {
        return false;
    }

    if (_open_h5_object(detectors_grp_id, H5O_GROUP, close_map, "Image", src_maps_grp_id, false, false))
    {
        herr_t err;
        ssize_t len;
        hsize_t nobj;
        int otype;
        char group_name[1024];
        char memb_name[1024];


        err = H5Gget_num_objs(detectors_grp_id, &nobj);
        if(nobj > 0)
        {
            len = H5Gget_objname_by_idx(detectors_grp_id, (hsize_t)0, memb_name, (size_t)1024);
            //printf("   %d ", len); fflush(stdout);
            //printf("  Member: %s ", memb_name); fflush(stdout);
            otype = H5Gget_objtype_by_idx(detectors_grp_id, (size_t)0);
            if (otype == H5G_GROUP)
            {
                if (_open_h5_object(hash_grp_id, H5O_GROUP, close_map, memb_name, detectors_grp_id, false, false))
                {
                    if (_open_h5_object(data_id, H5O_DATASET, close_map, "Data", hash_grp_id, false, false))
                    {
                        hid_t dspace_id = H5Dget_space(data_id);
                        close_map.push({ dspace_id, H5O_DATASPACE });
                    }
                }
            }
        }
    }

    if (_open_h5_object(spectrumstream_grp_id, H5O_GROUP, close_map, "SpectrumStream", src_maps_grp_id, false, false))
    {
        herr_t err;
        ssize_t len;
        hsize_t nobj;
        int otype;
        hsize_t names_cnt[1] = { 1 };
        char group_name[1024];
        char memb_name[1024];

        char acqui_data[10240];


        err = H5Gget_num_objs(spectrumstream_grp_id, &nobj);
        if (nobj > 0)
        {
            len = H5Gget_objname_by_idx(spectrumstream_grp_id, (hsize_t)0, memb_name, (size_t)1024);
            //printf("   %d ", len); fflush(stdout);
            //printf("  Member: %s ", memb_name); fflush(stdout);
            otype = H5Gget_objtype_by_idx(spectrumstream_grp_id, (size_t)0);
            if (otype == H5G_GROUP)
            {
                if (_open_h5_object(hash2_grp_id, H5O_GROUP, close_map, memb_name, spectrumstream_grp_id, false, false))
                {
                    if (_open_h5_object(data2_id, H5O_DATASET, close_map, "AcquisitionSettings", hash2_grp_id, false, false))
                    {
                        hid_t name_type = H5Tcopy(H5T_C_S1);
                        hid_t status = H5Tset_size(name_type, 255);

                        hid_t dspace_id = H5Dget_space(data2_id);
                        close_map.push({ dspace_id, H5O_DATASPACE });

                        hid_t mem_name_space = H5Screate_simple(1, &names_cnt[0], &names_cnt[0]);
                        close_map.push({ mem_name_space, H5O_DATASPACE });

                        hid_t scalers_type = H5Dget_type(data2_id);

                        err = H5Dread(data2_id, scalers_type, mem_name_space, dspace_id, H5P_DEFAULT, acqui_data);
                        if (err > -1)
                        {

                        }
                    }
                }
            }
        }
    }

    scan_info->meta_info.detectors.push_back(0);

    _close_h5_objects(close_map);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::get_scalers_and_metadata_confocal(std::string path, data_struct::Scan_Info* scan_info)
{
    std::lock_guard<std::mutex> lock(_mutex);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    logI << path << "\n";
    hid_t    file_id, src_maps_grp_id, detectors_grp_id, x_id, y_id;
    hid_t scaler_name_id, scaler_val_id, scaler_grp_id, attr_detector_names_id;
    hsize_t names_off[1] = { 0 };
    hsize_t names_cnt[1] = { 1 };
    hsize_t scalers_offset[3] = { 0,0,0 };
    hsize_t scalers_count[3] = { 1,1,1 };
    hsize_t value_offset[3] = { 0,0,0 };
    hsize_t value_count[3] = { 1,1,1 };
    hsize_t mem_count[2] = { 1,1 };
    hsize_t x_offset[3] = { 0,0,0 };
    hsize_t x_count[3] = { 1,1,1 };
    hsize_t y_offset[2] = { 0,0 };
    hsize_t y_count[2] = { 1,1 };
    char* detector_names[256];


    if (scan_info == nullptr)
    {
        return false;
    }

    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
    {
        return false;
    }
    
    if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "2D Scan", file_id))
    {
        return false;
    }


    //  post 2020 version
    if (_open_h5_object(detectors_grp_id, H5O_GROUP, close_map, "Detectors", src_maps_grp_id, false, false))
    {
        bool first_save = true;
        real_t* buffer = nullptr;
        hsize_t nobj = 0;
        H5Gget_num_objs(detectors_grp_id, &nobj);
        hid_t values_id;
        hid_t value_space;
        hid_t mem_space;
        hid_t names_id;
        hid_t name_space;
        hid_t mem_name_space = H5Screate_simple(1, &names_cnt[0], &names_cnt[0]);
        close_map.push({ mem_name_space, H5O_DATASPACE });
        hid_t name_type = H5Tcopy(H5T_C_S1);
        hid_t status = H5Tset_size(name_type, 255);


        for (hsize_t i = 0; i < nobj; i++)
        {
            char str_dset_name[2048] = { 0 };
            hid_t dsid;
            hsize_t len = H5Gget_objname_by_idx(detectors_grp_id, i, str_dset_name, 2048);
            if (_open_h5_object(dsid, H5O_DATASET, close_map, str_dset_name, detectors_grp_id))
            {
                hid_t scaler_space = H5Dget_space(dsid);
                close_map.push({ scaler_space, H5O_DATASPACE });
                hid_t scalers_type = H5Dget_type(dsid);

                if (first_save)
                {
                    H5Sget_simple_extent_dims(scaler_space, &scalers_count[0], NULL);
                    scan_info->meta_info.requested_cols = scalers_count[0];
                    scan_info->meta_info.requested_rows = scalers_count[1];
                    first_save = false;
                }
                Scaler_Map sm;
                sm.name = string(str_dset_name, len);
                sm.values.resize(scalers_count[0], scalers_count[1]);
                status = H5Dread(dsid, scalers_type, scaler_space, scaler_space, H5P_DEFAULT, sm.values.data());
                scan_info->scaler_maps.push_back(sm);
            }
        }
    }
    // pre 2020 version
    if (_open_h5_object(detectors_grp_id, H5O_DATASET, close_map, "Detectors", src_maps_grp_id, false, false))
    {
        hid_t scaler_space = H5Dget_space(detectors_grp_id);
        close_map.push({ scaler_space, H5O_DATASPACE });
        H5Sget_simple_extent_dims(scaler_space, &scalers_count[0], NULL);
        hid_t scalers_type = H5Dget_type(detectors_grp_id);
        hsize_t scaler_amt = scalers_count[2];
        scalers_count[2] = 1;
        mem_count[0] = scalers_count[0];
        mem_count[1] = scalers_count[1];
        scan_info->meta_info.requested_cols = scalers_count[0];
        scan_info->meta_info.requested_rows = scalers_count[1];
        hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);
        close_map.push({ mem_space, H5O_DATASPACE });

        for (hsize_t s = 0; s < scaler_amt; s++)
        {
            Scaler_Map sm;
            sm.values.resize(scalers_count[0], scalers_count[1]);
            scalers_offset[2] = s;
            H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, scalers_offset, nullptr, scalers_count, nullptr);
            H5Dread(detectors_grp_id, scalers_type, mem_space, scaler_space, H5P_DEFAULT, sm.values.data());
            scan_info->scaler_maps.push_back(sm);
        }

        //save the scalers names
        attr_detector_names_id = H5Aopen(detectors_grp_id, "Detector Names", H5P_DEFAULT);
        close_map.push({ attr_detector_names_id, H5O_ATTRIBUTE });

        hid_t ftype = H5Aget_type(attr_detector_names_id);
        hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
        hid_t error = H5Aread(attr_detector_names_id, type, detector_names);

        if (error == 0)
        {
            int i = 0;
            for (auto& itr : scan_info->scaler_maps)
            {
                itr.name = string(detector_names[i], strlen(detector_names[i]));
                i++;
            }
        }
    }

    if (_open_h5_object(x_id, H5O_DATASET, close_map, "MCA 1", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(0);
    }
    if (_open_h5_object(x_id, H5O_DATASET, close_map, "MCA 2", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(1);
    }
    if (_open_h5_object(x_id, H5O_DATASET, close_map, "MCA 3", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(2);
    }
    if (_open_h5_object(x_id, H5O_DATASET, close_map, "MCA 4", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(3);
    }

    /*
    if (_open_h5_object(x_id, H5O_DATASET, close_map, "X Positions", src_maps_grp_id, false, false))
    {
        
    }
    if (_open_h5_object(y_id, H5O_DATASET, close_map, "Y Positions", src_maps_grp_id, false, false))
    {
        
    }
    */
    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::get_scalers_and_metadata_gsecars(std::string path, data_struct::Scan_Info* scan_info)
{
    std::lock_guard<std::mutex> lock(_mutex);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    logI << path << "\n";
    hid_t    file_id, src_maps_grp_id, scalers_grp_id, pos_grp_id;
    hid_t config_grp_id, environ_grp_id, scaler_grp_id, name_id, value_id;

    hsize_t single_offset[1] = { 0 };
    hsize_t single_count[1] = { 1 };
    hsize_t* val_dims_in = nullptr;

    if (scan_info == nullptr)
    {
        return false;
    }

    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
    {
        return false;
    }
    
    if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "xrmmap", file_id, false, false))
    {
        if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "xrfmap", file_id))
        {
            return false;
        }
    }

    hid_t mem_single_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
    close_map.push({ mem_single_space, H5O_DATASPACE });

    if (_open_h5_object(config_grp_id, H5O_GROUP, close_map, "config", src_maps_grp_id, false, false))
    {
        if (_open_h5_object(environ_grp_id, H5O_GROUP, close_map, "environ", config_grp_id, false, false))
        {
            if (_open_h5_object(name_id, H5O_DATASET, close_map, "name", environ_grp_id, false, false))
            {
                if (_open_h5_object(value_id, H5O_DATASET, close_map, "value", environ_grp_id, false, false))
                {
                    hid_t name_type = H5Dget_type(name_id);
                    hid_t value_type = H5Dget_type(value_id);
                    hid_t value_space = H5Dget_space(value_id);
                    hid_t val_rank = H5Sget_simple_extent_ndims(value_space);
                    val_dims_in = new hsize_t[val_rank];
                    H5Sget_simple_extent_dims(value_space, &val_dims_in[0], NULL);

                    hsize_t amt = val_dims_in[0];
                    //for (hsize_t i = 0; i < amt; i++)
                    //{
                        //Extra_PV epv;
                        //scan_info->extra_pvs.push_back(epv);
                    //}
                }
            }
        }
    }

    //if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "positions", src_maps_grp_id, false, false))
    //{

    //}

    if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "mca1", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(0);
    }
    if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "mca2", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(1);
    }
    if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "mca3", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(2);
    }
    if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "mca4", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(3);
    }

    if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "det1", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(0);
    }
    if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "det2", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(1);
    }
    if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "det3", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(2);
    }
    if (_open_h5_object(pos_grp_id, H5O_GROUP, close_map, "det4", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(3);
    }

    if (_open_h5_object(scalers_grp_id, H5O_GROUP, close_map, "scalers", src_maps_grp_id, false, false))
    {
        /*
        hsize_t amt = val_dims_in[0];
        for (hsize_t i = 0; i < amt; i++)
        {
            Scaler_Map sm;
            scan_info->scaler_maps.push_back(sm);
        }
        */
    }

    _close_h5_objects(close_map);

    
    if (val_dims_in != nullptr)
        delete[] val_dims_in;

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::get_scalers_and_metadata_bnl(std::string path, data_struct::Scan_Info* scan_info)
{
    std::lock_guard<std::mutex> lock(_mutex);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    logI << path << "\n";
    char* adata[255];
    hid_t    file_id, src_maps_grp_id;
    hid_t scaler_name_id, scaler_val_id, scaler_grp_id;
    hid_t meta_data_id, tmp_id, status;
    hsize_t* val_dims_in = nullptr;
    hsize_t mem_offset[2] = { 0,0 };
    hsize_t mem_count[2] = { 1,1 };
    hsize_t single_offset[1] = { 0 };
    hsize_t single_count[1] = { 1 };
    hsize_t scaler_offset[3] = { 0,0,0 };

    if (scan_info == nullptr)
    {
        return false;
    }

    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
    {
        return false;
    }
    if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "xrfmap", file_id))
    {
        return false;
    }
    if (false == _open_h5_object(scaler_grp_id, H5O_GROUP, close_map, "scalers", src_maps_grp_id))
    {
        return false;
    }
    if (false == _open_h5_object(scaler_name_id, H5O_DATASET, close_map, "name", scaler_grp_id))
    {
        return false;
    }
    if (false == _open_h5_object(scaler_val_id, H5O_DATASET, close_map, "val", scaler_grp_id))
    {
        return false;
    }
    
    hid_t name_type = H5Dget_type(scaler_name_id);
    hid_t scaler_type = H5Dget_type(scaler_val_id);
    hid_t scaler_val_space = H5Dget_space(scaler_val_id);
    hid_t val_rank = H5Sget_simple_extent_ndims(scaler_val_space);
    val_dims_in = new hsize_t[val_rank];
    H5Sget_simple_extent_dims(scaler_val_space, &val_dims_in[0], NULL);

    mem_count[0] = val_dims_in[0];
    mem_count[1] = val_dims_in[1];

    scan_info->meta_info.requested_cols = val_dims_in[0];
    scan_info->meta_info.requested_rows = val_dims_in[1];
    //names
    hid_t mem_single_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
    close_map.push({ mem_single_space, H5O_DATASPACE });
    single_count[0] = val_dims_in[2];
    hid_t name_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
    close_map.push({ name_space, H5O_DATASPACE });
    single_count[0] = 1;

    char tmp_char[255] = { 0 };
    hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);
    close_map.push({ mem_space, H5O_DATASPACE });
    size_t scaler_cnt = val_dims_in[2];
    val_dims_in[2] = 1;
    for (hsize_t i = 0; i < scaler_cnt; i++)
    {
        data_struct::Scaler_Map scaler_map;
        scaler_map.values.resize(val_dims_in[1], val_dims_in[0]);
        scaler_map.unit = "cts";

        single_offset[0] = i;
        scaler_offset[2] = i;

        H5Sselect_hyperslab(scaler_val_space, H5S_SELECT_SET, scaler_offset, NULL, val_dims_in, NULL);
        H5Sselect_hyperslab(name_space, H5S_SELECT_SET, single_offset, nullptr, single_count, nullptr);

        status = H5Dread(scaler_name_id, name_type, mem_single_space, name_space, H5P_DEFAULT, (void*)tmp_char);
        if (status > -1)
        {
            scaler_map.name = std::string(tmp_char, 255);
            scaler_map.name.erase(std::remove_if(scaler_map.name.begin(), scaler_map.name.end(), ::isspace), scaler_map.name.end());
            scaler_map.name.erase(std::find(scaler_map.name.begin(), scaler_map.name.end(), '\0'), scaler_map.name.end());
        }
        status = H5Dread(scaler_val_id, H5T_NATIVE_REAL, mem_space, scaler_val_space, H5P_DEFAULT, scaler_map.values.data());

        scan_info->scaler_maps.push_back(scaler_map);
    }


    if (_open_h5_object(tmp_id, H5O_GROUP, close_map, "det1", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(0);
    }
    if (_open_h5_object(tmp_id, H5O_GROUP, close_map, "det2", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(1);
    }
    if (_open_h5_object(tmp_id, H5O_GROUP, close_map, "det3", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(2);
    }
    if (_open_h5_object(tmp_id, H5O_GROUP, close_map, "det4", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(3);
    }
    if (_open_h5_object(tmp_id, H5O_GROUP, close_map, "detsum", src_maps_grp_id, false, false))
    {
        scan_info->meta_info.detectors.push_back(-1);
    }
    
    if(_open_h5_object(tmp_id, H5O_GROUP, close_map, "scan_metadata", src_maps_grp_id))
    {
    
        // load attributes from this folder
        int na = H5Aget_num_attrs(tmp_id);

        for (int i = 0; i < na; i++) 
        {
            hid_t aid = H5Aopen_idx(tmp_id, (unsigned int)i);
            hid_t atype;
            hid_t aspace;
            char buf[1000];
            ssize_t len = H5Aget_name(aid, 1000, buf);
            data_struct::Extra_PV e_pv;
            e_pv.name = std::string(buf, len);
      
            atype = H5Aget_type(aid);
            if(H5Tis_variable_str(atype) > 0)
            {
                //size_t alen = H5Tget_size(atype);
                hid_t type = H5Tget_native_type(atype, H5T_DIR_ASCEND);
                if (H5Aread(aid, type, &adata) > -1)
                {
                    //e_pv.value = std::string(adata[0], alen);
                    e_pv.value = std::string(adata[0]);
                }
                
            }

            scan_info->extra_pvs.push_back(e_pv);

            H5Tclose(atype);
            H5Aclose(aid);
        }
        
    }

    _close_h5_objects(close_map);

    //free(adata[0]);

    if (val_dims_in != nullptr)
        delete[] val_dims_in;

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::start_save_seq(const std::string filename, bool force_new_file)
{

    if (_cur_file_id > -1)
    {
        logI<<" file already open, calling close() before opening new file. "<<"\n";
        end_save_seq();
    }

    if(false == force_new_file)
        _cur_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    if (_cur_file_id < 1)
    {
        logI<<"Creating file "<<filename<<"\n";
        _cur_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }
    if(_cur_file_id < 0)
    {
        logE<<"opening file "<<filename<<"\n";
        return false;
    }
	_cur_filename = filename;

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::end_save_seq(bool loginfo)
{

    if(_cur_file_id > 0)
    {
		if (loginfo)
		{
			logI << "closing file " << _cur_filename << "\n";
		}
        H5Fflush(_cur_file_id, H5F_SCOPE_LOCAL);

        ssize_t obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_DATASET | H5F_OBJ_LOCAL );
        if(obj_cnt > 0)
        {
            logI<<"closing forgotten datasets: "<<obj_cnt<<"\n";
            hid_t* objs = new hid_t[obj_cnt];
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_DATASET, -1, objs ) > -1)
            {
                for(int i=0; i<obj_cnt; i++)
                {
                    H5Dclose(objs[i]);
                }
            }
            delete [] objs;
        }
        obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_GROUP | H5F_OBJ_LOCAL );
        if(obj_cnt > 0)
        {
            logI<<"closing forgotten groups: "<<obj_cnt<<"\n";
            hid_t* objs = new hid_t[obj_cnt];
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_GROUP, -1, objs ) > -1)
            {
                for(int i=0; i<obj_cnt; i++)
                {
                    H5Gclose(objs[i]);
                }
            }
            delete [] objs;
        }
        obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_DATATYPE | H5F_OBJ_LOCAL );
        if(obj_cnt > 0)
        {
            logI<<"closing forgotten datatypes: "<<obj_cnt<<"\n";
            hid_t* objs = new hid_t[obj_cnt];
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_DATATYPE, -1, objs ) > -1)
            {
                for(int i=0; i<obj_cnt; i++)
                {
                    H5Tclose(objs[i]);
                }
            }
            delete [] objs;
        }
        obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_ATTR | H5F_OBJ_LOCAL );
        if(obj_cnt > 0)
        {
            logI<<"closing forgotten attributes: "<<obj_cnt<<"\n";
            hid_t* objs = new hid_t[obj_cnt];
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_ATTR, -1, objs ) > -1)
            {
                for(int i=0; i<obj_cnt; i++)
                {
                    H5Aclose(objs[i]);
                }
            }
            delete [] objs;
        }
        obj_cnt = H5Fget_obj_count( _cur_file_id, H5F_OBJ_ALL | H5F_OBJ_LOCAL );
        if(obj_cnt > 1) //file is still open
        {
            logW<<"**** did not close total objects "<<obj_cnt<<"\n";
        }


        H5Fclose(_cur_file_id);
        _cur_file_id = -1;
    }
    else
    {
        logW<<" could not close file because none is open"<<"\n";
        return false;
    }
    _cur_filename = "";
    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_spectra_volume(const std::string path,
                                  data_struct::Spectra_Volume * spectra_volume,
                                  real_t energy_offset,
                                  real_t energy_slope,
                                  real_t energy_quad,
                                  size_t row_idx_start,
                                  int row_idx_end,
                                  size_t col_idx_start,
                                  int col_idx_end)
{
    std::lock_guard<std::mutex> lock(_mutex);

    

    if(_cur_file_id < 0)
    {
        logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
        return false;
    }

//herr_t (*old_func)(void*);
//void *old_client_data;
//hid_t error_stack = H5Eget_current_stack();
//H5Eget_auto2(error_stack, &old_func, &old_client_data);

//H5Eset_auto2(error_stack, nullptr, nullptr);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    dset_id, spec_grp_id, int_spec_grp_id, dataspace_id, memoryspace_id, memoryspace_time_id, dataspace_time_id, file_time_id, maps_grp_id, dcpl_id;
    hid_t   dset_rt_id, dset_lt_id, incnt_dset_id, outcnt_dset_id;
    //herr_t   error;
    int status = 0;

    hsize_t chunk_dims[3];
    hsize_t dims_out[3];
	hsize_t maxdims[3] = { H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED };
    hsize_t offset[3];
    hsize_t count[3];
    hsize_t dims_time_out[2];
    hsize_t offset_time[2];
    hsize_t count_time[2];
	hsize_t tmp_dims[3];

    if (row_idx_end < (int)row_idx_start || (size_t)row_idx_end > spectra_volume->rows() -1)
    {
        row_idx_end = spectra_volume->rows();
    }
    if (col_idx_end < (int)col_idx_start || (size_t)col_idx_end > spectra_volume->cols() -1)
    {
        col_idx_end = spectra_volume->cols();
    }

    //get one element
    //data_struct::Fit_Element_Map* element;

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


    memoryspace_id = H5Screate_simple(3, count, nullptr);

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    memoryspace_time_id = H5Screate_simple(2, count_time, nullptr);
    file_time_id = H5Screate_simple(2, dims_time_out, nullptr);
    dataspace_time_id = H5Screate_simple (2, dims_time_out, nullptr);

    H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        logE<<"creating group MAPS"<<"\n";
        return false;
    }

    spec_grp_id = H5Gopen(maps_grp_id, "Spectra", H5P_DEFAULT);
    if(spec_grp_id < 0)
        spec_grp_id = H5Gcreate(maps_grp_id, "Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(spec_grp_id < 0)
    {
        logE<<"creating group MAPS/Spectra"<<"\n";
        return false;
    }
/*
    scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
    if(scalers_grp_id < 0)
        scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(scalers_grp_id < 0)
    {
        logE<<"creating group MAPS/Scalers"<<"\n";
        return false;
    }
*/
	// try to open mca dataset and expand before creating 
	dset_id = H5Dopen(spec_grp_id, path.c_str(), H5P_DEFAULT);
	if (dset_id < 0)
	{
		dataspace_id = H5Screate_simple(3, dims_out, maxdims);
		dset_id = H5Dcreate(spec_grp_id, path.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
	}
	else
	{
		dataspace_id = H5Dget_space(dset_id);
		int status_n = H5Sget_simple_extent_dims(dataspace_id, &tmp_dims[0], NULL);
		if (status_n > -1)
		{
			bool expand = false;
			for (int i = 0; i < 3; i++)
			{
				if (tmp_dims[i] < dims_out[i])
				{
					expand = true;
					break;
				}
			}
			if (expand)
			{
				herr_t err = H5Dset_extent(dset_id, dims_out);
			}
		}
	}

    dset_rt_id = H5Dcreate (spec_grp_id, "Elapsed_Realtime", H5T_INTEL_R, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dset_lt_id = H5Dcreate (spec_grp_id, "Elapsed_Livetime", H5T_INTEL_R, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    incnt_dset_id = H5Dcreate (spec_grp_id, "Input_Counts", H5T_INTEL_R, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    outcnt_dset_id = H5Dcreate (spec_grp_id, "Output_Counts", H5T_INTEL_R, dataspace_time_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Sselect_hyperslab (memoryspace_time_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

    real_t real_time;
    real_t life_time;
    real_t in_cnt;
    real_t out_cnt;
    for(size_t row=row_idx_start; row < (size_t)row_idx_end; row++)
    {
        offset[1] = row;
        offset_time[0] = row;
        for(size_t col=col_idx_start; col < (size_t)col_idx_end; col++)
        {
            const data_struct::Spectra *spectra = &((*spectra_volume)[row][col]);
            offset[2] = col;
            offset_time[1] = col;
            H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

            H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);

            H5Sselect_hyperslab (file_time_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

            real_time = spectra->elapsed_realtime();
            life_time = spectra->elapsed_livetime();
            in_cnt = spectra->input_counts();
            out_cnt = spectra->output_counts();
            H5Dwrite (dset_rt_id, H5T_NATIVE_REAL, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&real_time);
            H5Dwrite (dset_lt_id, H5T_NATIVE_REAL, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&life_time);
            H5Dwrite (incnt_dset_id, H5T_NATIVE_REAL, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&in_cnt);
            H5Dwrite (outcnt_dset_id, H5T_NATIVE_REAL, memoryspace_time_id, file_time_id, H5P_DEFAULT, (void*)&out_cnt);
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
    H5Sclose(dataspace_time_id);
    H5Sclose(dataspace_id);
    H5Pclose(dcpl_id);


    int_spec_grp_id = H5Gopen(spec_grp_id, STR_INT_SPEC.c_str(), H5P_DEFAULT);
    if(int_spec_grp_id < 0)
        int_spec_grp_id = H5Gcreate(spec_grp_id, STR_INT_SPEC.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //save quantification_standard integrated spectra
    data_struct::Spectra spectra = spectra_volume->integrate();
    count[0] = spectra.size();
    memoryspace_id = H5Screate_simple(1, count, nullptr);
    dataspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (int_spec_grp_id, "Spectra", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    offset[0] = 0;
    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&spectra[0]);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);


    //save real_time
    count[0] = 1;
    real_t save_val = spectra.elapsed_realtime();
    dataspace_id = H5Screate_simple (1, count, nullptr);
    memoryspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (int_spec_grp_id, "Elapsed_Realtime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save life_time
    save_val = spectra.elapsed_livetime();
    dataspace_id = H5Screate_simple (1, count, nullptr);
    memoryspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (int_spec_grp_id, "Elapsed_Livetime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save input_counts
    save_val = spectra.input_counts();
    dataspace_id = H5Screate_simple (1, count, nullptr);
    memoryspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (int_spec_grp_id, "Input_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save output_counts
    save_val = spectra.output_counts();
    dataspace_id = H5Screate_simple (1, count, nullptr);
    memoryspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (int_spec_grp_id, "Output_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save energy vector
    std::vector<real_t> out_vec;
    data_struct::gen_energy_vector(spectra.size(), energy_offset, energy_slope, &out_vec);
    count[0] = out_vec.size();
    memoryspace_id = H5Screate_simple(1, count, nullptr);
    dataspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (spec_grp_id, "Energy", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    offset[0] = 0;
    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&out_vec[0]);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save energy calibration
    count[0] = 3;
    dataspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (spec_grp_id, "Energy_Calibration", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    count[0] = 1;
    memoryspace_id = H5Screate_simple(1, count, nullptr);
    offset[0] = 0;
    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&energy_offset);
    offset[0] = 1;
    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&energy_slope);
    offset[0] = 2;
    H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&energy_quad);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(dataspace_id);

    //save file version
    save_val = HDF5_SAVE_VERSION;
    dataspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (maps_grp_id, "version", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    H5Dclose(dset_id);
    H5Sclose(dataspace_id);

    H5Gclose(int_spec_grp_id);
    H5Gclose(spec_grp_id);
    //H5Gclose(scalers_grp_id);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_element_fits(std::string path,
                                const data_struct::Fit_Count_Dict * const element_counts,
                                size_t row_idx_start,
                                int row_idx_end,
                                size_t col_idx_start,
                                int col_idx_end)
{
    std::lock_guard<std::mutex> lock(_mutex);

    

    if(_cur_file_id < 0)
    {
        logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
        return false;
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::string xrf_grp_name = "XRF_Analyzed";
    hid_t   dset_id, dset_ch_id, dset_un_id;
    hid_t   memoryspace, filespace, dataspace_id, dataspace_ch_id, dataspace_ch_off_id;
    hid_t   filetype, memtype, status;
    hid_t   dcpl_id;
    hid_t   xrf_grp_id, fit_grp_id, maps_grp_id;
    //herr_t   error;

    dset_id = -1;
    dset_ch_id = -1;
    hsize_t dims_out[3];
    hsize_t offset[1] = {0};
    hsize_t offset2[1] = {0};
    hsize_t offset_3d[3] = {0, 0, 0};
    hsize_t count[1] = {1};
    hsize_t count_3d[3] = {1, 1, 1};
    hsize_t chunk_dims[3];
	hsize_t tmp_dims[3];

    //fix this
    for(const auto& iter : *element_counts)
    {
		dims_out[1] = iter.second.rows();
		dims_out[2] = iter.second.cols();
        break;
    }
    //H5T_FLOAT

    dims_out[0] = element_counts->size();
    offset_3d[0] = 0;
    offset_3d[1] = 0;
    offset_3d[2] = 0;
    count_3d[0] = 1;
    count_3d[1] = dims_out[1];
    count_3d[2] = dims_out[2];
    chunk_dims[0] = 1;
    chunk_dims[1] = dims_out[1]; 
    chunk_dims[2] = dims_out[2];

	hsize_t      maxdims[3] = {H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED };

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    memoryspace = H5Screate_simple(3, count_3d, nullptr);
    filespace = H5Screate_simple(3, dims_out, nullptr);

    dataspace_ch_id = H5Screate_simple (1, dims_out, nullptr);
    dataspace_ch_off_id = H5Screate_simple (1, dims_out, nullptr);

    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        logE<<"creating group 'MAPS'"<<"\n";
        return false;
    }

    xrf_grp_id = H5Gopen(maps_grp_id, xrf_grp_name.c_str(), H5P_DEFAULT);
    if(xrf_grp_id < 0)
        xrf_grp_id = H5Gcreate(maps_grp_id, xrf_grp_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(xrf_grp_id < 0)
    {
        logE<<"creating group MAPS/"<<xrf_grp_name<<"\n";
        return false;
    }

    fit_grp_id = H5Gopen(xrf_grp_id, path.c_str(), H5P_DEFAULT);
    if(fit_grp_id < 0)
        fit_grp_id = H5Gcreate(xrf_grp_id, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(fit_grp_id < 0)
    {
        logE<<"creating group MAPS/"<<xrf_grp_name<<"/"<<path<<"\n";
        return false;
    }

    dset_id = H5Dopen (fit_grp_id, "Counts_Per_Sec", H5P_DEFAULT);
	if (dset_id < 0)
	{
		dataspace_id = H5Screate_simple(3, dims_out, maxdims);
		dset_id = H5Dcreate(fit_grp_id, "Counts_Per_Sec", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
	}
	else
	{
		//if we are opening, check the dim size to see if we have to expand to fit new counts size
		dataspace_id = H5Dget_space(dset_id);
		int status_n = H5Sget_simple_extent_dims(dataspace_id, &tmp_dims[0], NULL);
		if (status_n > -1)
		{
			bool expand = false;
			for (int i = 0; i < 3; i++)
			{
				if (tmp_dims[i] < dims_out[i])
				{
					expand = true;
					break;
				}
			}
			if (expand)
			{
				herr_t err = H5Dset_extent(dset_id, dims_out);
			}
		}
	}
    if(dset_id < 0)
    {
        logE<<"creating dataset MAPS/"<<xrf_grp_name<<"/"<<path<<"/Counts_Per_Sec"<<"\n";
        return false;
    }

    //filetype = H5Tcopy (H5T_FORTRAN_S1);
    filetype = H5Tcopy (H5T_C_S1);
    H5Tset_size (filetype, 256);
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, 255);

    dset_ch_id = H5Dopen (fit_grp_id, "Channel_Names", H5P_DEFAULT);
    if(dset_ch_id < 0)
        dset_ch_id = H5Dcreate (fit_grp_id, "Channel_Names", filetype, dataspace_ch_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dset_ch_id < 0)
    {
        logE<<"creating dataset MAPS/"<<xrf_grp_name<<"/"<<path<<"/Channel_Names"<<"\n";
        return false;
    }

	dset_un_id = H5Dopen(fit_grp_id, "Channel_Units", H5P_DEFAULT);
	if (dset_un_id < 0)
		dset_un_id = H5Dcreate(fit_grp_id, "Channel_Units", filetype, dataspace_ch_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset_un_id < 0)
	{
		logE << "creating dataset MAPS/" << xrf_grp_name << "/" << path << "/Channel_Units" << "\n";
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
    for (std::string el_name : data_struct::Element_Symbols)
    {
        element_lines.push_back(el_name);
    }
    for (std::string el_name : data_struct::Element_Symbols)
    {
        element_lines.push_back(el_name+"_L");
    }
    for (std::string el_name : data_struct::Element_Symbols)
    {
        element_lines.push_back(el_name+"_M");
    }

    //add pileups
    for (const auto& itr : *element_counts)
    {
        if( std::find(element_lines.begin(), element_lines.end(), itr.first) == element_lines.end() )
        {
            element_lines.push_back(itr.first);
        }

        /*
        int idx = itr.first.find("_");
        int len = itr.first.length();
        if(idx > -1)
        {
            std::string e1 = itr.first.substr(0, idx);
            std::string e2 = itr.first.substr(idx+1, len);
            auto it1 = std::find(element_lines.begin(), element_lines.end(), e1);
            auto it2 = std::find(element_lines.begin(), element_lines.end(), e2);
            if(it1 != element_lines.end() && it2 != element_lines.end())
            {
                element_lines.push_back(itr.first);
            }
        }
        */
    }
/*
    element_lines.push_back(STR_COHERENT_SCT_AMPLITUDE);
    element_lines.push_back(STR_COMPTON_AMPLITUDE);
	element_lines.push_back(STR_SUM_ELASTIC_INELASTIC_AMP);
    element_lines.push_back(STR_TOTAL_FLUORESCENCE_YIELD);
    element_lines.push_back(STR_NUM_ITR);
*/
    int i=0;
    //save by element Z order
    //for(const auto& iter : *element_counts)
	std::string units = "cts/s";
    for (std::string el_name : element_lines)
    {
        char tmp_char[256] = {0};
        if(element_counts->count(el_name) < 1 )
        {
            continue;
        }
        offset[0] = i;
        offset_3d[0] = i;

        H5Sselect_hyperslab (dataspace_ch_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        H5Sselect_hyperslab (dataspace_ch_off_id, H5S_SELECT_SET, offset2, nullptr, count, nullptr);

        el_name.copy(tmp_char, 254);

        status = H5Dwrite (dset_ch_id, memtype, dataspace_ch_off_id, dataspace_ch_id, H5P_DEFAULT, (void*)tmp_char);

		for (int z = 0; z < 256; z++)
		{
			tmp_char[z] = '\0';
		}
		if (el_name != STR_NUM_ITR && el_name != STR_RESIDUAL)
		{
			units.copy(tmp_char, 256);
		}
		status = H5Dwrite(dset_un_id, memtype, dataspace_ch_off_id, dataspace_ch_id, H5P_DEFAULT, (void*)tmp_char);

        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset_3d, nullptr, chunk_dims, nullptr);

        status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace, filespace, H5P_DEFAULT, (void*)element_counts->at(el_name).data());
        
        i++;
    }

    H5Dclose(dset_id);
    H5Dclose(dset_ch_id);
	H5Dclose(dset_un_id);
    H5Sclose(memoryspace);
    H5Sclose(filespace);
    H5Sclose(dataspace_ch_off_id);
    H5Sclose(dataspace_ch_id);
    H5Tclose(filetype);
    H5Tclose(memtype);
    H5Pclose(dcpl_id);    
    H5Sclose(dataspace_id);
    H5Gclose(fit_grp_id);
    H5Gclose(xrf_grp_id);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";


    return true;

}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_fitted_int_spectra(const std::string path,
                                     const data_struct::Spectra& spectra,
                                     const data_struct::Range& spectra_range,
                                     const data_struct::Spectra& background,
									 const size_t save_spectra_size)
{
    std::lock_guard<std::mutex> lock(_mutex);

    if(_cur_file_id < 0)
    {
        logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
        return false;
    }

    bool ret_val = true;
    hid_t   dset_id, dataspace_id, status;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    std::string dset_name = "/MAPS/XRF_Analyzed/" + path + "/"+ STR_FIT_INT_SPEC;
    std::string background_name = "/MAPS/XRF_Analyzed/" + path + "/"+ STR_FIT_INT_BACKGROUND;

    hsize_t offset[1] = {0};
    hsize_t count[1] = {1};
    count[0] = save_spectra_size;
    dataspace_id = H5Screate_simple(1, count, nullptr);

    // resize to the size of collected spectra
    data_struct::ArrayXr save_spectra;
    save_spectra.setZero(save_spectra_size);
    data_struct::ArrayXr save_background;
    save_background.setZero(save_spectra_size);

    int j = 0;
    for(int i= spectra_range.min; i <= spectra_range.max; i++)
    {
        if(std::isfinite(spectra[j]))
        {
            save_spectra[i] = spectra[j];
        }
        if (std::isfinite(background[j]))
        {
            save_background[i] = background[j];
        }
        j++;
    }

    // save spectra
    dset_id = H5Dopen (_cur_file_id, dset_name.c_str(), H5P_DEFAULT);
    if(dset_id < 0)
    {
        dset_id = H5Dcreate2 (_cur_file_id, dset_name.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    if(dset_id < 0)
    {
        logE<<"creating dataset "<<dset_name<<"\n";
        ret_val = false;
    }
    else
    {
        status = H5Dwrite(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)save_spectra.data());
        if (status < 0)
        {
            logW<<"Failed to save "<< dset_name <<"\n";
            ret_val = false;
        }

        H5Dclose(dset_id);
    }

    // save background
    dset_id = H5Dopen(_cur_file_id, background_name.c_str(), H5P_DEFAULT);
    if (dset_id < 0)
    {
        dset_id = H5Dcreate2(_cur_file_id, background_name.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    if (dset_id < 0)
    {
        logE << "creating dataset " << background_name << "\n";
        ret_val = false;
    }
    else
    {
        status = H5Dwrite(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)save_background.data());
        if (status < 0)
        {
            logW << "Failed to save " << background_name << "\n";
            ret_val = false;
        }

        H5Dclose(dset_id);
    }


    if(dataspace_id > -1)
    {
        H5Sclose(dataspace_id);
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";
    return ret_val;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_max_10_spectra(const std::string path,
	const data_struct::Range& spectra_range,
	const data_struct::Spectra& max_spectra,
	const data_struct::Spectra& max_10_spectra,
    const data_struct::Spectra& fit_int_background)
{
	std::lock_guard<std::mutex> lock(_mutex);

	if (_cur_file_id < 0)
	{
		logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
		return false;
	}

	bool ret_val = true;
	hid_t   dset_id, maps_grp_id, spec_grp_id, int_spec_grp_id, dataspace_id, status;
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	

	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };
	count[0] = max_spectra.size();
	dataspace_id = H5Screate_simple(1, count, nullptr);

	maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
	if (maps_grp_id < 0)
	{
		logE << "opening group MAPS" << "\n";
		ret_val = false;
	}
	spec_grp_id = H5Gopen(maps_grp_id, "Spectra", H5P_DEFAULT);
	if (spec_grp_id < 0)
	{
		spec_grp_id = H5Gcreate(maps_grp_id, "Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	if (spec_grp_id < 0)
	{
		logE << "creating group MAPS/Spectra" << "\n";
		ret_val = false;
	}

	int_spec_grp_id = H5Gopen(spec_grp_id, STR_INT_SPEC.c_str(), H5P_DEFAULT);
	if (int_spec_grp_id < 0)
	{
		int_spec_grp_id = H5Gcreate(spec_grp_id, STR_INT_SPEC.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	if (int_spec_grp_id < 0)
	{
		logE << "creating group MAPS/Spectra/"<< STR_INT_SPEC << "\n";
		ret_val = false;
	}

	dset_id = H5Dopen(int_spec_grp_id, STR_MAX_CHANNELS_INT_SPEC.c_str(), H5P_DEFAULT);
	if (dset_id < 0)
	{
		dset_id = H5Dcreate2(int_spec_grp_id, STR_MAX_CHANNELS_INT_SPEC.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	if (dset_id < 0)
	{
		logE << "creating dataset " << STR_MAX_CHANNELS_INT_SPEC << "\n";
		ret_val = false;
	}
	else
	{
		status = H5Dwrite(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)max_spectra.data());
		if (status < 0)
		{
			logW << "Failed to save " << STR_MAX_CHANNELS_INT_SPEC << "\n";
			ret_val = false;
		}
		H5Dclose(dset_id);
	}

	dset_id = H5Dopen(int_spec_grp_id, STR_MAX10_INT_SPEC.c_str(), H5P_DEFAULT);
	if (dset_id < 0)
	{
		dset_id = H5Dcreate2(int_spec_grp_id, STR_MAX10_INT_SPEC.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	}
	if (dset_id < 0)
	{
		logE << "creating dataset " << STR_MAX10_INT_SPEC << "\n";
		ret_val = false;
	}
	else
	{
		status = H5Dwrite(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)max_10_spectra.data());
		if (status < 0)
		{
			logW << "Failed to save " << STR_MAX10_INT_SPEC << "\n";
			ret_val = false;
		}
		H5Dclose(dset_id);
	}

	if (maps_grp_id > -1)
	{
		H5Dclose(maps_grp_id);
	}
	if (spec_grp_id > -1)
	{
		H5Dclose(spec_grp_id);
	}
	if (int_spec_grp_id > -1)
	{
		H5Dclose(int_spec_grp_id);
	}
	if (dataspace_id > -1)
	{
		H5Sclose(dataspace_id);
	}

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";
	return ret_val;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_quantification(data_struct::Detector* detector)
{

    std::lock_guard<std::mutex> lock(_mutex);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    dset_id, memoryspace_id, filespace_id, dataspace_id, filetype, dataspace_ch_id, memtype,  dset_ch_id, dset_un_id, q_int_spec_grp_id;
    hid_t   memtype_label, filetype_label, q_memoryspace_label_id, q_dataspace_label_id;
    hid_t   count_dataspace_id, count_dset_id, standard_grp_id, calib_grp_id;

    hid_t q_dataspace_id, q_memoryspace_id, q_filespace_id, q_dset_id, q_grp_id, q_fit_grp_id, maps_grp_id, status, scalers_grp_id, xrf_fits_grp_id;
    hid_t dset_labels_id;
    hsize_t offset[3];
    hsize_t count[3];

    char unit_char[255] = "cts/s";

    hsize_t q_dims_out[2];

    if(_cur_file_id < 0)
    {
        logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
        return false;
    }


    offset[0] = 0;
    offset[1] = 0;
    offset[2] = 0;
    count[0] = 1;
    count[1] = 0;
    count[2] = 0;

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if(maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(maps_grp_id < 0)
    {
        logE<<"creating group 'MAPS'"<<"\n";
        return false;
    }

    filetype = H5Tcopy (H5T_FORTRAN_S1);
    H5Tset_size (filetype, 256);
    memtype = H5Tcopy (H5T_C_S1);
    status = H5Tset_size (memtype, 255);


    filetype_label = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype_label, 10);
    memtype_label = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(memtype_label, 10);


    //--                        save calibration curve                  --

    if (detector != nullptr)
    {

        q_grp_id = H5Gopen(maps_grp_id, "Quantification", H5P_DEFAULT);
        if (q_grp_id < 0)
            q_grp_id = H5Gcreate(maps_grp_id, "Quantification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (q_grp_id < 0)
        {
            logE << "creating group MAPS/Quantification" << "\n";
            return false;
        }

        calib_grp_id = H5Gopen(q_grp_id, "Calibration", H5P_DEFAULT);
        if (calib_grp_id < 0)
            calib_grp_id = H5Gcreate(q_grp_id, "Calibration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (calib_grp_id < 0)
        {
            logE << "creating group MAPS/Quantification/Calibration" << "\n";
            return false;
        }

        //create dataset telling how many standards there are
        count[0] = 1;
        count_dataspace_id = H5Screate_simple(1, count, nullptr);
        count_dset_id = H5Dcreate(q_grp_id, "Number_Of_Standards", H5T_NATIVE_INT, count_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        int quant_size = detector->quantification_standards.size();
        status = H5Dwrite(count_dset_id, H5T_NATIVE_INT, count_dataspace_id, count_dataspace_id, H5P_DEFAULT, (void*)&quant_size);
        H5Sclose(count_dataspace_id);
        H5Dclose(count_dset_id);

        int standard_idx = 0;

        // ----------------------------------------- start per standard  ------------------------------------------------
        for (const auto& quant_itr : detector->quantification_standards)
        {
            //create group
            string standard_group_name = "Standard" + std::to_string(standard_idx);
            standard_grp_id = H5Gopen(q_grp_id, standard_group_name.c_str(), H5P_DEFAULT);
            if (standard_grp_id < 0)
                standard_grp_id = H5Gcreate(q_grp_id, standard_group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (standard_grp_id < 0)
            {
                logE << "creating group MAPS/Quantification/"<< standard_group_name << "\n";
                return false;
            }

            scalers_grp_id = H5Gopen(standard_grp_id, "Scalers", H5P_DEFAULT);
            if (scalers_grp_id < 0)
                scalers_grp_id = H5Gcreate(standard_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (scalers_grp_id < 0)
            {
                logE << "creating group MAPS/Quantification/"<< standard_group_name << "/"<< "Scalers" << "\n";
                return false;
            }

            xrf_fits_grp_id = H5Gopen(standard_grp_id, "XRF_Analyzed", H5P_DEFAULT);
            if (xrf_fits_grp_id < 0)
                xrf_fits_grp_id = H5Gcreate(standard_grp_id, "XRF_Analyzed", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (xrf_fits_grp_id < 0)
            {
                logE << "creating group MAPS/Quantification/"<< standard_group_name <<"/XRF_Analyzed" << "\n";
                return false;
            }

            //save quantification_standard element weights
            count[0] = quant_itr.second.element_standard_weights.size();
            memoryspace_id = H5Screate_simple(1, count, nullptr);
            filespace_id = H5Screate_simple(1, count, nullptr);
            dataspace_id = H5Screate_simple(1, count, nullptr);
            dataspace_ch_id = H5Screate_simple(1, count, nullptr);
            dset_id = H5Dcreate(standard_grp_id, "Element_Weights", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            dset_ch_id = H5Dcreate(standard_grp_id, "Element_Weights_Names", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            offset[0] = 0;
            count[0] = 1;
            H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            int offset_idx = 0;
            for (auto itr : quant_itr.second.element_standard_weights)
            {
                offset[0] = offset_idx;
                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

                status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(itr.second));
                char tmp_char[256] = { 0 };
                itr.first.copy(tmp_char, 254);
                status = H5Dwrite(dset_ch_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)tmp_char);

                offset_idx++;
            }
            H5Dclose(dset_id);
            H5Dclose(dset_ch_id);
            H5Sclose(memoryspace_id);
            H5Sclose(filespace_id);
            H5Sclose(dataspace_id);
            H5Sclose(dataspace_ch_id);


            q_int_spec_grp_id = H5Gopen(standard_grp_id, STR_INT_SPEC.c_str(), H5P_DEFAULT);
            if (q_int_spec_grp_id < 0)
                q_int_spec_grp_id = H5Gcreate(standard_grp_id, STR_INT_SPEC.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (q_int_spec_grp_id < 0)
            {
                logE << "creating group MAPS/Quantification/"<< standard_group_name<<"/"<< STR_INT_SPEC << "\n";
                return false;
            }

            //save quantification_standard integrated spectra
            data_struct::Spectra spectra = quant_itr.second.integrated_spectra;
            if (spectra.size() > 0)
            {
                count[0] = spectra.size();
                memoryspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dataspace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(q_int_spec_grp_id, "Spectra", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                offset[0] = 0;
                H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&spectra[0]);
                H5Dclose(dset_id);
                H5Sclose(memoryspace_id);
                H5Sclose(filespace_id);
                H5Sclose(dataspace_id);
            }

            //save standard name
            count[0] = 1;
            dataspace_id = H5Screate_simple(1, count, nullptr);
            memoryspace_id = H5Screate_simple(1, count, nullptr);
            dset_ch_id = H5Dcreate(standard_grp_id, "Standard_Name", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            char tmp_char[255] = { 0 };
            quant_itr.second.standard_filename.copy(tmp_char, 254);
            status = H5Dwrite(dset_ch_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
            H5Dclose(dset_ch_id);
            H5Sclose(dataspace_id);

            //save sr_current
            dataspace_id = H5Screate_simple(1, count, nullptr);
            dset_id = H5Dcreate(scalers_grp_id, STR_SR_CURRENT.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quant_itr.second.sr_current));
            H5Dclose(dset_id);
            H5Sclose(dataspace_id);

            //save us_ic
            dataspace_id = H5Screate_simple(1, count, nullptr);
            dset_id = H5Dcreate(scalers_grp_id, STR_US_IC.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quant_itr.second.US_IC));
            H5Dclose(dset_id);
            H5Sclose(dataspace_id);

            //save ds_ic
            dataspace_id = H5Screate_simple(1, count, nullptr);
            dset_id = H5Dcreate(scalers_grp_id, STR_DS_IC.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quant_itr.second.DS_IC));
            H5Dclose(dset_id);
            H5Sclose(dataspace_id);


            //save real_time
            real_t save_val = spectra.elapsed_realtime();
            dataspace_id = H5Screate_simple(1, count, nullptr);
            dset_id = H5Dcreate(q_int_spec_grp_id, "Elapsed_Realtime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
            H5Dclose(dset_id);
            H5Sclose(dataspace_id);

            //save life_time
            save_val = spectra.elapsed_livetime();
            dataspace_id = H5Screate_simple(1, count, nullptr);
            dset_id = H5Dcreate(q_int_spec_grp_id, "Elapsed_Livetime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
            H5Dclose(dset_id);
            H5Sclose(memoryspace_id);
            H5Sclose(dataspace_id);

            //save input counts
            save_val = spectra.input_counts();
            dataspace_id = H5Screate_simple(1, count, nullptr);
            dset_id = H5Dcreate(q_int_spec_grp_id, "Input_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
            H5Dclose(dset_id);
            H5Sclose(memoryspace_id);
            H5Sclose(dataspace_id);

            //save output counts
            save_val = spectra.output_counts();
            dataspace_id = H5Screate_simple(1, count, nullptr);
            dset_id = H5Dcreate(q_int_spec_grp_id, "Output_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
            H5Dclose(dset_id);
            H5Sclose(memoryspace_id);
            H5Sclose(dataspace_id);


            for (const auto& fit_itr : detector->fitting_quant_map)
            {

                q_fit_grp_id = H5Gopen(xrf_fits_grp_id, Fitting_Routine_To_Str.at(fit_itr.first).c_str(), H5P_DEFAULT);
                if (q_fit_grp_id < 0)
                    q_fit_grp_id = H5Gcreate(xrf_fits_grp_id, Fitting_Routine_To_Str.at(fit_itr.first).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (q_fit_grp_id < 0)
                {
                    logE << "creating group MAPS/Quantification/" << standard_group_name<<"/XRF_Analyzed/"<< Fitting_Routine_To_Str.at(fit_itr.first).c_str() << "\n";
                    return false;
                }

                //save quantification_standard counts
                std::unordered_map<std::string, real_t> element_counts = quant_itr.second.element_counts.at(fit_itr.first);
                count[0] = element_counts.size();
                dataspace_id = H5Screate_simple(1, count, nullptr);
                memoryspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(q_fit_grp_id, "Counts_Per_Sec", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_ch_id = H5Dcreate(q_fit_grp_id, "Channel_Names", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                dset_un_id = H5Dcreate(q_fit_grp_id, "Channel_Units", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                //create save ordered vector by element Z number with K , L, M lines
                std::vector<std::string> element_lines;
                for (std::string el_name : data_struct::Element_Symbols)
                {
                    element_lines.push_back(el_name);
                }
                for (std::string el_name : data_struct::Element_Symbols)
                {
                    element_lines.push_back(el_name + "_L");
                }
                for (std::string el_name : data_struct::Element_Symbols)
                {
                    element_lines.push_back(el_name + "_M");
                }

                element_lines.push_back(STR_NUM_ITR);
                element_lines.push_back(STR_RESIDUAL);

                offset_idx = 0;
                offset[0] = 0;
                count[0] = 1;
                H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                //save by element Z order
                for (std::string el_name : element_lines)
                {
                    if (element_counts.count(el_name) < 1)
                    {
                        continue;
                    }
                    offset[0] = offset_idx;
                    real_t val = element_counts.at(el_name);

                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                    std::memset(&tmp_char[0], 0, 255);
                    el_name.copy(tmp_char, 254);
                    status = H5Dwrite(dset_ch_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)tmp_char);
                    status = H5Dwrite(dset_un_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)unit_char);
                    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&val);

                    offset_idx++;
                }
                H5Gclose(q_fit_grp_id);
                H5Dclose(dset_id);
                H5Dclose(dset_ch_id);
                H5Dclose(dset_un_id);
                H5Sclose(memoryspace_id);
                H5Sclose(filespace_id);
                H5Sclose(dataspace_id);
            }
            H5Gclose(scalers_grp_id);
            H5Gclose(q_int_spec_grp_id);
            H5Gclose(xrf_fits_grp_id);
            H5Gclose(standard_grp_id);
            standard_idx++;
            
        }
        //save calibration curves
        //data_struct::Quantifiers quantifiers = quantification_standard->quantifier_map.at(path);
        //auto shell_itr = quantification_standard->_calibration_curves.begin();


        int element_cnt = CALIBRATION_CURVE_SIZE; // from element H to U
        q_dims_out[0] = 3;// shells K, L, and M
        q_dims_out[1] = element_cnt;

        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;

        count[0] = 1;
        count[1] = 1;
        count[2] = 0;

        q_memoryspace_label_id = H5Screate_simple(2, count, nullptr);

        count[0] = 1;
        count[1] = q_dims_out[1];
        count[2] = 0;

        q_memoryspace_id = H5Screate_simple(2, count, nullptr);
        q_filespace_id = H5Screate_simple(2, q_dims_out, nullptr);
        q_dataspace_label_id = H5Screate_simple(2, q_dims_out, nullptr);
        q_dataspace_id = H5Screate_simple(2, q_dims_out, nullptr);
        //q_dataspace_ch_id = H5Screate_simple (1, q_dims_out, nullptr);

		for (const auto& qitr : detector->fitting_quant_map)
		{
            q_fit_grp_id = H5Gopen(calib_grp_id, Fitting_Routine_To_Str.at(qitr.first).c_str(), H5P_DEFAULT);
            if (q_fit_grp_id < 0)
                q_fit_grp_id = H5Gcreate(calib_grp_id, Fitting_Routine_To_Str.at(qitr.first).c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (q_fit_grp_id < 0)
            {
                logE << "creating group MAPS/Quantification/Calibration/" << Fitting_Routine_To_Str.at(qitr.first).c_str() << "\n";
                continue;
            }

			dset_labels_id = H5Dopen(q_fit_grp_id, "Calibration_Curve_Labels", H5P_DEFAULT);
			if (dset_labels_id < 0)
				dset_labels_id = H5Dcreate(q_fit_grp_id, "Calibration_Curve_Labels", filetype_label, q_dataspace_label_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


			for (const auto & quant_scaler_itr : qitr.second.quant_scaler_map)
			{

				std::string q_dset_name = "Calibration_Curve_" + quant_scaler_itr.first;

				q_dset_id = H5Dopen(q_fit_grp_id, q_dset_name.c_str(), H5P_DEFAULT);
				if (q_dset_id < 0)
					q_dset_id = H5Dcreate(q_fit_grp_id, q_dset_name.c_str(), H5T_INTEL_R, q_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (q_dset_id < 0)
				{
					logE << "creating dataset MAPS/Quantification/Calibration/" << Fitting_Routine_To_Str.at(qitr.first).c_str() << "/" << q_dset_name  << "\n";
					continue;
				}

                int j = 0;
				//for(auto& shell_itr : quant_itr.second)
				for(const auto & calib_itr : quant_scaler_itr.second.curve_quant_map)
				{
					char label[10] = { 0 };
					//int element_offset = 0;
					//create dataset for different shell curves
                    std::vector<real_t> calibration_curve;
                    std::vector<std::string> calibration_curve_labels;

                    for (const auto& citr : calib_itr.second)
                    {
                        string name = citr.name;
                        if (calib_itr.first == Electron_Shell::L_SHELL)
                        {
                            name += "_L";
                        }
                        if (calib_itr.first == Electron_Shell::M_SHELL)
                        {
                            name += "_M";
                        }
                        calibration_curve.push_back(citr.calib_curve_val);
                        calibration_curve_labels.push_back(name);
                    }

					offset[0] = j;
					offset[1] = 0;

					count[0] = 1;
					count[1] = q_dims_out[1];
					//H5Sselect_hyperslab (q_dataspace_ch_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
					//H5Sselect_hyperslab (q_dataspace_ch_off_id, H5S_SELECT_SET, &offset[2], nullptr, count, nullptr);
					//status = H5Dwrite (q_dset_ch_id, memtype, q_dataspace_ch_off_id, q_dataspace_ch_id, H5P_DEFAULT, (void*)(el_name.c_str()));

					H5Sselect_hyperslab(q_filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
					status = H5Dwrite(q_dset_id, H5T_NATIVE_REAL, q_memoryspace_id, q_filespace_id, H5P_DEFAULT, (void*)&calibration_curve[0]);

					for (size_t k = 0; k < calibration_curve_labels.size(); k++)
					{
						memset(label, 0, 10);
						calibration_curve_labels[k].copy(&label[0], 9);
						offset[1] = k;
						count[1] = 1;
						status = H5Sselect_hyperslab(q_filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
						status = H5Dwrite(dset_labels_id, memtype_label, q_memoryspace_label_id, q_filespace_id, H5P_DEFAULT, (void*)&label[0]);
					}
                    j++;

				}
				H5Dclose(q_dset_id);
			}

			H5Dclose(dset_labels_id);
			H5Gclose(q_fit_grp_id);
		}
		
        //H5Dclose(q_dset_ch_id);
        H5Sclose(q_filespace_id);
        H5Sclose(q_memoryspace_id);
        H5Sclose(q_dataspace_label_id);
        H5Sclose(q_memoryspace_label_id);
        H5Sclose(q_dataspace_id);
        H5Gclose(calib_grp_id);
        H5Gclose(q_grp_id);
    }

    H5Tclose(filetype);
    H5Tclose(memtype);
    H5Tclose(filetype_label);
    H5Tclose(memtype_label);
    H5Gclose(maps_grp_id);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_save_params_override(hid_t group_id, data_struct::Params_Override * params_override)
{
	//TODO : save param override to hdf5 

	hid_t fit_params_grp_id;

	fit_params_grp_id = H5Gopen(group_id, "Fit_Parameters", H5P_DEFAULT);
	if (fit_params_grp_id < 0)
		fit_params_grp_id = H5Gcreate(group_id, "Fit_Parameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (fit_params_grp_id < 0)
	{
		logE << "creating group MAPS/Fit_Parameters_Override/Fit_Parameters" << "\n";
		return false;
	}

    H5Dclose(fit_params_grp_id);

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_save_scan_meta_data(hid_t scan_grp_id, data_struct::Scan_Meta_Info* meta_info)
{
    
    hid_t dataspace_id = -1, memoryspace_id = -1, filespace_id = -1;
	hid_t status;
	hid_t filetype, memtype;
    hid_t dset_id = -1;
	
	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };

	try
	{

        filetype = H5Tcopy(H5T_FORTRAN_S1);
		H5Tset_size(filetype, 256);
		memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

        //save y axis
        count[0] = meta_info->y_axis.size();
        memoryspace_id = H5Screate_simple(1, count, nullptr);
        dataspace_id = H5Screate_simple(1, count, nullptr);
        filespace_id = H5Screate_simple(1, count, nullptr);
        dset_id = H5Dcreate(scan_grp_id, "y_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id > 0)
		{
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)meta_info->y_axis.data());
            H5Dclose(dset_id);
            H5Sclose(dataspace_id);
            H5Sclose(filespace_id);
            H5Sclose(memoryspace_id);
		}

        
        count[0] = meta_info->x_axis.size();
        memoryspace_id = H5Screate_simple(1, count, nullptr);
        dataspace_id = H5Screate_simple(1, count, nullptr);
        filespace_id = H5Screate_simple(1, count, nullptr);
        dset_id = H5Dcreate(scan_grp_id, "x_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id > 0)
		{
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)meta_info->x_axis.data());
            H5Dclose(dset_id);
            H5Sclose(dataspace_id);
            H5Sclose(filespace_id);
            H5Sclose(memoryspace_id);
		}
		
        //save requested rows
        count[0] = 1;
        memoryspace_id = H5Screate_simple(1, count, nullptr);
        dataspace_id = H5Screate_simple(1, count, nullptr);
        filespace_id = H5Screate_simple(1, count, nullptr);
        dset_id = H5Dcreate(scan_grp_id, "requested_rows", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dset_id > 0)
        {
            status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(meta_info->requested_rows));
            H5Dclose(dset_id);
        }
        
        //save requested cols
        dset_id = H5Dcreate(scan_grp_id, "requested_cols", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (dset_id > 0)
        {
            status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(meta_info->requested_cols));
            H5Dclose(dset_id);
            H5Sclose(dataspace_id);
            H5Sclose(filespace_id);
        }
        
        //Save theta
        dset_id = H5Dcreate(scan_grp_id, "theta", H5T_NATIVE_REAL, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(dset_id > 0)
        {
            status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)&meta_info->theta);
            H5Dclose(dset_id);
        }

		dset_id = H5Dcreate(scan_grp_id, "scan_time_stamp", filetype, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id > 0)
		{
            char tmp_char[255] = { 0 };
            meta_info->scan_time_stamp.copy(tmp_char, 254);
            status = H5Dwrite(dset_id, memtype, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)tmp_char);
            H5Dclose(dset_id);
		}
		
		dset_id = H5Dcreate(scan_grp_id, "name", filetype, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id > 0)
		{
            char tmp_char[255] = { 0 };
            meta_info->name.copy(tmp_char, 254);
            status = H5Dwrite(dset_id, memtype, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)tmp_char);
            H5Dclose(dset_id);
        }

        H5Tclose(filetype);
        H5Tclose(memtype);

        if(memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
	}
	catch (...)
	{
        logE << "saving MAPS/Scan meta data" << "\n";
		return false;
	}

	return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_save_extras(hid_t scan_grp_id, std::vector<data_struct::Extra_PV>* extra_pvs)
{
    
    hid_t filespace_id = -1, filespace_id2 = -1, filespace_id3 = -1, filespace_id4 = -1;
    hid_t memoryspace_id = -1, memoryspace_id2 = -1, memoryspace_id3 = -1, memoryspace_id4 = -1;
    hid_t status = -1;
	hid_t filetype, memtype;
    hid_t dset_desc_id = -1, dset_unit_id = -1, dset_id = -1, dset_val_id = -1;
    hid_t extra_grp_id = -1;
	
	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };

    

	if (extra_pvs == nullptr)
	{
		return false;
	}

	try
	{

        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
		memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

		extra_grp_id = H5Gopen(scan_grp_id, "Extra_PVs", H5P_DEFAULT);
		if (extra_grp_id < 0)
			extra_grp_id = H5Gcreate(scan_grp_id, "Extra_PVs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (extra_grp_id < 0)
		{
            logE << "creating group MAPS/Scan/Extra_PVs" << "\n";
			return false;
		}

		//save extra pv's
		count[0] = extra_pvs->size();
        filespace_id = H5Screate_simple(1, count, nullptr);
        dset_id = H5Dcreate(extra_grp_id, "Names", filetype, filespace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
            logE << "creating dataset MAPS/Scan/Extra_PVs/Names" << "\n";
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

        filespace_id2 = H5Screate_simple(1, count, nullptr);
        dset_val_id = H5Dcreate(extra_grp_id, "Values", filetype, filespace_id2, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_val_id < 0)
		{
            logE << "creating dataset MAPS/Scan/Extra_PVs/Values" << "\n";
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

        filespace_id3 = H5Screate_simple(1, count, nullptr);
        dset_desc_id = H5Dcreate(extra_grp_id, "Description", filetype, filespace_id3, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_desc_id < 0)
		{
            logE << "creating dataset MAPS/Scan/Extra_PVs/Description" << "\n";
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

        filespace_id4 = H5Screate_simple(1, count, nullptr);
        dset_unit_id = H5Dcreate(extra_grp_id, "Unit", filetype, filespace_id4, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_unit_id < 0)
		{
            logE << "creating dataset MAPS/Scan/Extra_PVs/Unit" << "\n";
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

		count[0] = 1;
        memoryspace_id = H5Screate_simple(1, count, nullptr);
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        memoryspace_id2 = H5Screate_simple(1, count, nullptr);
        H5Sselect_hyperslab(memoryspace_id2, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        memoryspace_id3 = H5Screate_simple(1, count, nullptr);
        H5Sselect_hyperslab(memoryspace_id3, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        memoryspace_id4 = H5Screate_simple(1, count, nullptr);
        H5Sselect_hyperslab(memoryspace_id4, H5S_SELECT_SET, offset, nullptr, count, nullptr);

		for (int16_t i = 0; i < extra_pvs->size(); i++)
		{
			offset[0] = i;
			data_struct::Extra_PV pv = extra_pvs->at(i);
            
			H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id2, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id3, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(filespace_id4, H5S_SELECT_SET, offset, NULL, count, NULL);

            char tmp_char[255] = {0};
            pv.name.copy(tmp_char, 254);
            status = H5Dwrite(dset_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)tmp_char);
			
            char tmp_char2[255] = {0};
            pv.value.copy(tmp_char2, 254);
            status = H5Dwrite(dset_val_id, memtype, memoryspace_id2, filespace_id2, H5P_DEFAULT, (void*)tmp_char2);
			
            char tmp_char3[255] = {0};
            pv.description.copy(tmp_char, 254);
            status = H5Dwrite(dset_desc_id, memtype, memoryspace_id3, filespace_id3, H5P_DEFAULT, (void*)tmp_char3);

            char tmp_char4[255] = {0};
            pv.unit.copy(tmp_char, 254);
            status = H5Dwrite(dset_unit_id, memtype, memoryspace_id4, filespace_id4, H5P_DEFAULT, (void*)tmp_char4);
		}

        H5Tclose(filetype);
        H5Tclose(memtype);

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
        logE << "creating group MAPS/Scan/Extra_PVs" << "\n";
		return false;
	}
    
	return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_save_scalers(hid_t maps_grp_id, std::vector<data_struct::Scaler_Map>* scalers_map, data_struct::Params_Override* params_override, bool hasNetcdf)
{

    hid_t dataspace_id = -1, memoryspace_id = -1, filespace_id = -1, filespace_name_id = -1, memoryspace_str_id = -1;
    hid_t filetype, memtype;
    hid_t dset_names_id = -1;
    hid_t dset_units_id = -1;
    hid_t dset_values_id = -1;
    hid_t scalers_grp_id = -1;
    hid_t dcpl_id = -1, status;
    

    hsize_t offset[1] = { 0 };
    hsize_t count[1] = { 1 };

    hsize_t count_2d[2] = { 1, 1 };

    hsize_t offset_3d[3] = { 0, 0, 0 };
    hsize_t count_3d[3] = { 1, 1, 1 };

    real_t time_scaler_clock = 1.0;

    bool single_row_scan = false;

    try
    {
        memoryspace_str_id = H5Screate_simple(1, count, NULL);
        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
        memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

        count[0] = 1;

        real_t val;
        std::string units;

        scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
        if (scalers_grp_id < 0)
            scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (scalers_grp_id < 0)
        {
            logE << "creating group MAPS/Scalers" << "\n";
            return false;
        }

        _save_amps(scalers_grp_id, params_override);

        if (scalers_map != nullptr)
        {
            int cols = 0;
            int rows = 0;
            // rename from pv's to scaler names
            if (params_override != nullptr)
            {

                data_struct::ArrayXXr* time_map = nullptr;
                // try to find time scaler in scalers_map
                if (params_override->time_scaler_clock.length() > 0)
                {
                    time_scaler_clock = std::stod(params_override->time_scaler_clock);
                }

                for (auto& itr : *scalers_map)
                {
                    if (itr.name == params_override->time_scaler)
                    {
                        time_map = &(itr.values);
                    }
                }

                //update time map
                if (time_map != nullptr && time_scaler_clock > 0)
                {
                    (*time_map) /= time_scaler_clock;
                }

                // now iterate through scalers_map and update names and time normalized values
                for (auto& itr : *scalers_map)
                {
                    bool found_pv = false;

                    if (cols == 0 || rows == 0)
                    {
                        rows = itr.values.rows();
                        cols = itr.values.cols();
                    }

                    for (const auto& ts_itr : params_override->time_normalized_scalers)
                    {
                        if (ts_itr.second == itr.name)
                        {
                            itr.name = ts_itr.first;
                            if (time_map != nullptr)
                            {
                                itr.values /= (*time_map);
                            }
                            found_pv = true;
                            break;
                        }
                    }

                    if (found_pv == false)
                    {
                        for (const auto& s_itr : params_override->scaler_pvs)
                        {
                            if (s_itr.second == itr.name)
                            {
                                itr.name = s_itr.first;
                                break;
                            }
                        }
                    }
                }

                // add summed_scalers
                for (const auto& summed_scaler_itr : params_override->summed_scalers)
                {
                    data_struct::Scaler_Map s_map;
                    s_map.values.resize(rows, cols);
                    s_map.values.setZero(rows, cols);
                    s_map.name = summed_scaler_itr.scaler_name;
                    s_map.unit = " ";
                    // look for scaler names and add them up
                    for (const auto& scaler_names_itr : summed_scaler_itr.scalers_to_sum)
                    {
                        for (const auto& scaler : *scalers_map)
                        {
                            if(scaler.name == scaler_names_itr)
                            {
                                s_map.values += scaler.values;
                                break;
                            }
                        }
                    }
                    scalers_map->push_back(s_map);
                }
            }


            if( rows > 0 && cols > 0)
            {
                // create calculated scalers
                data_struct::Scaler_Map abs_ic_map, abs_cfg_map, H_dpc_cfg_map, V_dpc_cfg_map, dia1_dpc_cfg_map, dia2_dpc_cfg_map;
                abs_ic_map.name = "abs_ic";
                abs_ic_map.unit = " ";
                abs_ic_map.values.resize(rows, cols);

                abs_cfg_map.name = "abs_cfg";
                abs_cfg_map.unit = " ";
                abs_cfg_map.values.resize(rows, cols);

                H_dpc_cfg_map.name = "H_dpc_cfg";
                H_dpc_cfg_map.unit = " ";
                H_dpc_cfg_map.values.resize(rows, cols);

                V_dpc_cfg_map.name = "V_dpc_cfg";
                V_dpc_cfg_map.unit = " ";
                V_dpc_cfg_map.values.resize(rows, cols);

                dia1_dpc_cfg_map.name = "dia1_dpc_cfg";
                dia1_dpc_cfg_map.unit = " ";
                dia1_dpc_cfg_map.values.resize(rows, cols);

                dia2_dpc_cfg_map.name = "dia2_dpc_cfg";
                dia2_dpc_cfg_map.unit = " ";
                dia2_dpc_cfg_map.values.resize(rows, cols);

                // CFG_2 - 5
                data_struct::ArrayXXr* us_ic_map = nullptr;
                data_struct::ArrayXXr* ds_ic_map = nullptr;
                data_struct::ArrayXXr* cfg_2_map = nullptr;
                data_struct::ArrayXXr* cfg_3_map = nullptr;
                data_struct::ArrayXXr* cfg_4_map = nullptr;
                data_struct::ArrayXXr* cfg_5_map = nullptr;

                // search for scalers
                for (auto& scaler : *scalers_map)
                {
                    std::string upper_scaler_name = scaler.name;
                    std::transform(upper_scaler_name.begin(), upper_scaler_name.end(), upper_scaler_name.begin(), ::toupper);
                    if (upper_scaler_name == STR_US_IC)
                    {
                        us_ic_map = &(scaler.values);
                    }
                    if (upper_scaler_name == STR_DS_IC)
                    {
                        ds_ic_map = &(scaler.values);
                    }
                    if (upper_scaler_name == "CFG_2")
                    {
                        cfg_2_map = &(scaler.values);
                    }
                    if (upper_scaler_name == "CFG_3")
                    {
                        cfg_3_map = &(scaler.values);
                    }
                    if (upper_scaler_name == "CFG_4")
                    {
                        cfg_4_map = &(scaler.values);
                    }
                    if (upper_scaler_name == "CFG_5")
                    {
                        cfg_5_map = &(scaler.values);
                    }

                    if (us_ic_map != nullptr && ds_ic_map != nullptr && cfg_2_map != nullptr && cfg_3_map != nullptr && cfg_4_map != nullptr && cfg_5_map != nullptr)
                    {
                        break;
                    }
                }

                if (us_ic_map != nullptr && ds_ic_map != nullptr)
                {
                    abs_ic_map.values = (*ds_ic_map) / (*us_ic_map);
                }

                if (us_ic_map != nullptr && cfg_2_map != nullptr && cfg_3_map != nullptr && cfg_4_map != nullptr && cfg_5_map != nullptr)
                {
                    data_struct::ArrayXXr t_abs_map;
                    t_abs_map.resize(rows, cols);
                    t_abs_map = (*cfg_2_map) + (*cfg_3_map) + (*cfg_4_map) + (*cfg_5_map);
                    abs_cfg_map.values = t_abs_map / (*us_ic_map);
                    
                    if (t_abs_map.sum() != 0.0)
                    {
                        H_dpc_cfg_map.values = ((*cfg_2_map) - (*cfg_3_map) - (*cfg_4_map) + (*cfg_5_map)) / t_abs_map;
                        V_dpc_cfg_map.values = ((*cfg_2_map) + (*cfg_3_map) - (*cfg_4_map) - (*cfg_5_map)) / t_abs_map;
                        dia1_dpc_cfg_map.values = ((*cfg_2_map) - (*cfg_4_map)) / t_abs_map;
                        dia2_dpc_cfg_map.values = ((*cfg_3_map) - (*cfg_5_map)) / t_abs_map;
                    }

                }

                scalers_map->push_back(abs_ic_map);
                scalers_map->push_back(abs_cfg_map);
                scalers_map->push_back(H_dpc_cfg_map);
                scalers_map->push_back(V_dpc_cfg_map);
                scalers_map->push_back(dia1_dpc_cfg_map);
                scalers_map->push_back(dia2_dpc_cfg_map);
            }

            if (scalers_map->size() > 0)
            {
                count_3d[0] = scalers_map->size();
                for (const auto &itr : *scalers_map)
                {
                    count_3d[1] = itr.values.rows();
                    count_2d[0] = count_3d[1];
                    count_3d[2] = itr.values.cols();
                    count_2d[1] = count_3d[2];
                    break;
                }

                dataspace_id = H5Screate_simple(3, count_3d, NULL);
                filespace_id = H5Screate_simple(3, count_3d, NULL);

                count_3d[0] = 1;

                dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(dcpl_id, 3, count_3d);
                H5Pset_deflate(dcpl_id, 7);
               
                count_3d[0] = scalers_map->size();
                count[0] = count_3d[0];
                filespace_name_id = H5Screate_simple(1, count, NULL);

                dset_values_id = H5Dopen(scalers_grp_id, "Values", H5P_DEFAULT);
                if (dset_values_id < 0)
                {
                    dset_values_id = H5Dcreate(scalers_grp_id, "Values", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
                    if (dset_values_id < 0)
                    {
                        logE << " Could not open or create /MAPS/Scalers/Values\n";
                        return false;
                    }
                }
                dset_names_id = H5Dopen(scalers_grp_id, "Names", H5P_DEFAULT);
                if (dset_names_id < 0)
                {
                    dset_names_id = H5Dcreate(scalers_grp_id, "Names", filetype, filespace_name_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    if (dset_names_id < 0)
                    {
                        logE << " Could not open or create /MAPS/Scalers/Names\n";
                        return false;
                    }
                }
                dset_units_id = H5Dopen(scalers_grp_id, "Units", H5P_DEFAULT);
                if (dset_units_id < 0)
                {
                    dset_units_id = H5Dcreate(scalers_grp_id, "Units", filetype, filespace_name_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    if (dset_units_id < 0)
                    {
                        logE << " Could not open or create /MAPS/Scalers/Units\n";
                        return false;
                    }
                }
                count_3d[0] = 1;
                count[0] = 1;

                memoryspace_id = H5Screate_simple(2, count_2d, NULL);

                int idx = 0;
                for (auto &itr : *scalers_map)
                {
                    offset[0] = idx;
                    offset_3d[0] = idx;
                    idx++;
                    char tmp_char[255] = { 0 };
                    char tmp_char_units[255] = { 0 };
                    itr.name.copy(tmp_char, 254);
                    itr.unit.copy(tmp_char_units, 254);
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);
                    status = H5Dwrite(dset_units_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char_units);
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    itr.values = itr.values.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });
                    status = H5Dwrite(dset_values_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)itr.values.data());

                }
            }
        }

        H5Tclose(filetype);
        H5Tclose(memtype);

        if (dset_names_id > -1)
        {
            H5Dclose(dset_names_id);
            dset_names_id = -1;
        }
        if (dset_units_id > -1)
        {
            H5Dclose(dset_units_id);
            dset_units_id = -1;
        }
        if (dset_values_id > -1)
        {
            H5Dclose(dset_values_id);
            dset_values_id = -1;
        }
        if (dataspace_id > -1)
        {
            H5Sclose(dataspace_id);
            dataspace_id = -1;
        }
        if (filespace_id > -1)
        {
            H5Sclose(filespace_id);
            filespace_id = -1;
        }
        if (filespace_name_id > -1)
        {
            H5Sclose(filespace_name_id);
            filespace_name_id = -1;
        }
        if (memoryspace_str_id > -1)
        {
            H5Sclose(memoryspace_str_id);
            memoryspace_str_id = -1;
        }
        if (memoryspace_id > -1)
        {
            H5Sclose(memoryspace_id);
            memoryspace_id = -1;
        }
        if (dcpl_id > -1)
        {
            H5Pclose(dcpl_id);
            dcpl_id = -1;
        }
        if (scalers_grp_id > -1)
        {
            H5Gclose(scalers_grp_id);
            scalers_grp_id = -1;
        }
    }
    catch (...)
    {
        logE << "creating group MAPS/Scalers" << "\n";
        return false;
    }
    logI << "Done" << "\n";

    return true;
}

//-----------------------------------------------------------------------------

void HDF5_IO::_save_amps(hid_t scalers_grp_id, data_struct::Params_Override * params_override)
{
    
    hid_t dataspace_id = -1, memoryspace_id = -1;
    hid_t dset_id = -1;
    hid_t status;
    hid_t filetype, memtype;

    char tmp_char[255] = {0};
    hsize_t offset[1] = { 0 };
    hsize_t count[1] = { 3 };

    filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
    memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(memtype, 255);
	std::string units;
    real_t us_amp_sens_num_val = params_override->us_amp_sens_num;
    real_t us_amp_sens_unit_val = params_override->us_amp_sens_unit;

    real_t ds_amp_sens_num_val = params_override->ds_amp_sens_num;
    real_t ds_amp_sens_unit_val = params_override->ds_amp_sens_unit;

    real_t trans_us_amp_sens_num_val;
    std::string trans_us_amp_sens_unit;
    real_t trans_ds_amp_sens_num_val;
    std::string trans_ds_amp_sens_unit;

    switch((int)us_amp_sens_num_val)
    {
        case 0:
            trans_us_amp_sens_num_val = 1;
            break;
        case 1:
            trans_us_amp_sens_num_val = 2;
            break;
        case 2:
            trans_us_amp_sens_num_val = 5;
            break;
        case 3:
            trans_us_amp_sens_num_val = 10;
            break;
        case 4:
            trans_us_amp_sens_num_val = 20;
            break;
        case 5:
            trans_us_amp_sens_num_val = 50;
            break;
        case 6:
            trans_us_amp_sens_num_val = 100;
            break;
        case 7:
            trans_us_amp_sens_num_val = 200;
            break;
        case 8:
            trans_us_amp_sens_num_val = 500;
            break;
    }
    switch((int)ds_amp_sens_num_val)
    {
        case 0:
            trans_ds_amp_sens_num_val = 1;
            break;
        case 1:
            trans_ds_amp_sens_num_val = 2;
            break;
        case 2:
            trans_ds_amp_sens_num_val = 5;
            break;
        case 3:
            trans_ds_amp_sens_num_val = 10;
            break;
        case 4:
            trans_ds_amp_sens_num_val = 20;
            break;
        case 5:
            trans_ds_amp_sens_num_val = 50;
            break;
        case 6:
            trans_ds_amp_sens_num_val = 100;
            break;
        case 7:
            trans_ds_amp_sens_num_val = 200;
            break;
        case 8:
            trans_ds_amp_sens_num_val = 500;
            break;
    }

    switch((int)us_amp_sens_unit_val)
    {
        case 0:
            trans_us_amp_sens_unit= "pA/V";
            break;
        case 1:
            trans_us_amp_sens_unit= "nA/V";
            break;
        case 2:
            trans_us_amp_sens_unit= "uA/V";
            break;
        case 3:
            trans_us_amp_sens_unit= "mA/V";
            break;
    }

    switch((int)ds_amp_sens_unit_val)
    {
        case 0:
            trans_ds_amp_sens_unit= "pA/V";
            break;
        case 1:
            trans_ds_amp_sens_unit= "nA/V";
            break;
        case 2:
            trans_ds_amp_sens_unit= "uA/V";
            break;
        case 3:
            trans_ds_amp_sens_unit= "mA/V";
            break;
    }

    dataspace_id = H5Screate_simple(1, count, NULL);

    count[0] = 1;
    memoryspace_id = H5Screate_simple(1, count, NULL);

    dset_id = H5Dcreate(scalers_grp_id, "us_amp", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&us_amp_sens_num_val);

    offset[0] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&us_amp_sens_unit_val);

    offset[0] = 2;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&trans_us_amp_sens_num_val);

    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }


    offset[0] = 0;
    dset_id = H5Dcreate(scalers_grp_id, "ds_amp", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&ds_amp_sens_num_val);

    offset[0] = 1;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&ds_amp_sens_unit_val);

    offset[0] = 2;
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&trans_ds_amp_sens_num_val);
    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }
    if(dataspace_id > -1)
    {
        H5Sclose(dataspace_id);
        dataspace_id = -1;
    }

    offset[0] = 0;
    count[0] = 1;
    dataspace_id = H5Screate_simple(1, count, NULL);

    dset_id = H5Dcreate(scalers_grp_id, "us_amp_num", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&trans_us_amp_sens_num_val);
    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }

    dset_id = H5Dcreate(scalers_grp_id, "ds_amp_num", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&trans_ds_amp_sens_num_val);
    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }


    trans_us_amp_sens_unit.copy(tmp_char, 254);
    dset_id = H5Dcreate(scalers_grp_id, "us_amp_unit", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }

    trans_ds_amp_sens_unit.copy(tmp_char, 254);
    dset_id = H5Dcreate(scalers_grp_id, "ds_amp_unit", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);


    if(dset_id > -1)
    {
        H5Dclose(dset_id);
        dset_id = -1;
    }
    if(dataspace_id > -1)
    {
        H5Sclose(dataspace_id);
        dataspace_id = -1;
    }
    if(memoryspace_id > -1)
    {
        H5Sclose(memoryspace_id);
        memoryspace_id = -1;
    }


    H5Tclose(filetype);
    H5Tclose(memtype);
    
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_scan_scalers(size_t detector_num,
                                data_struct::Scan_Info *scan_info,
                                data_struct::Params_Override * params_override,
                                bool hasNetcdf,
                                size_t row_idx_start,
                                int row_idx_end,
                                size_t col_idx_start,
                                int col_idx_end)
{

    std::lock_guard<std::mutex> lock(_mutex);
	std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t scan_grp_id, maps_grp_id, po_grp_id;

	if (scan_info == nullptr)
    {
        logW << "scalers_map == nullptr. Not returning from save_scan_scalers" << "\n";
		return false;
    }

    if(_cur_file_id < 0)
    {
        logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
        return false;
    }

    logI << "Saving scalers to hdf5"<< "\n";

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if (maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (maps_grp_id < 0)
    {
        logE << "creating group 'MAPS'" << "\n";
        return false;
    }

    scan_grp_id = H5Gopen(maps_grp_id, "Scan", H5P_DEFAULT);
    if (scan_grp_id < 0)
        scan_grp_id = H5Gcreate(maps_grp_id, "Scan", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scan_grp_id < 0)
    {
        logE << "creating group MAPS/Scan" << "\n";
        return false;
    }

	po_grp_id = H5Gopen(maps_grp_id, "Fit_Parameters_Override", H5P_DEFAULT);
	if (po_grp_id < 0)
		po_grp_id = H5Gcreate(maps_grp_id, "Fit_Parameters_Override", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (po_grp_id < 0)
	{
		logE << "creating group MAPS/Fit_Parameters_Override" << "\n";
		return false;
	}

	_save_params_override(po_grp_id, params_override);

    _save_scan_meta_data(scan_grp_id, &(scan_info->meta_info));
	
    _save_extras(scan_grp_id, &(scan_info->extra_pvs));
	
    _save_scalers(maps_grp_id, &(scan_info->scaler_maps), params_override, hasNetcdf);

	H5Gclose(po_grp_id);
    H5Gclose(scan_grp_id);
    H5Gclose(maps_grp_id);

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_scan_scalers_confocal(std::string path,
                                size_t detector_num,
                                size_t row_idx_start,
                                int row_idx_end,
                                size_t col_idx_start,
                                int col_idx_end)
{

    std::lock_guard<std::mutex> lock(_mutex);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    hid_t scan_grp_id, maps_grp_id, scalers_grp_id, status, error;
    hid_t    file_id, src_maps_grp_id;
    hid_t    dataspace_detectors_id, dset_detectors_id, attr_detector_names_id;
    hid_t   xpos_dataspace_id, xpos_id, ypos_dataspace_id, ypos_id;
    hid_t dataspace_id, dataset_id;
    char* detector_names[256];
    int det_rank;
    hsize_t* det_dims_in;
    hsize_t names_off[1] = { 0 };
    hsize_t names_cnt[1] = { 1 };
    hsize_t scalers_offset[3] = { 0,0,0 };
	hsize_t scalers_count[3] = { 1,1,1 };
	hsize_t value_offset[3] = { 0,0,0 };
	hsize_t value_count[3] = { 1,1,1 };
	hsize_t mem_count[2] = { 1,1 };
    hsize_t x_offset[3] = {0,0,0};
    hsize_t x_count[3] = {1,1,1};
    hsize_t y_offset[2] = {0,0};
    hsize_t y_count[2] = {1,1};
    void *f_data;
    bool confocal_ver_2020 = false;

    if(_cur_file_id < 0)
    {
        logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
        return false;
    }

    logI << "Saving scalers to hdf5"<< "\n";

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if (maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (maps_grp_id < 0)
    {
        logE << "creating group 'MAPS'" << "\n";
        return false;
    }
	close_map.push({ maps_grp_id, H5O_GROUP });

    scan_grp_id = H5Gopen(maps_grp_id, "Scan", H5P_DEFAULT);
    if (scan_grp_id < 0)
        scan_grp_id = H5Gcreate(maps_grp_id, "Scan", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scan_grp_id < 0)
    {
		_close_h5_objects(close_map);
        logE << "creating group MAPS/Scan" << "\n";
        return false;
    }
	close_map.push({ scan_grp_id, H5O_GROUP });

    scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
    if (scalers_grp_id < 0)
        scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scalers_grp_id < 0)
    {
		_close_h5_objects(close_map);
        logE << "creating group MAPS/Scalers" << "\n";
        return false;
    }
	close_map.push({ scalers_grp_id, H5O_GROUP });

    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
    {
        _close_h5_objects(close_map);
        return false;
    }

    if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "2D Scan", file_id))
    {
        _close_h5_objects(close_map);
        return false;
    }

    if (false == _open_h5_object(xpos_id, H5O_DATASET, close_map, "X Positions", src_maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    xpos_dataspace_id = H5Dget_space(xpos_id);
    close_map.push({xpos_dataspace_id, H5O_DATASPACE});

    if (false == _open_h5_object(ypos_id, H5O_DATASET, close_map, "Y Positions", src_maps_grp_id))
    {
        _close_h5_objects(close_map);
        return false;
    }
    ypos_dataspace_id = H5Dget_space(ypos_id);
    close_map.push({ypos_dataspace_id, H5O_DATASPACE});

    hid_t xy_type = H5Dget_type(xpos_id);
    det_rank = H5Sget_simple_extent_ndims(xpos_dataspace_id);
    det_dims_in = new hsize_t[det_rank];
    H5Sget_simple_extent_dims(xpos_dataspace_id, &det_dims_in[0], NULL);

    dataspace_id = H5Screate_simple (1, &det_dims_in[1], NULL);
    dataset_id = H5Dcreate(scan_grp_id, "x_axis", xy_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    x_count[1] = det_dims_in[1];
    f_data = malloc(sizeof(float) * det_dims_in[1]);
    H5Sselect_hyperslab (xpos_dataspace_id, H5S_SELECT_SET, x_offset, NULL, x_count, NULL);
    status = H5Dread (xpos_id, xy_type, dataspace_id, xpos_dataspace_id, H5P_DEFAULT, f_data);
    status = H5Dwrite (dataset_id, xy_type, H5S_ALL, dataspace_id, H5P_DEFAULT, f_data);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    free(f_data);


    dataspace_id = H5Screate_simple (1, &det_dims_in[0], NULL);
    dataset_id = H5Dcreate(scan_grp_id, "y_axis", xy_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    y_count[1] = det_dims_in[0];
    f_data = malloc(sizeof(float) * det_dims_in[0]);
    H5Sselect_hyperslab (ypos_dataspace_id, H5S_SELECT_SET, y_offset, NULL, y_count, NULL);
    status = H5Dread (ypos_id, xy_type, dataspace_id, ypos_dataspace_id, H5P_DEFAULT, f_data);
    status = H5Dwrite (dataset_id, xy_type, H5S_ALL, dataspace_id, H5P_DEFAULT, f_data);
    H5Sclose(dataspace_id);
    H5Dclose(dataset_id);
    free(f_data);
    delete [] det_dims_in;

    //Save scalers
    if (false == _open_h5_object(dset_detectors_id, H5O_DATASET, close_map, "Detectors", src_maps_grp_id, false, false))
    {
        // try Detectors as group since it changed in 2020
        if (false == _open_h5_object(dset_detectors_id, H5O_GROUP, close_map, "Detectors", src_maps_grp_id))
        {
            _close_h5_objects(close_map);
            return false;
        }
        else
        {
            confocal_ver_2020 = true;
        }
    }
    if (confocal_ver_2020)
    {
        bool first_save = true;
        real_t* buffer = nullptr;
        hsize_t nobj = 0;
        H5Gget_num_objs(dset_detectors_id, &nobj);
        hid_t values_id;
        hid_t value_space;
        hid_t mem_space;
        hid_t names_id;
        hid_t name_space;
        hid_t mem_name_space = H5Screate_simple(1, &names_cnt[0], &names_cnt[0]);
        close_map.push({ mem_name_space, H5O_DATASPACE });
        hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
        hid_t name_type = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(name_type, 255);

        for (hsize_t i = 0; i < nobj; i++)
        {
            char str_dset_name[2048] = { 0 };
            hid_t dsid;
            hsize_t len = H5Gget_objname_by_idx(dset_detectors_id, i, str_dset_name, 2048);
            if (_open_h5_object(dsid, H5O_DATASET, close_map, str_dset_name, dset_detectors_id))
            {
                hid_t scaler_space = H5Dget_space(dsid);
                close_map.push({ scaler_space, H5O_DATASPACE });
                hid_t scalers_type = H5Dget_type(dsid);
                names_off[0] = i;
                if (first_save)
                {
                    H5Sget_simple_extent_dims(scaler_space, &scalers_count[0], NULL);
                    value_count[0] = nobj;
                    value_count[1] = scalers_count[0];
                    value_count[2] = scalers_count[1];
                    value_space = H5Screate_simple(3, &value_count[0], &value_count[0]);
                    close_map.push({ value_space, H5O_DATASPACE });
                    // create values
                    values_id = H5Dcreate(scalers_grp_id, "Values", scalers_type, value_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    close_map.push({ values_id, H5O_DATASET });
                    // create names 
                    value_count[0] = 1;
                    names_cnt[0] = nobj;
                    name_space = H5Screate_simple(1, &names_cnt[0], NULL);
                    names_id = H5Dcreate(scalers_grp_id, "Names", name_type, name_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    close_map.push({ names_id, H5O_DATASET });
                    // reset name_cnt to 1 for writing 
                    names_cnt[0] = 1;
                    buffer = new real_t[scalers_count[0] * scalers_count[1]];
                    mem_count[0] = scalers_count[0];
                    mem_count[1] = scalers_count[1];
                    mem_space = H5Screate_simple(2, mem_count, mem_count);
                    close_map.push({ mem_space, H5O_DATASPACE });
                    first_save = false;
                }
                value_offset[0] = i;
                H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, scalers_offset, nullptr, scalers_count, nullptr);
                H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
                H5Sselect_hyperslab(name_space, H5S_SELECT_SET, names_off, nullptr, names_cnt, nullptr);
                status = H5Dread(dsid, scalers_type, mem_space, scaler_space, H5P_DEFAULT, &buffer[0]);
                status = H5Dwrite(values_id, scalers_type, mem_space, value_space, H5P_DEFAULT, &buffer[0]);
                status = H5Dwrite(names_id, name_type, mem_name_space, name_space, H5P_DEFAULT, str_dset_name);
            }
        }

        if (buffer != nullptr)
        {
            delete[] buffer;
        }
    }
    else
    {
        hid_t scaler_space = H5Dget_space(dset_detectors_id);
        close_map.push({ scaler_space, H5O_DATASPACE });
        H5Sget_simple_extent_dims(scaler_space, &scalers_count[0], NULL);
        value_count[0] = scalers_count[2];
        value_count[1] = scalers_count[0];
        value_count[2] = scalers_count[1];
        hid_t value_space = H5Screate_simple(3, &value_count[0], &value_count[0]);
        close_map.push({ value_space, H5O_DATASPACE });
        hid_t scalers_type = H5Dget_type(dset_detectors_id);
        hid_t values_id = H5Dcreate(scalers_grp_id, "Values", scalers_type, value_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        close_map.push({ values_id, H5O_DATASET });
        real_t* buffer = new real_t[scalers_count[0] * scalers_count[1]];
        hsize_t scaler_amt = scalers_count[2];
        scalers_count[2] = 1;
        value_count[0] = 1;
        mem_count[0] = scalers_count[0];
        mem_count[1] = scalers_count[1];
        hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);
        close_map.push({ mem_space, H5O_DATASPACE });

        for (hsize_t s = 0; s < scaler_amt; s++)
        {
            value_offset[0] = s;
            scalers_offset[2] = s;
            H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, scalers_offset, nullptr, scalers_count, nullptr);
            H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
            status = H5Dread(dset_detectors_id, scalers_type, mem_space, scaler_space, H5P_DEFAULT, &buffer[0]);
            status = H5Dwrite(values_id, scalers_type, mem_space, value_space, H5P_DEFAULT, &buffer[0]);
        }
        delete[] buffer;

        //save the scalers names
        dataspace_detectors_id = H5Dget_space(dset_detectors_id);
        close_map.push({ dataspace_detectors_id, H5O_DATASPACE });
        attr_detector_names_id = H5Aopen(dset_detectors_id, "Detector Names", H5P_DEFAULT);
        close_map.push({ dataspace_detectors_id, H5O_ATTRIBUTE });

        det_rank = H5Sget_simple_extent_ndims(dataspace_detectors_id);
        det_dims_in = new hsize_t[det_rank];
        H5Sget_simple_extent_dims(dataspace_detectors_id, &det_dims_in[0], NULL);

        hid_t ftype = H5Aget_type(attr_detector_names_id);
        hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
        error = H5Aread(attr_detector_names_id, type, detector_names);

        if (error == 0)
        {
            dataspace_id = H5Screate_simple(1, &det_dims_in[2], NULL);
            dataset_id = H5Dcreate(scalers_grp_id, "Names", type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            status = H5Dwrite(dataset_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, detector_names);
            H5Sclose(dataspace_id);
            H5Dclose(dataset_id);
            //for (size_t z = 0; z < det_dims_in[2]; z++)
            //{
                //TODO: look into why this is causing exception in windows
                
            //}
            //free(detector_names);
        }
        delete[] det_dims_in;
    }

    _close_h5_objects(close_map);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_scan_scalers_gsecars(std::string path,
									size_t detector_num,
									size_t row_idx_start,
									int row_idx_end,
									size_t col_idx_start,
									int col_idx_end)
{

	std::lock_guard<std::mutex> lock(_mutex);
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

	hid_t scan_grp_id, maps_grp_id, scalers_grp_id, status, error;
	hid_t    file_id, src_maps_grp_id;
	hid_t    dataspace_detectors_id, dset_detectors_id;
	hid_t   xypos_dataspace_id, xypos_id;
	hid_t x_dataspace_id, y_dataspace_id, x_dataset_id, y_dataset_id;
	char* detector_names[256];
	int det_rank;
	hsize_t* det_dims_in;
	hsize_t single_offset[1] = { 0 };
	hsize_t single_count[1] = { 1 };
	hsize_t value_offset[3] = { 0,0,0 };
	hsize_t value_count[3] = { 1,1,1 };
	hsize_t mem_offset[2] = { 0,0 };
	hsize_t mem_count[2] = { 1,1 };
	hsize_t xy_offset[3] = { 0,0,0 };
	hsize_t xy_count[3] = { 1,1,1 };
	GSE_CARS_SAVE_VER version = GSE_CARS_SAVE_VER::UNKNOWN;
	float *x_data;
	float *y_data;

	if (_cur_file_id < 0)
	{
		logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
		return false;
	}
	
	logI << "Saving scalers to hdf5" << "\n";

	maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
	if (maps_grp_id < 0)
		maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (maps_grp_id < 0)
	{
		logE << "creating group 'MAPS'" << "\n";
		return false;
	}

	close_map.push({ maps_grp_id, H5O_GROUP });

	scan_grp_id = H5Gopen(maps_grp_id, "Scan", H5P_DEFAULT);
	if (scan_grp_id < 0)
		scan_grp_id = H5Gcreate(maps_grp_id, "Scan", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (scan_grp_id < 0)
	{
		_close_h5_objects(close_map);
		logE << "creating group MAPS/Scan" << "\n";
		return false;
	}

	close_map.push({ scan_grp_id, H5O_GROUP });

	scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
	if (scalers_grp_id < 0)
		scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (scalers_grp_id < 0)
	{
		_close_h5_objects(close_map);
		logE << "creating group MAPS/Scalers" << "\n";
		return false;
	}
	close_map.push({ scalers_grp_id, H5O_GROUP });

	if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
		return false;

	if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "xrmmap", file_id, true, false))
	{
		if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "xrfmap", file_id))
		{
			return false;
		}
		version = GSE_CARS_SAVE_VER::XRFMAP;
	}
	else
	{
		version = GSE_CARS_SAVE_VER::XRMMAP;
	}

	if (false == _open_h5_object(xypos_id, H5O_DATASET, close_map, "positions/pos", src_maps_grp_id))
		return false;
	xypos_dataspace_id = H5Dget_space(xypos_id);
	close_map.push({ xypos_dataspace_id, H5O_DATASPACE });

	hid_t xy_type = H5Dget_type(xypos_id);
	det_rank = H5Sget_simple_extent_ndims(xypos_dataspace_id);
	det_dims_in = new hsize_t[det_rank];
	H5Sget_simple_extent_dims(xypos_dataspace_id, &det_dims_in[0], NULL);

	x_dataspace_id = H5Screate_simple(1, &det_dims_in[1], NULL);
	y_dataspace_id = H5Screate_simple(1, &det_dims_in[0], NULL);
	x_dataset_id = H5Dcreate(scan_grp_id, "x_axis", xy_type, x_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	y_dataset_id = H5Dcreate(scan_grp_id, "y_axis", xy_type, y_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	close_map.push({ x_dataspace_id, H5O_DATASPACE });
	close_map.push({ y_dataspace_id, H5O_DATASPACE });
	close_map.push({ x_dataset_id, H5O_DATASET });
	close_map.push({ y_dataset_id, H5O_DATASET });

	xy_offset[2] = 0;
	xy_count[1] = det_dims_in[1];
	x_data = new float[det_dims_in[1]];
	y_data = new float[det_dims_in[0]];

	mem_count[0] = det_dims_in[0];
	mem_count[1] = det_dims_in[1];

	status = H5Sselect_hyperslab(xypos_dataspace_id, H5S_SELECT_SET, xy_offset, NULL, xy_count, NULL);
	status = H5Dread(xypos_id, xy_type, x_dataspace_id, xypos_dataspace_id, H5P_DEFAULT, &x_data[0]);
	if (status > -1)
	{
		status = H5Dwrite(x_dataset_id, xy_type, H5S_ALL, x_dataspace_id, H5P_DEFAULT, &x_data[0]);
	}

	xy_offset[2] = 1;
	xy_count[1] = 1;
	xy_count[0] = det_dims_in[0];

	status = H5Sselect_hyperslab(xypos_dataspace_id, H5S_SELECT_SET, xy_offset, NULL, xy_count, NULL);
	status = H5Dread(xypos_id, xy_type, y_dataspace_id, xypos_dataspace_id, H5P_DEFAULT, &y_data[0]);
	if (status > -1)
	{
		status = H5Dwrite(y_dataset_id, xy_type, H5S_ALL, y_dataspace_id, H5P_DEFAULT, &y_data[0]);
	}

	//Save scalers
	//names
	std::vector<std::string> names_array = { STR_US_IC, STR_DS_IC, "I2", "TSCALER" };
    std::vector<std::string> units_array = { "cts/us", "cts/us", "cts/us", "us" };
	hid_t mem_single_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
	single_count[0] = 4;
	hid_t name_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
	single_count[0] = 1;
	hid_t dtype = H5Tcopy(H5T_C_S1);
	H5Tset_size(dtype, 255);
	hid_t names_id = H5Dcreate(scalers_grp_id, "Names", dtype, name_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	close_map.push({ names_id, H5O_DATASET });

    hid_t units_id = H5Dcreate(scalers_grp_id, "Units", dtype, name_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    close_map.push({ units_id, H5O_DATASET });

	for (int i = 0; i < 4; i++)
	{
		char tmp_char[255] = { 0 };
        char unit_char[255] = { 0 };
		single_offset[0] = i;
		names_array[i].copy(tmp_char, 254);
        units_array[i].copy(unit_char, 254);
		H5Sselect_hyperslab(name_space, H5S_SELECT_SET, single_offset, nullptr, single_count, nullptr);
		H5Dwrite(names_id, dtype, mem_single_space, name_space, H5P_DEFAULT, (void*)tmp_char);
        H5Dwrite(units_id, dtype, mem_single_space, name_space, H5P_DEFAULT, (void*)unit_char);
	}
	//values

	value_count[0] = 4;
	value_count[1] = det_dims_in[0];
	value_count[2] = det_dims_in[1];
	hid_t value_space = H5Screate_simple(3, &value_count[0], &value_count[0]);
	hid_t values_id = H5Dcreate(scalers_grp_id, "Values", H5T_NATIVE_REAL, value_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	close_map.push({ values_id, H5O_DATASET });
	real_t *buffer = new real_t[det_dims_in[0] * det_dims_in[1]];
	value_count[0] = 1;


	if (version == GSE_CARS_SAVE_VER::XRMMAP)
	{
        hid_t upstream_ic_id = H5Dopen(src_maps_grp_id, "scalars/I0_raw", H5P_DEFAULT);
        hid_t downstream_ic_id = H5Dopen(src_maps_grp_id, "scalars/I1_raw", H5P_DEFAULT);
        hid_t i2_id = H5Dopen(src_maps_grp_id, "scalars/I2_raw", H5P_DEFAULT);
        hid_t tscaler_id = H5Dopen(src_maps_grp_id, "scalars/TSCALER_raw", H5P_DEFAULT);

		close_map.push({ upstream_ic_id, H5O_DATASET });
		close_map.push({ downstream_ic_id, H5O_DATASET });
		close_map.push({ i2_id, H5O_DATASET });
		close_map.push({ tscaler_id, H5O_DATASET });

		if (upstream_ic_id > -1)
		{
			hid_t scaler_space = H5Dget_space(upstream_ic_id);

			status = H5Dread(upstream_ic_id, H5T_NATIVE_REAL, scaler_space, scaler_space, H5P_DEFAULT, &buffer[0]);
			if (status > -1)
			{
				value_offset[0] = 0;
				H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
				status = H5Dwrite(values_id, H5T_NATIVE_REAL, scaler_space, value_space, H5P_DEFAULT, &buffer[0]);
			}

			status = H5Dread(downstream_ic_id, H5T_NATIVE_REAL, scaler_space, scaler_space, H5P_DEFAULT, &buffer[0]);
			if (status > -1)
			{
				value_offset[0] = 1;
				H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
				status = H5Dwrite(values_id, H5T_NATIVE_REAL, scaler_space, value_space, H5P_DEFAULT, &buffer[0]);
			}

			status = H5Dread(i2_id, H5T_NATIVE_REAL, scaler_space, scaler_space, H5P_DEFAULT, &buffer[0]);
			if (status > -1)
			{
				value_offset[0] = 2;
				H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
				status = H5Dwrite(values_id, H5T_NATIVE_REAL, scaler_space, value_space, H5P_DEFAULT, &buffer[0]);
			}

			status = H5Dread(tscaler_id, H5T_NATIVE_REAL, scaler_space, scaler_space, H5P_DEFAULT, &buffer[0]);
			if (status > -1)
			{
				value_offset[0] = 3;
				H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
				status = H5Dwrite(values_id, H5T_NATIVE_REAL, scaler_space, value_space, H5P_DEFAULT, &buffer[0]);
			}
		}
	}
	else if (version == GSE_CARS_SAVE_VER::XRFMAP)
	{
		std::map<std::string, int> save_name_idx = { {"I0", 0}, {"I1", 1}, {"I2", 2}, {"TSCALER", 3} };
        hid_t scalers_ds_id = H5Dopen(src_maps_grp_id, "roimap/det_raw", H5P_DEFAULT);
		hid_t scalers_names_id = H5Dopen(src_maps_grp_id, "roimap/det_name", H5P_DEFAULT);
		close_map.push({ scalers_ds_id, H5O_DATASET });
		close_map.push({ scalers_names_id, H5O_DATASET });
		
		hid_t scaler_name_space = H5Dget_space(scalers_names_id);
		hid_t scaler_space = H5Dget_space(scalers_ds_id);

		hid_t buffer_space = H5Screate_simple(2, &mem_count[0], &mem_count[0]);
		close_map.push({ buffer_space, H5O_DATASPACE });

		H5Sget_simple_extent_dims(scaler_name_space, &single_count[0], NULL);
		int amt = single_count[0];
		single_count[0] = 1;

		for (int i = 0; i < amt; i++)
		{
			single_offset[0] = i;
			char tmp_char[256] = { 0 };
			H5Sselect_hyperslab(scaler_name_space, H5S_SELECT_SET, single_offset, nullptr, single_count, nullptr);
			status = H5Dread(scalers_names_id, dtype, mem_single_space, scaler_name_space, H5P_DEFAULT, (void*)tmp_char);
			if (status > -1)
			{
				std::string sname(tmp_char);
				if (save_name_idx.count(sname) > 0)
				{
					xy_offset[0] = 0;
					xy_offset[1] = 0;
					xy_offset[2] = i;
					xy_count[0] = det_dims_in[0];
					xy_count[1] = det_dims_in[1];
					xy_count[2] = 1;

					H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, xy_offset, nullptr, xy_count, nullptr);
					status = H5Dread(scalers_ds_id, H5T_NATIVE_REAL, buffer_space, scaler_space, H5P_DEFAULT, &buffer[0]);
					if (status > -1)
					{
						value_offset[0] = save_name_idx[sname];
						H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
						status = H5Dwrite(values_id, H5T_NATIVE_REAL, buffer_space, value_space, H5P_DEFAULT, &buffer[0]);
					}
				}
			}
		}
	}
	
	hid_t ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
	status = H5Ocopy(src_maps_grp_id, "config", scan_grp_id, "config", ocpypl_id, H5P_DEFAULT);

	_close_h5_objects(close_map);

	delete[] x_data;
	delete[] y_data;
	delete[] det_dims_in;
	delete[] buffer;

	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

	return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_scan_scalers_bnl(std::string path,
    size_t detector_num,
    data_struct::Params_Override* params_override,
    size_t row_idx_start,
    int row_idx_end,
    size_t col_idx_start,
    int col_idx_end)
{

    std::lock_guard<std::mutex> lock(_mutex);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

    hid_t scan_grp_id, maps_grp_id, scalers_grp_id, status, error;
    hid_t    file_id, src_maps_grp_id;
    hid_t    dataspace_detectors_id, dset_detectors_id;
    hid_t   xypos_dataspace_id, xypos_id;
    hid_t x_dataspace_id, y_dataspace_id, x_dataset_id, y_dataset_id;
    char* detector_names[256];
    int det_rank;
    hsize_t *det_dims_in = nullptr;
    hsize_t* val_dims_in = nullptr;;
    hsize_t scaler_offset[3] = { 0,0,0 };
    hsize_t single_offset[1] = { 0 };
    hsize_t single_count[1] = { 1 };
    hsize_t value_offset[3] = { 0,0,0 };
    hsize_t value_count[3] = { 1,1,1 };
    hsize_t mem_offset[2] = { 0,0 };
    hsize_t mem_count[2] = { 1,1 };
    hsize_t xy_offset[3] = { 0,0,0 };
    hsize_t xy_count[3] = { 1,1,1 };
    hid_t ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
    double* x_data = nullptr;
    double* y_data = nullptr;
    double* buffer = nullptr;

    if (_cur_file_id < 0)
    {
        logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
        return false;
    }

    logI << "Saving scalers to hdf5" << "\n";

    maps_grp_id = H5Gopen(_cur_file_id, "MAPS", H5P_DEFAULT);
    if (maps_grp_id < 0)
        maps_grp_id = H5Gcreate(_cur_file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (maps_grp_id < 0)
    {
        logE << "creating group 'MAPS'" << "\n";
        return false;
    }

    close_map.push({ maps_grp_id, H5O_GROUP });

    scan_grp_id = H5Gopen(maps_grp_id, "Scan", H5P_DEFAULT);
    if (scan_grp_id < 0)
        scan_grp_id = H5Gcreate(maps_grp_id, "Scan", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scan_grp_id < 0)
    {
        _close_h5_objects(close_map);
        logE << "creating group MAPS/Scan" << "\n";
        return false;
    }

    close_map.push({ scan_grp_id, H5O_GROUP });

    scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
    if (scalers_grp_id < 0)
        scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scalers_grp_id < 0)
    {
        _close_h5_objects(close_map);
        logE << "creating group MAPS/Scalers" << "\n";
        return false;
    }
    close_map.push({ scalers_grp_id, H5O_GROUP });

    if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
        return false;

    
    if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "xrfmap", file_id))
    {
        return false;
    }

    if (false == _open_h5_object(xypos_id, H5O_DATASET, close_map, "positions/pos", src_maps_grp_id))
        return false;
    xypos_dataspace_id = H5Dget_space(xypos_id);
    close_map.push({ xypos_dataspace_id, H5O_DATASPACE });

    hid_t xy_type = H5Dget_type(xypos_id);
    det_rank = H5Sget_simple_extent_ndims(xypos_dataspace_id);
    det_dims_in = new hsize_t[det_rank];
    H5Sget_simple_extent_dims(xypos_dataspace_id, &det_dims_in[0], NULL);

    x_dataspace_id = H5Screate_simple(1, &det_dims_in[1], NULL);
    y_dataspace_id = H5Screate_simple(1, &det_dims_in[2], NULL);
    x_dataset_id = H5Dcreate(scan_grp_id, "x_axis", xy_type, x_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    y_dataset_id = H5Dcreate(scan_grp_id, "y_axis", xy_type, y_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    close_map.push({ x_dataspace_id, H5O_DATASPACE });
    close_map.push({ y_dataspace_id, H5O_DATASPACE });
    close_map.push({ x_dataset_id, H5O_DATASET });
    close_map.push({ y_dataset_id, H5O_DATASET });

    xy_offset[0] = 0;
    xy_count[1] = det_dims_in[1];
    x_data = new double[det_dims_in[1]];
    y_data = new double[det_dims_in[2]];

    status = H5Sselect_hyperslab(xypos_dataspace_id, H5S_SELECT_SET, xy_offset, NULL, xy_count, NULL);
    status = H5Dread(xypos_id, xy_type, x_dataspace_id, xypos_dataspace_id, H5P_DEFAULT, &x_data[0]);
    if (status > -1)
    {
        status = H5Dwrite(x_dataset_id, xy_type, x_dataspace_id, x_dataspace_id, H5P_DEFAULT, &x_data[0]);
    }

    xy_offset[0] = 1;
    xy_count[1] = 1;
    xy_count[2] = det_dims_in[2];

    status = H5Sselect_hyperslab(xypos_dataspace_id, H5S_SELECT_SET, xy_offset, NULL, xy_count, NULL);
    status = H5Dread(xypos_id, xy_type, y_dataspace_id, xypos_dataspace_id, H5P_DEFAULT, &y_data[0]);
    if (status > -1)
    {
        status = H5Dwrite(y_dataset_id, xy_type, y_dataspace_id, y_dataspace_id, H5P_DEFAULT, &y_data[0]);
    }

    //Save scalers
    // name, val
    hid_t scaler_name_id, scaler_val_id, scaler_grp_id;
    if (false == _open_h5_object(scaler_grp_id, H5O_GROUP, close_map, "scalers", src_maps_grp_id))
    {
        return false;
    }
    if (false == _open_h5_object(scaler_name_id, H5O_DATASET, close_map, "name", scaler_grp_id))
    {
        return false;
    }
    if (false == _open_h5_object(scaler_val_id, H5O_DATASET, close_map, "val", scaler_grp_id))
    {
        return false;
    }

    hid_t name_type = H5Dget_type(scaler_name_id);
    hid_t scaler_type = H5Dget_type(scaler_val_id);
    hid_t scaler_val_space = H5Dget_space(scaler_val_id);
    hid_t val_rank = H5Sget_simple_extent_ndims(scaler_val_space);
    val_dims_in = new hsize_t[val_rank];
    H5Sget_simple_extent_dims(scaler_val_space, &val_dims_in[0], NULL);

    value_count[0] = val_dims_in[2];
    value_count[1] = val_dims_in[0];
    value_count[2] = val_dims_in[1];
    mem_count[0] = val_dims_in[0];
    mem_count[1] = val_dims_in[1];
    //names
    hid_t mem_single_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
    single_count[0] = val_dims_in[2];
    hid_t name_space = H5Screate_simple(1, &single_count[0], &single_count[0]);
    single_count[0] = 1;
    hid_t dtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(dtype, 255);
    hid_t names_id = H5Dcreate(scalers_grp_id, "Names", dtype, name_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    close_map.push({ names_id, H5O_DATASET });

    hid_t units_id = H5Dcreate(scalers_grp_id, "Units", dtype, name_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    close_map.push({ units_id, H5O_DATASET });

    hid_t value_space = H5Screate_simple(3, &value_count[0], &value_count[0]);
    close_map.push({ value_space, H5O_DATASPACE });
    hid_t values_id = H5Dcreate(scalers_grp_id, "Values", scaler_type, value_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    close_map.push({ values_id, H5O_DATASET });
    

    char tmp_char[255] = { 0 };
    char unit_char[255] = "cts";
    buffer = new double[val_dims_in[0] * val_dims_in[1]];
    hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);
    close_map.push({ mem_space, H5O_DATASPACE });
    size_t scaler_cnt = val_dims_in[2];
    val_dims_in[2] = 1;
    value_count[0] = 1;
    for (hsize_t i = 0; i < scaler_cnt; i++)
    {
        single_offset[0] = i;
        value_offset[0] = i;
        scaler_offset[2] = i;

        H5Sselect_hyperslab(scaler_val_space, H5S_SELECT_SET, scaler_offset, NULL, val_dims_in, NULL);
        H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, NULL, value_count, NULL);
        H5Sselect_hyperslab(name_space, H5S_SELECT_SET, single_offset, nullptr, single_count, nullptr);

        status = H5Dread(scaler_name_id, name_type, mem_single_space, name_space, H5P_DEFAULT, (void*)tmp_char);
        if (status > -1)
        {
            std::string read_name = std::string(tmp_char, 255);
            read_name.erase(std::remove_if(read_name.begin(), read_name.end(), ::isspace), read_name.end());
            read_name.erase(std::find(read_name.begin(), read_name.end(), '\0'), read_name.end());
            for (const auto& s_itr : params_override->scaler_pvs)
            {
                if (read_name == s_itr.second)
                {
                    for (int j = 0; j < 255; j++)
                    {
                        tmp_char[j] = '\0';
                    }
                    s_itr.first.copy(tmp_char, 255);
                    break;
                }
            }
            H5Dwrite(names_id, dtype, mem_single_space, name_space, H5P_DEFAULT, (void*)tmp_char);
        }
        H5Dwrite(units_id, dtype, mem_single_space, name_space, H5P_DEFAULT, (void*)unit_char);
        status = H5Dread(scaler_val_id, scaler_type, mem_space, scaler_val_space, H5P_DEFAULT, (void*)buffer);
        if (status > -1)
        {
            H5Dwrite(values_id, scaler_type, mem_space, value_space, H5P_DEFAULT, (void*)buffer);
        }
    }

    // todo: save scan meta data

    _close_h5_objects(close_map);

    if (x_data != nullptr)
        delete[] x_data;
    if (y_data != nullptr)
        delete[] y_data;
    if(det_dims_in != nullptr)
        delete[] det_dims_in;
    if (val_dims_in != nullptr)
        delete[] val_dims_in;
    if (buffer != nullptr)
        delete[] buffer;

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::generate_avg(std::string avg_filename, std::vector<std::string> files_to_avg)
{
    std::lock_guard<std::mutex> lock(_mutex);
    logI  << avg_filename << "\n";

    hid_t ocpypl_id, status, src_maps_grp_id, src_analyzed_grp_id, dst_fit_grp_id, src_quant_grp_id, dst_quant_grp_id;
    std::vector<hid_t> hdf5_file_ids;
    std::string group_name = "";

    //save to restore later
    hid_t saved_file_id = _cur_file_id;

    for(auto& filename : files_to_avg)
    {
        hid_t    rd_file_id;

        rd_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (rd_file_id > -1)
        {
            hdf5_file_ids.push_back(rd_file_id);
        }
    }

    if(hdf5_file_ids.size() > 1)
    {
        hid_t file_id = H5Fcreate(avg_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

        hid_t dst_maps_grp_id = H5Gcreate(file_id, "MAPS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
        status = H5Pset_copy_object(ocpypl_id, H5O_COPY_MERGE_COMMITTED_DTYPE_FLAG);

        src_maps_grp_id = H5Gopen(hdf5_file_ids[0], "MAPS", H5P_DEFAULT);
        //Copy Scan and Scalers
        if(src_maps_grp_id > -1)
        {
            status = H5Ocopy(src_maps_grp_id, "Scan", dst_maps_grp_id, "Scan", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_maps_grp_id, "Scalers", dst_maps_grp_id, "Scalers", ocpypl_id, H5P_DEFAULT);
        }

        //Average XRF_Analyzed
        _generate_avg_analysis(src_maps_grp_id, dst_maps_grp_id, "MAPS", ocpypl_id, hdf5_file_ids);

        //Average Spectra
        group_name = "Spectra";
        src_analyzed_grp_id = H5Gopen(src_maps_grp_id, group_name.c_str(), H5P_DEFAULT);
        if(src_analyzed_grp_id > -1)
        {
            dst_fit_grp_id = H5Gcreate(dst_maps_grp_id, group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            _gen_average("/MAPS/"+group_name+"/Elapsed_Livetime", "Elapsed_Livetime", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
            _gen_average("/MAPS/"+group_name+"/Elapsed_Realtime", "Elapsed_Realtime", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
            _gen_average("/MAPS/"+group_name+"/Input_Counts", "Input_Counts", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
            _gen_average("/MAPS/"+group_name+"/Output_Counts", "Output_Counts", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
			
            //for now copy the first detector energy and energy_calib, in the future would be nice to normalize over 4 detectors
            status = H5Ocopy(src_analyzed_grp_id, "Energy", dst_fit_grp_id, "Energy", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_analyzed_grp_id, "Energy_Calibration", dst_fit_grp_id, "Energy_Calibration", ocpypl_id, H5P_DEFAULT);
            //sum mca_arr
            _gen_average("/MAPS/"+group_name+"/mca_arr", "mca_arr", src_analyzed_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids, false);


            _generate_avg_integrated_spectra(src_analyzed_grp_id, dst_fit_grp_id, "MAPS/Spectra", ocpypl_id, hdf5_file_ids);

            H5Gclose(dst_fit_grp_id);
            H5Gclose(src_analyzed_grp_id);
        }

        //Average Quantification
        src_quant_grp_id = H5Gopen(src_maps_grp_id, "Quantification", H5P_DEFAULT);
        if(src_quant_grp_id > -1)
        {
            dst_quant_grp_id = H5Gopen(dst_maps_grp_id, "Quantification", H5P_DEFAULT);
            if(dst_quant_grp_id < 0)
                dst_quant_grp_id = H5Gcreate(dst_maps_grp_id, "Quantification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            hid_t dst_calib_fit_grp_id = H5Gopen(dst_quant_grp_id, "Calibration", H5P_DEFAULT);
            if(dst_calib_fit_grp_id < 0)
                dst_calib_fit_grp_id = H5Gcreate(dst_quant_grp_id, "Calibration", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            hid_t cablib_grp_id = H5Gopen(src_quant_grp_id, "Calibration", H5P_DEFAULT);
            if (cablib_grp_id > -1)
            {
                for (std::string analysis_grp_name : xrf_analysis_save_names)
                {
                    hid_t src_fit_grp_id = H5Gopen(cablib_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT);
                    if (src_fit_grp_id > -1)
                    {
                        hid_t dst_fit_grp_id = H5Gcreate(dst_calib_fit_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/Calibration_Curve_SR_Current", "Calibration_Curve_SR_Current", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/Calibration_Curve_DS_IC", "Calibration_Curve_DS_IC", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/Calibration_Curve_US_IC", "Calibration_Curve_US_IC", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        H5Gclose(src_fit_grp_id);
                        H5Gclose(dst_fit_grp_id);
                    }
                }
                H5Gclose(cablib_grp_id);
            }
            H5Gclose(dst_calib_fit_grp_id);
            
            status = H5Ocopy(src_quant_grp_id, "Number_Of_Standards", dst_quant_grp_id, "Number_Of_Standards", ocpypl_id, H5P_DEFAULT);

            hid_t num_stand_id = H5Dopen(src_quant_grp_id, "Number_Of_Standards", H5P_DEFAULT);
            if (num_stand_id > -1)
            {
                int num_standards = 0;
                hid_t dspace = H5Dget_space(num_stand_id);
                if (H5Dread(num_stand_id, H5T_NATIVE_INT, dspace, dspace, H5P_DEFAULT, &num_standards) > -1)
                {
                    for (int i = 0; i < num_standards; i++)
                    {
                        string standard_group_name = "Standard" + std::to_string(i);
                        string whole_standard_name = "MAPS/Quantification/" + standard_group_name;
                        hid_t standard_grp_id = H5Gopen(src_quant_grp_id, standard_group_name.c_str(), H5P_DEFAULT);
                        hid_t dst_standard_id = H5Gcreate(dst_quant_grp_id, standard_group_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                        status = H5Ocopy(standard_grp_id, "Element_Weights", dst_standard_id, "Element_Weights", ocpypl_id, H5P_DEFAULT);
                        status = H5Ocopy(standard_grp_id, "Element_Weights_Names", dst_standard_id, "Element_Weights_Names", ocpypl_id, H5P_DEFAULT);
                        status = H5Ocopy(standard_grp_id, "Scalers", dst_standard_id, "Scalers", ocpypl_id, H5P_DEFAULT);
                        status = H5Ocopy(standard_grp_id, "Standard_Name", dst_standard_id, "Standard_Name", ocpypl_id, H5P_DEFAULT);
                        _generate_avg_integrated_spectra(standard_grp_id, dst_standard_id, whole_standard_name.c_str(), ocpypl_id, hdf5_file_ids);
                        _generate_avg_analysis(standard_grp_id, dst_standard_id, whole_standard_name, ocpypl_id, hdf5_file_ids);
                        H5Gclose(standard_grp_id);
                        H5Gclose(dst_standard_id);
                    }
                }
                H5Sclose(dspace);
                H5Dclose(num_stand_id);
            }
            H5Gclose(dst_quant_grp_id);
            H5Gclose(src_quant_grp_id);
        }

        status = H5Ocopy(src_maps_grp_id, "version", dst_maps_grp_id, "version", ocpypl_id, H5P_DEFAULT);

        if(src_maps_grp_id > -1)
            H5Gclose(src_maps_grp_id);

        H5Gclose(dst_maps_grp_id);
        _cur_file_id = file_id;
        end_save_seq();
    }
    else
    {
        logE<<"Could not file any files to average for dataset "<<avg_filename<<"\n";
    }

    for(auto& f_id : hdf5_file_ids)
    {
        _cur_file_id = f_id;
        end_save_seq(false);
    }

    logI<<"closing file"<<"\n";

    _cur_file_id = saved_file_id;

    return true;
}

//-----------------------------------------------------------------------------

void HDF5_IO::_gen_average(std::string full_hdf5_path, std::string dataset_name, hid_t src_fit_grp_id, hid_t dst_fit_grp_id, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids, bool avg)
{
    std::vector<hid_t> analysis_ids;
	hid_t error;
	//hid_t status;

//    status = H5Ocopy(src_fit_grp_id, dataset_name.c_str(), dst_fit_grp_id, dataset_name.c_str(), ocpypl_id, H5P_DEFAULT);

//    if (status == 0)
    {
        hid_t dset_id = H5Dopen2(src_fit_grp_id, dataset_name.c_str(), H5P_DEFAULT);
        hid_t dataspace_id = H5Dget_space(dset_id);
        hid_t file_type = H5Dget_type(dset_id);
        hid_t props = H5Dget_create_plist(dset_id);
        hid_t dst_dset_id = H5Dcreate(dst_fit_grp_id, dataset_name.c_str(), file_type, dataspace_id, H5P_DEFAULT, props, H5P_DEFAULT);
        if(dst_dset_id < 1)
        {
            logE<<""<<full_hdf5_path<< " " <<dataset_name<<"\n";
            return;
        }
        int rank = H5Sget_simple_extent_ndims(dataspace_id);

        hsize_t* tmp_dims = new hsize_t[rank];
        hsize_t* dims_in = new hsize_t[rank];
        int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
        if(status_n < 0)
        {
            logE<<"could not get dataset dimensions for "<<full_hdf5_path<< " " <<dataset_name<<"\n";
            return;
        }

        //get all the other files dataset ids
        for(long unsigned int k=1; k<hdf5_file_ids.size(); k++)
        {
            hid_t det_analysis_dset_id = H5Dopen2(hdf5_file_ids[k], full_hdf5_path.c_str(), H5P_DEFAULT);
            if( det_analysis_dset_id > -1 )
            {
                // ran into a bug where detectors had different dims for counts per sec. need to check the min and use that to generate avg
                hid_t tdataspace_id = H5Dget_space(det_analysis_dset_id);
                int status_n = H5Sget_simple_extent_dims(tdataspace_id, &tmp_dims[0], NULL);
                if(status_n > -1)
                {
                    for(int i=0; i< rank; i++)
                    {
                        if(tmp_dims[i] < dims_in[i] )
                        {
                            dims_in[i] = tmp_dims[i];
                        }
                    }
                    analysis_ids.push_back(det_analysis_dset_id);
                }
            }
        }

        long long total = 1;
        hsize_t* offset_rank = new hsize_t[rank];
        for(int i=0; i< rank; i++)
        {
            offset_rank[i] = 0;
            total *= dims_in[i];
        }

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_rank, nullptr, dims_in, nullptr);

        delete []offset_rank;

		//long long avail_mem = get_total_mem() * .95;
		long long avail_mem = get_available_mem();
		if ((total * sizeof(real_t) * 2) > (avail_mem) && rank == 3) //mca_arr
		{
			//read in partial dataset at a time by chunks.
			hsize_t* offset = new hsize_t[rank];
			hsize_t* chunk_dims = new hsize_t[rank];
			unsigned long chunk_total = 1;
			error = H5Pget_chunk(dataspace_id, rank, chunk_dims);
			if (error < 0)
			{
				for (int i = 0; i < rank; i++)
				{
					offset[i] = 0;
				}
				chunk_total *= dims_in[0];
				chunk_dims[0] = dims_in[0];
				chunk_dims[1] = 1;
				chunk_dims[2] = 1;
			}
			else
			{
				for (int i = 0; i < rank; i++)
				{
					offset[i] = 0;
					chunk_total *= chunk_dims[i];
				}
			}
            data_struct::ArrayXr buffer1(chunk_total); //don't need to zero because we are reading in full buffer
            data_struct::ArrayXr buffer2(chunk_total); //don't need to zero because we are reading in full buffer
			
			hid_t memoryspace_id = H5Screate_simple(rank, chunk_dims, nullptr);
			for (int w = 0; w < dims_in[1]; w++)
			{
				for (int h = 0; h < dims_in[2]; h++)
				{
					offset[1] = w;
					offset[2] = h;
					float divisor = 1.0;
					H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
					error = H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
					if (error > -1)
					{
                        buffer1 = buffer1.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });
					}
					for (long unsigned int k = 0; k < analysis_ids.size(); k++)
					{
						//hid_t dspace_id = H5Dget_space(analysis_ids[k]);
						//H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
						//error = H5Dread(analysis_ids[k], H5T_NATIVE_REAL, memoryspace_id, dspace_id, H5P_DEFAULT, buffer2.data());
						error = H5Dread(analysis_ids[k], H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer2.data());
						if (error > -1)
						{
							buffer2 = buffer2.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });
							buffer1 += buffer2;
							divisor += 1.0;
						}
						else
						{
							logE << "reading " << full_hdf5_path << " dataset " << "\n";
						}
					}

					if (avg)
					{
						buffer1 /= divisor;
					}

					error = H5Dwrite(dst_dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
				}
			}
			delete[] offset;
			delete[] chunk_dims;
			H5Sclose(memoryspace_id);
		}
		else
		{
			//read in the whole dataset
			data_struct::ArrayXr buffer1(total);
			data_struct::ArrayXr buffer2(total);
			float divisor = 1.0;
            error = H5Dread(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
			if (error > -1)
			{
				buffer1 = buffer1.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });
			}
			for (long unsigned int k = 0; k < analysis_ids.size(); k++)
			{
                error = H5Dread(analysis_ids[k], H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, buffer2.data());
				if (error > -1)
				{
					buffer2 = buffer2.unaryExpr([](real_t v) { return std::isfinite(v) ? v : (real_t)0.0; });
					buffer1 += buffer2;
					divisor += 1.0;
				}
				else
				{
					logE << "reading " << full_hdf5_path << " dataset " << "\n";
				}
			}

			if (avg)
			{
				buffer1 /= divisor;
			}

			//hid_t dst_dset_id = H5Dopen2(dst_fit_grp_id, dataset_name.c_str(), H5P_DEFAULT);
            error = H5Dwrite(dst_dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
		}
        for(long unsigned int k=0; k<analysis_ids.size(); k++)
        {
            H5Dclose(analysis_ids[k]);
        }

        //clean up
        delete [] dims_in;
        delete [] tmp_dims;
        H5Dclose(dset_id);
        H5Dclose(dst_dset_id);
    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_generate_avg_analysis(hid_t src_maps_grp_id, hid_t dst_maps_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids)
{
    hid_t src_analyzed_grp_id = H5Gopen(src_maps_grp_id, "XRF_Analyzed", H5P_DEFAULT);
    if (src_analyzed_grp_id > -1)
    {
        hid_t dst_analyzed_grp_id = H5Gcreate(dst_maps_grp_id, "XRF_Analyzed", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for (std::string analysis_grp_name : xrf_analysis_save_names)
        {
            hid_t src_fit_grp_id = H5Gopen(src_analyzed_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT);
            if (src_fit_grp_id > -1)
            {
                //std::vector<hid_t> analysis_ids;

                hid_t dst_fit_grp_id = H5Gcreate(dst_analyzed_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                //copy channel names
                hid_t status = H5Ocopy(src_fit_grp_id, "Channel_Names", dst_fit_grp_id, "Channel_Names", ocpypl_id, H5P_DEFAULT);
                if (status < 0)
                {
                    logE << "coping Channel_Names!\n";
                }
                status = H5Ocopy(src_fit_grp_id, "Channel_Units", dst_fit_grp_id, "Channel_Units", ocpypl_id, H5P_DEFAULT);
                if (status < 0)
                {
                    logE << "coping Channel_Units!\n";
                }
                status = H5Ocopy(src_fit_grp_id, "Calibration_Curve_Labels", dst_fit_grp_id, "Calibration_Curve_Labels", ocpypl_id, H5P_DEFAULT);

                _gen_average(group_name + "/XRF_Analyzed/" + analysis_grp_name + "/Counts_Per_Sec", "Counts_Per_Sec", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);

                if (analysis_grp_name == STR_FIT_NNLS || analysis_grp_name == STR_FIT_GAUSS_MATRIX)
                {
                    _gen_average(group_name + "/XRF_Analyzed/" + analysis_grp_name + "/"+ STR_FIT_INT_SPEC, STR_FIT_INT_SPEC, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids, false);
                    _gen_average(group_name + "/XRF_Analyzed/" + analysis_grp_name + "/"+ STR_FIT_INT_BACKGROUND, STR_FIT_INT_BACKGROUND, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids, false);
                }
                H5Gclose(dst_fit_grp_id);
                H5Gclose(src_fit_grp_id);

            }
        }
        H5Gclose(dst_analyzed_grp_id);
+        H5Gclose(src_analyzed_grp_id);
    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_generate_avg_integrated_spectra(hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids)
{
    hid_t src_inner_grp_id = H5Gopen(src_analyzed_grp_id, STR_INT_SPEC.c_str(), H5P_DEFAULT);
    if(src_inner_grp_id > -1)
    {
        hid_t dst_inner_grp_id = H5Gcreate(dst_fit_grp_id, STR_INT_SPEC.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        _gen_average(group_name+"/Integrated_Spectra/Elapsed_Livetime", "Elapsed_Livetime", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Elapsed_Realtime", "Elapsed_Realtime", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Input_Counts", "Input_Counts", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Output_Counts", "Output_Counts", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
		_gen_average(group_name + "/Integrated_Spectra/"+ STR_MAX_CHANNELS_INT_SPEC, STR_MAX_CHANNELS_INT_SPEC, src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
		_gen_average(group_name + "/Integrated_Spectra/"+ STR_MAX10_INT_SPEC, STR_MAX10_INT_SPEC, src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);

        //don't average the integrated spectra, just sum it
        _gen_average(group_name+"/Integrated_Spectra/Spectra", "Spectra", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids, false);

        H5Gclose(dst_inner_grp_id);
        H5Gclose(src_inner_grp_id);
    }
}

//-----------------------------------------------------------------------------

bool HDF5_IO::generate_stream_dataset(std::string dataset_directory,
                                      std::string dataset_name,
                                      int detector_num,
                                      size_t height,
                                      size_t width)
{

    std::string str_detector_num = std::to_string(detector_num);
    std::string full_save_path = dataset_directory+ DIR_END_CHAR+"img.dat"+ DIR_END_CHAR +dataset_name+".h5"+str_detector_num;
    //io::file::HDF5_IO::inst()->set_filename(full_save_path);


    return false;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_stream_row(size_t d_hash,
                             size_t detector_num,
                             size_t row,
                             std::vector< data_struct::Spectra* >  *spectra_row)
{
    return false;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_itegrade_spectra(data_struct::Spectra * spectra)
{
    return false;
}

//-----------------------------------------------------------------------------


bool HDF5_IO::close_dataset(size_t d_hash)
{
    return false;
}

//-----------------------------------------------------------------------------

void HDF5_IO::update_theta(std::string dataset_file, std::string theta_pv_str)
{
    std::lock_guard<std::mutex> lock(_mutex);
	hid_t file_id, theta_id, extra_names, extra_values;
	std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
	char tmp_char[256] = { 0 };
	hsize_t dims_in[1] = { 0 };
	hsize_t offset_1d[1] = { 0 };
	hsize_t count_1d[1] = { 1 };
	hid_t rerror = 0;
	real_t theta_value = 0;

	file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
		return;
	close_map.push({ file_id, H5O_FILE });
	//if (false == _open_h5_object(file_id, H5O_FILE, close_map, dataset_file, -1))
	//	return;

	if (false == _open_h5_object(theta_id, H5O_DATASET, close_map, "/MAPS/Scan/theta", file_id))
		return;

	if (false == _open_h5_object(extra_names, H5O_DATASET, close_map, "/MAPS/Scan/Extra_PVs/Names", file_id))
		return;

	if (false == _open_h5_object(extra_values, H5O_DATASET, close_map, "/MAPS/Scan/Extra_PVs/Values", file_id))
		return;
	

	hid_t theta_space = H5Dget_space(theta_id);
	//hid_t theta_type = H5Dget_type(theta_id);
	hid_t name_space = H5Dget_space(extra_names);
	hid_t name_type = H5Dget_type(extra_names);
	H5Sget_simple_extent_dims(name_space, &dims_in[0], nullptr);
	hid_t memoryspace_id = H5Screate_simple(1, count_1d, nullptr);
	close_map.push({ memoryspace_id, H5O_DATASPACE });

	for (hsize_t i = 0; i < dims_in[0]; i++)
	{
		for (int z = 0; z < 256; z++)
			tmp_char[z] = 0;

		offset_1d[0] = i;
		H5Sselect_hyperslab(name_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
		rerror = H5Dread(extra_names, name_type, memoryspace_id, name_space, H5P_DEFAULT, (void*)tmp_char);

		std::string value(tmp_char, 255);
		value.erase(std::remove(value.begin(), value.end(), ' '), value.end());
		if (theta_pv_str == value)
		{
			for (int z = 0; z < 256; z++)
				tmp_char[z] = 0;
			rerror = H5Dread(extra_values, name_type, memoryspace_id, name_space, H5P_DEFAULT, (void*)tmp_char);
			theta_value = atof(tmp_char);
			rerror = H5Dwrite(theta_id, H5T_NATIVE_REAL, memoryspace_id, theta_space, H5P_DEFAULT, (void*)&theta_value);
			break;
		}

	}

	_close_h5_objects(close_map);

}

//-----------------------------------------------------------------------------

void HDF5_IO::update_scalers(std::string dataset_file, data_struct::Params_Override* params_override)
{
    std::lock_guard<std::mutex> lock(_mutex);
    if (params_override == nullptr)
    {
        return;
    }

    hid_t maps_grp_id;
    real_t time_scaler_clock = 1;
    try
    {
        hid_t scaler_dset_id = -1;
        hid_t scaler_names_id = -1;
        int time_idx = -1;
        std::map<std::string, int> scaler_name_idx_map;
        data_struct::ArrayXXr time_map;
        data_struct::ArrayXXr val_map;
        hid_t status;
        hid_t file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (file_id > -1)
        {
            scaler_dset_id = H5Dopen(file_id, "/MAPS/Scalers/Values", H5P_DEFAULT);
            scaler_names_id = H5Dopen(file_id, "/MAPS/Scalers/Names", H5P_DEFAULT);
        }
        if (scaler_dset_id > 0 && scaler_names_id > 0)
        {
            hsize_t offset_1d[1] = { 0 };
            hsize_t count_1d[1] = { 1 };
            hsize_t offset_2d[2] = { 0,0 };
            hsize_t count_2d[2] = { 1,1 };
            hsize_t offset_3d[3] = { 0,0,0 };
            hsize_t count_3d[3] = { 1,1,1 };

            hsize_t scaler_dims[3] = { 1,1,1 };

            hid_t scalername_space = H5Dget_space(scaler_names_id);
            hid_t scalername_type = H5Dget_type(scaler_names_id);
            hid_t scaler_space = H5Dget_space(scaler_dset_id);
            hid_t scaler_type = H5Dget_type(scaler_dset_id);

            H5Sget_simple_extent_dims(scaler_space, &scaler_dims[0], nullptr);

            time_map.resize(scaler_dims[1], scaler_dims[2]);
            val_map.resize(scaler_dims[1], scaler_dims[2]);
            count_2d[0] = count_3d[1] = scaler_dims[1];
            count_2d[1] = count_3d[2] = scaler_dims[2];

            hid_t space_2d = H5Screate_simple(2, &count_2d[0], nullptr);

            std::string scaler_name_str;
            char char_data[256] = { 0 };

            // scaler name , hdf offset
            hid_t single_space = H5Screate_simple(1, &count_1d[0], nullptr);
            hid_t status;
            // read in all scaler names and generate offset map
            for (hsize_t i = 0; i < scaler_dims[0]; i++)
            {
                offset_1d[0] = i;

                H5Sselect_hyperslab(scalername_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
                status = H5Dread(scaler_names_id, scalername_type, single_space, scalername_space, H5P_DEFAULT, (void*)&char_data[0]);
                if (status > -1)
                {
                    scaler_name_str = std::string(char_data, 255);
                    scaler_name_str.erase(std::remove(scaler_name_str.begin(), scaler_name_str.end(), ' '), scaler_name_str.end());
                    scaler_name_idx_map[scaler_name_str] = i;
                }
            }


            // search for time scaler pv
            if (params_override->time_scaler_clock.length() > 0)
            {
                time_scaler_clock = std::stod(params_override->time_scaler_clock);
            }

            for (const auto& itr : scaler_name_idx_map)
            {
                if (itr.first == params_override->time_scaler)
                {
                    offset_3d[0] = itr.second;
                    H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                    status = H5Dread(scaler_dset_id, scaler_type, space_2d, scaler_space, H5P_DEFAULT, (void*)time_map.data());
                    if (status > -1)
                    {
                        time_map /= time_scaler_clock;
                        break;
                    }
                }
            }


            for (const auto& ts_itr : params_override->time_normalized_scalers)
            {
                for (auto& scaler_name_itr : scaler_name_idx_map)
                {
                    if (ts_itr.second == scaler_name_itr.first)
                    {
                        memset(char_data, 0, 256);
                        ts_itr.first.copy(char_data, 255);
                        offset_1d[0] = scaler_name_itr.second;
                        H5Sselect_hyperslab(scalername_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
                        status = H5Dwrite(scaler_names_id, scalername_type, single_space, scalername_space, H5P_DEFAULT, (void*)&char_data[0]);

                        if (status > -1)
                        {
                            offset_3d[0] = scaler_name_itr.second;
                            H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                            status = H5Dread(scaler_dset_id, scaler_type, space_2d, scaler_space, H5P_DEFAULT, (void*)val_map.data());
                            if (status > -1)
                            {
                                val_map /= time_map;
                                status = H5Dwrite(scaler_dset_id, scaler_type, space_2d, scaler_space, H5P_DEFAULT, (void*)val_map.data());
                                if (status < 0)
                                {
                                    logW << "Failed to write normalize values for scaler " << char_data << "!\n";
                                }
                            }
                        }
                        break;
                    }
                }
            }

            // similar as above but don't have to normalize by time, just update name
            for (const auto& s_itr : params_override->scaler_pvs)
            {
                for (auto& scaler_name_itr : scaler_name_idx_map)
                {
                    if (s_itr.second == scaler_name_itr.first)
                    {
                        memset(char_data, 0, 256);
                        s_itr.first.copy(char_data, 255);
                        offset_1d[0] = scaler_name_itr.second;
                        H5Sselect_hyperslab(scalername_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
                        status = H5Dwrite(scaler_names_id, scalername_type, single_space, scalername_space, H5P_DEFAULT, (void*)&char_data[0]);
                        break;
                    }
                }
            }

            //  summed_scalers
            for (const auto& summed_scaler_itr : params_override->summed_scalers)
            {
                for (const auto& scaler : scaler_name_idx_map)
                {
                    // look for scaler to update
                    if (summed_scaler_itr.scaler_name == scaler.first)
                    {
                        //now we have to search through pv's and add them up
                        data_struct::ArrayXXr summed_map;
                        summed_map.resize(val_map.rows(), val_map.cols());
                        for (const auto& scaler_names_itr : summed_scaler_itr.scalers_to_sum)
                        {
                            for (const auto& scaler2 : scaler_name_idx_map)
                            {
                                if (scaler2.first == scaler_names_itr)
                                {
                                    //read and add
                                    offset_3d[0] = scaler2.second;
                                    H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                                    status = H5Dread(scaler_dset_id, scaler_type, space_2d, scaler_space, H5P_DEFAULT, (void*)val_map.data());
                                    if (status > -1)
                                    {
                                        summed_map += val_map;
                                    }
                                    break;
                                }
                            }
                        }
                        // pre normalized from above
                        //if (summed_scaler_itr.normalize_by_time == true)
                        //{
                        //    summed_map /= time_map;
                        //}

                        //write back
                        status = H5Dwrite(scaler_dset_id, scaler_type, space_2d, scaler_space, H5P_DEFAULT, (void*)val_map.data());
                        if (status < 0)
                        {
                            logW << "Failed to write normalize values for summed scaler "<<char_data<<"!\n";
                        }
                    }
                }
            }
        }
    }
    catch (...)
    {

    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_add_v9_quant(hid_t file_id, 
							hid_t quant_space,
							hid_t chan_names,
							hid_t chan_space,
							int chan_amt,
							std::string quant_str,
							std::string new_loc)
{
    hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
    hid_t memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, 255);
	char tmp_char1[256] = { 0 };

    //create quantification dataset. In v9 the array starts at element Z 10 insead of element Z 1
    std::string currnt_quant_str = "/MAPS/Quantification/Calibration/" + quant_str + "/Calibration_Curve_SR_Current";
    std::string us_quant_str = "/MAPS/Quantification/Calibration/" + quant_str + "/Calibration_Curve_US_IC";
    std::string ds_quant_str = "/MAPS/Quantification/Calibration/" + quant_str + "/Calibration_Curve_DS_IC";
	hid_t cc_current = H5Dopen(file_id, currnt_quant_str.c_str(), H5P_DEFAULT);
	hid_t cc_us_ic = H5Dopen(file_id, us_quant_str.c_str(), H5P_DEFAULT);
	hid_t cc_ds_ic = H5Dopen(file_id, ds_quant_str.c_str(), H5P_DEFAULT);
    hid_t quant_dset = H5Dcreate1(file_id, new_loc.c_str(), H5T_NATIVE_REAL, quant_space, H5P_DEFAULT);
	if (quant_dset < 0)
	{
		quant_dset = H5Dopen(file_id, new_loc.c_str(), H5P_DEFAULT);
	}
    if(quant_dset > -1 && cc_current > -1 && cc_us_ic > -1 && cc_ds_ic > -1)
    {
		hid_t chan_units = H5Dopen(file_id, "/MAPS/channel_units", H5P_DEFAULT);
        //share the space between sr_current, us_ic, and ds_ic
        hid_t cc_space = H5Dget_space(cc_current);
        hid_t us_space = H5Dget_space(cc_us_ic);
        hid_t ds_space = H5Dget_space(cc_ds_ic);
        
        hsize_t count_1d[1] = {1};
        hsize_t count_2d[2] = {1,1};
        hsize_t count_3d[3] = {1,1,1};
        hsize_t offset_1d[1] = {0};
        hsize_t offset_2d[2] = {0,0};
        hsize_t offset_3d[3] = {0,0,0};
        hid_t memoryspace_id = H5Screate_simple(1, count_1d, nullptr);
        char sr_current_carr[255] = "SRCURRENT";
        char us_ic_carr[255];
        char ds_ic_carr[255];
        real_t real_val = 0.0;
        hid_t err;

        STR_US_IC.copy(us_ic_carr, 254);
        STR_DS_IC.copy(ds_ic_carr, 254);

        count_1d[0] = 3;
        hid_t quant_name_space = H5Screate_simple(1, count_1d, nullptr);
        count_1d[0] = 1;
        //save quant_names to know what each index is
        new_loc += "_names";
        hid_t quant_names_dset = H5Dcreate1(file_id, new_loc.c_str(), filetype, quant_name_space, H5P_DEFAULT);
        if(quant_names_dset > -1)
        {
            offset_1d[0] = 0;
            H5Sselect_hyperslab(quant_name_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
            H5Dwrite(quant_names_dset, filetype, memoryspace_id, quant_name_space, H5P_DEFAULT, (void*)sr_current_carr);

            offset_1d[0] = 1;
            H5Sselect_hyperslab(quant_name_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
            H5Dwrite(quant_names_dset, filetype, memoryspace_id, quant_name_space, H5P_DEFAULT, (void*)us_ic_carr);

            offset_1d[0] = 2;
            H5Sselect_hyperslab(quant_name_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
            H5Dwrite(quant_names_dset, filetype, memoryspace_id, quant_name_space, H5P_DEFAULT, (void*)ds_ic_carr);


            H5Sclose(quant_name_space);
            H5Dclose(quant_names_dset);
        }
        //reset back
        offset_1d[0] = 0;


        for(int chan_idx=0; chan_idx<chan_amt; chan_idx++)
        {
			offset_1d[0] = chan_idx;
			
            H5Sselect_hyperslab(chan_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
            char tmp_char[256] = {0};
            if (H5Dread(chan_names, memtype, memoryspace_id, chan_space, H5P_DEFAULT, (void*)tmp_char) > -1)
            {
                std::string el_name_str = std::string(tmp_char);
                unsigned long underscore_idx = el_name_str.find("_");
                //can check if > 0 instead of -1 since it shouldn't start with an '_'
                if (underscore_idx > 0 && underscore_idx < el_name_str.length())
                {
                    if(el_name_str[underscore_idx+1] == 'L')
                        offset_2d[0] = 1;
                    if(el_name_str[underscore_idx+1] == 'M')
                        offset_2d[0] = 2;
                }
                else
                {
                    offset_2d[0] = 0;
                }
                auto element = data_struct::Element_Info_Map::inst()->get_element(el_name_str);
                if (element != nullptr)
                {
                    offset_2d[1] = element->number - 1;
                    H5Sselect_hyperslab(cc_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);
                    H5Sselect_hyperslab(us_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);
                    H5Sselect_hyperslab(ds_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);

                    offset_3d[2] = chan_idx;
                    offset_3d[0] = 0;
                    real_val = 0.0;
                    H5Sselect_hyperslab(quant_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                    if (H5Dread(cc_current, H5T_NATIVE_REAL, memoryspace_id, cc_space, H5P_DEFAULT, (void*)&real_val) > -1)
                    {
                        err = H5Dwrite(quant_dset, H5T_NATIVE_REAL, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);
                    }
                    offset_3d[0] = 1;
                    real_val = 0.0;
                    H5Sselect_hyperslab(quant_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                    if (H5Dread(cc_us_ic, H5T_NATIVE_REAL, memoryspace_id, us_space, H5P_DEFAULT, (void*)&real_val) > -1)
                    {
                        err = H5Dwrite(quant_dset, H5T_NATIVE_REAL, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);
                    }
                    offset_3d[0] = 2;
                    real_val = 0.0;
                    H5Sselect_hyperslab(quant_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                    if (H5Dread(cc_ds_ic, H5T_NATIVE_REAL, memoryspace_id, ds_space, H5P_DEFAULT, (void*)&real_val) > -1)
                    {
                        err = H5Dwrite(quant_dset, H5T_NATIVE_REAL, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);
                    }
                    //change /MAPS/channel_units from cts/s to ug/cm2 for the first 3 of 4 in dim[0]
                    if (chan_units > -1)
                    {
                        hid_t unit_type = H5Dget_type(chan_units);
                        hid_t unit_space = H5Dget_space(chan_units);
                        std::string update = "ug/cm2";
                        for (int z = 0; z < 256; z++)
                        {
                            tmp_char1[z] = 0;
                        }
                        update.copy(tmp_char1, 256);
                        for (int i = 0; i < 3; i++)
                        {
                            for (int j = 0; j < chan_amt; j++)
                            {
                                offset_2d[0] = i;
                                offset_2d[1] = j;
                                H5Sselect_hyperslab(unit_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);
                                H5Dwrite(chan_units, unit_type, memoryspace_id, unit_space, H5P_DEFAULT, (void*)tmp_char1);
                            }
                        }
                    }
                }
            }
        }
        if (chan_units > -1)
        {
            H5Dclose(chan_units);
        }
    }
    if (quant_dset > -1)
    {
        H5Dclose(quant_dset);
    }
    if (cc_current > -1)
    {
        H5Dclose(cc_current);
    }
    if (cc_us_ic > -1)
    {
        H5Dclose(cc_us_ic);
    }
    if (cc_ds_ic > -1)
    {
        H5Dclose(cc_ds_ic);
    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_add_extra_pvs(hid_t file_id, std::string group_name)
{
	char tmp_char[256] = { 0 };
    //open scan extras and create 4xN array
    hid_t extra_names = H5Dopen(file_id, "/MAPS/Scan/Extra_PVs/Names", H5P_DEFAULT);
    hid_t extra_units = H5Dopen(file_id, "/MAPS/Scan/Extra_PVs/Unit", H5P_DEFAULT);
    hid_t extra_values = H5Dopen(file_id, "/MAPS/Scan/Extra_PVs/Values", H5P_DEFAULT);
    hid_t extra_desc = H5Dopen(file_id, "/MAPS/Scan/Extra_PVs/Description", H5P_DEFAULT);
    if(extra_names > -1 && extra_units > -1 && extra_values > -1 && extra_desc > -1)
    {
        hid_t name_space = H5Dget_space(extra_names);
        int rank = H5Sget_simple_extent_ndims(name_space);
        hsize_t* dims_in = new hsize_t[rank];
        H5Sget_simple_extent_dims(name_space, &dims_in[0], nullptr);
        hsize_t extra_pv_dims[2];
        extra_pv_dims[0] = 4;
        extra_pv_dims[1] = dims_in[0];

        hid_t file_space = H5Screate_simple(2, &extra_pv_dims[0], &extra_pv_dims[0]);

        hid_t name_type = H5Dget_type(extra_names);
        std::string extra_pvs_str = group_name + "/extra_pvs";
        std::string extra_pvs_as_csv_str = group_name + "/extra_pvs_as_csv";
        hid_t extra_pvs = H5Dcreate1(file_id, extra_pvs_str.c_str(), name_type, file_space, H5P_DEFAULT);
        hid_t extra_pvs_as_csv = H5Dcreate1(file_id, extra_pvs_as_csv_str.c_str(), name_type, name_space, H5P_DEFAULT);

        hsize_t offset_1d[1] = {0};
        hsize_t count_1d[1] = {1};
        hsize_t offset_2d[2] = {0,0};
        hsize_t count_2d[2] = {1,1};

        hid_t memoryspace_id = H5Screate_simple(1, count_1d, nullptr);
        for(hsize_t i=0; i<dims_in[0]; i++)
        {
			for (int z = 0; z < 256; z++)
				tmp_char[z] = 0;
            std::string as_csv = "";
            offset_1d[0] = i;
            offset_2d[1] = i;
            //name
            offset_2d[0] = 0;
            H5Sselect_hyperslab(name_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);
            H5Dread(extra_names, name_type, memoryspace_id, name_space, H5P_DEFAULT, (void*)tmp_char);
            H5Dwrite(extra_pvs, name_type, memoryspace_id, file_space, H5P_DEFAULT, (void*)tmp_char);
            std::string val(tmp_char);
            val = val.substr(0,val.find("\x20"));
            as_csv += val;
            as_csv += ",";
            //value
            offset_2d[0] = 1;
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);
            H5Dread(extra_values, name_type, memoryspace_id, name_space, H5P_DEFAULT, (void*)tmp_char);
            H5Dwrite(extra_pvs, name_type, memoryspace_id, file_space, H5P_DEFAULT, (void*)tmp_char);
            std::string val2(tmp_char);
            val2 = val2.substr(0,val2.find("\x20"));
            as_csv += val2;
            //description
            offset_2d[0] = 2;
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);
            H5Dread(extra_desc, name_type, memoryspace_id, name_space, H5P_DEFAULT, (void*)tmp_char);
            H5Dwrite(extra_pvs, name_type, memoryspace_id, file_space, H5P_DEFAULT, (void*)tmp_char);
            //units
            offset_2d[0] = 3;
            H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);
            H5Dread(extra_units, name_type, memoryspace_id, name_space, H5P_DEFAULT, (void*)tmp_char);
            H5Dwrite(extra_pvs, name_type, memoryspace_id, file_space, H5P_DEFAULT, (void*)tmp_char);

            as_csv.copy(tmp_char, 254);
            H5Dwrite(extra_pvs_as_csv, name_type, memoryspace_id, name_space, H5P_DEFAULT, (void*)tmp_char);
        }
        delete [] dims_in;

        if (extra_pvs > -1)
            H5Dclose(extra_pvs);
        if (extra_pvs_as_csv > -1)
            H5Dclose(extra_pvs_as_csv);
    }

    if (extra_names > -1)
        H5Dclose(extra_names);
    if (extra_units > -1)
        H5Dclose(extra_units);
    if (extra_values > -1)
        H5Dclose(extra_values);
    if (extra_desc > -1)
        H5Dclose(extra_desc);
}

//-----------------------------------------------------------------------------

void HDF5_IO::add_v9_layout(std::string dataset_file)
{
    std::lock_guard<std::mutex> lock(_mutex);

    logI  << dataset_file << "\n";
    hid_t saved_file_id = _cur_file_id;

    hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
    hid_t memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, 255);
    char tmp_char1[256] = { 0 };
	hsize_t offset2d[2] = { 0,0 };
	hsize_t count2d[2] = { 1,1 };

    hid_t file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
        _cur_file_id = saved_file_id;
		return;
	}
    //Scan
    if (H5Gget_objinfo(file_id, "/MAPS/x_axis", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Scan/x_axis", H5L_SAME_LOC, "/MAPS/x_axis", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for x_axis" << "\n";
        }
    }

    if (H5Gget_objinfo(file_id, "/MAPS/y_axis", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Scan/y_axis", H5L_SAME_LOC, "/MAPS/y_axis", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for y_axis" << "\n";
        }
    }

    if (H5Gget_objinfo(file_id, "/MAPS/scan_time_stamp", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Scan/scan_time_stamp", H5L_SAME_LOC, "/MAPS/scan_time_stamp", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for scan_time_stamp" << "\n";
        }
    }
    //create extra_pvs, extra_pvs_as_csv, extra_strings

    //Scalers
    if (H5Gget_objinfo(file_id, "/MAPS/ds_amp", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Scalers/ds_amp", H5L_SAME_LOC, "/MAPS/ds_amp", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for ds_amp" << "\n";
        }
    }

    if (H5Gget_objinfo(file_id, "/MAPS/us_amp", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Scalers/us_amp", H5L_SAME_LOC, "/MAPS/us_amp", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for us_amp" << "\n";
        }
    }

    //need to only add scalers that were in maps_fit_parameter_override since old GUI doesn't support having a lot of scalers
    _add_v9_scalers(file_id);
    /*
    if( H5Lcreate_hard(file_id, "/MAPS/Scalers/Names", H5L_SAME_LOC, "/MAPS/scaler_names", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for scaler_names"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Scalers/Values", H5L_SAME_LOC, "/MAPS/scalers", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for scalers"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Scalers/Units", H5L_SAME_LOC, "/MAPS/scaler_units", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for scaler_units"<<  "\n";
    }
	*/


    //Spectra
    if (H5Gget_objinfo(file_id, "/MAPS/energy", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Spectra/Energy", H5L_SAME_LOC, "/MAPS/energy", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for energy" << "\n";
        }
    }
    if (H5Gget_objinfo(file_id, "/MAPS/energy_calib", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Spectra/Energy_Calibration", H5L_SAME_LOC, "/MAPS/energy_calib", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for energy_calib" << "\n";
        }
    }
    if (H5Gget_objinfo(file_id, "/MAPS/int_spec", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Spectra/Integrated_Spectra/Spectra", H5L_SAME_LOC, "/MAPS/int_spec", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for int_spec" << "\n";
        }
    }
    if (H5Gget_objinfo(file_id, "/MAPS/mca_arr", 0, NULL) < 0)
    {
        if (H5Lcreate_hard(file_id, "/MAPS/Spectra/mca_arr", H5L_SAME_LOC, "/MAPS/mca_arr", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW << " Couldn't create soft link for mca_arr" << "\n";
        }
    }


	std::string fit_int_name = "/MAPS/XRF_Analyzed/Fitted/"+ STR_FIT_INT_SPEC;
	std::string max_name = "/MAPS/Spectra/Integrated_Spectra/"+ STR_MAX_CHANNELS_INT_SPEC;
	std::string max10_name = "/MAPS/Spectra/Integrated_Spectra/"+ STR_MAX10_INT_SPEC;
	std::string v9_max_name = "/MAPS/max_chan_spec";
	hid_t fit_int_id, max_id, max_10_id, max_space, max_type, v9_max_id, v9_space;
	
	max_id = H5Dopen(file_id, max_name.c_str(), H5P_DEFAULT);
	if (max_id > -1)
	{
		fit_int_id = H5Dopen(file_id, fit_int_name.c_str(), H5P_DEFAULT);
		max_10_id = H5Dopen(file_id, max10_name.c_str(), H5P_DEFAULT);

		max_space = H5Dget_space(max_id);
		max_type = H5Dget_type(max_id);
		H5Sget_simple_extent_dims(max_space, &count2d[0], nullptr);
		//make 5 x spectra_size matrix
		count2d[1] = count2d[0];
		count2d[0] = 5;
		v9_space = H5Screate_simple(2, &count2d[0], &count2d[0]);

		v9_max_id = H5Dopen(file_id, v9_max_name.c_str(), H5P_DEFAULT);
		if (v9_max_id < 0)
		{
			v9_max_id = H5Dcreate(file_id, v9_max_name.c_str(), max_type, v9_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		}
		if (v9_max_id > -1)
		{
			real_t *buf = new real_t[count2d[1]];
			count2d[0] = 1;
			
			if (H5Dread(max_id, max_type, max_space, max_space, H5P_DEFAULT, (void*)buf) > -1)
			{
				offset2d[0] = 0;
				H5Sselect_hyperslab(v9_space, H5S_SELECT_SET, offset2d, nullptr, count2d, nullptr);
				H5Dwrite(v9_max_id, max_type, max_space, v9_space, H5P_DEFAULT, (void*)buf);
			}
			if (max_10_id > -1)
			{
				if (H5Dread(max_10_id, max_type, max_space, max_space, H5P_DEFAULT, (void*)buf) > -1)
				{
					offset2d[0] = 1;
					H5Sselect_hyperslab(v9_space, H5S_SELECT_SET, offset2d, nullptr, count2d, nullptr);
					H5Dwrite(v9_max_id, max_type, max_space, v9_space, H5P_DEFAULT, (void*)buf);
				}
				H5Dclose(max_10_id);
			}
			if (fit_int_id > -1)
			{
				if (H5Dread(fit_int_id, max_type, max_space, max_space, H5P_DEFAULT, (void*)buf) > -1)
				{
					offset2d[0] = 2;
					H5Sselect_hyperslab(v9_space, H5S_SELECT_SET, offset2d, nullptr, count2d, nullptr);
					H5Dwrite(v9_max_id, max_type, max_space, v9_space, H5P_DEFAULT, (void*)buf);
				}
				H5Dclose(fit_int_id);
			}
			delete[]buf;
			H5Dclose(v9_max_id);
		}
		H5Sclose(v9_space);
		H5Dclose(max_id);
	}
	

    hsize_t quant_dims[3];
    quant_dims[0] = 3;
    quant_dims[1] = 1;
    quant_dims[2] = 1; //num channel names


    //XRF_Analyzed
    if (H5Gget_objinfo(file_id, "/MAPS/channel_names", 0, NULL) < 0)
    {
        if (H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/ROI/Channel_Names", 0, NULL) >= 0)
        {
            H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/ROI/Channel_Names", H5L_SAME_LOC, "/MAPS/channel_names", H5P_DEFAULT, H5P_DEFAULT);
        }
        else if (H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/Fitted/Channel_Names", 0, NULL) >= 0)
        {
            H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/Fitted/Channel_Names", H5L_SAME_LOC, "/MAPS/channel_names", H5P_DEFAULT, H5P_DEFAULT);
        }
        else if (H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/NNLS/Channel_Names", 0, NULL) >= 0)
        {
            H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/NNLS/Channel_Names", H5L_SAME_LOC, "/MAPS/channel_names", H5P_DEFAULT, H5P_DEFAULT);
        }
    }

    hid_t chan_names = H5Dopen(file_id, "/MAPS/channel_names", H5P_DEFAULT);
    if (chan_names < 0)
    {
        _cur_file_id = saved_file_id;
        return;
    }
        
    hid_t chan_space = H5Dget_space(chan_names);
    if(chan_names > -1)
    {
        hsize_t chan_size = 1;
        H5Sget_simple_extent_dims(chan_space, &chan_size, nullptr);
        quant_dims[2] = chan_size; //num channel names
    }
    hid_t quant_space = H5Screate_simple(3, &quant_dims[0], &quant_dims[0]);

    //Channel Units are a 4 x channels so we can't do a hardlink
    // the 4 are SR_current, US_IC, DS_IC, and cts/s
    if (H5Gget_objinfo(file_id, "/MAPS/channel_units", 0, NULL) < 0)
    {
        hsize_t unit_dims[2];
        hsize_t offset_dims[2] = { 0,0 };
        unit_dims[0] = 4;
        unit_dims[1] = quant_dims[2];
        hid_t units_space = H5Screate_simple(2, &unit_dims[0], &unit_dims[0]);
        hid_t ch_unit_id = H5Dcreate1(file_id, "/MAPS/channel_units", filetype, units_space, H5P_DEFAULT);
        if (ch_unit_id > 0)
        {
            hsize_t mem_dims[1] = { 1 };
            hsize_t count_2d[2] = { 1,1 };
            hid_t mem_space = H5Screate_simple(1, &mem_dims[0], &mem_dims[0]);

            std::string str_val = "cts/s";
            for (int z = 0; z < 256; z++)
            {
                tmp_char1[z] = 0;
            }
            str_val.copy(tmp_char1, 256);
            for (hsize_t i = 0; i < unit_dims[0]; i++)
            {
                for (hsize_t j = 0; j < unit_dims[1]; j++)
                {
                    offset_dims[0] = i;
                    offset_dims[1] = j;
                    H5Sselect_hyperslab(units_space, H5S_SELECT_SET, offset_dims, nullptr, count_2d, nullptr);
                    H5Dwrite(ch_unit_id, filetype, mem_space, units_space, H5P_DEFAULT, (void*)&str_val[0]);
                }
            }
            H5Dclose(ch_unit_id);
        }
        else
        {
            logW << "Couldn't create /MAPS/channel_units" << "\n";
        }
    }

    if (H5Gget_objinfo(file_id, "/MAPS/XRF_roi", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/ROI/Counts_Per_Sec", 0, NULL) >= 0)
    {
        H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/ROI/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_roi", H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Gget_objinfo(file_id, "/MAPS/XRF_roi_quant", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_roi", 0, NULL) >= 0)
    {
        _add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], "ROI", "/MAPS/XRF_roi_quant");
    }

    if (H5Gget_objinfo(file_id, "/MAPS/XRF_roi_plus", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec", 0, NULL) >= 0)
    {
        H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_roi_plus", H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Gget_objinfo(file_id, "/MAPS/XRF_roi_plus_quant", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_roi_plus", 0, NULL) >= 0)
    {
        _add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], "NNLS", "/MAPS/XRF_roi_plus_quant");
    }

    if (H5Gget_objinfo(file_id, "/MAPS/XRF_fits", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec", 0, NULL) >= 0)
    {
        H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_fits", H5P_DEFAULT, H5P_DEFAULT);
    }
    if (H5Gget_objinfo(file_id, "/MAPS/XRF_fits_quant", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_fits", 0, NULL) >= 0)
    {
        _add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], STR_FIT_GAUSS_MATRIX, "/MAPS/XRF_fits_quant");
    }

    _add_extra_pvs(file_id, "/MAPS");

    if (chan_names > -1)
    {
        H5Dclose(chan_names);
    }
    if (quant_space > -1)
    {
        H5Sclose(quant_space);
    }

    //change version to 9
    real_t version = 9;
    hid_t version_id = H5Dopen(file_id, "/MAPS/version", H5P_DEFAULT);
    hid_t ver_space = H5Dget_space(version_id);
    hid_t ver_type = H5Dget_type(version_id);
    H5Dwrite(version_id, ver_type, ver_space, ver_space, H5P_DEFAULT, (void*)&version);
    H5Dclose(version_id);
    if (H5Gget_objinfo(file_id, "/version", 0, NULL) < 0)
    {
        H5Lcreate_hard(file_id, "/MAPS/version", H5L_SAME_LOC, "/version", H5P_DEFAULT, H5P_DEFAULT);
    }

    _cur_file_id = file_id;
    end_save_seq();
    _cur_file_id = saved_file_id;

}

//-----------------------------------------------------------------------------

void HDF5_IO::_add_v9_scalers(hid_t file_id)
{

    std::map<std::string, int> scaler_map;
    hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
    hid_t memtype = H5Tcopy(H5T_C_S1);

    hid_t names_id = H5Dopen(file_id, "/MAPS/Scalers/Names", H5P_DEFAULT);
    hid_t values_id = H5Dopen(file_id, "/MAPS/Scalers/Values", H5P_DEFAULT);
    hid_t units_id = H5Dopen(file_id, "/MAPS/Scalers/Units", H5P_DEFAULT);
    if (names_id < 0 || values_id < 0 || units_id < 0)
    {
        logW << "Could not open /MAPS/Scalers/Names Values, or Units dataset to add v9 layout\n";
        return;
    }

    hid_t name_space = H5Dget_space(names_id);
    hid_t value_space = H5Dget_space(values_id);
    hid_t unit_space = H5Dget_space(units_id);


    hsize_t name_size = 1;
    hsize_t value_size[3] = { 1,1,1 };
    H5Sget_simple_extent_dims(name_space, &name_size, nullptr);
    H5Sget_simple_extent_dims(value_space, &value_size[0], nullptr);

    hsize_t offset_1d[1] = { 0 };
    hsize_t count_1d[1] = { 1 };
    hsize_t offset_3d[3] = { 0,0,0 };
    hsize_t count_3d[3] = { 1,1,1 };

    hsize_t new_offset_1d[1] = { 0 };
    hsize_t new_offset_3d[3] = { 0,0,0 };

    hsize_t new_max_1d[1] = { name_size };
    hsize_t new_max_3d[3] = { name_size, value_size[1], value_size[2] };

    count_3d[1] = value_size[1];
    count_3d[2] = value_size[2];


    hid_t mem_space_1d = H5Screate_simple(1, &count_1d[0], &count_1d[0]);

    for (hsize_t i = 0; i < name_size; i++)
    {
        offset_1d[0] = i;
        H5Sselect_hyperslab(name_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
        char tmp_char[256] = { 0 };
        if (H5Dread(names_id, filetype, mem_space_1d, name_space, H5P_DEFAULT, (void*)tmp_char) > -1)
        {
            std::string scaler_name_str = std::string(tmp_char, 255);
            scaler_name_str.erase(std::remove(scaler_name_str.begin(), scaler_name_str.end(), ' '), scaler_name_str.end());
            int c_idx = scaler_name_str.find(':');

            if (c_idx < 0 && scaler_name_str.length() > 0)
            {
                scaler_map[scaler_name_str] = i;
            }
        }
    }

    count_1d[0] = { scaler_map.size() };
    hid_t new_name_space = H5Screate_simple(1, &count_1d[0], &new_max_1d[0]);

    count_3d[0] = { scaler_map.size() };
    hid_t new_value_space = H5Screate_simple(3, &count_3d[0], &new_max_3d[0]);

    hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 1, count_1d);
    H5Pset_deflate(dcpl_id, 7);

    hid_t dcpl_id2 = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id2, 3, count_3d);
    H5Pset_deflate(dcpl_id, 7);

    hid_t new_names_id = H5Dopen(file_id, "/MAPS/scaler_names", H5P_DEFAULT);
    if (new_names_id < 0)
    {
        new_names_id = H5Dcreate(file_id, "/MAPS/scaler_names", filetype, new_name_space, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    }
    else
    {
        H5Pclose(new_name_space);
        herr_t err = H5Dset_extent(new_names_id, &count_1d[0]);
        if (err < 0)
        {
            logW << "Error extending dataset /MAPS/scaler_names";
        }
        new_name_space = H5Dget_space(new_names_id);
    }

    hid_t new_units_id = H5Dopen(file_id, "/MAPS/scaler_units", H5P_DEFAULT);
    if (new_units_id < 0)
    {
        new_units_id = H5Dcreate(file_id, "/MAPS/scaler_units", filetype, new_name_space, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
    }
    else
    {
        herr_t err = H5Dset_extent(new_units_id, &count_1d[0]);
        if (err < 0)
        {
            logW << "Error extending dataset /MAPS/scaler_units";
        }
    }

    hid_t new_values_id = H5Dopen(file_id, "/MAPS/scalers", H5P_DEFAULT);
    if (new_values_id < 0)
    {
        new_values_id = H5Dcreate(file_id, "/MAPS/scalers", H5T_NATIVE_REAL, new_value_space, H5P_DEFAULT, dcpl_id2, H5P_DEFAULT);
    }
    else
    {
        H5Pclose(new_value_space);
        herr_t err = H5Dset_extent(new_values_id, &count_3d[0]);
        if (err < 0)
        {
            logW << "Error extending dataset /MAPS/scalers";
        }
        new_value_space = H5Dget_space(new_values_id);
    }


    if (new_names_id < 0 || new_units_id < 0 || new_values_id < 0)
    {
        logW << "Could not open /MAPS/scalers _names, or _units dataset to add v9 layout\n";
        return;
    }

    count_1d[0] = 1;
    count_3d[0] = 1;

    hid_t value_mem_space = H5Screate_simple(3, &count_3d[0], &count_3d[0]);

    data_struct::ArrayXXr tmp_values;
    tmp_values.resize(count_3d[1], count_3d[2]);

    hsize_t i = 0;
    for (const auto& itr : scaler_map)
    {
        offset_1d[0] = itr.second;
        offset_3d[0] = itr.second;

        new_offset_1d[0] = i;
        new_offset_3d[0] = i;

        H5Sselect_hyperslab(name_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
        H5Sselect_hyperslab(new_name_space, H5S_SELECT_SET, new_offset_1d, nullptr, count_1d, nullptr);
        

        char tmp_char_name[256] = { 0 };
        itr.first.copy(tmp_char_name, 254);
        H5Dwrite(new_names_id, filetype, mem_space_1d, new_name_space, H5P_DEFAULT, (void*)&tmp_char_name[0]);

        char tmp_char[256] = { 0 };
        if (H5Dread(units_id, filetype, mem_space_1d, name_space, H5P_DEFAULT, (void*)tmp_char) > -1)
        {
            H5Dwrite(new_units_id, filetype, mem_space_1d, new_name_space, H5P_DEFAULT, (void*)tmp_char);
        }

        H5Sselect_hyperslab(value_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
        if (H5Dread(values_id, H5T_NATIVE_REAL, value_mem_space, value_space, H5P_DEFAULT, (void*)tmp_values.data()) > -1)
        {
            H5Sselect_hyperslab(new_value_space, H5S_SELECT_SET, new_offset_3d, nullptr, count_3d, nullptr);
            H5Dwrite(new_values_id, H5T_NATIVE_REAL, value_mem_space, new_value_space, H5P_DEFAULT, (void*)tmp_values.data());
        }

        i++;
    }

    H5Pclose(dcpl_id);
    H5Pclose(dcpl_id2);

    H5Sclose(value_mem_space);
    H5Sclose(new_name_space);
    H5Sclose(new_value_space);
    H5Sclose(mem_space_1d);
    H5Sclose(unit_space);
    H5Sclose(name_space);
    H5Sclose(value_space);

    H5Dclose(names_id);
    H5Dclose(units_id);
    H5Dclose(values_id);
    H5Dclose(new_names_id);
    H5Dclose(new_units_id);
    H5Dclose(new_values_id);

}

//-----------------------------------------------------------------------------

bool HDF5_IO::_add_exchange_meta(hid_t file_id, std::string exchange_idx, std::string fits_link, std::string normalize_scaler)
{
    char desc[256] = {0};
    std::string exhange_str = "/exchange_" + exchange_idx;
    std::string xaxis_str = exhange_str + "/x_axis" ;
    std::string yaxis_str = exhange_str + "/y_axis" ;
    std::string str_scalers = "/MAPS/Scalers/Values";
    std::string str_scaler_names = "/MAPS/Scalers/Names";
    std::string str_scaler_units = "/MAPS/Scalers/Units";

    std::string str_scan_theta = "/MAPS/Scan/theta";

    std::string exchange_scalers = exhange_str + "/scalers";
    std::string exchange_scaler_names = exhange_str + "/scaler_names";
    std::string exchange_scaler_units = exhange_str + "/scaler_units";

    std::string str_desc = exhange_str + "/description";
    std::string str_version = exhange_str + "/version";
    std::string str_theta = exhange_str + "/theta";

    std::string exchange_images = exhange_str + "/data";
    std::string exchange_image_names = exhange_str + "/data_names";
    std::string exchange_image_units = exhange_str + "/data_units";



    std::string full_fit_link_path = "/MAPS/XRF_Analyzed/" + fits_link;
    hid_t fits_grp = H5Gopen(file_id, full_fit_link_path.c_str(), H5P_DEFAULT);
    if(fits_grp < 0)
    {
        //if we don't find the analysis, don't make an exchange for it
        return false;
    }



    hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
    hsize_t count [1] = {1};
    hid_t dataspace_id = H5Screate_simple (1, count, nullptr);

    hid_t exchange_id = H5Dopen(file_id, exhange_str.c_str(), H5P_DEFAULT);
    if(exchange_id < 0)
    {
        exchange_id = H5Gcreate(file_id, exhange_str.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }

    //Scan
    if( H5Lcreate_hard(file_id, "/MAPS/Scan/x_axis", H5L_SAME_LOC, xaxis_str.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for x_axis"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Scan/y_axis", H5L_SAME_LOC, yaxis_str.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for y_axis"<<  "\n";
    }

    _add_extra_pvs(file_id, exhange_str);

    if(normalize_scaler.size() > 0)
    {
        //Save description
        std::string norm_desc = fits_link + " normalized by " + normalize_scaler;
        hid_t dset_id = H5Dcreate(file_id, str_desc.c_str(), filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(dset_id > 0)
        {
            norm_desc.copy(desc, 256);
            H5Dwrite(dset_id, filetype, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&desc);
            H5Dclose(dset_id);
        }

        dset_id = H5Dopen(file_id, "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec", H5P_DEFAULT);
        //hid_t chan_units_id = H5Dopen(file_id, "/MAPS/XRF_Analyzed/Fitted/Channel_Units", H5P_DEFAULT);
        hid_t chan_names_id = H5Dopen(file_id, "/MAPS/XRF_Analyzed/Fitted/Channel_Names", H5P_DEFAULT);
        hid_t scaler_dset_id = H5Dopen(file_id, "/MAPS/Scalers/Values", H5P_DEFAULT);
        hid_t scaler_units_id = H5Dopen(file_id, "/MAPS/Scalers/Units", H5P_DEFAULT);
        hid_t scaler_names_id = H5Dopen(file_id, "/MAPS/Scalers/Names", H5P_DEFAULT);
        hid_t ds_ic_quant_id = H5Dopen(file_id, "/MAPS/Quantification/Calibration/Fitted/Calibration_Curve_DS_IC", H5P_DEFAULT);
        hid_t quant_space = H5Dget_space(ds_ic_quant_id);
        hid_t quant_type = H5Dget_type(ds_ic_quant_id);
        if(dset_id > 0 && chan_names_id > 0 && scaler_dset_id > 0 && scaler_units_id > 0 &&  scaler_names_id > 0)
        {
            hsize_t chan_dims[3] = {1,1,1};
            hsize_t scaler_dims[3] = {1,1,1};
            hsize_t image_dims[3] = {1,1,1};
            hsize_t image_dims_single[1] = {1};
            hsize_t offset[3] = {0,0,0};
            hsize_t offset_image[3] = {0,0,0};
            hsize_t offset_single[3] = {0};

            hsize_t offset_quant[2] = {0,0};
            hsize_t count_quant[2] = {1,1};

            hid_t chan_type = H5Dget_type(dset_id);
            hid_t scalername_type = H5Dget_type(scaler_names_id);

            hid_t chan_space = H5Dget_space(dset_id);
            hid_t chan_name_space = H5Dget_space(chan_names_id);
            hid_t scaler_space = H5Dget_space(scaler_dset_id);

            H5Sget_simple_extent_dims(chan_space, &chan_dims[0], nullptr);
            H5Sget_simple_extent_dims(scaler_space, &scaler_dims[0], nullptr);

            image_dims_single[0] = chan_dims[0] + scaler_dims[0];
            image_dims[0] = chan_dims[0] + scaler_dims[0];
            image_dims[1] = chan_dims[1];
            image_dims[2] = chan_dims[2];
            hid_t image_space = H5Screate_simple(3, &image_dims[0], &image_dims[0]);
            hid_t image_single_space = H5Screate_simple(1, &image_dims_single[0], &image_dims_single[0]);

            image_dims[0] = 1;
            hid_t readwrite_space = H5Screate_simple(3, &image_dims[0], &image_dims[0]);

            image_dims_single[0] = {1};
            hid_t readwrite_single_space = H5Screate_simple(1, &image_dims_single[0], &image_dims_single[0]);


            hid_t image_dset_id = H5Dcreate(file_id, exchange_images.c_str(), chan_type, image_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            hid_t image_names_dset_id = H5Dcreate(file_id, exchange_image_names.c_str(), filetype, image_single_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            hid_t image_units_dset_id = H5Dcreate(file_id, exchange_image_units.c_str(), filetype, image_single_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            real_t *data = new real_t[chan_dims[1] * chan_dims[2]];
            real_t *ds_ic_data = new real_t[chan_dims[1] * chan_dims[2]];
            std::string scaler_name_str;
            char char_data[256]={0};
            char char_ug_data[256]="ug/cm2";
            int k =0;

            real_t quant_value = 1.0;

            for (std::string::size_type x=0; x<normalize_scaler.length(); ++x)
            {
                normalize_scaler[x] = std::tolower(normalize_scaler[x]);
            }
            // save scalers first
            for(hsize_t i=0; i < scaler_dims[0]; i++)
            {
                offset[0] = i;
                offset_single[0] = i;
                offset_image[0] = k;
                k++;
                H5Sselect_hyperslab (image_space, H5S_SELECT_SET, offset, nullptr, image_dims, nullptr);
                //read write values
                hid_t status = H5Dread(scaler_dset_id, chan_type, readwrite_space, image_space, H5P_DEFAULT, (void*)&data[0]);
                if(status > -1)
                {
                    H5Dwrite(image_dset_id, chan_type, readwrite_space, image_space, H5P_DEFAULT, (void*)&data[0]);
                }

                //read write names
                H5Sselect_hyperslab (image_single_space, H5S_SELECT_SET, offset_single, nullptr, image_dims_single, nullptr);
                status = H5Dread(scaler_names_id, scalername_type, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                if(status > -1)
                {
                    H5Dwrite(image_names_dset_id, filetype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                }


                scaler_name_str = std::string(char_data, 256);
                scaler_name_str.erase(std::remove(scaler_name_str.begin(), scaler_name_str.end(), ' '), scaler_name_str.end());
                //to lower
                for (std::string::size_type x=0; x<scaler_name_str.length(); ++x)
                {
                    scaler_name_str[x] = std::tolower(scaler_name_str[x]);
                }
                if(scaler_name_str == normalize_scaler)
                {
                    for(hsize_t z=0; z < (chan_dims[1] * chan_dims[2]); z++)
                    {
                        ds_ic_data[z] = data[z];
                    }
                }

                //read write units
                status = H5Dread(scaler_units_id, scalername_type, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                if(status > -1)
                {
                    H5Dwrite(image_units_dset_id, filetype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                }
            }

            for(hsize_t i=0; i < chan_dims[0]; i++)
            {
                offset[0] = i;
                offset_image[0] = k;
                offset_single[0] = i;
                k++;

                // read write names
                H5Sselect_hyperslab (chan_name_space, H5S_SELECT_SET, offset_single, nullptr, image_dims_single, nullptr);
                H5Sselect_hyperslab (image_single_space, H5S_SELECT_SET, offset_image, nullptr, image_dims_single, nullptr);
                hid_t status = H5Dread(chan_names_id, scalername_type, readwrite_single_space, chan_name_space, H5P_DEFAULT, (void*)&char_data[0]);
                if(status > -1)
                {
                    H5Dwrite(image_names_dset_id, filetype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                }

                // get quantification for ds_ic and store in quant_value
                if(ds_ic_quant_id > -1)
                {
                    std::string chan_name_str = std::string(char_data, 256);
                    chan_name_str.erase(std::remove(chan_name_str.begin(), chan_name_str.end(), ' '), chan_name_str.end());
                    data_struct::Element_Info* element = data_struct::Element_Info_Map::inst()->get_element(chan_name_str);
                    if(element != nullptr)
                    {
                        offset_quant[1] = element->number - 1;

                        H5Sselect_hyperslab (quant_space, H5S_SELECT_SET, offset_quant, nullptr, count_quant, nullptr);
                        hid_t status = H5Dread(ds_ic_quant_id, quant_type, readwrite_single_space, quant_space, H5P_DEFAULT, (void*)&quant_value);
                        if(status < 0)
                        {
                            quant_value = 1.0;
                        }
                    }
                }
                else
                {
                    quant_value = 1.0;
                }


                H5Sselect_hyperslab (chan_space, H5S_SELECT_SET, offset, nullptr, image_dims, nullptr);
                H5Sselect_hyperslab (image_space, H5S_SELECT_SET, offset_image, nullptr, image_dims, nullptr);
                //read write values
                status = H5Dread(dset_id, chan_type, readwrite_space, chan_space, H5P_DEFAULT, (void*)&data[0]);
                if(status > -1)
                {
                    for(hsize_t z=0; z < (chan_dims[1] * chan_dims[2]); z++)
                    {
                        data[z] = data[z] / quant_value / ds_ic_data[z];
                    }
                    H5Dwrite(image_dset_id, chan_type, readwrite_space, image_space, H5P_DEFAULT, (void*)&data[0]);
                }

                //read write units
                //status = H5Dread(chan_units_id, scalername_type, readwrite_single_space, chan_name_space, H5P_DEFAULT, (void*)&char_data[0]);
                if(status > -1)
                {
                    H5Dwrite(image_units_dset_id, filetype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_ug_data[0]);
                }
            }

            delete [] data;
            delete [] ds_ic_data;
            H5Dclose(ds_ic_quant_id);
            H5Dclose(image_dset_id);
            H5Dclose(image_names_dset_id);
            H5Dclose(image_units_dset_id);
            //H5Dclose(chan_units_id);
            H5Dclose(chan_names_id);
            H5Dclose(scaler_units_id);
            H5Dclose(scaler_names_id);
            H5Dclose(scaler_dset_id);
            H5Dclose(dset_id);
        }

    }
    else
    {
        if( H5Lcreate_hard(file_id, str_scalers.c_str(), H5L_SAME_LOC, exchange_scalers.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW  << " Couldn't create soft link for"<< exchange_scalers <<  "\n";
        }

        if( H5Lcreate_hard(file_id, str_scaler_names.c_str(), H5L_SAME_LOC, exchange_scaler_names.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW  << " Couldn't create soft link for"<< exchange_scaler_names <<  "\n";
        }

        if( H5Lcreate_hard(file_id, str_scaler_units.c_str(), H5L_SAME_LOC, exchange_scaler_units.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW  << " Couldn't create soft link for"<< exchange_scaler_units <<  "\n";
        }

        if(fits_grp > 0)
        {
            std::string str_images = "/MAPS/XRF_Analyzed/" + fits_link + "/Counts_Per_Sec";
            std::string str_image_names = "/MAPS/XRF_Analyzed/" + fits_link + "/Channel_Names";
            std::string str_image_units = "/MAPS/XRF_Analyzed/" + fits_link + "/Channel_Units";


            if( H5Lcreate_hard(file_id, str_images.c_str(), H5L_SAME_LOC, exchange_images.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
            {
                logW  << " Couldn't create soft link for"<< exchange_images <<  "\n";
            }

            if( H5Lcreate_hard(file_id, str_image_names.c_str(), H5L_SAME_LOC, exchange_image_names.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
            {
                logW  << " Couldn't create soft link for"<< exchange_image_names <<  "\n";
            }

            if( H5Lcreate_hard(file_id, str_image_units.c_str(), H5L_SAME_LOC, exchange_image_units.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
            {
                logW  << " Couldn't create soft link for"<< exchange_image_units <<  "\n";
            }

            // Add quantification
            hsize_t quant_dims[3] = {3,1,1};
            hid_t chan_names = H5Dopen(file_id, str_image_names.c_str(), H5P_DEFAULT);
            hid_t chan_space = H5Dget_space(chan_names);
            if(chan_names > -1)
            {
                hsize_t chan_size = 1;
                H5Sget_simple_extent_dims(chan_space, &chan_size, nullptr);
                quant_dims[2] = chan_size; //num channel names
            }
            hid_t quant_space = H5Screate_simple(3, &quant_dims[0], &quant_dims[0]);
            _add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], fits_link, exhange_str+"/quant");
            H5Sclose(quant_space);
            H5Dclose(chan_names);
        }

        //Save description
        hid_t dset_id = H5Dcreate(file_id, str_desc.c_str(), filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(dset_id > 0)
        {
            fits_link.copy(desc, 256);
            H5Dwrite(dset_id, filetype, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&desc);
            H5Dclose(dset_id);
        }
    }




    //Add version dataset
    real_t save_val = HDF5_EXCHANGE_VERSION;
    hid_t dset_id = H5Dcreate(file_id, str_version.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if(dset_id > 0)
    {
        H5Dwrite(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
    }

    //Add theta
    if( H5Lcreate_hard(file_id, str_scan_theta.c_str(), H5L_SAME_LOC, str_theta.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << "Couldn't create soft link for"<< str_theta <<  "\n";
    }


    H5Sclose(dataspace_id);
    H5Gclose(exchange_id);

    return true;
}

//-----------------------------------------------------------------------------

void HDF5_IO::add_exchange_layout(std::string dataset_file)
{
    std::lock_guard<std::mutex> lock(_mutex);

    logI  << dataset_file << "\n";
    hid_t saved_file_id = _cur_file_id;
    hid_t file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    std::string exchange_fits[4] = { STR_FIT_GAUSS_MATRIX, STR_FIT_ROI, STR_FIT_NNLS, STR_FIT_GAUSS_MATRIX };
    std::string exchange_normalize[4] = {STR_DS_IC, "", "", ""};

    int ex_idx = 0;
    for(int i=0; i < 4; i++)
    {
        if(true == _add_exchange_meta(file_id, std::to_string(ex_idx), exchange_fits[i], exchange_normalize[i]))
        {
            ex_idx++;
        }
    }
    if (H5Gget_objinfo(file_id, "/version", 0, NULL) < 0)
    {
        H5Lcreate_hard(file_id, "/MAPS/version", H5L_SAME_LOC, "/version", H5P_DEFAULT, H5P_DEFAULT);
    }

    _cur_file_id = file_id;
    end_save_seq();
    logI<<"closing file"<<"\n";

    _cur_file_id = saved_file_id;
}

//-----------------------------------------------------------------------------

void HDF5_IO::export_int_fitted_to_csv(std::string dataset_file)
{
    std::lock_guard<std::mutex> lock(_mutex);

    logI << dataset_file << "\n";

    hid_t file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    std::string exchange_fits[4] = { STR_FIT_GAUSS_MATRIX, STR_FIT_ROI, STR_FIT_NNLS };

    if (file_id > 0)
    {
        ArrayXr energy_array;
        ArrayXr int_spectra;
        ArrayXr model_spectra;
        ArrayXr background_array;
        std::string csv_path;


        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        hid_t dset_id, dataspace_id;
        hsize_t dims_in[1] = { 0 };

        if (_open_h5_object(dset_id, H5O_DATASET, close_map, "/MAPS/Spectra/Integrated_Spectra/Spectra", file_id))
        {
            dataspace_id = H5Dget_space(dset_id);
            close_map.push({ dataspace_id, H5O_DATASPACE });

            int rank = H5Sget_simple_extent_ndims(dataspace_id);
            if (rank == 1)
            {
                int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
                if (status_n > -1)
                {
                    int_spectra.resize(dims_in[0]);
                    H5Dread(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, int_spectra.data());
                }
            }
        }

        if (_open_h5_object(dset_id, H5O_DATASET, close_map, "/MAPS/Spectra/Energy", file_id))
        {
            if (dims_in[0] > 0)
            {
                dataspace_id = H5Dget_space(dset_id);
                close_map.push({ dataspace_id, H5O_DATASPACE });
                energy_array.resize(dims_in[0]);
                H5Dread(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, energy_array.data());
            }
        }

        for (const auto& fit_itr : exchange_fits)
        {

            int s_idx = dataset_file.find("img.dat");
            if (s_idx > 0)
            {
                csv_path = dataset_file.substr(0, s_idx);
                csv_path += "output";
                csv_path += DIR_END_CHAR;
                csv_path += dataset_file.substr(s_idx+8);
            }
            else
            {
                csv_path = dataset_file;
            }
            csv_path += "_";
            csv_path += fit_itr;
            csv_path += ".csv";

            std::string h5_path = "/MAPS/XRF_Analyzed/" + fit_itr + "/"+ STR_FIT_INT_BACKGROUND;
            if (_open_h5_object(dset_id, H5O_DATASET, close_map, h5_path.c_str(), file_id))
            {
                if (dims_in[0] > 0)
                {
                    background_array.resize(dims_in[0]);
                    dataspace_id = H5Dget_space(dset_id);
                    close_map.push({ dataspace_id, H5O_DATASPACE });
                    H5Dread(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, background_array.data());
                }
            }

            h5_path = "/MAPS/XRF_Analyzed/" + fit_itr + "/"+ STR_FIT_INT_SPEC;
            if (_open_h5_object(dset_id, H5O_DATASET, close_map, h5_path.c_str(), file_id))
            {
                if (dims_in[0] > 0)
                {
                    model_spectra.resize(dims_in[0]);
                    dataspace_id = H5Dget_space(dset_id);
                    close_map.push({ dataspace_id, H5O_DATASPACE });
                    H5Dread(dset_id, H5T_NATIVE_REAL, dataspace_id, dataspace_id, H5P_DEFAULT, model_spectra.data());
                    csv::save_fit_and_int_spectra(csv_path, &energy_array, &int_spectra, &model_spectra, &background_array);
                }
            }
        }
        _close_h5_objects(close_map);
    }


    hid_t saved_file_id = _cur_file_id;
    _cur_file_id = file_id;
    end_save_seq();
    logI << "closing file" << "\n";
    _cur_file_id = saved_file_id;
}

} //end namespace file
}// end namespace io
