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

#define HDF5_SAVE_VERSION 10.0

#define HDF5_EXCHANGE_VERSION 1.0

const std::vector<std::string> hdf5_copy_dset_names = {"Element_Weights",
                                                       "Element_Weights_Names",
                                                       "DS_IC",
                                                       "US_IC",
                                                       "Standard_Name",
                                                       "Channel_Names",
													   "Channel_Units",
                                                       "Names",
                                                       "version"
                                                      };


const std::vector<std::string> hdf5_copy_grp_names = {"Scan"
                                                      };

const std::vector<std::string> xrf_analysis_save_names = {"ROI",
                                                          "Params",
                                                          "Fitted",
                                                          "SVD",
                                                          "NNLS"
                                                         };





namespace io
{
namespace file
{

std::mutex HDF5_IO::_mutex;

HDF5_IO* HDF5_IO::_this_inst(0);


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

bool HDF5_IO::_open_h5_object(hid_t &id, H5_OBJECTS obj, std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map, std::string s1, hid_t id2, bool log_error)
{
    if (obj == H5O_FILE)
    {
        id = H5Fopen(s1.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if(id < 0)
        {
            _close_h5_objects(close_map);
			if (log_error)
			{
				logE << "Failed to open file " << s1 << "\n";
			}
            return false;
        }
    }
    else if (obj == H5O_GROUP)
    {
        id = H5Gopen(id2, s1.c_str(), H5P_DEFAULT);
        if(id < 0)
        {
            _close_h5_objects(close_map);
			if (log_error)
			{
				logE << "Failed to open group " << s1 << "\n";
			}
            return false;
        }
    }
    else if (obj ==  H5O_DATASET)
    {
        id = H5Dopen2(id2, s1.c_str(), H5P_DEFAULT);
        if(id < 0)
        {
            _close_h5_objects(close_map);
			if (log_error)
			{
				logE << "Failed to open dataset " << s1 << "\n";
			}
            return false;
        }
    }
    else
    {
        _close_h5_objects(close_map);
		if (log_error)
		{
			logE << "Unknown H5_OBJECT type " << obj << "\n";
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

	std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
	if (log_error)
	{
		logI << path << " detector : " << detector_num << "\n";
	}
   hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_detectors_id;
   //hid_t    dset_xpos_id, dset_ypos_id, dataspace_xpos_id, dataspace_ypos_id;
   hid_t    dataspace_detectors_id;
   hid_t    attr_detector_names_id, attr_timebase_id;
   herr_t   error;
   std::string detector_path;
   char* detector_names[256];
   real_t time_base = 1.0f;
   real_t el_time = 1.0f;
   real_t * buffer;
   hsize_t offset_row[2] = {0,0};
   hsize_t count_row[2] = {0,0};
   hsize_t offset_meta[3] = {0,0,0};
   hsize_t count_meta[3] = {1,1,1};
   std::unordered_map<std::string, int> detector_lookup;
   std::string incnt_str = "ICR Ch ";
   std::string outcnt_str = "OCR Ch ";

   switch(detector_num)
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

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1, log_error) )
        return false;

    if ( false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "2D Scan", file_id) )
       return false;

    if ( false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id) )
       return false;
    dataspace_id = H5Dget_space(dset_id);
    close_map.push({dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(dset_detectors_id, H5O_DATASET, close_map, "Detectors", maps_grp_id) )
       return false;
    dataspace_detectors_id = H5Dget_space(dset_detectors_id);
    close_map.push({dataspace_detectors_id, H5O_DATASPACE});
    attr_detector_names_id=H5Aopen(dset_detectors_id, "Detector Names", H5P_DEFAULT);
    close_map.push({dataspace_detectors_id, H5O_ATTRIBUTE});
    attr_timebase_id=H5Aopen(dset_detectors_id, "TIMEBASE", H5P_DEFAULT);
    close_map.push({dataspace_detectors_id, H5O_ATTRIBUTE});

//    if ( false == _open_h5_object(dset_xpos_id, H5O_DATASET, close_map, "X Positions", maps_grp_id) )
//       return false;
//    dataspace_xpos_id = H5Dget_space(dset_xpos_id);
//    close_map.push({dataspace_xpos_id, H5O_DATASPACE});

//    if ( false == _open_h5_object(dset_ypos_id, H5O_DATASET, close_map, "Y Positions", maps_grp_id) )
//       return false;
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
		return false;
	}

    if(spec_vol->rows() < dims_in[0] || spec_vol->cols() < dims_in[1] || spec_vol->samples_size() < dims_in[2])
    {
        spec_vol->resize_and_zero(dims_in[0], dims_in[1], dims_in[2]);
    }



    int det_rank = H5Sget_simple_extent_ndims(dataspace_detectors_id);
    hsize_t* det_dims_in = new hsize_t[det_rank];
    H5Sget_simple_extent_dims(dataspace_detectors_id, &det_dims_in[0], nullptr);

    hid_t ftype = H5Aget_type(attr_detector_names_id);
    hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
    error = H5Aread(attr_detector_names_id, type, detector_names);
    if(error == 0)
    {
        //generate lookup map
        for(size_t z = 0; z < det_dims_in[2]; z++)
        {
            detector_lookup[std::string(detector_names[z])] = z;
            //free(detector_names[z]);
        }
    }
    delete [] det_dims_in;


    error = H5Aread(attr_timebase_id, H5T_NATIVE_REAL, &time_base);

    count[0] = 1; //1 row

    memoryspace_id = H5Screate_simple(2, count_row, nullptr);
    close_map.push({memoryspace_id, H5O_DATASPACE});
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

                 offset_meta[2] = detector_lookup["Timer"];
                 H5Sselect_hyperslab (dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                 error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &live_time);
                 el_time = live_time/time_base;
                 spectra->elapsed_livetime(el_time);

//                 H5Sselect_hyperslab (dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
//                 error = H5Dread(dset_rt_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
//                 spectra->elapsed_realtime(real_time);

                 offset_meta[2] = detector_lookup[incnt_str];
                 H5Sselect_hyperslab (dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                 error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &in_cnt);
                 spectra->input_counts(in_cnt*1000.0);

                 offset_meta[2] = detector_lookup[outcnt_str];
                 H5Sselect_hyperslab (dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                 error = H5Dread(dset_detectors_id, H5T_NATIVE_REAL, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &out_cnt);
                 spectra->output_counts(out_cnt*1000.0);

//                 spectra->recalc_elapsed_livetime();

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
    data_struct::Spectra_Volume *spec_vol)
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
   herr_t   error;
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

             logI<<"read row "<<row<<"\n";
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

    hid_t    dset_id, spec_grp_id, int_spec_grp_id, dataspace_id, memoryspace_id, memoryspace_time_id, dataspace_time_id, file_time_id, filespace_id, maps_grp_id, dcpl_id;
    hid_t   dset_rt_id, dset_lt_id, incnt_dset_id, outcnt_dset_id;
    //herr_t   error;
    int status = 0;

    hsize_t chunk_dims[3];
    hsize_t dims_out[3];
    hsize_t offset[3];
    hsize_t count[3];
    hsize_t dims_time_out[2];
    hsize_t offset_time[2];
    hsize_t count_time[2];

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
    filespace_id = H5Screate_simple(3, dims_out, nullptr);

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    dataspace_id = H5Screate_simple (3, dims_out, nullptr);

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
    dset_id = H5Dcreate (spec_grp_id, path.c_str(), H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);

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
            H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

            H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);

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
    H5Sclose(filespace_id);
    H5Sclose(dataspace_time_id);
    H5Sclose(dataspace_id);
    H5Pclose(dcpl_id);


    int_spec_grp_id = H5Gopen(spec_grp_id, "Integrated_Spectra", H5P_DEFAULT);
    if(int_spec_grp_id < 0)
        int_spec_grp_id = H5Gcreate(spec_grp_id, "Integrated_Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    //save quantification_standard integrated spectra
    data_struct::Spectra spectra = spectra_volume->integrate();
    count[0] = spectra.size();
    memoryspace_id = H5Screate_simple(1, count, nullptr);
    filespace_id = H5Screate_simple(1, count, nullptr);
    dataspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (int_spec_grp_id, "Spectra", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&spectra[0]);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(filespace_id);
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
    filespace_id = H5Screate_simple(1, count, nullptr);
    dataspace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (spec_grp_id, "Energy", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&out_vec[0]);
    H5Dclose(dset_id);
    H5Sclose(memoryspace_id);
    H5Sclose(filespace_id);
    H5Sclose(dataspace_id);

    //save energy calibration
    count[0] = 3;
    dataspace_id = H5Screate_simple (1, count, nullptr);
    filespace_id = H5Screate_simple (1, count, nullptr);
    dset_id = H5Dcreate (spec_grp_id, "Energy_Calibration", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    count[0] = 1;
    memoryspace_id = H5Screate_simple(1, count, nullptr);
    offset[0] = 0;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_offset);
    offset[0] = 1;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_slope);
    offset[0] = 2;
    H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
    status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&energy_quad);
    H5Dclose(dset_id);
    H5Sclose(filespace_id);
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

    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl_id, 3, chunk_dims);
    H5Pset_deflate (dcpl_id, 7);

    memoryspace = H5Screate_simple(3, count_3d, nullptr);
    filespace = H5Screate_simple(3, dims_out, nullptr);

    dataspace_id = H5Screate_simple (3, dims_out, nullptr);

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
    if(dset_id < 0)
        dset_id = H5Dcreate (fit_grp_id, "Counts_Per_Sec", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
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

    element_lines.push_back(STR_COHERENT_SCT_AMPLITUDE);
    element_lines.push_back(STR_COMPTON_AMPLITUDE);
	element_lines.push_back(STR_SUM_ELASTIC_INELASTIC_AMP);
    element_lines.push_back(STR_TOTAL_FLUORESCENCE_YIELD);
    element_lines.push_back(STR_NUM_ITR);

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
		if (el_name != STR_NUM_ITR)
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

void HDF5_IO::save_quantifications(std::map<string, data_struct::Quantification_Standard*> &quants)
{
    for(auto& itr: quants)
    {
        HDF5_IO::save_quantification(itr.second);
    }
}

//-----------------------------------------------------------------------------

bool HDF5_IO::save_quantification(data_struct::Quantification_Standard * quantification_standard)
{
    std::lock_guard<std::mutex> lock(_mutex);

    

//hid_t error_stack = H5Eget_current_stack();
//H5Eset_auto2(error_stack, nullptr, nullptr);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    hid_t    dset_id, memoryspace_id, filespace_id, dataspace_id, filetype, dataspace_ch_id, memtype,  dset_ch_id, q_int_spec_grp_id;
    hid_t   memtype_label, filetype_label, q_memoryspace_label_id, q_dataspace_label_id;
    //herr_t   error;

    hid_t q_dataspace_id, q_memoryspace_id, q_filespace_id, q_dset_id, q_grp_id, q_fit_grp_id, maps_grp_id, status, scalers_grp_id, xrf_fits_grp_id;
    hid_t dset_labels_id;
    hsize_t offset[3];
    hsize_t count[3];

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
//    chunk_dims[0] = 1;
//    chunk_dims[1] = 1;
//    chunk_dims[2] = dims_out[2];

//    dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
//    H5Pset_chunk(dcpl_id, 3, chunk_dims);
//    H5Pset_deflate (dcpl_id, 7);

//    memoryspace = H5Screate_simple(3, count, nullptr);
//    filespace = H5Screate_simple(3, dims_out, nullptr);

//    dataspace_id = H5Screate_simple (3, dims_out, nullptr);

//    dataspace_ch_id = H5Screate_simple (1, dims_out, nullptr);


//    H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset, nullptr, count, nullptr);

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

    if( quantification_standard != nullptr && quantification_standard->processed() == true)
    {

        q_grp_id = H5Gopen(maps_grp_id, "Quantification", H5P_DEFAULT);
        if(q_grp_id < 0)
            q_grp_id = H5Gcreate(maps_grp_id, "Quantification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(q_grp_id < 0)
        {
            logE<<"creating group MAPS/Quantification"<<"\n";
            return false;
        }

        scalers_grp_id = H5Gopen(q_grp_id, "Scalers", H5P_DEFAULT);
        if(scalers_grp_id < 0)
            scalers_grp_id = H5Gcreate(q_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(scalers_grp_id < 0)
        {
            logE<<"creating group MAPS/Quantification/Scalers"<<"\n";
            return false;
        }

        xrf_fits_grp_id = H5Gopen(q_grp_id, "XRF_Analyzed", H5P_DEFAULT);
        if(xrf_fits_grp_id < 0)
            xrf_fits_grp_id = H5Gcreate(q_grp_id, "XRF_Analyzed", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(xrf_fits_grp_id < 0)
        {
            logE<<"creating group MAPS/Quantification/XRF_Analyzed"<<"\n";
            return false;
        }

        //save quantification_standard element weights
        const std::unordered_map<std::string, data_struct::Element_Quant> e_weights = quantification_standard->element_quants;
        count[0] = e_weights.size();
        memoryspace_id = H5Screate_simple(1, count, nullptr);
        filespace_id = H5Screate_simple(1, count, nullptr);
        dataspace_id = H5Screate_simple (1, count, nullptr);
        dataspace_ch_id = H5Screate_simple (1, count, nullptr);
        dset_id = H5Dcreate (q_grp_id, "Element_Weights", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        dset_ch_id = H5Dcreate (q_grp_id, "Element_Weights_Names", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        offset[0] = 0;
        count[0] = 1;
        H5Sselect_hyperslab (memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        int offset_idx = 0;
        for(auto itr: e_weights)
        {
            offset[0] = offset_idx;
            H5Sselect_hyperslab (filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

            status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&(itr.second.weight));
            char tmp_char[256] = {0};
            itr.first.copy(tmp_char, 254);
            status = H5Dwrite (dset_ch_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)tmp_char);

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
            logE<<"creating group MAPS/Quantification/Integrated_Spectra"<<"\n";
            return false;
        }

        //save quantification_standard integrated spectra
        data_struct::Spectra spectra = quantification_standard->integrated_spectra;
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
        dataspace_id = H5Screate_simple (1, count, nullptr);
        memoryspace_id = H5Screate_simple (1, count, nullptr);
        dset_ch_id = H5Dcreate (q_grp_id, "Standard_Name", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        char tmp_char[255] = {0};
        quantification_standard->standard_filename.copy(tmp_char, 254);
        status = H5Dwrite (dset_ch_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
        H5Dclose(dset_ch_id);
        H5Sclose(dataspace_id);

        //save sr_current
        dataspace_id = H5Screate_simple (1, count, nullptr);
        dset_id = H5Dcreate (scalers_grp_id, "SR_Current", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->sr_current));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save us_ic
        dataspace_id = H5Screate_simple (1, count, nullptr);
        dset_id = H5Dcreate (scalers_grp_id, "US_IC", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->US_IC));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save ds_ic
        dataspace_id = H5Screate_simple (1, count, nullptr);
        dset_id = H5Dcreate (scalers_grp_id, "DS_IC", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quantification_standard->DS_IC));
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);


        //save real_time
        real_t save_val = spectra.elapsed_realtime();
        dataspace_id = H5Screate_simple (1, count, nullptr);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Elapsed_Realtime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(dataspace_id);

        //save life_time
        save_val = spectra.elapsed_livetime();
        dataspace_id = H5Screate_simple (1, count, nullptr);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Elapsed_Livetime", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(dataspace_id);

        //save input counts
        save_val = spectra.input_counts();
        dataspace_id = H5Screate_simple (1, count, nullptr);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Input_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(dataspace_id);

        //save output counts
        save_val = spectra.output_counts();
        dataspace_id = H5Screate_simple (1, count, nullptr);
        dset_id = H5Dcreate (q_int_spec_grp_id, "Output_Counts", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        H5Dclose(dset_id);
        H5Sclose(memoryspace_id);
        H5Sclose(dataspace_id);

        //save calibration curves
        //data_struct::Quantifiers quantifiers = quantification_standard->quantifier_map.at(path);
        //auto shell_itr = quantification_standard->_calibration_curves.begin();

        q_dims_out[0] = 3;// shells K, L, and M
        q_dims_out[1] = 93; // elements 1 - 92

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
        q_dataspace_label_id = H5Screate_simple (2, q_dims_out, nullptr);
        q_dataspace_id = H5Screate_simple (2, q_dims_out, nullptr);
        //q_dataspace_ch_id = H5Screate_simple (1, q_dims_out, nullptr);

		if (quantification_standard->element_counts.size() > 0)
		{
			for (auto& qitr : quantification_standard->quantifier_map)
			{

				q_fit_grp_id = H5Gopen(xrf_fits_grp_id, qitr.first.c_str(), H5P_DEFAULT);
				if (q_fit_grp_id < 0)
					q_fit_grp_id = H5Gcreate(xrf_fits_grp_id, qitr.first.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (q_fit_grp_id < 0)
				{
					logE << "creating group MAPS/Quantification/" << qitr.first << "\n";
					return false;
				}

				//save quantification_standard counts
				const std::unordered_map<std::string, std::unordered_map<std::string, real_t> > all_element_counts = quantification_standard->element_counts;
				std::unordered_map<std::string, real_t> element_counts = all_element_counts.at(qitr.first);
				count[0] = element_counts.size();
				dataspace_id = H5Screate_simple(1, count, nullptr);
				memoryspace_id = H5Screate_simple(1, count, nullptr);
				filespace_id = H5Screate_simple(1, count, nullptr);
				dset_id = H5Dcreate(q_fit_grp_id, "Counts_Per_Sec", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				dset_ch_id = H5Dcreate(q_fit_grp_id, "Channel_Names", filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
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
					status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&val);

					offset_idx++;
				}
				H5Dclose(dset_id);
				H5Dclose(dset_ch_id);
				H5Sclose(memoryspace_id);
				H5Sclose(filespace_id);
				H5Sclose(dataspace_id);


				dset_labels_id = H5Dopen(q_fit_grp_id, "Calibration_Curve_Labels", H5P_DEFAULT);
				if (dset_labels_id < 0)
					dset_labels_id = H5Dcreate(q_fit_grp_id, "Calibration_Curve_Labels", filetype_label, q_dataspace_label_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);


				for (size_t i = 0; i < qitr.second.calib_curves.size(); i++)
				{

					std::string q_dset_name = "Calibration_Curve_" + qitr.second.calib_curves[i].quantifier_name;

					q_dset_id = H5Dopen(q_fit_grp_id, q_dset_name.c_str(), H5P_DEFAULT);
					if (q_dset_id < 0)
						q_dset_id = H5Dcreate(q_fit_grp_id, q_dset_name.c_str(), H5T_INTEL_R, q_dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
					if (q_dset_id < 0)
					{
						logE << "creating dataset MAPS/Quantification/" << q_dset_name << "/" << q_dset_name << "\n";
						continue;
					}


					//for(auto& shell_itr : quant_itr.second)
					for (size_t j = 0; j < 3; j++)
					{
						char label[10] = { 0 };
						//int element_offset = 0;
						//create dataset for different shell curves
						std::vector<real_t> calibration_curve = qitr.second.calib_curves[i].shell_curves[j];
						std::vector<std::string> calibration_curve_labels = qitr.second.calib_curves[i].shell_curves_labels[j];
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

					}
					H5Dclose(q_dset_id);
				}

				H5Dclose(dset_labels_id);
				H5Gclose(q_fit_grp_id);
			}
		}
        //H5Dclose(q_dset_ch_id);
        H5Sclose(q_filespace_id);
        H5Sclose(q_memoryspace_id);
        H5Sclose(q_dataspace_label_id);
        H5Sclose(q_memoryspace_label_id);
        H5Sclose(q_dataspace_id);
        //H5Sclose(q_dataspace_ch_id);
        H5Gclose(xrf_fits_grp_id);
        H5Gclose(scalers_grp_id);
        H5Gclose(q_int_spec_grp_id);
        H5Gclose(q_grp_id);



    }

    H5Tclose(filetype);
    H5Tclose(memtype);

    H5Tclose(filetype_label);
    H5Tclose(memtype_label);



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

bool HDF5_IO::_save_scan_meta_data(hid_t scan_grp_id, struct mda_file *mda_scalers, data_struct::Params_Override * params_override)
{
    hid_t dataspace_id = -1, memoryspace_id = -1, filespace_id = -1;
	hid_t status;
	hid_t filetype, memtype;
    hid_t dset_id = -1;
	
	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };

	bool single_row_scan = false;

    

	try
	{

        filetype = H5Tcopy(H5T_FORTRAN_S1);
		H5Tset_size(filetype, 256);
		memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);


        memoryspace_id = H5Screate_simple(1, count, nullptr);

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
                dataspace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(scan_grp_id, "y_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
                    logE << "creating dataset 'y_axis'" << "\n";
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
					return false;
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);

                if(mda_scalers->scan->last_point == 0)
                    count[0] = 1;
                else
                    count[0] = mda_scalers->scan->last_point;
                dataspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(scan_grp_id, "x_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
                    logE << "creating dataset 'x_axis'" << "\n";
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->positioners_data[0][i];
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);


                //save requested rows
                count[0] = 1;
                dataspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(scan_grp_id, "requested_rows", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset_id < 0)
                {
                    logE << "creating dataset 'requested_rows'" << "\n";
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
                    return false;
                }
                int value = 1;
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&value);
                H5Dclose(dset_id);
                H5Sclose(dataspace_id);
                H5Sclose(filespace_id);

                //save requested cols
                count[0] = 1;
                dataspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(scan_grp_id, "requested_cols", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset_id < 0)
                {
                    logE << "creating dataset 'requested_cols'" << "\n";
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
                    return false;
                }
                value = mda_scalers->header->dimensions[0];
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&value);
                H5Dclose(dset_id);
                H5Sclose(dataspace_id);
                H5Sclose(filespace_id);

			}
			else
			{
                //save y axis
                if(mda_scalers->scan->last_point == 0)
                    count[0] = 1;
                else
                    count[0] = mda_scalers->scan->last_point;
                dataspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(scan_grp_id, "y_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
                    logE << "creating dataset 'y_axis'" << "\n";
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->positioners_data[0][i];
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);

                //save x axis
                if(mda_scalers->scan->sub_scans[0]->last_point == 0)
                    count[0] = 1;
                else
                    count[0] = mda_scalers->scan->sub_scans[0]->last_point;
                dataspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(scan_grp_id, "x_axis", H5T_INTEL_F64, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				if (dset_id < 0)
				{
                    logE << "creating dataset 'x_axis'" << "\n";
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
					return false;
				}
				count[0] = 1;
				for (int32_t i = 0; i < mda_scalers->scan->sub_scans[0]->last_point; i++)
				{
					offset[0] = i;
					double pos = mda_scalers->scan->sub_scans[0]->positioners_data[0][i];
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&pos);
				}
				H5Dclose(dset_id);
				H5Sclose(dataspace_id);
				H5Sclose(filespace_id);


                //save requested rows
                count[0] = 1;
                dataspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(scan_grp_id, "requested_rows", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset_id < 0)
                {
                    logE << "creating dataset 'requested_rows'" << "\n";
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
                    return false;
                }
                int value = mda_scalers->header->dimensions[0];
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&value);
                H5Dclose(dset_id);
                H5Sclose(dataspace_id);
                H5Sclose(filespace_id);

                //save requested cols
                count[0] = 1;
                dataspace_id = H5Screate_simple(1, count, nullptr);
                filespace_id = H5Screate_simple(1, count, nullptr);
                dset_id = H5Dcreate(scan_grp_id, "requested_cols", H5T_INTEL_I32, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                if (dset_id < 0)
                {
                    logE << "creating dataset 'requested_cols'" << "\n";
                    H5Sclose(dataspace_id);
                    H5Sclose(filespace_id);
                    return false;
                }
                value = mda_scalers->header->dimensions[1];
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)&value);
                H5Dclose(dset_id);
                H5Sclose(dataspace_id);
                H5Sclose(filespace_id);

			}

		}

		//save write date
		count[0] = 1;

        //Save theta
        dset_id = H5Dcreate(scan_grp_id, "theta", H5T_NATIVE_REAL, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(dset_id > 0)
        {
            real_t theta = 0.0;
            struct mda_pv * pv = nullptr;
            //find theta by param_override->theta_pv in extra names
            for (int16_t i = 0; i < mda_scalers->extra->number_pvs; i++)
            {
                pv = mda_scalers->extra->pvs[i];
                if(pv == nullptr)
                {
                    continue;
                }

                if (pv->name != nullptr && params_override->theta_pv.compare(pv->name) == 0)
                {
                    break;
                }
            }
            if( pv != nullptr)
            {
                switch (pv->type)
                {
                case EXTRA_PV_STRING:
                    //str_val = std::string(pv->values);
                    break;
                case EXTRA_PV_INT16:
                    //s_val = (short*)pv->values;
                    //str_val = std::to_string(*s_val);
                    break;
                case EXTRA_PV_INT32:
                    //i_val = (int*)pv->values;
                    //str_val = std::to_string(*i_val);
                    break;
                case EXTRA_PV_FLOAT:
                    theta = *((float*)pv->values);
                    //str_val = std::to_string(*f_val);
                    break;
                case EXTRA_PV_DOUBLE:
                    theta = *((double*)pv->values);
                    //str_val = std::to_string(*d_val);
                    break;
                }

                status = H5Dwrite(dset_id, H5T_NATIVE_REAL, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)&theta);
            }
            H5Dclose(dset_id);
            dset_id = -1;
        }

		dset_id = H5Dcreate(scan_grp_id, "scan_time_stamp", filetype, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
            logE << "creating dataset 'scan_time_stamp'" << "\n";
			return false;
		}
		if (mda_scalers->scan->time != nullptr)
		{
			std::string str_time = std::string(mda_scalers->scan->time);
            char tmp_char[255] = {0};
            str_time.copy(tmp_char, 254);
            status = H5Dwrite(dset_id, memtype, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)tmp_char);
		}
        if(dset_id > -1)
        {
            H5Dclose(dset_id);
            dset_id = -1;
        }
		dset_id = H5Dcreate(scan_grp_id, "name", filetype, memoryspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (dset_id < 0)
		{
            logE << "creating dataset 'name'" << "\n";
			return false;
		}
		if(mda_scalers->scan->name != nullptr)
		{
			std::string str_name = std::string(mda_scalers->scan->name);
            char tmp_char[255] = {0};
            str_name.copy(tmp_char, 254);
            status = H5Dwrite(dset_id, memtype, memoryspace_id, memoryspace_id, H5P_DEFAULT, (void*)tmp_char);
		}

        H5Tclose(filetype);
        H5Tclose(memtype);

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
        logE << "saving MAPS/Scan meta data" << "\n";
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

    

	if (mda_scalers->extra == nullptr)
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
		count[0] = (size_t)mda_scalers->extra->number_pvs;
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
                char tmp_char[255] = {0};
                name.copy(tmp_char, 254);
                status = H5Dwrite(dset_id, memtype, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)tmp_char);
			}

            char tmp_char2[255] = {0};
            str_val.copy(tmp_char2, 254);
            status = H5Dwrite(dset_val_id, memtype, memoryspace_id2, filespace_id2, H5P_DEFAULT, (void*)tmp_char2);
			if (pv->description != nullptr)
			{
                std::string description = std::string(pv->description);
                char tmp_char[255] = {0};
                description.copy(tmp_char, 254);
                status = H5Dwrite(dset_desc_id, memtype, memoryspace_id3, filespace_id3, H5P_DEFAULT, (void*)tmp_char);
			}

			if (pv->unit != nullptr)
			{
                std::string unit = std::string(pv->unit);
                char tmp_char[255] = {0};
                unit.copy(tmp_char, 254);
                status = H5Dwrite(dset_unit_id, memtype, memoryspace_id4, filespace_id4, H5P_DEFAULT, (void*)tmp_char);
			}
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

bool HDF5_IO::_save_scalers(hid_t maps_grp_id, struct mda_file *mda_scalers, data_struct::Spectra_Volume * spectra_volume, data_struct::Params_Override * params_override, bool hasNetcdf)
{
    hid_t dataspace_id = -1, memoryspace_id = -1, filespace_id = -1, filespace_name_id = -1, memoryspace_str_id = -1;
    hid_t filetype, memtype;
    hid_t dset_names_id = -1;
	hid_t dset_units_id = -1;
    hid_t scalers_grp_id = -1;
    hid_t dcpl_id = -1, status;

    hsize_t offset[1] = { 0 };
    hsize_t count[1] = { 1 };

    //hsize_t offset_2d[2] = { 0, 0 };
	hsize_t count_2d[2] = { 1, 1 };

    hsize_t offset_3d[3] = { 0, 0, 0 };
    hsize_t count_3d[3] = { 1, 1, 1 };

    int mda_time_scaler_idx = -1;
    real_t time_scaler_clock = 1.0;

    MDA_IO mda_io;

    bool single_row_scan = false;

    Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> scaler_mat, abs_cfg_mat, H_dpc_cfg_mat, V_dpc_cfg_mat, dia1_dpc_cfg_mat, dia2_dpc_cfg_mat;

    

    //don't save these scalers
    std::list<std::string> ignore_scaler_strings = { "ELT1", "ERT1", "ICR1", "OCR1" };
    std::list<struct scaler_struct> scalers;
    std::list<data_struct::Summed_Scaler> summed_scalers;
    hid_t dset_cps_id = -1;

	memoryspace_str_id = H5Screate_simple(1, count, NULL);

    try
    {

        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
        memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

        count[0] = 1;

        real_t val;
		std::string units;
        bool save_cfg_abs = false;
        if (params_override != nullptr && mda_scalers->scan->scan_rank > 1)
        {
            int us_ic_idx = -1;
            int ds_ic_idx = -1;
            int cfg_2_idx = -1;
            int cfg_3_idx = -1;
            int cfg_4_idx = -1;
            int cfg_5_idx = -1;

            int hdf_idx = 0;

            if(params_override->time_scaler_clock.length() > 0)
            {
                time_scaler_clock = std::stod(params_override->time_scaler_clock);
            }

            for (auto itr : params_override->scaler_pvs)
            {
                //don't save ELT1, ERT1, ICR1, OCR1. these are saved elsewhere
                std::list<std::string>::iterator s_itr = std::find(ignore_scaler_strings.begin(), ignore_scaler_strings.end(), itr.first);
                if (s_itr != ignore_scaler_strings.end())
                    continue;

                int mda_idx = mda_io.find_scaler_index(mda_scalers, itr.second, val, units);
                scalers.push_back(scaler_struct(itr.first, units, mda_idx, hdf_idx, false));
                hdf_idx++;
                std::string scaler_name = itr.first;
                std::transform(scaler_name.begin(), scaler_name.end(), scaler_name.begin(), ::toupper);
                if (mda_idx > -1)
                {
                    if (scaler_name == "US_IC")
                        us_ic_idx = mda_idx;
                    else if (scaler_name == "DS_IC")
                        ds_ic_idx = mda_idx;
                    else if (scaler_name == "CFG_2")
                        cfg_2_idx = mda_idx;
                    else if (scaler_name == "CFG_3")
                        cfg_3_idx = mda_idx;
                    else if (scaler_name == "CFG_4")
                        cfg_4_idx = mda_idx;
                    else if (scaler_name == "CFG_5")
                        cfg_5_idx = mda_idx;
                }
            }
            for (auto itr : params_override->time_normalized_scalers)
            {
                //don't save ELT1, ERT1, ICR1, OCR1. these are saved elsewhere
                std::list<std::string>::iterator s_itr = std::find(ignore_scaler_strings.begin(), ignore_scaler_strings.end(), itr.first);
                if (s_itr != ignore_scaler_strings.end())
                    continue;

                std::string scaler_name = itr.first;
                std::transform(scaler_name.begin(), scaler_name.end(), scaler_name.begin(), ::toupper);


                int mda_idx = mda_io.find_scaler_index(mda_scalers, itr.second, val, units);
                if (mda_idx > -1)
                {
                    bool found_scaler = false;
                    for(auto& subitr : scalers)
                    {
                        if(subitr.hdf_name == itr.first)
                        {
                            subitr.mda_idx = mda_idx;
                            subitr.normalize_by_time = true;
                            subitr.hdf_units = units;
                            found_scaler = true;
                            break;
                        }
                    }
                    if(found_scaler == false)
                    {
                        scalers.push_back(scaler_struct(itr.first, units, mda_idx, hdf_idx, true));
                        hdf_idx++;
                    }
                    if (scaler_name == "US_IC")
                        us_ic_idx = mda_idx;
                    else if (scaler_name == "DS_IC")
                        ds_ic_idx = mda_idx;
                    else if (scaler_name == "CFG_2")
                        cfg_2_idx = mda_idx;
                    else if (scaler_name == "CFG_3")
                        cfg_3_idx = mda_idx;
                    else if (scaler_name == "CFG_4")
                        cfg_4_idx = mda_idx;
                    else if (scaler_name == "CFG_5")
                        cfg_5_idx = mda_idx;
                }
            }

            // Maps summed scaler name to scaler mda index
            for (auto& itr : params_override->summed_scalers)
            {
                bool found = false;
                for(auto &scaler_itr : itr.scalers_to_sum)
                {
                    for(auto &found_scalers_itr : scalers)
                    {
                        if(scaler_itr.first == found_scalers_itr.hdf_name && itr.normalize_by_time == found_scalers_itr.normalize_by_time)
                        {
                            scaler_itr.second = found_scalers_itr.mda_idx;
                            found = true;
                            break;
                        }
                    }
                }
                if(found)
                {
                    summed_scalers.push_back(itr);
                }
            }

            //search for time scaler index
            mda_time_scaler_idx = mda_io.find_scaler_index(mda_scalers, params_override->time_scaler, val, units);

            scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
            if (scalers_grp_id < 0)
                scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if (scalers_grp_id < 0)
            {
                logE << "creating group MAPS/Scalers" << "\n";
                return false;
            }

            _save_amps(scalers_grp_id, mda_scalers, params_override);

            if (scalers.size() > 0)
            {

                if (mda_scalers->header->data_rank == 2)
                {
                    if(hasNetcdf)
                    {
                        count_3d[0] = 1;
                        if(mda_scalers->scan->last_point == 0)
                            count_3d[1] = 1;
                        else
                            count_3d[1] = mda_scalers->scan->last_point;
                        if(mda_scalers->scan->sub_scans[0]->last_point == 0)
                            count_3d[2] = 1;
                        else
                            count_3d[2] = mda_scalers->scan->sub_scans[0]->last_point;
                    }
                    else
                    {
                        if(mda_scalers->header->dimensions[1] == 2000)
                        {
                            count_3d[0] = 1;
                            count_3d[1] = 1;
                            if(mda_scalers->scan->last_point == 0)
                                count_3d[2] = 1;
                            else
                                count_3d[2] = mda_scalers->scan->last_point;
                            single_row_scan = true;
                        }
                        else
                        {
                            logE<<"Unknown or bad mda file"<<"\n";
                        }
                    }
                }
                else if (mda_scalers->header->data_rank == 3)
                {
                    count_3d[0] = 1;
                    if(mda_scalers->scan->last_point == 0)
                        count_3d[1] = 1;
                    else
                        count_3d[1] = mda_scalers->scan->last_point;
                    if(mda_scalers->scan->sub_scans[0]->last_point == 0)
                        count_3d[2] = 1;
                    else
                        count_3d[2] = mda_scalers->scan->sub_scans[0]->last_point;
                }
                else
                {
                    logE << "Unsupported rank " << mda_scalers->header->data_rank << " . Skipping scalers" << "\n";
                    return false;
                }
                dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(dcpl_id, 3, count_3d);
                H5Pset_deflate(dcpl_id, 7);

                if (us_ic_idx > -1 && ds_ic_idx > -1 && cfg_2_idx > -1 && cfg_3_idx > -1 && cfg_4_idx > -1 && cfg_5_idx > -1)
                {
                    count_3d[0] = scalers.size() + summed_scalers.size() + 6; //abs_ic, abs_cfg, H_dpc_cfg, V_dpc_cfg, dia1_dpc_cfg, dia2_dpc_cfg
                    save_cfg_abs = true;
                }
                else
                {
                    count_3d[0] = scalers.size() + summed_scalers.size();
                }

                if(spectra_volume != nullptr)
                {
                    // if we have spectra volume loaded, save elt, ert, in_cnt, and out_cnt scalers
                    count_3d[0] += 4;
                }

                dataspace_id = H5Screate_simple(3, count_3d, NULL);
                filespace_id = H5Screate_simple(3, count_3d, NULL);

                count[0] = count_3d[0];
                filespace_name_id = H5Screate_simple(1, count, NULL);

                dset_cps_id = H5Dcreate(scalers_grp_id, "Values", H5T_INTEL_R, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
                dset_names_id = H5Dcreate(scalers_grp_id, "Names", filetype, filespace_name_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
				dset_units_id = H5Dcreate(scalers_grp_id, "Units", filetype, filespace_name_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                count_3d[0] = 1;
                count[0] = 1;

				//save scalers
				if (single_row_scan)
				{
					count_2d[0] = 1;
                    if(mda_scalers->scan->last_point == 0)
                        count_2d[1] = 1;
                    else
                        count_2d[1] = mda_scalers->scan->last_point;
				}
				else
				{
                    if(mda_scalers->scan->last_point == 0)
                        count_2d[0] = 1;
                    else
                        count_2d[0] = mda_scalers->scan->last_point;
                    if(mda_scalers->scan->sub_scans[0]->last_point == 0)
                        count_2d[1] = 1;
                    else
                        count_2d[1] = mda_scalers->scan->sub_scans[0]->last_point;
				}
				count_3d[1] = count_2d[0];
				count_3d[2] = count_2d[1];

				offset_3d[1] = 0;
				offset_3d[2] = 0;

				memoryspace_id = H5Screate_simple(2, count_2d, NULL);

				scaler_mat.resize(count_2d[0], count_2d[1]);
                scaler_mat.setZero(count_2d[0], count_2d[1]);

                for (auto &itr : scalers)
                {
                    scaler_mat.Zero(count_2d[0], count_2d[1]);
                    offset[0] = itr.hdf_idx;
                    char tmp_char[255] = {0};
					char tmp_char_units[255] = { 0 };
                    itr.hdf_name.copy(tmp_char, 254);
					itr.hdf_units.copy(tmp_char_units, 254);
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);
					status = H5Dwrite(dset_units_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char_units);

                    if(itr.mda_idx < 0)
                    {
                        continue;
                    }

                    if (single_row_scan)
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
							val = mda_scalers->scan->detectors_data[itr.mda_idx][i];
                            if(itr.normalize_by_time && mda_time_scaler_idx > -1)
                            {
                                real_t scaler_time_normalizer = 1.0;
                                real_t det_time = mda_scalers->scan->detectors_data[mda_time_scaler_idx][i];
                                scaler_time_normalizer = det_time / time_scaler_clock;
                                val /= scaler_time_normalizer;
                            }
							scaler_mat(0, i) = val;
                        }
                    }
                    else
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            for (int32_t j = 0; j < mda_scalers->scan->sub_scans[0]->last_point; j++)
                            {
								val = mda_scalers->scan->sub_scans[i]->detectors_data[itr.mda_idx][j];
                                if(itr.normalize_by_time && mda_time_scaler_idx > -1)
                                {
                                    real_t scaler_time_normalizer = 1.0;
                                    real_t det_time = mda_scalers->scan->sub_scans[i]->detectors_data[mda_time_scaler_idx][j];
                                    scaler_time_normalizer = det_time / time_scaler_clock;
									val /= scaler_time_normalizer;
                                }
								scaler_mat(i, j) = val;
                            }
                        }
                    }

                    offset_3d[0] = itr.hdf_idx;
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)scaler_mat.data());

                }

                for (auto &itr : summed_scalers)
                {
                    scaler_mat.Zero(count_2d[0], count_2d[1]);
                    // sum values before saving. If time normalized then divide by time val
                    offset[0] = hdf_idx;
                    char tmp_char[255] = {0};
                    //char tmp_char_units[255] = { 0 };
                    itr.scaler_name.copy(tmp_char, 254);
                    //itr.hdf_units.copy(tmp_char_units, 254);
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);
                    //status = H5Dwrite(dset_units_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char_units);

                    if (single_row_scan)
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            scaler_mat(0, i) = 0.0;
                            for(auto &scaler_itr : itr.scalers_to_sum)
                            {
                                if(scaler_itr.second < 0)
                                {
                                    continue;
                                }
                                val = mda_scalers->scan->detectors_data[scaler_itr.second][i];
                                if(itr.normalize_by_time && mda_time_scaler_idx > -1)
                                {
                                    real_t scaler_time_normalizer = 1.0;
                                    real_t det_time = mda_scalers->scan->detectors_data[mda_time_scaler_idx][i];
                                    scaler_time_normalizer = det_time / time_scaler_clock;
                                    val /= scaler_time_normalizer;
                                }
                                scaler_mat(0, i) += val;
                            }
                        }
                    }
                    else
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            for (int32_t j = 0; j < mda_scalers->scan->sub_scans[0]->last_point; j++)
                            {
                                scaler_mat(i, j) = 0.0;
                                for(auto &scaler_itr : itr.scalers_to_sum)
                                {
                                    if(scaler_itr.second < 0)
                                    {
                                        continue;
                                    }
                                    val = mda_scalers->scan->sub_scans[i]->detectors_data[scaler_itr.second][j];
                                    if(itr.normalize_by_time && mda_time_scaler_idx > -1)
                                    {
                                        real_t scaler_time_normalizer = 1.0;
                                        real_t det_time = mda_scalers->scan->sub_scans[i]->detectors_data[mda_time_scaler_idx][j];
                                        scaler_time_normalizer = det_time / time_scaler_clock;
                                        val /= scaler_time_normalizer;
                                    }
                                    scaler_mat(i, j) += val;
                                }
                            }
                        }
                    }

                    offset_3d[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)scaler_mat.data());
                    hdf_idx++;
                }

                if(spectra_volume != nullptr)
                {
                    Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> elt_map, ert_map, in_cnt_map, out_cnt_map;
                    spectra_volume->generate_scaler_maps(&elt_map, &ert_map, &in_cnt_map, &out_cnt_map);
                    char elt_char[255] = "Elapsed Live Time";
                    char ert_char[255] = "Elapsed Real Time";
                    char in_char[255] = "Input Counts";
                    char out_char[255] = "Output Counts";

                    offset[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)elt_char);
                    offset_3d[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)elt_map.data());
                    hdf_idx++;

                    offset[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)ert_char);
                    offset_3d[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)ert_map.data());
                    hdf_idx++;

                    offset[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)in_char);
                    offset_3d[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)in_cnt_map.data());
                    hdf_idx++;

                    offset[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)out_char);
                    offset_3d[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)out_cnt_map.data());
                    hdf_idx++;
                }

                if (save_cfg_abs)
                {
                    scaler_mat.Zero(count_2d[0], count_2d[1]);
					//save calculated names
                    std::string tmp_name = "abs_ic";
                    offset[0] = hdf_idx;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    char tmp_char[255] = {0};
                    tmp_name.copy(tmp_char, 254);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "abs_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 1;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "H_dpc_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 2;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "V_dpc_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 3;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "dia1_dpc_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 4;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

                    tmp_name = "dia2_dpc_cfg";
                    tmp_name.copy(tmp_char, 254);
                    offset[0] = hdf_idx + 5;
                    H5Sselect_hyperslab(filespace_name_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, filespace_name_id, H5P_DEFAULT, (void*)tmp_char);

					abs_cfg_mat.resize(count_2d[0], count_2d[1]);
					H_dpc_cfg_mat.resize(count_2d[0], count_2d[1]);
					V_dpc_cfg_mat.resize(count_2d[0], count_2d[1]);
					dia1_dpc_cfg_mat.resize(count_2d[0], count_2d[1]);
					dia2_dpc_cfg_mat.resize(count_2d[0], count_2d[1]);

                    abs_cfg_mat.setZero(count_2d[0], count_2d[1]);
                    H_dpc_cfg_mat.setZero(count_2d[0], count_2d[1]);
                    V_dpc_cfg_mat.setZero(count_2d[0], count_2d[1]);
                    dia1_dpc_cfg_mat.setZero(count_2d[0], count_2d[1]);
                    dia2_dpc_cfg_mat.setZero(count_2d[0], count_2d[1]);

                    if (single_row_scan)
                    {

                        for (int32_t j = 0; j < mda_scalers->scan->last_point; j++)
                        {
                            real_t us_ic = mda_scalers->scan->detectors_data[us_ic_idx][j];
                            real_t ds_ic = mda_scalers->scan->detectors_data[ds_ic_idx][j];
                            real_t t_2 = mda_scalers->scan->detectors_data[cfg_2_idx][j];
                            real_t t_3 = mda_scalers->scan->detectors_data[cfg_3_idx][j];
                            real_t t_4 = mda_scalers->scan->detectors_data[cfg_4_idx][j];
                            real_t t_5 = mda_scalers->scan->detectors_data[cfg_5_idx][j];

                            real_t t_abs = t_2 + t_3 + t_4 + t_5;
                            scaler_mat(0,j) = ds_ic / us_ic;
                            abs_cfg_mat(0,j) = t_abs / us_ic;
                            H_dpc_cfg_mat(0, j) = (t_2 - t_3 - t_4 + t_5) / t_abs;
                            V_dpc_cfg_mat(0, j) = (t_2 + t_3 - t_4 - t_5) / t_abs;
                            dia1_dpc_cfg_mat(0, j) = (t_2 - t_4) / t_abs;
                            dia2_dpc_cfg_mat(0, j) = (t_3 - t_5) / t_abs;
							/*
							if (itr.normalize_by_time)
							{
								real_t scaler_time_normalizer = 1.0;
								if (mda_time_scaler_idx > -1)
								{
                                    real_t det_time = mda_scalers->scan->sub_scans[i]->detectors_data[mda_time_scaler_idx][i];
									scaler_time_normalizer = det_time / time_scaler_clock;
								}
								val /= scaler_time_normalizer;
							}
                            */
                        }
                    }
                    else
                    {
                        for (int32_t i = 0; i < mda_scalers->scan->last_point; i++)
                        {
                            for (int32_t j = 0; j < mda_scalers->scan->sub_scans[0]->last_point; j++)
                            {
                                real_t us_ic = mda_scalers->scan->sub_scans[i]->detectors_data[us_ic_idx][j];
                                real_t ds_ic = mda_scalers->scan->sub_scans[i]->detectors_data[ds_ic_idx][j];
                                real_t t_2 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_2_idx][j];
                                real_t t_3 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_3_idx][j];
                                real_t t_4 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_4_idx][j];
                                real_t t_5 = mda_scalers->scan->sub_scans[i]->detectors_data[cfg_5_idx][j];

                                real_t t_abs = t_2 + t_3 + t_4 + t_5;
								scaler_mat(i, j) = ds_ic / us_ic;
								abs_cfg_mat(i, j) = t_abs / us_ic;
								H_dpc_cfg_mat(i, j) = (t_2 - t_3 - t_4 + t_5) / t_abs;
								V_dpc_cfg_mat(i, j) = (t_2 + t_3 - t_4 - t_5) / t_abs;
								dia1_dpc_cfg_mat(i, j) = (t_2 - t_4) / t_abs;
								dia2_dpc_cfg_mat(i, j) = (t_3 - t_5) / t_abs;
                            }
                        }
                    }
					offset_3d[0] = hdf_idx;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)scaler_mat.data());
					offset_3d[0] = hdf_idx + 1;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)abs_cfg_mat.data());
					offset_3d[0] = hdf_idx + 2;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)H_dpc_cfg_mat.data());
					offset_3d[0] = hdf_idx + 3;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)V_dpc_cfg_mat.data());
					offset_3d[0] = hdf_idx + 4;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)dia1_dpc_cfg_mat.data());
					offset_3d[0] = hdf_idx + 5;
					H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    status = H5Dwrite(dset_cps_id, H5T_NATIVE_REAL, memoryspace_id, filespace_id, H5P_DEFAULT, (void*)dia2_dpc_cfg_mat.data());
                }
            }
        }

        H5Tclose(filetype);
        H5Tclose(memtype);

        if(dset_names_id > -1)
        {
            H5Dclose(dset_names_id);
            dset_names_id = -1;
        }
		if (dset_units_id > -1)
		{
			H5Dclose(dset_units_id);
			dset_units_id = -1;
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
		if (memoryspace_str_id > -1)
		{
			H5Sclose(memoryspace_str_id);
			memoryspace_str_id = -1;
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
		
		if (memoryspace_str_id > -1)
			H5Sclose(memoryspace_str_id);
        if(memoryspace_id > -1)
            H5Sclose(memoryspace_id);
        if(dset_names_id > -1)
            H5Dclose(dset_names_id);
		if (dset_units_id > -1)
			H5Dclose(dset_units_id);
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

        logE << "creating group MAPS/Scalers" << "\n";
        return false;
    }
    logI << "Done" << "\n";
    return true;
}

//-----------------------------------------------------------------------------

void HDF5_IO::_save_amps(hid_t scalers_grp_id, struct mda_file *mda_scalers, data_struct::Params_Override * params_override)
{
    hid_t dataspace_id = -1, memoryspace_id = -1;
    hid_t dset_id = -1;
    hid_t status;
    hid_t filetype, memtype;

    char tmp_char[255] = {0};
    hsize_t offset[1] = { 0 };
    hsize_t count[1] = { 3 };
    io::file::MDA_IO mda_io;

    filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
    memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(memtype, 255);
	std::string units;
    real_t us_amp_sens_num_val = params_override->us_amp_sens_num;
    mda_io.find_scaler_index(mda_scalers, params_override->us_amp_sens_num_pv, us_amp_sens_num_val, units);
    real_t us_amp_sens_unit_val = params_override->us_amp_sens_unit;
    mda_io.find_scaler_index(mda_scalers, params_override->us_amp_sens_unit_pv, us_amp_sens_unit_val, units);

    real_t ds_amp_sens_num_val = params_override->ds_amp_sens_num;
    mda_io.find_scaler_index(mda_scalers, params_override->ds_amp_sens_num_pv, ds_amp_sens_num_val, units);
    real_t ds_amp_sens_unit_val = params_override->ds_amp_sens_unit;
    mda_io.find_scaler_index(mda_scalers, params_override->ds_amp_sens_unit_pv, ds_amp_sens_unit_val, units);

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
                                struct mda_file *mda_scalers,
                                data_struct::Spectra_Volume * spectra_volume,
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

	if (mda_scalers == nullptr)
    {
        logW << "mda_scalers == nullptr. Not returning from save_scan_scalers" << "\n";
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

    _save_scan_meta_data(scan_grp_id, mda_scalers, params_override);
	
    _save_extras(scan_grp_id, mda_scalers);
	
    _save_scalers(maps_grp_id, mda_scalers, spectra_volume, params_override, hasNetcdf);

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
    //hsize_t n_offset[1] = {0};
    //hsize_t n_count[1] = {1};
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
        H5Gclose(maps_grp_id);
        logE << "creating group MAPS/Scan" << "\n";
        return false;
    }

    scalers_grp_id = H5Gopen(maps_grp_id, "Scalers", H5P_DEFAULT);
    if (scalers_grp_id < 0)
        scalers_grp_id = H5Gcreate(maps_grp_id, "Scalers", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (scalers_grp_id < 0)
    {
        H5Gclose(scan_grp_id);
        H5Gclose(scalers_grp_id);
        H5Gclose(maps_grp_id);
        logE << "creating group MAPS/Scalers" << "\n";
        return false;
    }

    if ( false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1) )
        return false;

    if ( false == _open_h5_object(src_maps_grp_id, H5O_GROUP, close_map, "2D Scan", file_id) )
       return false;

    //_save_scan_meta_data(scan_grp_id, mda_scalers);
    if ( false == _open_h5_object(xpos_id, H5O_DATASET, close_map, "X Positions", src_maps_grp_id) )
       return false;
    xpos_dataspace_id = H5Dget_space(xpos_id);
    close_map.push({xpos_dataspace_id, H5O_DATASPACE});

    if ( false == _open_h5_object(ypos_id, H5O_DATASET, close_map, "Y Positions", src_maps_grp_id) )
       return false;
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
    //ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
    //status = H5Ocopy(src_maps_grp_id, "Detectors", scalers_grp_id, "Values", ocpypl_id, H5P_DEFAULT);
	hid_t scaler_id = H5Dopen(src_maps_grp_id, "Detectors", H5P_DEFAULT);
	if (scaler_id > -1)
	{
		hid_t scaler_space = H5Dget_space(scaler_id);
		H5Sget_simple_extent_dims(scaler_space, &scalers_count[0], NULL);
		value_count[0] = scalers_count[2];
		value_count[1] = scalers_count[0];
		value_count[2] = scalers_count[1];
		hid_t value_space = H5Screate_simple(3, &value_count[0], &value_count[0]);
		hid_t scalers_type = H5Dget_type(scaler_id);
		hid_t values_id = H5Dcreate(scalers_grp_id, "Values", scalers_type, value_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		real_t *buffer = new real_t[scalers_count[0] * scalers_count[1]];
		hsize_t scaler_amt = scalers_count[2];
		scalers_count[2] = 1;
		value_count[0] = 1;
		mem_count[0] = scalers_count[0];
		mem_count[1] = scalers_count[1];
		hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);

		for (hsize_t s = 0; s < scaler_amt; s++)
		{
			value_offset[0] = s;
			scalers_offset[2] = s;
			H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, scalers_offset, nullptr, scalers_count, nullptr);
			H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
			status = H5Dread(scaler_id, scalers_type, mem_space, scaler_space, H5P_DEFAULT, &buffer[0]);
			status = H5Dwrite(values_id, scalers_type, mem_space, value_space, H5P_DEFAULT, &buffer[0]);
		}
		delete[] buffer;

		H5Sclose(scaler_space);
		H5Sclose(mem_space);
		H5Sclose(value_space);
		H5Dclose(values_id);
		H5Dclose(scaler_id);
	}


    //save the scalers names
    if ( false == _open_h5_object(dset_detectors_id, H5O_DATASET, close_map, "Detectors", src_maps_grp_id) )
       return false;
    dataspace_detectors_id = H5Dget_space(dset_detectors_id);
    close_map.push({dataspace_detectors_id, H5O_DATASPACE});
    attr_detector_names_id=H5Aopen(dset_detectors_id, "Detector Names", H5P_DEFAULT);
    close_map.push({dataspace_detectors_id, H5O_ATTRIBUTE});

    det_rank = H5Sget_simple_extent_ndims(dataspace_detectors_id);
    det_dims_in = new hsize_t[det_rank];
    H5Sget_simple_extent_dims(dataspace_detectors_id, &det_dims_in[0], NULL);

    hid_t ftype = H5Aget_type(attr_detector_names_id);
    hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
    error = H5Aread(attr_detector_names_id, type, detector_names);

    if(error == 0)
    {
        dataspace_id = H5Screate_simple (1, &det_dims_in[2], NULL);
        dataset_id = H5Dcreate(scalers_grp_id, "Names", type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (dataset_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, detector_names);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        for(size_t z = 0; z < det_dims_in[2]; z++)
        {
			//TODO: look into why this is causing exception in windows
            //free(detector_names[z]);
        }
    }
    delete [] det_dims_in;


    _close_h5_objects(close_map);

    H5Gclose(scan_grp_id);
    H5Gclose(scalers_grp_id);
    H5Gclose(maps_grp_id);


    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::generate_avg(std::string avg_filename, std::vector<std::string> files_to_avg)
{
    std::lock_guard<std::mutex> lock(_mutex);
    logI  << avg_filename << "\n";

    hid_t ocpypl_id, status, src_maps_grp_id, src_analyzed_grp_id, dst_fit_grp_id;
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
        src_analyzed_grp_id = H5Gopen(src_maps_grp_id, "Quantification", H5P_DEFAULT);
        if(src_analyzed_grp_id > -1)
        {
            dst_fit_grp_id = H5Gcreate(dst_maps_grp_id, "Quantification", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

            _generate_avg_analysis(src_analyzed_grp_id, dst_fit_grp_id, "MAPS/Quantification", ocpypl_id, hdf5_file_ids);

            status = H5Ocopy(src_analyzed_grp_id, "Element_Weights", dst_fit_grp_id, "Element_Weights", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_analyzed_grp_id, "Element_Weights_Names", dst_fit_grp_id, "Element_Weights_Names", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_analyzed_grp_id, "Scalers", dst_fit_grp_id, "Scalers", ocpypl_id, H5P_DEFAULT);
            status = H5Ocopy(src_analyzed_grp_id, "Standard_Name", dst_fit_grp_id, "Standard_Name", ocpypl_id, H5P_DEFAULT);


            _generate_avg_integrated_spectra(src_analyzed_grp_id, dst_fit_grp_id, "MAPS/Quantification", ocpypl_id, hdf5_file_ids);

            H5Gclose(dst_fit_grp_id);
            H5Gclose(src_analyzed_grp_id);
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

        hsize_t* dims_in = new hsize_t[rank];
        int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
        if(status_n < 0)
        {
            logE<<"could not get dataset dimensions for "<<full_hdf5_path<< " " <<dataset_name<<"\n";
            return;
        }
        unsigned long total = 1;
        for(int i=0; i< rank; i++)
        {
            total *= dims_in[i];
        }
        //get all the other files dataset ids
        for(long unsigned int k=1; k<hdf5_file_ids.size(); k++)
        {
            hid_t det_analysis_dset_id = H5Dopen2(hdf5_file_ids[k], full_hdf5_path.c_str(), H5P_DEFAULT);
            if( det_analysis_dset_id > -1 )
                analysis_ids.push_back(det_analysis_dset_id);
        }

        data_struct::ArrayXr buffer1(total);
		data_struct::ArrayXr buffer2(total);
        float divisor = 1.0;
        error = H5Dread(dset_id, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer1.data());
        for(long unsigned int k=0; k<analysis_ids.size(); k++)
        {
            error = H5Dread(analysis_ids[k], H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer2.data());
            if(error > -1)
            {
                buffer1 += buffer2;
                divisor += 1.0;
            }
            else
            {
            logE<<"reading "<<full_hdf5_path<<" dataset "<<"\n";
            }
        }

        if(avg)
        {
            buffer1 /= divisor;
        }

        //hid_t dst_dset_id = H5Dopen2(dst_fit_grp_id, dataset_name.c_str(), H5P_DEFAULT);
        error = H5Dwrite(dst_dset_id, H5T_NATIVE_REAL, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer1.data());

        for(long unsigned int k=0; k<analysis_ids.size(); k++)
        {
            H5Dclose(analysis_ids[k]);
        }

        //clean up
        delete [] dims_in;
        H5Dclose(dset_id);
        H5Dclose(dst_dset_id);
    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_generate_avg_analysis(hid_t src_maps_grp_id, hid_t dst_maps_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids)
{
    hid_t src_analyzed_grp_id = H5Gopen(src_maps_grp_id, "XRF_Analyzed", H5P_DEFAULT);
    if(src_analyzed_grp_id > -1)
    {
        hid_t dst_analyzed_grp_id = H5Gcreate(dst_maps_grp_id, "XRF_Analyzed", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        for(std::string analysis_grp_name : xrf_analysis_save_names)
        {
            hid_t src_fit_grp_id = H5Gopen(src_analyzed_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT);
            if(src_fit_grp_id > -1)
            {
                //std::vector<hid_t> analysis_ids;

                hid_t dst_fit_grp_id = H5Gcreate(dst_analyzed_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                //copy channel names
                hid_t status = H5Ocopy(src_fit_grp_id, "Channel_Names", dst_fit_grp_id, "Channel_Names", ocpypl_id, H5P_DEFAULT);
                if(status < 0)
                {
                    logE<<"coping Channel_Names!\n";
                }
                status = H5Ocopy(src_fit_grp_id, "Channel_Units", dst_fit_grp_id, "Channel_Units", ocpypl_id, H5P_DEFAULT);
                if(status < 0)
                {
                    logE<<"coping Channel_Units!\n";
                }
                status = H5Ocopy(src_fit_grp_id, "Calibration_Curve_Labels", dst_fit_grp_id, "Calibration_Curve_Labels", ocpypl_id, H5P_DEFAULT);
                if(status < 0)
                {
                    logE<<"coping Calibration_Curve_Labels!\n";
                }
                _gen_average(group_name+"/XRF_Analyzed/"+analysis_grp_name+"/Counts_Per_Sec", "Counts_Per_Sec", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);

                _gen_average(group_name+"/XRF_Analyzed/"+analysis_grp_name+"/Calibration_Curve_SR_Current", "Calibration_Curve_SR_Current", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                _gen_average(group_name+"/XRF_Analyzed/"+analysis_grp_name+"/Calibration_Curve_DS_IC", "Calibration_Curve_DS_IC", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                _gen_average(group_name+"/XRF_Analyzed/"+analysis_grp_name+"/Calibration_Curve_US_IC", "Calibration_Curve_US_IC", src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);

                H5Gclose(dst_fit_grp_id);
                H5Gclose(src_fit_grp_id);

            }
        }
        H5Gclose(dst_analyzed_grp_id);
        H5Gclose(src_analyzed_grp_id);
    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_generate_avg_integrated_spectra(hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids)
{
    hid_t src_inner_grp_id = H5Gopen(src_analyzed_grp_id, "Integrated_Spectra", H5P_DEFAULT);
    if(src_inner_grp_id > -1)
    {
        hid_t dst_inner_grp_id = H5Gcreate(dst_fit_grp_id, "Integrated_Spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        _gen_average(group_name+"/Integrated_Spectra/Elapsed_Livetime", "Elapsed_Livetime", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Elapsed_Realtime", "Elapsed_Realtime", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Input_Counts", "Input_Counts", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);
        _gen_average(group_name+"/Integrated_Spectra/Output_Counts", "Output_Counts", src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids);

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

void HDF5_IO::add_v9_layout(std::string dataset_directory,
                            std::string dataset_file,
                            const std::vector<size_t>& detector_num_arr)
{
    for(size_t detector_num : detector_num_arr)
    {
        _add_v9_layout(dataset_directory+"img.dat"+ DIR_END_CHAR +dataset_file+".h5"+std::to_string(detector_num));
    }
    _add_v9_layout(dataset_directory+"img.dat"+ DIR_END_CHAR +dataset_file+".h5");
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
    std::string currnt_quant_str = "/MAPS/Quantification/XRF_Analyzed/" + quant_str + "/Calibration_Curve_SR_Current";
    std::string us_quant_str = "/MAPS/Quantification/XRF_Analyzed/" + quant_str + "/Calibration_Curve_US_IC";
    std::string ds_quant_str = "/MAPS/Quantification/XRF_Analyzed/" + quant_str + "/Calibration_Curve_DS_IC";
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
        
        hsize_t count_1d[1] = {1};
        hsize_t count_2d[2] = {1,1};
        hsize_t count_3d[3] = {1,1,1};
        hsize_t offset_1d[1] = {0};
        hsize_t offset_2d[2] = {0,0};
        hsize_t offset_3d[3] = {0,0,0};
        hid_t memoryspace_id = H5Screate_simple(1, count_1d, nullptr);
        char sr_current_carr[255] = "SRCURRENT";
        char us_ic_carr[255] = "US_IC";
        char ds_ic_carr[255] = "DS_IC";
        real_t real_val = 0.0;

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
            H5Dread(chan_names, memtype, memoryspace_id, chan_space, H5P_DEFAULT, (void*)tmp_char);
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
            if(element != nullptr)
            {
                offset_2d[1] = element->number - 1;
                H5Sselect_hyperslab(cc_space, H5S_SELECT_SET, offset_2d, nullptr, count_2d, nullptr);

                offset_3d[2] = chan_idx;
                offset_3d[0] = 0;
                H5Sselect_hyperslab(quant_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                H5Dread(cc_current, H5T_NATIVE_REAL, memoryspace_id, cc_space, H5P_DEFAULT, (void*)&real_val);
                H5Dwrite(quant_dset, H5T_NATIVE_REAL, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);
                offset_3d[0] = 1;
                H5Sselect_hyperslab(quant_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                H5Dread(cc_us_ic, H5T_NATIVE_REAL, memoryspace_id, cc_space, H5P_DEFAULT, (void*)&real_val);
                H5Dwrite(quant_dset, H5T_NATIVE_REAL, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);
                offset_3d[0] = 2;
                H5Sselect_hyperslab(quant_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                H5Dread(cc_ds_ic, H5T_NATIVE_REAL, memoryspace_id, cc_space, H5P_DEFAULT, (void*)&real_val);
                H5Dwrite(quant_dset, H5T_NATIVE_REAL, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);

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
                    for(int i=0; i < 3; i++)
                    {
                        for(int j=0; j < chan_amt; j++)
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
		if (quant_dset > -1)
		{
			H5Dclose(quant_dset);
		}
		if (chan_units > -1)
		{
			H5Dclose(chan_units);
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
    }
}

//-----------------------------------------------------------------------------

void HDF5_IO::_add_v9_layout(std::string dataset_file)
{
    std::lock_guard<std::mutex> lock(_mutex);

    logI  << dataset_file << "\n";
    hid_t saved_file_id = _cur_file_id;

    hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
    hid_t memtype = H5Tcopy(H5T_C_S1);
    H5Tset_size(memtype, 255);
    char tmp_char1[256] = { 0 };

    hid_t file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    //Scan
    if( H5Lcreate_hard(file_id, "/MAPS/Scan/x_axis", H5L_SAME_LOC, "/MAPS/x_axis", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for x_axis"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Scan/y_axis", H5L_SAME_LOC, "/MAPS/y_axis", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for y_axis"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Scan/scan_time_stamp", H5L_SAME_LOC, "/MAPS/scan_time_stamp", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for scan_time_stamp"<<  "\n";
    }
    //create extra_pvs, extra_pvs_as_csv, extra_strings

    //Scalers
    if( H5Lcreate_hard(file_id, "/MAPS/Scalers/ds_amp", H5L_SAME_LOC, "/MAPS/ds_amp", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for ds_amp"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Scalers/us_amp", H5L_SAME_LOC, "/MAPS/us_amp", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for us_amp"<<  "\n";
    }
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

    //Spectra
    if( H5Lcreate_hard(file_id, "/MAPS/Spectra/Energy", H5L_SAME_LOC, "/MAPS/energy", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for energy"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Spectra/Energy_Calibration", H5L_SAME_LOC, "/MAPS/energy_calib", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for energy_calib"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Spectra/Integrated_Spectra/Spectra", H5L_SAME_LOC, "/MAPS/int_spec", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for int_spec"<<  "\n";
    }
    if( H5Lcreate_hard(file_id, "/MAPS/Spectra/mca_arr", H5L_SAME_LOC, "/MAPS/mca_arr", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << " Couldn't create soft link for mca_arr"<<  "\n";
    }


    hsize_t quant_dims[3];
    quant_dims[0] = 3;
    quant_dims[1] = 1;
    quant_dims[2] = 1; //num channel names


    //XRF_Analyzed
    if( H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/ROI/Channel_Names", H5L_SAME_LOC, "/MAPS/channel_names", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        if( H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/Fitted/Channel_Names", H5L_SAME_LOC, "/MAPS/channel_names", H5P_DEFAULT, H5P_DEFAULT) < 0)
        {
            logW  << " Couldn't create soft link for channel_names"<<  "\n";
        }
    }

    hid_t chan_names = H5Dopen(file_id, "/MAPS/channel_names", H5P_DEFAULT);
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
    hsize_t unit_dims[2];
    hsize_t offset_dims[2] = {0,0};
    unit_dims[0] = 4;
    unit_dims[1] = quant_dims[2];
    hid_t units_space = H5Screate_simple(2, &unit_dims[0], &unit_dims[0]);
    hid_t ch_unit_id = H5Dcreate1(file_id, "/MAPS/channel_units", filetype, units_space, H5P_DEFAULT);
    if(ch_unit_id > 0)
    {
        hsize_t mem_dims[1] = {1};
        hsize_t count_2d[2] = {1,1};
        hid_t mem_space = H5Screate_simple(1, &mem_dims[0], &mem_dims[0]);

        std::string str_val = "cts/s";
        for (int z = 0; z < 256; z++)
        {
            tmp_char1[z] = 0;
        }
        str_val.copy(tmp_char1, 256);
        for(hsize_t i=0; i < unit_dims[0] ; i++)
        {
            for(hsize_t j=0; j < unit_dims[1]; j++)
            {
                offset_dims[0] = i;
                offset_dims[1] = j;
                H5Sselect_hyperslab (units_space, H5S_SELECT_SET, offset_dims, nullptr, count_2d, nullptr);
                H5Dwrite(ch_unit_id, filetype, mem_space, units_space, H5P_DEFAULT, (void*)&str_val[0]);
            }
        }
    }
    else
    {
        logW << "Couldn't create /MAPS/channel_units" << "\n";
    }

    if( H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/ROI/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_roi", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
		hid_t ana_id = H5Dopen(file_id, "/MAPS/XRF_roi", H5P_DEFAULT);
		if (ana_id > -1)
		{
			_add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], "ROI", "/MAPS/XRF_roi_quant");
			H5Dclose(ana_id);
		}
		else
		{
			logW << "Couldn't create soft link for XRF_roi" << "\n";
		}
    }
    else
    {
        _add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], "ROI", "/MAPS/XRF_roi_quant");
    }
    if( H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_roi_plus", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
		hid_t ana_id = H5Dopen(file_id, "/MAPS/XRF_roi_plus", H5P_DEFAULT);
		if (ana_id > -1)
		{
			_add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], "ROI", "/MAPS/XRF_roi_quant");
			H5Dclose(ana_id);
		}
        logW  << " Couldn't create soft link for XRF_roi_plus"<<  "\n";
    }
    else
    {
        _add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], "NNLS", "/MAPS/XRF_roi_plus_quant");
    }
    if( H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_fits", H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
		hid_t ana_id = H5Dopen(file_id, "/MAPS/XRF_fits", H5P_DEFAULT);
		if (ana_id > -1)
		{
			_add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], "ROI", "/MAPS/XRF_roi_quant");
			H5Dclose(ana_id);
		}
        logW  << " Couldn't create soft link for XRF_fits"<<  "\n";
    }
    else
    {
        _add_v9_quant(file_id, quant_space, chan_names, chan_space, quant_dims[2], "Fitted", "/MAPS/XRF_fits_quant");
    }

    _add_extra_pvs(file_id, "/MAPS");

    //change version to 9
    real_t version = 9;
    hid_t version_id = H5Dopen(file_id, "/MAPS/version", H5P_DEFAULT);
    hid_t ver_space = H5Dget_space(version_id);
    hid_t ver_type = H5Dget_type(version_id);
    H5Dwrite(version_id, ver_type, ver_space, ver_space, H5P_DEFAULT, (void*)&version);
    H5Dclose(version_id);

    _cur_file_id = file_id;
    end_save_seq();
    logI<<"closing file"<<"\n";

    _cur_file_id = saved_file_id;

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
        hid_t ds_ic_quant_id = H5Dopen(file_id, "/MAPS/Quantification/XRF_Analyzed/Fitted/Calibration_Curve_DS_IC", H5P_DEFAULT);
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

void HDF5_IO::_add_exchange_layout(std::string dataset_file)
{
    std::lock_guard<std::mutex> lock(_mutex);

    logI  << dataset_file << "\n";
    hid_t saved_file_id = _cur_file_id;
    hid_t file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    std::string exchange_fits[4] = {"Fitted", "ROI", "NNLS", "Fitted"};
    std::string exchange_normalize[4] = {"DS_IC", "", "", ""};

    int ex_idx = 0;
    for(int i=0; i < 4; i++)
    {
        if(true == _add_exchange_meta(file_id, std::to_string(ex_idx), exchange_fits[i], exchange_normalize[i]))
        {
            ex_idx++;
        }
    }

    _cur_file_id = file_id;
    end_save_seq();
    logI<<"closing file"<<"\n";

    _cur_file_id = saved_file_id;
}

//-----------------------------------------------------------------------------

void HDF5_IO::add_exchange_layout(std::string dataset_directory,
								std::string dataset_file,
								const std::vector<size_t>& detector_num_arr)
{
	
	for (size_t detector_num : detector_num_arr)
	{
		_add_exchange_layout(dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file + ".h5" + std::to_string(detector_num));
	}
	_add_exchange_layout(dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file + ".h5");

}

//-----------------------------------------------------------------------------

} //end namespace file
}// end namespace io
