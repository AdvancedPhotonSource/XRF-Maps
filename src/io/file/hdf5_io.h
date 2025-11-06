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


#ifndef HDF5_IO_H
#define HDF5_IO_H

#include <codecvt>
#include <cstdint>
#include <locale>
#include <list>
#include <mutex>
#include <queue>
#include <future>
#include <map>
#include <stack>
#include <type_traits>
#include "hdf5.h"
#include "data_struct/spectra_volume.h"
#include "data_struct/fit_element_map.h"
#include "data_struct/detector.h"
#include "data_struct/params_override.h"
#include "data_struct/scan_info.h"

#include "core/mem_info.h"
#include "core/defines.h"

#include "data_struct/element_info.h"
#include "data_struct/scaler_lookup.h"

#include "csv_io.h"
namespace io
{
namespace file
{

#define HDF5_SAVE_VERSION 10.0

#define HDF5_EXCHANGE_VERSION 1.0

enum H5_OBJECTS{H5O_FILE, H5O_GROUP, H5O_DATASPACE, H5O_DATASET, H5O_ATTRIBUTE, H5O_PROPERTY, H5O_DATATYPE};

enum H5_SPECTRA_LAYOUTS {MAPS_RAW, MAPS_V9, MAPS_V10, XSPRESS, APS_SEC20};

enum GSE_CARS_SAVE_VER {UNKNOWN, XRFMAP, XRMMAP};

using ROI_Vec = std::vector<std::pair<int, int>>;

DLL_EXPORT herr_t h5_ext_file_info(hid_t loc_id, const char *name, const H5L_info2_t *linfo, void *opdata);

//-----------------------------------------------------------------------------

template<typename T_real>
int parse_str_val_to_int(std::string start_delim, std::string end_delim, std::string lookup_str)
{
    size_t find_idx = (int)lookup_str.find(start_delim);
    if (find_idx == std::string::npos)
    {
        logW << "Could not find string property start delimeter " << start_delim << "\n";
        return -1;
    }
    std::string leftover = lookup_str.substr(find_idx + start_delim.length());
    find_idx = leftover.find(end_delim); // find end double quote
    if (find_idx == std::string::npos)
    {
        logW << "Could not find spectra sample count end delimeter \"\n";
        return -1;
    }
    std::string str_val = leftover.substr(0, find_idx);
    return std::stoi(str_val.c_str());
}

//-----------------------------------------------------------------------------

template<typename T_real>
T_real translate_back_sens_num(int value)
{
    if (value == 1)
    {
        return T_real(0.);
    }
    else if (value == 2)
    {
        return T_real(1.);
    }
    else if (value == 5)
    {
        return T_real(2.);
    }
    else if (value == 10)
    {
        return T_real(3.);
    }
    else if (value == 20)
    {
        return T_real(4.);
    }
    else if (value == 50)
    {
        return T_real(5.);
    }
    else if (value == 100)
    {
        return T_real(6.);
    }
    else if (value == 200)
    {
        return T_real(7.);
    }
    else if (value == 500)
    {
        return T_real(8.);
    }
    return T_real(-1.);
}

//-----------------------------------------------------------------------------

template<typename T_real>
T_real translate_back_sens_unit(std::string value)
{
    if (value == "pA/V")
    {
        return T_real(0);
    }
    else if (value == "nA/V")
    {
        return T_real(1);
    }
    else if (value == "uA/V")
    {
        return T_real(2);
    }
    else if (value == "mA/V")
    {
        return T_real(3);
    }
    return T_real (-1);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<typename T_real>
struct Detector_HDF5_Struct
{
    hid_t    dset_id;
    hid_t    dataspace_id;
    T_real* buffer;
};

class DLL_EXPORT HDF5_IO
{
public:

    static HDF5_IO* inst();

    ~HDF5_IO();

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_spectra_volume(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        //_is_loaded = ERROR_LOADING;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        logI << path << " detector : " << detector_num << "\n";

        hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
        hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
        herr_t   error;
        std::string detector_path;
        T_real* buffer;
        hsize_t offset_row[2] = { 0,0 };
        hsize_t count_row[2] = { 0,0 };
        hsize_t offset_meta[3] = { 0,0,0 };
        hsize_t count_meta[3] = { 1,1,1 };


        switch (detector_num)
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

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
            return false;

        if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id))
            return false;

        if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id))
            return false;
        dataspace_id = H5Dget_space(dset_id);
        close_map.push({ dataspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id))
            return false;
        dataspace_lt_id = H5Dget_space(dset_lt_id);
        close_map.push({ dataspace_lt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id))
            return false;
        dataspace_rt_id = H5Dget_space(dset_rt_id);
        close_map.push({ dataspace_rt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id))
            return false;
        dataspace_inct_id = H5Dget_space(dset_incnt_id);
        close_map.push({ dataspace_inct_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id))
            return false;
        dataspace_outct_id = H5Dget_space(dset_outcnt_id);
        close_map.push({ dataspace_outct_id, H5O_DATASPACE });

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

        buffer = new T_real[dims_in[0] * dims_in[2]]; // spectra_size x cols
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
        close_map.push({ memoryspace_id, H5O_DATASPACE });
        memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
        close_map.push({ memoryspace_meta_id, H5O_DATASPACE });
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
        H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;

        offset_meta[0] = detector_num;
        for (size_t row = 0; row < spec_vol->rows(); row++)
        {
            offset[1] = row;
            offset_meta[1] = row;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

            if (error > -1)
            {
                for (size_t col = 0; col < spec_vol->cols(); col++)
                {
                    offset_meta[2] = col;
                    data_struct::Spectra<T_real>* spectra = &((*spec_vol)[row][col]);

                    H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_lt_id, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                    spectra->elapsed_livetime(live_time);

                    H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_rt_id, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                    spectra->elapsed_realtime(real_time);

                    H5Sselect_hyperslab(dataspace_inct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_incnt_id, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                    spectra->input_counts(in_cnt);

                    H5Sselect_hyperslab(dataspace_outct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_outcnt_id, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                    spectra->output_counts(out_cnt);

                    spectra->recalc_elapsed_livetime();

                    for (size_t s = 0; s < count_row[0]; s++)
                    {
                        (*spectra)[s] = buffer[(count_row[1] * s) + col];
                    }
                    //logD<<"saved col "<<col<<"\n";
                }

                //logD<<"read row "<<row<<"\n";
            }
            else
            {
                logE << "reading row " << row << "\n";
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

    template<typename T_real>
    bool load_spectra_volume_with_callback(std::string path, const std::vector<size_t>& detector_num_arr, data_struct::IO_Callback_Func_Def<T_real> callback_func, void* user_data)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        //_is_loaded = ERROR_LOADING;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        std::map<size_t, struct Detector_HDF5_Struct<T_real>> detector_hid_map;

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        hid_t    file_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
        hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
        herr_t   error = -1;
        std::string detector_path;
        hsize_t offset_row[2] = { 0,0 };
        hsize_t count_row[2] = { 0,0 };
        hsize_t offset_meta[3] = { 0,0,0 };
        hsize_t count_meta[3] = { 1,1,1 };

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
            return false;

        if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id))
            return false;
        for (size_t detector_num : detector_num_arr)
        {
            detector_hid_map.insert({ detector_num, Detector_HDF5_Struct<T_real>() });

            switch (detector_num)
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


            if (false == _open_h5_object(detector_hid_map[detector_num].dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id))
                return false;
            detector_hid_map[detector_num].dataspace_id = H5Dget_space(detector_hid_map[detector_num].dset_id);
            close_map.push({ detector_hid_map[detector_num].dataspace_id, H5O_DATASPACE });
        }

        if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id))
            return false;
        dataspace_lt_id = H5Dget_space(dset_lt_id);
        close_map.push({ dataspace_lt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id))
            return false;
        dataspace_rt_id = H5Dget_space(dset_rt_id);
        close_map.push({ dataspace_rt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id))
            return false;
        dataspace_inct_id = H5Dget_space(dset_incnt_id);
        close_map.push({ dataspace_inct_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id))
            return false;
        dataspace_outct_id = H5Dget_space(dset_outcnt_id);
        close_map.push({ dataspace_outct_id, H5O_DATASPACE });


        int rank = H5Sget_simple_extent_ndims(detector_hid_map[detector_num_arr[0]].dataspace_id);
        if (rank != 3)
        {
            _close_h5_objects(close_map);
            logE << "Dataset /MAPS_RAW/" << detector_path << " rank != 3. rank = " << rank << ". Can't load dataset. returning" << "\n";
            return false;
            //throw exception ("Dataset is not a volume");
        }
        hsize_t* dims_in = new hsize_t[rank];
        hsize_t* offset = new hsize_t[rank];
        hsize_t* count = new hsize_t[rank];

        int status_n = H5Sget_simple_extent_dims(detector_hid_map[detector_num_arr[0]].dataspace_id, &dims_in[0], nullptr);
        if (status_n < 0)
        {
            _close_h5_objects(close_map);
            logE << "getting dataset rank for MAPS_RAW/" << detector_path << "\n";
            return false;
        }

        for (int i = 0; i < rank; i++)
        {

            offset[i] = 0;
            count[i] = dims_in[i];
        }

        for (size_t detector_num : detector_num_arr)
        {
            detector_hid_map[detector_num].buffer = new T_real[dims_in[0] * dims_in[2]]; // spectra_size x cols
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
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
        H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;

        for (size_t row = 0; row < dims_in[1]; row++)
        {
            offset[1] = row;
            offset_meta[1] = row;

            for (size_t detector_num : detector_num_arr)
            {
                H5Sselect_hyperslab(detector_hid_map[detector_num].dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                error = _read_h5d<T_real>(detector_hid_map[detector_num].dset_id, memoryspace_id, detector_hid_map[detector_num].dataspace_id, H5P_DEFAULT, detector_hid_map[detector_num].buffer);
            }

            if (error > -1)
            {
                for (size_t col = 0; col < count_row[1]; col++)
                {
                    offset_meta[2] = col;

                    for (size_t detector_num : detector_num_arr)
                    {
                        offset_meta[0] = detector_num;
                        data_struct::Spectra<T_real>* spectra = new data_struct::Spectra<T_real>(dims_in[0]);

                        H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_lt_id, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                        spectra->elapsed_livetime(live_time);

                        H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_rt_id, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                        spectra->elapsed_realtime(real_time);

                        H5Sselect_hyperslab(dataspace_inct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_incnt_id, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                        spectra->input_counts(in_cnt);

                        H5Sselect_hyperslab(dataspace_outct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_outcnt_id, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                        spectra->output_counts(out_cnt);

                        for (size_t s = 0; s < count_row[0]; s++)
                        {
                            (*spectra)[s] = detector_hid_map[detector_num].buffer[(count_row[1] * s) + col];
                        }
                        callback_func(row, col, dims_in[1], count_row[1], detector_num, spectra, user_data);
                    }
                }

            }
            else
            {
                logE << "reading row " << row << "\n";
            }
        }

        delete[] dims_in;
        delete[] offset;
        delete[] count;

        for (size_t detector_num : detector_num_arr)
        {
            delete[] detector_hid_map[detector_num].buffer;
        }

        detector_hid_map.clear();

        _close_h5_objects(close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
	bool load_spectra_volume_emd_with_callback(std::string path, const std::vector<size_t>& detector_num_arr, data_struct::IO_Callback_Func_Def<T_real> callback_func, void* user_data)
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
        hsize_t image_dims[3] = { 0,0,0 };
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
        if (error > -1)
        {
            width = image_dims[0];
            height = image_dims[1];
        }
        else
        {
            logE << "Could not get image dimensions from /Data/Image/" << str_grp_name << "/Data\n ";
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

        if (frame > 0)
        {
            if (false == _open_h5_object(frame_id, H5O_DATASET, close_map, "FrameLocationTable", spectra_grp_id))
                return false;
            dataspace_frame_id = H5Dget_space(frame_id);
            close_map.push({ dataspace_frame_id, H5O_DATASPACE });

            hsize_t frame_dims[2] = { 0,0 };
            error = H5Sget_simple_extent_dims(dataspace_frame_id, &frame_dims[0], nullptr);

            if (frame < frame_dims[0])
            {
                int* frames = new int[frame_dims[0]];
                //read frames
                hid_t frame_memoryspace_id = H5Screate_simple(2, frame_dims, nullptr);
                error = H5Dread(frame_id, H5T_NATIVE_INT, frame_memoryspace_id, dataspace_frame_id, H5P_DEFAULT, &frames[0]);
                start_offset = frames[frame];
                delete[] frames;
                close_map.push({ frame_memoryspace_id, H5O_DATASPACE });
            }
            else
            {
                logW << " detector_num_arr[0] " << frame << " > Max Frames" << frame_dims[0] << ". Setting frame_num_start = 0\n";
            }

        }
        if (false == _open_h5_object(meta_id, H5O_DATASET, close_map, "Metadata", spectra_grp_id))
            return false;
        dataspace_meta_id = H5Dget_space(meta_id);
        close_map.push({ dataspace_meta_id, H5O_DATASPACE });

        char* acquisition[1];
        hid_t ftype = H5Dget_type(acqui_id);
        close_map.push({ ftype, H5O_DATATYPE });
        hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
        error = H5Dread(acqui_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, acquisition);

        //prase bincount/samples width and height
        std::string str_acqui = std::string(acquisition[0]);
        int samples = parse_str_val_to_int<T_real>(STR_BINCOUNT, "\"", str_acqui);
        //int height = parse_str_val_to_int(STR_HEIGHT, "\"", str_acqui);; //wrong place, this data is incorrect
        //int width = parse_str_val_to_int(STR_WIDTH, "\"", str_acqui);; // need to read image size in from /Data/Image/...

        if (samples < 0 || height < 0 || width < 0)
        {
            logE << "Unknown spectra volume size Width: " << width << " Height: " << height << " Samples: " << samples << "\n";
            _close_h5_objects(close_map);
            return false;
        }


        hsize_t rank;
        rank = H5Sget_simple_extent_ndims(dataspace_id);
        unsigned short* buffer;
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

        if (start_offset > count[0])
        {
            logW << " frame start offset: " << start_offset << " > dataset count: " << count[0] << ". Setting start offset = 0\n.";
            start_offset = 0;
        }

        memoryspace_id = H5Screate_simple(2, chunk_dims, nullptr);
        close_map.push({ memoryspace_id, H5O_DATASPACE });
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);

        size_t read_leftover = (count[0] - start_offset) % chunk_dims[0];
        size_t read_amt = (count[0] - start_offset) / chunk_dims[0];

        int row = 0;
        int col = 0;

        T_real incnt = 0.0;
        T_real outcnt = 0.0;
        T_real elt = 1.0;
        offset[0] = start_offset; // start of frame offset
        data_struct::Spectra<T_real>* spectra = new data_struct::Spectra<T_real>(samples);
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
                        if (col >= width)
                        {
                            col = 0;
                            row++;
                            logI << "Reading row " << row << " of " << height << "\n";
                        }
                        if (row >= height)
                        {
                            col = 0;
                            row = 0;
                            frame_idx++;
                            frame = detector_num_arr[frame_idx];
                            if (frame_idx > detector_num_arr.size())
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
                        spectra = new data_struct::Spectra<T_real>(samples);
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
                logE << "reading chunk " << i << " of " << read_amt << "\n";
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
                    if (col >= width)
                    {
                        col = 0;
                        row++;
                        logI << "Reading row " << row << " of " << height << "\n";
                    }
                    if (row >= height)
                    {
                        col = 0;
                        row = 0;
                        frame++;
                    }
                    spectra = new data_struct::Spectra<T_real>(samples);
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

    template<typename T_real>
    bool load_spectra_vol_apsu(std::string path, std::string filename, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol, [[maybe_unused]] data_struct::Scan_Info<T_real> &scan_info, [[maybe_unused]] bool logerr = true)
    {
        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
        hid_t    file_id, xspres_grp_id;
        {
            std::lock_guard<std::mutex> lock(_mutex);
            if (false == _open_h5_object(file_id, H5O_FILE, close_map, path+filename, -1))
            {
                return false;
            }
            /*
            if (false == _open_h5_object(pos_grp_id, H5O_GROUP, close_map, "/positions", file_id))
            {
                return false;
            }
            if (false == _open_h5_object(xspres_grp_id, H5O_GROUP, close_map, "/flyXRF", file_id))
            */
            if (false == _open_h5_object(xspres_grp_id, H5O_GROUP, close_map, "/detectors/XRF_ME7", file_id, false, false))
            {
                if (false == _open_h5_object(xspres_grp_id, H5O_GROUP, close_map, "/detectors/XRF_RS", file_id))
                {
                    return false;
                }
                return false;
            }
        
            size_t idx = filename.find_last_of("_master.h5") - 9;
            if(idx != std::string::npos)
            {
                std::string base_name = filename.substr(0, idx);
                std::map<std::string, hid_t> ext_links;
                H5Literate2(xspres_grp_id, H5_INDEX_NAME, H5_ITER_INC, NULL, h5_ext_file_info, &ext_links);
                spec_vol->resize_and_zero(ext_links.size(),1,4096);
                for(unsigned long i = 0; i<ext_links.size(); i++)
                {
                    std::string search_name = base_name + "_" + std::to_string(i) + ".hdf5"; 
                    logI<< "searching "<< search_name << "\n";
                    if(ext_links.count(search_name) > 0)
                    {
                        logI<< "loading "<< search_name << "\n";
                        _load_spectra_line_xspress3(ext_links.at(search_name), detector_num, &(*spec_vol)[i]);
                        H5Gclose(ext_links.at(search_name));
                    }
                }
            }
/*
            // load spectra link
            err = H5Lget_name_by_idx(xspres_grp_id, ".", H5_INDEX_NAME, H5_ITER_NATIVE, 0, &xspress3_link[0], 2048, H5P_DEFAULT);
            if (err < 0) 
            {
                logE << "Error retrieving flyXRF link value\n";
                _close_h5_objects(close_map);
                return false;
            }
            err = H5Lget_name_by_idx(pos_grp_id, ".", H5_INDEX_NAME, H5_ITER_NATIVE, 0, &positions_link[0], 2048, H5P_DEFAULT);
            if (err < 0) 
            {
                logE << "Error retrieving positions link name\n";
                _close_h5_objects(close_map);
                return false;
            }

            // load positions
            if (false == _open_h5_object(pos_file_id, H5O_FILE, close_map, path + "positions/" + positions_link, -1))
            {
                return false;
            }
            if (false == _open_h5_object(enc1_id, H5O_DATASET, close_map, "/stream/Encoder1", pos_file_id, -1))
            {
                return false;
            }
            if (false == _open_h5_object(enc2_id, H5O_DATASET, close_map, "/stream/Encoder2", pos_file_id, -1))
            {
                return false;
            }
            space_id = H5Dget_space(enc1_id);
            close_map.push({ space_id, H5O_DATASPACE });

            err = H5Sget_simple_extent_dims(space_id, &dims[0], nullptr);
            if (err < 0)
            {
                logE<< "Could not read encoder dims\n";
                _close_h5_objects(close_map);
                return false;
            }
            scan_info.meta_info.x_axis.resize(dims[0]);
            scan_info.meta_info.y_axis.resize(dims[0]);
            err = H5Dread(enc1_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, scan_info.meta_info.x_axis.data());
            if (err < 0)
            {
                logE<< "Could not read encoder 1 values\n";
                _close_h5_objects(close_map);
                return false;
            }

            err = H5Dread(enc2_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, scan_info.meta_info.y_axis.data());
            if (err < 0)
            {
                logE<< "Could not read encoder 2 values\n";
                _close_h5_objects(close_map);
                return false;
            }

            // read in scalers from TertraMM

*/
            
            _close_h5_objects(close_map);
        }
        /*
        spec_vol->resize_and_zero(1,1,1);
        if(false == load_spectra_line_xspress3(path + "flyXRF/" + xspress3_link, detector_num, &(*spec_vol)[0]) )
        if(false == load_spectra_line_xspress3(path + "XRF_ME7/" + xspress3_link, detector_num, &(*spec_vol)[0]) )
        {
            return false;
        }
*/
        return true;
    }

    //-----------------------------------------------------------------------------
    // load spectra volume as float (T_real) but scan_info as double (T_real2)
    template<typename T_real, typename T_real2>
    bool load_spectra_vol_polar_energy_scan(std::string path, std::string filename, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol, data_struct::Scan_Info<T_real2> &scan_info, data_struct::Params_Override<T_real>* params_override,  bool logerr = true)
    {
        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
        hid_t    file_id, dset_id, space_id;
        hid_t    pr1_order_id, pr2_order_id;
        hid_t    ert_dset_id, ert_space_id;
        hid_t    elt_dset_id, elt_space_id;
        hid_t    icr_dset_id, icr_space_id;
        hid_t    ocr_dset_id, ocr_space_id;
        hid_t    energy_dset_id, energy_space_id;
        hid_t    i0_dset_id, i0_space_id;
        hid_t    title_id, start_time_id;
        hsize_t dims3[3] = { 0, 0, 0 };
        data_struct::ArrayTr<T_real> pr1_array;
        data_struct::ArrayTr<T_real> pr2_array;
        hsize_t offset3[3] = { 0, 0, 0 };
        hsize_t count3[3] = { 1, 1, 1 };

        bool use_pr1 = true;
        hsize_t dims1[1] = { 0 };
        hsize_t count[1] = { 1 };

        herr_t err;
        std::string str_elt_path = "/entry/externals/vortex/NDAttributes/LiveTime_"+std::to_string(detector_num);
        std::string str_ert_path = "/entry/externals/vortex/NDAttributes/RealTime_"+std::to_string(detector_num);
        std::string str_icr_path = "/entry/externals/vortex/NDAttributes/ICR_"+std::to_string(detector_num);
        std::string str_ocr_path = "/entry/externals/vortex/NDAttributes/OCR_"+std::to_string(detector_num);

        std::string str_energy_path = "/entry/data/energy";
        std::string str_i0_path = "/entry/data/4idgI0";

        std::string str_title_path = "/entry/title";
        std::string str_start_time_path = "/entry/start_time";

        std::string str_pr1_path = "/entry/data/pr1_pzt_localdc";
        std::string str_pr2_path = "/entry/data/pr2_pzt_localdc";

        offset3[1] = detector_num;
        {
            std::lock_guard<std::mutex> lock(_mutex);
            if (false == _open_h5_object(file_id, H5O_FILE, close_map, path+filename, -1))
            {
                return false;
            }
            
            if(false == _open_h5_object(dset_id, H5O_DATASET, close_map, "/entry/externals/vortex/detector/data", file_id, false))
            {
                logW << "Tried to open /entry/externals/vortex/detector/data but failed.\n";
                _close_h5_objects(close_map);
                return false;
            }

            if(false == _open_h5_object(elt_dset_id, H5O_DATASET, close_map, str_elt_path, file_id, false))
            {
                logW << "Tried to open "<<str_elt_path<<" but failed.\n";
                _close_h5_objects(close_map);
                return false;
            }
            if(false == _open_h5_object(ert_dset_id, H5O_DATASET, close_map, str_ert_path, file_id, false))
            {
                logW << "Tried to open "<<str_ert_path<<" but failed.\n";
                _close_h5_objects(close_map);
                return false;
            }
            if(false == _open_h5_object(icr_dset_id, H5O_DATASET, close_map, str_icr_path, file_id, false))
            {
                logW << "Tried to open "<<str_icr_path<<" but failed.\n";
                _close_h5_objects(close_map);
                return false;
            }
            if(false == _open_h5_object(ocr_dset_id, H5O_DATASET, close_map, str_ocr_path, file_id, false))
            {
                logW << "Tried to open "<<str_ocr_path<<" but failed.\n";
                _close_h5_objects(close_map);
                return false;
            }
            if(false == _open_h5_object(energy_dset_id, H5O_DATASET, close_map, str_energy_path, file_id, false))
            {
                logW << "Tried to open "<<str_energy_path<<" but failed.\n";
                _close_h5_objects(close_map);
                return false;
            }
            if(false == _open_h5_object(i0_dset_id, H5O_DATASET, close_map, str_i0_path, file_id, false))
            {
                logW << "Tried to open "<<str_i0_path<<" but failed.\n";
                _close_h5_objects(close_map);
                return false;
            }

            if(_open_h5_object(pr1_order_id, H5O_DATASET, close_map, str_pr1_path, file_id, false, false))
            {
                hid_t pr_space_id = H5Dget_space(pr1_order_id);
                close_map.push({ pr_space_id, H5O_DATASPACE });

                err = H5Sget_simple_extent_dims(pr_space_id, &dims1[0], nullptr);
                if (err >= 0)
                {
                    pr1_array.resize(dims1[0]);
                    err = _read_h5d<T_real>(pr1_order_id, H5S_ALL, pr_space_id, H5P_DEFAULT, pr1_array.data());
                    if(err < 0)
                    {
                        pr1_array.setOnes();
                        use_pr1 = false;
                    }
                    else
                    {
                        if(pr1_array.size() > 2)
                        {
                            if(pr1_array[0] != pr1_array[1])
                            {
                                use_pr1 = true;
                            }
                        }
                    }
                }
            }

            if(_open_h5_object(pr2_order_id, H5O_DATASET, close_map, str_pr2_path, file_id, false, false))
            {
                hid_t pr_space_id = H5Dget_space(pr2_order_id);
                close_map.push({ pr_space_id, H5O_DATASPACE });

                err = H5Sget_simple_extent_dims(pr_space_id, &dims1[0], nullptr);
                if (err >= 0)
                {
                    pr2_array.resize(dims1[0]);
                    err = _read_h5d<T_real>(pr2_order_id, H5S_ALL, pr_space_id, H5P_DEFAULT, pr2_array.data());
                    if(err < 0)
                    {
                        pr2_array.setOnes();
                    }
                    else
                    {
                        if(pr2_array.size() > 2)
                        {
                            if(pr2_array[0] != pr2_array[1])
                            {
                                use_pr1 = false;
                            }
                        }
                    }
                }
            }


            space_id = H5Dget_space(dset_id);
            close_map.push({ space_id, H5O_DATASPACE });

            elt_space_id = H5Dget_space(elt_dset_id);
            close_map.push({ elt_space_id, H5O_DATASPACE });

            ert_space_id = H5Dget_space(ert_dset_id);
            close_map.push({ ert_space_id, H5O_DATASPACE });

            icr_space_id = H5Dget_space(icr_dset_id);
            close_map.push({ icr_space_id, H5O_DATASPACE });

            ocr_space_id = H5Dget_space(ocr_dset_id);
            close_map.push({ ocr_space_id, H5O_DATASPACE });

            energy_space_id = H5Dget_space(energy_dset_id);
            close_map.push({ energy_space_id, H5O_DATASPACE });

            i0_space_id = H5Dget_space(i0_dset_id);
            close_map.push({ i0_space_id, H5O_DATASPACE });

            err = H5Sget_simple_extent_dims(space_id, &dims3[0], nullptr);
            if (err < 0)
            {
                logE<< "Could not read spectra dims\n";
                _close_h5_objects(close_map);
                return false;
            }

            count3[2] = dims3[2];
            count[0] = dims3[2];

            hid_t memoryspace_id = H5Screate_simple(1, count, nullptr);
            close_map.push({ memoryspace_id, H5O_DATASPACE });

            hid_t tid1 = H5Tcopy (H5T_C_S1);
            H5Tset_size (tid1, H5T_VARIABLE);
            H5Tset_cset(tid1, H5T_CSET_UTF8);

            if (_open_h5_object(title_id, H5O_DATASET, close_map, str_title_path, file_id, false, false))   
            {
                char* tmp_char_arr[255];
                err = H5Dread(title_id, tid1, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)tmp_char_arr);
                if (err == 0)
                {
                    scan_info.meta_info.name = std::string(*tmp_char_arr);
                }
            }

            if (_open_h5_object(start_time_id, H5O_DATASET, close_map, str_start_time_path,  file_id, false, false))
            {
                char* tmp_char_arr[255];
                err = H5Dread(start_time_id, tid1, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)tmp_char_arr);

                if (err == 0)
                {
                    scan_info.meta_info.scan_time_stamp = std::string(*tmp_char_arr);
                }
            }
            
            H5Tclose(tid1);

            scan_info.meta_info.theta = 0;
            scan_info.meta_info.scan_type = STR_SCAN_TYPE_POLAR_XANES;
            // Take 1 spectra per magnetic polarization, divide them by 2 and store in each row
            scan_info.meta_info.requested_cols = dims3[0] / 2;
            scan_info.meta_info.requested_rows = 2;  // two polarities
            scan_info.meta_info.x_axis.resize(scan_info.meta_info.requested_cols);
            scan_info.meta_info.y_axis.resize(scan_info.meta_info.requested_rows);
            scan_info.meta_info.polarity_pattern = params_override->polarity_pattern;
            spec_vol->resize_and_zero(scan_info.meta_info.requested_rows, scan_info.meta_info.requested_cols, dims3[2]);

            struct data_struct::Scaler_Map<T_real2> energy_map;
            energy_map.name = "Energy";
            energy_map.unit = "";
            energy_map.time_normalized = false;
            energy_map.values.resize(scan_info.meta_info.requested_rows, scan_info.meta_info.requested_cols);

            struct data_struct::Scaler_Map<T_real2> i0_map;
            i0_map.name = "I0";
            i0_map.unit = "counts";
            i0_map.time_normalized = false;
            i0_map.values.resize(scan_info.meta_info.requested_rows, scan_info.meta_info.requested_cols);

            data_struct::ArrayTr<T_real> elt_array(dims3[0]);
            data_struct::ArrayTr<T_real> ert_array(dims3[0]);
            data_struct::ArrayTr<T_real> icr_array(dims3[0]);
            data_struct::ArrayTr<T_real> ocr_array(dims3[0]);

            data_struct::ArrayTr<T_real2> i0_array(dims3[0]);
            data_struct::ArrayTr<T_real2> energy_array(dims3[0]);

            err = _read_h5d<T_real>(elt_dset_id, H5S_ALL, elt_space_id, H5P_DEFAULT, elt_array.data());
            if(err < 0)
            {
                logE<< "Could not read live time data, setting it all to 1's\n";
                elt_array.setOnes();
            }
            err = _read_h5d<T_real>(ert_dset_id, H5S_ALL, ert_space_id, H5P_DEFAULT, ert_array.data());
            if(err < 0)
            {
                logE<< "Could not read real time data, setting it all to 1's\n";
                ert_array.setOnes();
            }
            err = _read_h5d<T_real>(icr_dset_id, H5S_ALL, icr_space_id, H5P_DEFAULT, icr_array.data());
            if(err < 0)
            {
                logE<< "Could not read icr data, setting it all to 1's\n";
                icr_array.setOnes();
            }
            err = _read_h5d<T_real>(ocr_dset_id, H5S_ALL, ocr_space_id, H5P_DEFAULT, ocr_array.data());
            if(err < 0)
            {
                logE<< "Could not read ocr data, setting it all to 1's\n";
                ocr_array.setOnes();
            }

            err = _read_h5d<T_real2>(energy_dset_id, H5S_ALL, energy_space_id, H5P_DEFAULT, energy_array.data());
            if(err < 0)
            {
                logE<< "Could not read energy data, setting it all to 1's\n";
                energy_array.setOnes();
            }
            err = _read_h5d<T_real2>(i0_dset_id, H5S_ALL, i0_space_id, H5P_DEFAULT, i0_array.data());
            if(err < 0)
            {
                logE<< "Could not read i0 data, setting it all to 1's\n";
                i0_array.setOnes();
            }


            int polarity = 0;
            int idx = 0;

            if (params_override != nullptr && params_override->polarity_pattern.size() == 0)
            {
                if (use_pr1 && (pr1_array.size() != dims3[0]))
                {
                    logE << str_pr1_path << " is not the same size as spectra data. Can not load properly.\n";
                    _close_h5_objects(close_map);
                    return false;
                }
                else if (use_pr1 == false && (pr2_array.size() != dims3[0]))
                {
                    logE << str_pr2_path << " is not the same size as spectra data. Can not load properly.\n";
                    _close_h5_objects(close_map);
                    return false;
                }
            }
            int polarity_str_idx = 0;
            int left_pol_idx = 0;
            int right_pol_idx = 0;

            for(size_t i = 0; i < dims3[0]; i++)
            {
                offset3[0] = i;
                H5Sselect_hyperslab(space_id, H5S_SELECT_SET, offset3, nullptr, count3, nullptr);

                if (params_override != nullptr && params_override->polarity_pattern.size() > 0)
                {
                    if (polarity_str_idx >= params_override->polarity_pattern.size())
                    {
                        polarity_str_idx = 0;
                    }
                    if (params_override->polarity_pattern.at(polarity_str_idx) == 'L')
                    {
                        polarity = 0;
                        idx = left_pol_idx;
                        left_pol_idx++;
                    }
                    else if (params_override->polarity_pattern.at(polarity_str_idx) == 'R')
                    {
                        polarity = 1;
                        idx = right_pol_idx;
                        right_pol_idx++;
                    }
                    else
                    {
                        logW << "Unknown polarity override " << params_override->polarity_pattern[polarity_str_idx] << " in override string : " << params_override->polarity_pattern << "\n";
                    }
                }
                else
                {
                    if (i > 0)
                    {
                        if (use_pr1)
                        {
                            if (pr1_array[i] != pr1_array[i - 1])
                            {
                                if (polarity == 1)
                                {
                                    polarity = 0;
                                    idx++;
                                }
                                else
                                {
                                    polarity = 1;
                                }
                            }
                        }
                        else
                        {
                            if (pr2_array[i] != pr2_array[i - 1])
                            {
                                if (polarity == 1)
                                {
                                    polarity = 0;
                                }
                                else
                                {
                                    polarity = 1;
                                }
                            }
                            else
                            {
                                idx++;
                            }
                        }
                    }
                }

                err = _read_h5d<T_real>(dset_id, memoryspace_id, space_id, H5P_DEFAULT, &(*spec_vol)[polarity][idx][0]);
                if(err < 0)
                {
                    logE<< "Could not read spectra data at "<<i<<". Setting this to all 0's\n";
                    (*spec_vol)[idx][polarity].setZero();
                }
                (*spec_vol)[polarity][idx].elapsed_livetime(elt_array[i]);
                (*spec_vol)[polarity][idx].elapsed_realtime(ert_array[i]);
                (*spec_vol)[polarity][idx].input_counts(icr_array[i]);
                (*spec_vol)[polarity][idx].output_counts(ocr_array[i]);

                energy_map.values(polarity, idx) = energy_array[i];
                i0_map.values(polarity, idx) = i0_array[i];
            }

            scan_info.scaler_maps.push_back(energy_map);
            scan_info.scaler_maps.push_back(i0_map);

            /*
            // x, y, and z dataset are all rows x cols size 
            // read x motor scaler from /entry/measurement/pseudo/x
            // read y motor scaler from /entry/measurement/pseudo/y
            // read z motor scaler from /entry/measurement/pseudo/z
            bool has_x = _open_h5_object(x_dset_id, H5O_DATASET, close_map, "/entry/measurement/pseudo/x", file_id, false, false);
            x_space_id = H5Dget_space(x_dset_id);
            close_map.push({ x_space_id, H5O_DATASPACE });
            bool has_y = _open_h5_object(y_dset_id, H5O_DATASET, close_map, "/entry/measurement/pseudo/y", file_id, false, false);
            y_space_id = H5Dget_space(y_dset_id);
            close_map.push({ y_space_id, H5O_DATASPACE });
            hsize_t frame_dims[1] = {1};
            hsize_t c_dims[1] = {1};
            hid_t frame_memoryspace_id = H5Screate_simple(1, frame_dims, nullptr);
                
            if(has_y)
            {
                hsize_t i =0;
                for (hsize_t y = 0; y < dims[1]; y++)
                {
                    frame_dims[0] = i;
                    H5Sselect_hyperslab(y_space_id, H5S_SELECT_SET, frame_dims, nullptr, c_dims, nullptr);    
                    herr_t err = _read_h5d<T_real>(y_dset_id, frame_memoryspace_id, y_space_id, H5P_DEFAULT, &scan_info.meta_info.y_axis[y]);
                    if(err > 0)
                    {
                        logE<< "can not read y pos "<<y<<" : "<<i<<"\n";
                    }
                    i++;
                }
            }
            if(has_x)
            {
                hsize_t i =0;
                for (hsize_t x = 0; x < dims[0]; x++)
                {
                    frame_dims[0] = i;
                    H5Sselect_hyperslab(x_space_id, H5S_SELECT_SET, frame_dims, nullptr, c_dims, nullptr);
                    herr_t err = _read_h5d<T_real>(x_dset_id, frame_memoryspace_id, x_space_id, H5P_DEFAULT, &scan_info.meta_info.x_axis[x]);
                    if(err > 0)
                    {
                        logE<< "can not read x pos "<<x<<" : "<<i<<"\n";
                    }
                    i += dims[1];
                }
            }
                */
            // read scaler /entry/snapshots/pre_scan/energy or /entry/snapshots/post_scan/energy for incident energy

            // Unpack the external link value
            //err = H5Lunpack_elink_val((const void*)buff, 1024, &flags, (const char**)&xspress3_link_name[0], (const char**)&xspress3_link[0]);
            /*
            err = H5Lunpack_elink_val((const void*)buff, 1024, &flags, &fname, &objname);
            if (err < 0) 
            {
                logE<< "\n\nError unpacking external link value\n";
                _close_h5_objects(close_map);
                return false;
            }
            link_name = std::string(fname);
            logI<<link_name<<" : "<< objname<<"\n\n";
            //logI<<fname<<" : "<< objname<<"\n\n";
            */
            _close_h5_objects(close_map);
            return true;
        }
        //return load_spectra_volume_xspress3(path + DIR_END_CHAR + link_name, detector_num, spec_vol);
        return false;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_spectra_vol_esrf(std::string path, std::string &title, data_struct::Spectra_Volume<T_real>* spec_vol, data_struct::Scan_Info<T_real> &scan_info_edf, [[maybe_unused]] bool logerr = true)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
        hid_t    file_id, root_grp_id, fluo_grp_id, title_id, scanDim1_id, scanDim2_id, start_time_id;
        char root_group_name[2048] = { 0 };

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
            return false;

        H5Gget_objname_by_idx(file_id, 0, &root_group_name[0], 2047);

        if (false == _open_h5_object(root_grp_id, H5O_GROUP, close_map, root_group_name, file_id))
            return false;

        if (false == _open_h5_object(fluo_grp_id, H5O_GROUP, close_map, "FLUO", root_grp_id))
            return false;

        if (false == _open_h5_object(scanDim1_id, H5O_DATASET, close_map, "scanDim_1", fluo_grp_id))
            return false;

        if (false == _open_h5_object(scanDim2_id, H5O_DATASET, close_map, "scanDim_2", fluo_grp_id))
            return false;

        if (false == _open_h5_object(title_id, H5O_DATASET, close_map, "title", root_grp_id))
            return false;

        if (false == _open_h5_object(start_time_id, H5O_DATASET, close_map, "start_time", root_grp_id))
            return false;
        
        char tmp_name[256] = { 0 };
        hid_t ftype = H5Dget_type(title_id);
        close_map.push({ ftype, H5O_DATATYPE });
        hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
        herr_t error = H5Dread(title_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)&tmp_name[0]);

        if (error == 0)
        {
            title = std::string(tmp_name);
        }

        char tmp_time[256] = { 0 };
        hid_t ftype2 = H5Dget_type(start_time_id);
        close_map.push({ ftype2, H5O_DATATYPE });
        hid_t type2 = H5Tget_native_type(ftype2, H5T_DIR_ASCEND);
        error = H5Dread(start_time_id, type2, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)&tmp_time[0]);

        if (error == 0)
        {
            scan_info_edf.meta_info.scan_time_stamp = std::string(tmp_time);
        }
        else
        {
            scan_info_edf.meta_info.scan_time_stamp = "";
        }

        int rows = 0;
        int cols = 0;

        error = H5Dread(scanDim1_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cols);
        error = H5Dread(scanDim2_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &rows);


        scan_info_edf.meta_info.x_axis.resize(cols);
        scan_info_edf.meta_info.y_axis.resize(rows);
        scan_info_edf.meta_info.name = title;

        spec_vol->resize_and_zero(rows, cols, 2048);

        _close_h5_objects(close_map);

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_spectra_volume_emd(std::string path, size_t frame_num, data_struct::Spectra_Volume<T_real> *spec_vol, bool logerr = true)
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
        hsize_t image_dims[3] = { 0,0,0 };
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
        if (error > -1)
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

        if (frame_num > 0)
        {
            if (false == _open_h5_object(frame_id, H5O_DATASET, close_map, "FrameLocationTable", spectra_grp_id))
                return false;
            dataspace_frame_id = H5Dget_space(frame_id);
            close_map.push({ dataspace_frame_id, H5O_DATASPACE });

            hsize_t frame_dims[2] = { 0,0 };
            error = H5Sget_simple_extent_dims(dataspace_frame_id, &frame_dims[0], nullptr);

            if (frame_num < frame_dims[0])
            {
                int* frames = new int[frame_dims[0]];
                //read frames
                hid_t frame_memoryspace_id = H5Screate_simple(2, frame_dims, nullptr);
                error = H5Dread(frame_id, H5T_NATIVE_INT, frame_memoryspace_id, dataspace_frame_id, H5P_DEFAULT, &frames[0]);
                start_offset = frames[frame_num];
                delete[] frames;
                close_map.push({ frame_memoryspace_id, H5O_DATASPACE });
            }
            else
            {
                logW << "frame_num " << frame_num << " > Max Frames" << frame_dims[0] << ". Setting frame_num = 0\n";
            }

        }
        if (false == _open_h5_object(meta_id, H5O_DATASET, close_map, "Metadata", spectra_grp_id))
            return false;
        dataspace_meta_id = H5Dget_space(meta_id);
        close_map.push({ dataspace_meta_id, H5O_DATASPACE });

        char* acquisition[1];
        hid_t ftype = H5Dget_type(acqui_id);
        close_map.push({ ftype, H5O_DATATYPE });
        hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
        error = H5Dread(acqui_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, acquisition);

        //prase bincount/samples width and height
        std::string str_acqui = std::string(acquisition[0]);
        int samples = parse_str_val_to_int<T_real>(STR_BINCOUNT, "\"", str_acqui);

        if (samples < 0 || height < 0 || width < 0)
        {
            logE << "Unknown spectra volume size Width: " << width << " Height: " << height << " Samples: " << samples << "\n";
            _close_h5_objects(close_map);
            return false;
        }

        spec_vol->resize_and_zero(height, width, samples);

        hsize_t rank;
        rank = H5Sget_simple_extent_ndims(dataspace_id);
        unsigned short* buffer;
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

        if (start_offset > count[0])
        {
            logW << "frame start offset: " << start_offset << " > dataset count: " << count[0] << ". Setting start offset = 0\n.";
            start_offset = 0;
        }

        memoryspace_id = H5Screate_simple(2, chunk_dims, nullptr);
        close_map.push({ memoryspace_id, H5O_DATASPACE });
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);

        size_t read_leftover = (count[0] - start_offset) % chunk_dims[0];
        size_t read_amt = (count[0] - start_offset) / chunk_dims[0];

        int row = 0;
        int col = 0;

        offset[0] = start_offset; // start of frame offset
        T_real incnt = 0.0;
        T_real outcnt = 0.0;
        T_real elt = 1.0; // TODO: read from metadata
        data_struct::Spectra<T_real>* spectra = &((*spec_vol)[row][col]);
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
                        if (col >= width)
                        {
                            col = 0;
                            row++;
                            logI << "Reading row " << row << " of " << height << "\n";
                        }
                        if (row >= height)
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
                logE << "reading chunk " << i << " of " << read_amt << "\n";
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
                    if (col >= width)
                    {
                        col = 0;
                        row++;
                        logI << "Reading row " << row << " of " << height << "\n";
                    }
                    if (row >= height)
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

    template<typename T_real>
    bool load_spectra_line_xspress3(std::string path, size_t detector_num, data_struct::Spectra_Line<T_real>* spec_row)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        //_is_loaded = ERROR_LOADING;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        //logI << path << " detector : " << detector_num << "\n";
        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
        hid_t    file_id;

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
        {
            return false;
        }
        bool ret = _load_spectra_line_xspress3(file_id, detector_num, spec_row);
        _close_h5_objects(close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        //logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";
        return ret;
    }

        //-----------------------------------------------------------------------------

    template<typename T_real>
    bool _load_spectra_line_xspress3(hid_t file_id, size_t detector_num, data_struct::Spectra_Line<T_real>* spec_row)
    {
        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
        hid_t  dset_id = -1, dataspace_id = -1, maps_grp_id = -1, scaler_grp_id = -1, scaler2_grp_id = -1, memoryspace_id = -1, memoryspace_meta_id = -1;
        hid_t dset_lt_id = -1, dset_rt_id = -1, dset_outcnt_id = -1;
        hid_t    dataspace_lt_id = -1, dataspace_rt_id = -1, dataspace_outcnt_id = -1;
      
        herr_t   error = -1;
        std::string detector_path;
        T_real* buffer = nullptr;
        hsize_t offset_row[3] = { 0,0,0 };
        hsize_t count_row[3] = { 0,0,0 };
        hsize_t offset_meta[1] = { 0 };
        hsize_t count_meta[1] = { 1 };
        int bLoadElt = 1;
        bool bLoadErt = true;
        bool bLoadOutCnt = true;

        std::string live_time_dataset_name = "CHAN" + std::to_string(detector_num + 1) + "SCA0";
        std::string live_time_dataset_name2 = "TriggerLiveTime_" + std::to_string(detector_num);
        
        std::string real_time_dataset_name2 = "RealTime_" + std::to_string(detector_num);

        std::string outcounts_dataset_name2 = "OutputCounts_" + std::to_string(detector_num);

        if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/entry/data", file_id, false, false))
        {
            if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/data", file_id, false, false))
            {
                if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/entry/instrument/detector", file_id, false, false))
                {
                    if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/entry/instrument/xspress3", file_id, false, true))
                    {
                        return false;
                    }
                }
            }
        }

        _open_h5_object(scaler_grp_id, H5O_GROUP, close_map, "/entry/instrument/NDAttributes", file_id, false, false);
            

        _open_h5_object(scaler2_grp_id, H5O_GROUP, close_map, "/entry/instrument/detector/NDAttributes", file_id, false, false);
            

        if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, "data", maps_grp_id))
        {
            return false;
        }
        dataspace_id = H5Dget_space(dset_id);
        close_map.push({ dataspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, live_time_dataset_name, scaler_grp_id, false, false))
        {
            if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, live_time_dataset_name, scaler2_grp_id, false, false))
            {
                if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, live_time_dataset_name2, scaler_grp_id, false, false))
                {
                    bLoadElt = 0;
                }
                else
                {
                    bLoadElt = 2;
                }
            }
        }
        
        if (false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, real_time_dataset_name2, scaler_grp_id, false, false))
        {
            bLoadErt = false;
        }

        if (false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, outcounts_dataset_name2, scaler_grp_id, false, false))
        {
            bLoadOutCnt = false;
        }

        if (bLoadElt > 0)
        {
            dataspace_lt_id = H5Dget_space(dset_lt_id);
            close_map.push({ dataspace_lt_id, H5O_DATASPACE });
        }
        if (bLoadErt)
        {
            dataspace_rt_id = H5Dget_space(dset_rt_id);
            close_map.push({ dataspace_rt_id, H5O_DATASPACE });
        }
        if (bLoadOutCnt)
        {
            dataspace_outcnt_id = H5Dget_space(bLoadOutCnt);
            close_map.push({ dataspace_outcnt_id, H5O_DATASPACE });
        }


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
            offset[i] = 0;
            count[i] = dims_in[i];
        }

        buffer = new T_real[dims_in[0] * dims_in[2]]; // cols x spectra_size
        count_row[0] = dims_in[0];
        count_row[1] = 1;
        count_row[2] = dims_in[2];

        size_t greater_cols = std::max(spec_row->size(), (size_t)dims_in[0]);
        size_t greater_channels = std::max(spec_row[0].size(), (size_t)dims_in[2]);

        if (spec_row->size() < dims_in[0] || spec_row[0].size() < dims_in[2])
        {
            spec_row->resize_and_zero(greater_cols, greater_channels);
        }

        memoryspace_id = H5Screate_simple(3, count_row, nullptr);
        close_map.push({ memoryspace_id, H5O_DATASPACE });
        memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
        close_map.push({ memoryspace_meta_id, H5O_DATASPACE });
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
        //H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        //T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;

        //offset[1] = row;

        offset_row[1] = detector_num;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
        error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);


        if (error > -1)
        {

            for (size_t col = 0; col < dims_in[0]; col++)
            {
                offset_meta[0] = col;
                data_struct::Spectra<T_real>* spectra = &((*spec_row)[col]);

                if (bLoadElt > 0)
                {
                    H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_lt_id, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                    if (bLoadElt == 1)
                    {
                        live_time *= 0.000000125;
                    }
                    spectra->elapsed_livetime(live_time);
                }

                if (bLoadErt)
                {
                    H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_rt_id, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                    spectra->elapsed_realtime(real_time);
                }

                if (bLoadOutCnt)
                {
                    H5Sselect_hyperslab(dataspace_outcnt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_outcnt_id, memoryspace_meta_id, dataspace_outcnt_id, H5P_DEFAULT, &out_cnt);
                    spectra->output_counts(out_cnt);
                }
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



                for (size_t s = 0; s < count_row[2]; s++) // num samples
                {
                    //(*spectra)[s] = buffer[(count_row[1] * s) + col];
                    (*spectra)[s] = buffer[(count_row[2] * col) + s];
                }
                //logD<<"saved col "<<col<<"\n";
            }
        }

        delete[] dims_in;
        delete[] offset;
        delete[] count;
        delete[] buffer;

        _close_h5_objects(close_map);

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_spectra_volume_xspress3(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        //_is_loaded = ERROR_LOADING;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        logI << path << " detector : " << detector_num << "\n";

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
        hid_t    file_id, dset_id, dataspace_id, maps_grp_id, scaler_grp_id, scaler2_grp_id, memoryspace_id, memoryspace_meta_id = -1;
        hid_t dset_lt_id, dset_rt_id, dset_outcnt_id = -1;
        hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_outcnt_id = -1;
      
        herr_t   error = -1;
        std::string detector_path;
        T_real* buffer = nullptr;
        hsize_t offset_row[3] = { 0,0,0 };
        hsize_t count_row[3] = { 0,0,0 };
        hsize_t offset_meta[1] = { 0 };
        hsize_t count_meta[1] = { 1 };
        int bLoadElt = 1;
        bool bLoadErt = true;
        bool bLoadOutCnt = true;

        std::string live_time_dataset_name = "CHAN" + std::to_string(detector_num + 1) + "SCA0";
        std::string live_time_dataset_name2 = "TriggerLiveTime_" + std::to_string(detector_num);
        
        std::string real_time_dataset_name2 = "RealTime_" + std::to_string(detector_num);

        std::string outcounts_dataset_name2 = "OutputCounts_" + std::to_string(detector_num);


        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
        {
            return false;
        }
        if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/entry/data", file_id, false, false))
        {
            if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/entry/instrument/detector", file_id, false, false))
            {
                if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "/entry/instrument/xspress3", file_id, false, true))
                {
                    return false;
                }
            }
        }

        _open_h5_object(scaler_grp_id, H5O_GROUP, close_map, "/entry/instrument/NDAttributes", file_id, false, false);
            

        _open_h5_object(scaler2_grp_id, H5O_GROUP, close_map, "/entry/instrument/detector/NDAttributes", file_id, false, false);
            

        if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, "data", maps_grp_id))
        {
            return false;
        }
        dataspace_id = H5Dget_space(dset_id);
        close_map.push({ dataspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, live_time_dataset_name, scaler_grp_id, false, false))
        {
            if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, live_time_dataset_name, scaler2_grp_id, false, false))
            {
                if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, live_time_dataset_name2, scaler_grp_id, false, false))
                {
                    bLoadElt = 0;
                }
                else
                {
                    bLoadElt = 2;
                }
            }
        }
        
        if (false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, real_time_dataset_name2, scaler_grp_id, false, false))
        {
            bLoadErt = false;
        }

        if (false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, outcounts_dataset_name2, scaler_grp_id, false, false))
        {
            bLoadOutCnt = false;
        }

        if (bLoadElt > 0)
        {
            dataspace_lt_id = H5Dget_space(dset_lt_id);
            close_map.push({ dataspace_lt_id, H5O_DATASPACE });
        }
        if (bLoadErt)
        {
            dataspace_rt_id = H5Dget_space(dset_rt_id);
            close_map.push({ dataspace_rt_id, H5O_DATASPACE });
        }
        if (bLoadOutCnt)
        {
            dataspace_outcnt_id = H5Dget_space(bLoadOutCnt);
            close_map.push({ dataspace_outcnt_id, H5O_DATASPACE });
        }


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
            offset[i] = 0;
            count[i] = dims_in[i];
        }

        buffer = new T_real[dims_in[2]]; // cols x samples_size
        count_row[0] = 1;
        count_row[1] = 1;  // detector amount
        count_row[2] = dims_in[2]; // samples_size

        // can't do this because it is width * height and will always be greater
        //size_t greater_cols = std::max(spec_row->size(), (size_t)dims_in[0]);
        size_t greater_channels = std::max(spec_vol->samples_size(), (size_t)dims_in[2]);

        //if (spec_row->size() < dims_in[0] || spec_row[0].size() < dims_in[2])
        if ((*spec_vol)[0][0].size() < dims_in[2])
        {
            spec_vol->resize_samples(greater_channels);            
        }

        memoryspace_id = H5Screate_simple(3, count_row, nullptr);
        close_map.push({ memoryspace_id, H5O_DATASPACE });
        memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
        close_map.push({ memoryspace_meta_id, H5O_DATASPACE });
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
        //H5Sselect_hyperslab (memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        //T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;

        //offset[1] = row;

        offset_row[1] = detector_num;
        
        for (size_t col = 0; col < spec_vol->cols(); col++)
        {
            for (size_t row = 0; row < spec_vol->rows(); row++)
            {
                offset_row[0] = (col * spec_vol->rows()) + row;
                H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
                error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
                if (error > -1)
                {
                    data_struct::Spectra<T_real>* spectra = &((*spec_vol)[row][col]);
                    /*
                    offset_meta[0] = col;
                    if (bLoadElt > 0)
                    {
                        H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_lt_id, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                        if (bLoadElt == 1)
                        {
                            live_time *= 0.000000125;
                        }
                        spectra->elapsed_livetime(live_time);
                    }

                    if (bLoadErt)
                    {
                        H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_rt_id, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                        spectra->elapsed_realtime(real_time);
                    }

                    if (bLoadOutCnt)
                    {
                        H5Sselect_hyperslab(dataspace_outcnt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_outcnt_id, memoryspace_meta_id, dataspace_outcnt_id, H5P_DEFAULT, &out_cnt);
                        spectra->output_counts(out_cnt);
                    }
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
                    */


                    for (size_t s = 0; s < count_row[2]; s++) // num samples
                    {
                        (*spectra)[s] = buffer[s];
                    }
                    //logD<<"saved col "<<col<<"\n";
                }
                else
                {
                    logE<<" Could not read row "<< row <<" col "<<col<<" \n";
                }
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

    template<typename T_real>
    bool load_spectra_volume_confocal(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol, bool log_error=true)
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
                logI << path << "\n";
            }
            else
            {
                logI << path << " detector : " << detector_num << "\n";
            }
        }
        hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id, dset_detectors_id, memoryspace_id2 = -1;
        //hid_t    dset_xpos_id, dset_ypos_id, dataspace_xpos_id, dataspace_ypos_id;
        hid_t    dataspace_detectors_id = -1;
        hid_t    attr_detector_names_id, attr_timebase_id = -1;
        hid_t   elt_id = -1;
        hid_t   incnt_id = -1;
        hid_t   outcnt_id = -1;

        herr_t   error;
        std::string detector_path;
        char* detector_names[256];
        T_real time_base = 1.0f;
        T_real el_time = 1.0f;
        T_real* buffer;
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

        std::string det_num_p1 = std::to_string(detector_num + 1);
        detector_path = "MCA " + det_num_p1;
        incnt_str += det_num_p1;
        outcnt_str += det_num_p1;

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
            delete[] dims_in;
            delete[] offset;
            delete[] count;
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
        buffer = new T_real[dims_in[1] * dims_in[2]]; //  cols x spectra_size
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

        if (spec_vol->rows() < dims_in[0] || spec_vol->cols() < dims_in[1] || spec_vol->samples_size() < dims_in[2])
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
                    detector_lookup[std::string(detector_names[z])] = (int)z;
                    //free(detector_names[z]);
                }
            }
            delete[] det_dims_in;
        }


        if (std::is_same<T_real, float>::value)
        {
            error = H5Aread(attr_timebase_id, H5T_NATIVE_FLOAT, &time_base);
        }
        else if (std::is_same<T_real, double>::value)
        {
            error = H5Aread(attr_timebase_id, H5T_NATIVE_DOUBLE, &time_base);
        }
        else
        {
            error = 1;
        }

        count[0] = 1; //1 row

        memoryspace_id = H5Screate_simple(2, count_row, nullptr);
        close_map.push({ memoryspace_id, H5O_DATASPACE });

        memoryspace_id2 = H5Screate_simple(2, count2, nullptr);
        close_map.push({ memoryspace_id2, H5O_DATASPACE });

        memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
        close_map.push({ memoryspace_meta_id, H5O_DATASPACE });
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
        H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

        T_real live_time = 1.0;
        //T_real real_time = 1.0;
        T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;


        for (size_t row = 0; row < dims_in[0]; row++)
        {
            offset[0] = row;
            offset_meta[0] = row;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

            if (error > -1)
            {
                for (size_t col = 0; col < dims_in[1]; col++)
                {
                    offset_meta[1] = col;

                    data_struct::Spectra<T_real>* spectra = &((*spec_vol)[row][col]);



                    if (confocal_ver_2020)
                    {
                        offset2[0] = row;
                        offset2[1] = col;
                        H5Sselect_hyperslab(memoryspace_id2, H5S_SELECT_SET, offset2, nullptr, count2, nullptr);
                        if (elt_id > -1)
                        {
                            error = _read_h5d<T_real>(elt_id, memoryspace_id2, memoryspace_id2, H5P_DEFAULT, &live_time);
                            el_time = live_time / time_base;
                            spectra->elapsed_livetime(el_time);
                        }
                        if (incnt_id > -1)
                        {
                            error = _read_h5d<T_real>(incnt_id, memoryspace_id2, memoryspace_id2, H5P_DEFAULT, &in_cnt);
                            in_cnt *= 1000.0;
                            spectra->input_counts(in_cnt);
                        }
                        if (outcnt_id > -1)
                        {
                            error = _read_h5d<T_real>(outcnt_id, memoryspace_id2, memoryspace_id2, H5P_DEFAULT, &out_cnt);
                            out_cnt *= 1000.0;
                            spectra->output_counts(out_cnt);
                        }
                    }
                    else
                    {
                        offset_meta[2] = detector_lookup[elt_str];
                        H5Sselect_hyperslab(dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_detectors_id, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &live_time);
                        el_time = live_time / time_base;
                        spectra->elapsed_livetime(el_time);

                        offset_meta[2] = detector_lookup[incnt_str];
                        H5Sselect_hyperslab(dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_detectors_id, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &in_cnt);
                        spectra->input_counts(in_cnt * 1000.0);

                        offset_meta[2] = detector_lookup[outcnt_str];
                        H5Sselect_hyperslab(dataspace_detectors_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                        error = _read_h5d<T_real>(dset_detectors_id, memoryspace_meta_id, dataspace_detectors_id, H5P_DEFAULT, &out_cnt);
                        spectra->output_counts(out_cnt * 1000.0);
                    }

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

    template<typename T_real>
	bool load_spectra_volume_gsecars(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol, bool log_error = true)
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
        hid_t    file_id, dset_id, dataspace_id, maps_grp_id, memoryspace_id, memoryspace_meta_id;
        hid_t	 livetime_id, realtime_id, inpcounts_id, outcounts_id;
        hid_t    livetime_dataspace_id, realtime_dataspace_id, inpcounts_dataspace_id, outcounts_dataspace_id;
        herr_t   error;
        std::string detector_path;
        T_real* buffer;
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
        buffer = new T_real[dims_in[1] * dims_in[2]]; //  cols x spectra_size
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

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;


        for (size_t row = 0; row < dims_in[0]; row++)
        {
            offset[0] = row;
            offset_meta[0] = row;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

            if (error > -1) //no error
            {
                for (size_t col = 0; col < dims_in[1]; col++)
                {
                    offset_meta[1] = col;

                    data_struct::Spectra<T_real>* spectra = &((*spec_vol)[row][col]);

                    H5Sselect_hyperslab(livetime_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    H5Sselect_hyperslab(realtime_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    H5Sselect_hyperslab(inpcounts_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    H5Sselect_hyperslab(outcounts_dataspace_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

                    error = _read_h5d<T_real>(realtime_id, memoryspace_meta_id, realtime_dataspace_id, H5P_DEFAULT, &real_time);
                    if (error > -1)
                    {
                        spectra->elapsed_realtime(real_time);
                    }
                    error = _read_h5d<T_real>(livetime_id, memoryspace_meta_id, livetime_dataspace_id, H5P_DEFAULT, &live_time);
                    if (error > -1)
                    {
                        spectra->elapsed_livetime(live_time);
                    }
                    error = _read_h5d<T_real>(inpcounts_id, memoryspace_meta_id, inpcounts_dataspace_id, H5P_DEFAULT, &in_cnt);
                    if (error > -1)
                    {
                        spectra->input_counts(in_cnt);
                    }
                    error = _read_h5d<T_real>(outcounts_id, memoryspace_meta_id, outcounts_dataspace_id, H5P_DEFAULT, &out_cnt);
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

    template<typename T_real>
    bool load_spectra_volume_bnl(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol, bool log_error = true)
    {
        data_struct::Scan_Info<T_real> scan_info;
        T_real dwell = 1.0;
        bool has_scaninfo = get_scalers_and_metadata_bnl(path, &scan_info);
        if (has_scaninfo)
        {
            // search for dwell time
            for (auto &itr : scan_info.extra_pvs)
            {
                if (itr.name == STR_BNL_PARAM_DWELL)
                {
                    dwell = std::stof(itr.value);
                    break;
                }
            }
        }

        std::lock_guard<std::mutex> lock(_mutex);

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
        if (log_error)
        {
            logI << path << " detector : " << detector_num << "\n";
        }
        hid_t    file_id = -1, dset_id = -1, dataspace_id = -1, maps_grp_id = -1, memoryspace_id = -1, memoryspace_meta_id = -1;
        herr_t   error = -1;
        std::string detector_path;
        T_real* buffer;
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

        if (detector_num == (size_t)-1)
        {
            detector_path = "detsum";
        }
        else
        {
            detector_path = "det" + std::to_string(detector_num+1);
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
        buffer = new T_real[dims_in[1] * dims_in[2]]; //  cols x spectra_size
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



        for (size_t row = 0; row < dims_in[0]; row++)
        {
            offset[0] = row;
            offset_meta[0] = row;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

            if (error > -1) //no error
            {
                for (size_t col = 0; col < dims_in[1]; col++)
                {
                    offset_meta[1] = col;

                    data_struct::Spectra<T_real>* spectra = &((*spec_vol)[row][col]);
                    spectra->elapsed_livetime(dwell);
                    spectra->elapsed_realtime(dwell);
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
                    // The data is float 32bit . I think it is dead time corrected because it is fractional.
                    // Initially thought it was normalized by dwell time.
                    //*spectra *= dwell; // counts can be saved in counts/sec and we fit in counts
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

    template<typename T_real>
    bool load_integrated_spectra_bnl(std::string path, size_t detector_num, data_struct::Spectra<T_real>* spec, [[maybe_unused]] bool log_error)
    {
        data_struct::Spectra_Volume<T_real> spec_vol;
        bool retval = load_spectra_volume_bnl(path, detector_num, &spec_vol);
        if (retval)
        {
            *spec = spec_vol.integrate();
        }
        return retval;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_and_integrate_spectra_volume(std::string path, size_t detector_num, data_struct::Spectra<T_real>* spectra)
    {
        std::lock_guard<std::mutex> lock(_mutex);


        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        logI << path << " detector : " << detector_num << "\n";

        hid_t    file_id = -1, dset_id = -1, dataspace_id = -1, maps_grp_id = -1, memoryspace_id = -1, memoryspace_meta_id = -1, dset_incnt_id = -1, dset_outcnt_id = -1, dset_rt_id = -1, dset_lt_id = -1;
        hid_t    dataspace_lt_id = -1, dataspace_rt_id = -1, dataspace_inct_id = -1, dataspace_outct_id = -1;
        herr_t   error = -1;
        std::string detector_path;
        T_real* buffer;
        hsize_t offset_row[2] = { 0,0 };
        hsize_t count_row[2] = { 0,0 };
        hsize_t offset_meta[3] = { 0,0,0 };
        hsize_t count_meta[3] = { 1,1,1 };


        switch (detector_num)
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

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
            return false;

        if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS_RAW", file_id))
            return false;

        if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, detector_path, maps_grp_id))
            return false;
        dataspace_id = H5Dget_space(dset_id);
        close_map.push({ dataspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "livetime", maps_grp_id))
            return false;
        dataspace_lt_id = H5Dget_space(dset_lt_id);
        close_map.push({ dataspace_lt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "realtime", maps_grp_id))
            return false;
        dataspace_rt_id = H5Dget_space(dset_rt_id);
        close_map.push({ dataspace_rt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "inputcounts", maps_grp_id))
            return false;
        dataspace_inct_id = H5Dget_space(dset_incnt_id);
        close_map.push({ dataspace_inct_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "ouputcounts", maps_grp_id))
            return false;
        dataspace_outct_id = H5Dget_space(dset_outcnt_id);
        close_map.push({ dataspace_outct_id, H5O_DATASPACE });


        int rank = H5Sget_simple_extent_ndims(dataspace_id);
        if (rank != 3)
        {
            _close_h5_objects(close_map);
            logE << "Dataset /MAPS_RAW/" << detector_path << " rank != 3. rank = " << rank << ". Can't load dataset. returning" << "\n";
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
            logE << "getting dataset rank for MAPS_RAW/" << detector_path << "\n";
            return false;
        }

        for (int i = 0; i < rank; i++)
        {

            offset[i] = 0;
            count[i] = dims_in[i];
        }

        buffer = new T_real[dims_in[0] * dims_in[2]]; // spectra_size x cols
        count_row[0] = dims_in[0];
        count_row[1] = dims_in[2];

        count[1] = 1; //1 row

        if ((hsize_t)spectra->size() != dims_in[0])
        {
            spectra->resize(dims_in[0]);
            spectra->setZero(dims_in[0]);
        }

        memoryspace_id = H5Screate_simple(2, count_row, nullptr);
        memoryspace_meta_id = H5Screate_simple(1, count_meta, nullptr);
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset_row, nullptr, count_row, nullptr);
        H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;

        T_real live_time_total = 0.0;
        T_real real_time_total = 0.0;
        T_real in_cnt_total = 0.0;
        T_real out_cnt_total = 0.0;

        offset_meta[0] = detector_num;
        for (size_t row = 0; row < dims_in[1]; row++)
        {
            offset[1] = row;
            offset_meta[1] = row;

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);

            if (error > -1)
            {
                for (size_t col = 0; col < count_row[1]; col++)
                {
                    offset_meta[2] = col;

                    H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_lt_id, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, &live_time);
                    live_time_total += live_time;

                    H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_rt_id, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, &real_time);
                    real_time_total += real_time;

                    H5Sselect_hyperslab(dataspace_inct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_incnt_id, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, &in_cnt);
                    in_cnt_total += in_cnt;

                    H5Sselect_hyperslab(dataspace_outct_id, H5S_SELECT_SET, offset_meta, nullptr, count_meta, nullptr);
                    error = _read_h5d<T_real>(dset_outcnt_id, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, &out_cnt);
                    out_cnt_total += out_cnt;

                    for (size_t s = 0; s < count_row[0]; s++)
                    {
                        (*spectra)[s] += buffer[(count_row[1] * s) + col];
                    }
                }

                //logI<<"read row "<<row<<"\n";
            }
            else
            {
                logE << "reading row " << row << "\n";
            }
        }

        spectra->elapsed_livetime(live_time_total);
        spectra->elapsed_realtime(real_time_total);
        spectra->input_counts(in_cnt_total);
        spectra->output_counts(out_cnt_total);

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

    template<typename T_real>
    bool load_spectra_vol_analyzed_h5(std::string path,
                                      data_struct::Spectra_Volume<T_real>* spectra_volume,
                                      int row_idx_start = 0,
                                      int row_idx_end = -1,
                                      int col_idx_start = 0,
                                      int col_idx_end = -1)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        //_is_loaded = ERROR_LOADING;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        logI << path << "\n";

        hid_t    file_id = -1, dset_id = -1, dataspace_id = -1, spec_grp_id = -1, memoryspace_id = -1, memoryspace_meta_id = -1, dset_incnt_id = -1, dset_outcnt_id = -1, dset_rt_id = -1, dset_lt_id = -1;
        hid_t    dataspace_lt_id = -1, dataspace_rt_id = -1, dataspace_inct_id = -1, dataspace_outct_id = -1;
        herr_t   error = -1;
        hsize_t dims_in[3] = { 0,0,0 };
        hsize_t offset[3] = { 0,0,0 };
        hsize_t count[3] = { 1,1,1 };
        hsize_t offset_time[2] = { 0,0 };
        hsize_t count_time[2] = { 1,1 };

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1, false))
            return false;

        if (false == _open_h5_object(spec_grp_id, H5O_GROUP, close_map, "/MAPS/Spectra", file_id))
            return false;

        if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, "mca_arr", spec_grp_id))
            return false;
        dataspace_id = H5Dget_space(dset_id);
        close_map.push({ dataspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "Elapsed_Livetime", spec_grp_id))
            return false;
        dataspace_lt_id = H5Dget_space(dset_lt_id);
        close_map.push({ dataspace_lt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "Elapsed_Realtime", spec_grp_id))
            return false;
        dataspace_rt_id = H5Dget_space(dset_rt_id);
        close_map.push({ dataspace_rt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "Input_Counts", spec_grp_id))
            return false;
        dataspace_inct_id = H5Dget_space(dset_incnt_id);
        close_map.push({ dataspace_inct_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "Output_Counts", spec_grp_id))
            return false;
        dataspace_outct_id = H5Dget_space(dset_outcnt_id);
        close_map.push({ dataspace_outct_id, H5O_DATASPACE });


        int rank = H5Sget_simple_extent_ndims(dataspace_id);
        if (rank != 3)
        {
            _close_h5_objects(close_map);
            logE << "Dataset /MAPS/Spectra/mca_arr  rank != 3. rank = " << rank << ". Can't load dataset. returning" << "\n";
            return false;
            //throw exception ("Dataset is not a volume");
        }

        int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
        if (status_n < 0)
        {
            _close_h5_objects(close_map);
            logE << "getting dataset dims for /MAPS/Spectra/mca_arr" << "\n";
            return false;
        }

        for (int i = 0; i < rank; i++)
        {

            offset[i] = 0;
            count[i] = dims_in[i];
        }

        if (row_idx_end < row_idx_start)
        {
            row_idx_end = dims_in[1];
        }

        if (col_idx_end < col_idx_start)
        {
            col_idx_end = dims_in[2];
        }

        spectra_volume->resize_and_zero(dims_in[1], dims_in[2], dims_in[0]);

        //buffer = new T_real [dims_in[0] * dims_in[2]]; // cols x spectra_size
        count[0] = dims_in[0];
        count[1] = 1;
        count[2] = 1;

        memoryspace_id = H5Screate_simple(3, count, nullptr);
        memoryspace_meta_id = H5Screate_simple(2, count_time, nullptr);
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;

        //offset[1] = row;

        for (size_t row = (size_t)row_idx_start; row < (size_t)row_idx_end; row++)
        {
            offset[1] = row;
            offset_time[0] = row;
            for (size_t col = (size_t)col_idx_start; col < (size_t)col_idx_end; col++)
            {
                data_struct::Spectra<T_real>* spectra = &((*spectra_volume)[row][col]);
                offset[2] = col;
                offset_time[1] = col;
                H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

                //error = H5Dread (dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);
                error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)spectra->data());
                if (error > 0)
                {
                    logW << "Counld not read row " << row << " col " << col << "\n";
                }

                H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_inct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_outct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

                error = _read_h5d<T_real>(dset_rt_id, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, (void*)&real_time);
                error = _read_h5d<T_real>(dset_lt_id, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, (void*)&live_time);
                error = _read_h5d<T_real>(dset_incnt_id, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, (void*)&in_cnt);
                error = _read_h5d<T_real>(dset_outcnt_id, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, (void*)&out_cnt);

                spectra->elapsed_livetime(live_time);
                spectra->elapsed_realtime(real_time);
                spectra->input_counts(in_cnt);
                spectra->output_counts(out_cnt);
            }
        }

        _close_h5_objects(close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

        logI << "elapsed time: " << elapsed_seconds.count() << "s\n";

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_integrated_spectra_analyzed_h5(std::string path, data_struct::Spectra<T_real>* spectra, [[maybe_unused]] ROI_Vec* roi = nullptr, bool log_error=true)
    {
        std::lock_guard<std::mutex> lock(_mutex);
        hid_t    file_id = -1;

        file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file_id < 0)
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

    template<typename T_real>
    bool load_integrated_spectra_analyzed_h5_roi(std::string path, data_struct::Spectra<T_real>* int_spectra, ROI_Vec& roi)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        bool is_v9 = false;
        //_is_loaded = ERROR_LOADING;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        logI << path << "\n";

        hid_t    file_id = -1, dset_id = -1, dataspace_id = -1, spec_grp_id = -1, memoryspace_id = -1, memoryspace_meta_id = -1, dset_incnt_id = -1, dset_outcnt_id = -1;
        hid_t   memoryspace_1 = -1;
        hid_t    dset_rt_id = -1, dset_lt_id = -1, dset_scalers = -1, dset_scaler_names = -1;
        hid_t    dataspace_lt_id = -1, dataspace_rt_id = -1, dataspace_inct_id = -1, dataspace_outct_id = -1, dataspace_scalers = -1, dataspace_scaler_names = -1;
        herr_t   error = -1;
        hsize_t dims_in[3] = { 0,0,0 };
        hsize_t offset[3] = { 0,0,0 };
        hsize_t count[3] = { 1,1,1 };
        hsize_t offset_time[2] = { 0,0 };
        hsize_t count_time[2] = { 1,1 };
        hsize_t offset_1[1] = { 0 };
        hsize_t count_1[1] = { 1 };

        hsize_t elt_off = (hsize_t)-1, ert_off = (hsize_t)-1, in_off = (hsize_t)-1, out_off = (hsize_t)-1;

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1, false))
            return false;

        if (false == _open_h5_object(spec_grp_id, H5O_GROUP, close_map, "/MAPS/Spectra", file_id, false, false))
        {
            if (false == _open_h5_object(spec_grp_id, H5O_GROUP, close_map, "/MAPS", file_id))
            {
                return false;
            }
        }
        
        if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, "mca_arr", spec_grp_id))
        {
            return false;
        }
        
        dataspace_id = H5Dget_space(dset_id);
        close_map.push({ dataspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "Elapsed_Livetime", spec_grp_id, false, false))
        {
            if (false == _open_h5_object(dset_scalers, H5O_DATASET, close_map, "scalers", spec_grp_id))
            {
                return false;
            }
            if (false == _open_h5_object(dset_scaler_names, H5O_DATASET, close_map, "scaler_names", spec_grp_id))
            {
                return false;
            }
            dataspace_scalers = H5Dget_space(dset_scalers);
            close_map.push({ dataspace_scalers, H5O_DATASPACE });

            dataspace_scaler_names = H5Dget_space(dset_scaler_names);
            close_map.push({ dataspace_scaler_names, H5O_DATASPACE });
            hid_t memtype = H5Dget_type(dset_scaler_names);
            close_map.push({ memtype, H5O_DATATYPE });

            //  read scaler names and search for elt1, ert1, incnt1, outcnt1
            hsize_t dims_out[1];
            H5Sget_simple_extent_dims(dataspace_scaler_names, &dims_out[0], nullptr);
            char tmp_name[256] = { 0 };
            memoryspace_1 = H5Screate_simple(1, count, nullptr);
            for (hsize_t idx = 0; idx < dims_out[0]; idx++)
            {
                offset_1[0] = idx;
                memset(&tmp_name[0], 0, 254);
                H5Sselect_hyperslab(dataspace_scaler_names, H5S_SELECT_SET, offset_1, nullptr, count_1, nullptr);
                error = H5Dread(dset_scaler_names, memtype, memoryspace_1, dataspace_scaler_names, H5P_DEFAULT, (void*)&tmp_name[0]);
                if (error == 0)
                {
                    std::string name = std::string(tmp_name);

                    if (name == STR_ELT + "1")
                    {
                        elt_off = idx;
                    }
                    else if (name == STR_ERT + "1")
                    {
                        ert_off = idx;
                    }
                    else if (name == STR_ICR + "1")
                    {
                        in_off = idx;
                    }
                    else if (name == STR_OCR + "1")
                    {
                        out_off = idx;
                    }
                }
            }
            is_v9 = true;
            H5Tclose(memtype);
        }
        else
        {
            dataspace_lt_id = H5Dget_space(dset_lt_id);
            close_map.push({ dataspace_lt_id, H5O_DATASPACE });

            if (false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "Elapsed_Realtime", spec_grp_id))
                return false;
            dataspace_rt_id = H5Dget_space(dset_rt_id);
            close_map.push({ dataspace_rt_id, H5O_DATASPACE });

            if (false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "Input_Counts", spec_grp_id))
                return false;
            dataspace_inct_id = H5Dget_space(dset_incnt_id);
            close_map.push({ dataspace_inct_id, H5O_DATASPACE });

            if (false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "Output_Counts", spec_grp_id))
                return false;
            dataspace_outct_id = H5Dget_space(dset_outcnt_id);
            close_map.push({ dataspace_outct_id, H5O_DATASPACE });
        }

        int rank = H5Sget_simple_extent_ndims(dataspace_id);
        if (rank != 3)
        {
            _close_h5_objects(close_map);
            logE << "Dataset /MAPS/Spectra/mca_arr  rank != 3. rank = " << rank << ". Can't load dataset. returning" << "\n";
            return false;
            //throw exception ("Dataset is not a volume");
        }

        int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
        if (status_n < 0)
        {
            _close_h5_objects(close_map);
            logE << "getting dataset dims for /MAPS/Spectra/mca_arr" << "\n";
            return false;
        }

        for (int i = 0; i < rank; i++)
        {

            offset[i] = 0;
            count[i] = dims_in[i];
        }


        int_spectra->resize(dims_in[0]);
        int_spectra->setZero(dims_in[0]);

        count[0] = dims_in[0];
        count[1] = 1;
        count[2] = 1;

        memoryspace_id = H5Screate_simple(3, count, nullptr);
        memoryspace_meta_id = H5Screate_simple(2, count_time, nullptr);
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;


        for (auto& itr : roi)
        {
            hsize_t xoffset = itr.first;
            hsize_t yoffset = itr.second;
            
            offset[0] = 0;
            offset[1] = yoffset;
            offset_time[0] = yoffset;
            offset[2] = xoffset;
            offset_time[1] = xoffset;

            count[0] = dims_in[0];

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

            data_struct::Spectra<T_real> spectra;
            spectra.resize(dims_in[0]);

            error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)spectra.data());
            if (error > 0)
            {
                logW << "Counld not read row " << yoffset << " col " << xoffset << "\n";
            }

            if (is_v9)
            {
                count[0] = 1;

                offset[0] = elt_off;
                H5Sselect_hyperslab(dataspace_scalers, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                error = _read_h5d<T_real>(dset_scalers, memoryspace_meta_id, dataspace_scalers, H5P_DEFAULT, (void*)&live_time);
                if (error > 0)
                    live_time = 0;

                offset[0] = ert_off;
                H5Sselect_hyperslab(dataspace_scalers, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                error = _read_h5d<T_real>(dset_scalers, memoryspace_meta_id, dataspace_scalers, H5P_DEFAULT, (void*)&real_time);
                if (error > 0)
                    real_time = 0;

                offset[0] = in_off;
                H5Sselect_hyperslab(dataspace_scalers, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                error = _read_h5d<T_real>(dset_scalers, memoryspace_meta_id, dataspace_scalers, H5P_DEFAULT, (void*)&in_cnt);
                if (error > 0)
                    in_cnt = 0;

                offset[0] = out_off;
                H5Sselect_hyperslab(dataspace_scalers, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                error = _read_h5d<T_real>(dset_scalers, memoryspace_meta_id, dataspace_scalers, H5P_DEFAULT, (void*)&out_cnt);
                if (error > 0)
                    out_cnt = 0;
            }
            else
            {
                H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_inct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_outct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

                error = _read_h5d<T_real>(dset_rt_id, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, (void*)&real_time);
                error = _read_h5d<T_real>(dset_lt_id, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, (void*)&live_time);
                error = _read_h5d<T_real>(dset_incnt_id, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, (void*)&in_cnt);
                error = _read_h5d<T_real>(dset_outcnt_id, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, (void*)&out_cnt);
            }

            spectra.elapsed_livetime(live_time);
            spectra.elapsed_realtime(real_time);
            spectra.input_counts(in_cnt);
            spectra.output_counts(out_cnt);

            int_spectra->add(spectra);
        }

        int_spectra->recalc_elapsed_livetime();

        _close_h5_objects(close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

        logI << "elapsed time: " << elapsed_seconds.count() << "s\n";

        return true;
    }
 
    //-----------------------------------------------------------------------------

    /**
    * Loads only Upstream/Downstream Ion chambers and SR_Current
    */

    template<typename T_real>
    bool load_quantification_scalers_analyzed_h5(std::string path, data_struct::Params_Override<T_real> *override_values)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        //_is_loaded = ERROR_LOADING;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        logI << path << "\n";

        std::map<std::string, data_struct::ArrayXXr<T_real>> scalers_map;
        T_real srcurrent, us_ic, us_fm, ds_ic;
        hid_t    file_id = -1, ds_ic_id = -1, us_ic_id = -1, us_fm_id = -1, srcurrent_id = -1;
        hid_t d_space;
        hid_t d_type;
        hid_t status;
        //herr_t   error;
        //hsize_t offset[1] = {0};
        hsize_t count[1] = { 1 };
        hid_t readwrite_space = H5Screate_simple(1, &count[0], &count[0]);

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
            return false;

        if (false == _open_h5_object(ds_ic_id, H5O_DATASET, close_map, "/MAPS/Quantification/Standard0/Scalers/DS_IC", file_id, false, false))
        {
            if (false == _open_h5_object(ds_ic_id, H5O_DATASET, close_map, "/MAPS/Quantification/Standard1/Scalers/DS_IC", file_id, false, false))
            {
                ds_ic_id = -1;
                override_values->DS_IC = 1.0;
            }
        }
            
        if (false == _open_h5_object(us_ic_id, H5O_DATASET, close_map, "/MAPS/Quantification/Standard0/Scalers/US_IC", file_id, false, false))
        {
            if (false == _open_h5_object(us_ic_id, H5O_DATASET, close_map, "/MAPS/Quantification/Standard1/Scalers/US_IC", file_id, false, false))
            {
                us_ic_id = -1;
                override_values->US_IC = 1.0;
            }
        }

        if (false == _open_h5_object(us_fm_id, H5O_DATASET, close_map, "/MAPS/Quantification/Standard0/Scalers/US_FM", file_id, false, false))
        {
            if (false == _open_h5_object(us_fm_id, H5O_DATASET, close_map, "/MAPS/Quantification/Standard1/Scalers/US_FM", file_id, false, false))
            {
                us_fm_id = -1;
                override_values->US_FM = 1.0;
            }
        }

        if (false == _open_h5_object(srcurrent_id, H5O_DATASET, close_map, "/MAPS/Quantification/Standard0/Scalers/SR_Current", file_id, false, false))
        {
            if (false == _open_h5_object(srcurrent_id, H5O_DATASET, close_map, "/MAPS/Quantification/Standard1/Scalers/SR_Current", file_id, false, false))
            {
                srcurrent_id = -1;
                override_values->sr_current = 1.0;
            }
        }

        if (srcurrent_id == -1 || us_ic_id == -1 || us_fm_id == -1 || ds_ic_id == -1)
        {
            load_scalers_analyzed_h5(path, scalers_map);
        }


        if (srcurrent_id > -1)
        {
            //read in scaler
            d_space = H5Dget_space(srcurrent_id);
            close_map.push({ d_space, H5O_DATASPACE });
            d_type = H5Dget_type(srcurrent_id);
            close_map.push({ d_type, H5O_DATATYPE });
            status = H5Dread(srcurrent_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)&srcurrent);
            if (status > -1)
            {
                override_values->sr_current = (srcurrent);
            }
        }
        else if(scalers_map.count(STR_SR_CURRENT) > 0)
        {
            override_values->sr_current = scalers_map.at(STR_SR_CURRENT).sum() / scalers_map.at(STR_SR_CURRENT).size();
        }

        if (us_ic_id > -1)
        {
            d_space = H5Dget_space(us_ic_id);
            close_map.push({ d_space, H5O_DATASPACE });
            d_type = H5Dget_type(us_ic_id);
            close_map.push({ d_type, H5O_DATATYPE });
            status = H5Dread(us_ic_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)&us_ic);
            if (status > -1)
            {
                override_values->US_IC = (us_ic);
            }
        }
        else if (scalers_map.count(STR_US_IC) > 0)
        {
            override_values->US_IC = scalers_map.at(STR_US_IC).sum() / scalers_map.at(STR_US_IC).size();
        }

        if (us_fm_id > -1)
        {
            d_space = H5Dget_space(us_fm_id);
            close_map.push({ d_space, H5O_DATASPACE });
            d_type = H5Dget_type(us_fm_id);
            close_map.push({ d_type, H5O_DATATYPE });
            status = H5Dread(us_fm_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)&us_fm);
            if (status > -1)
            {
                override_values->US_FM = (us_fm);
            }
        }
        else if (scalers_map.count(STR_US_FM) > 0)
        {
            override_values->US_FM = scalers_map.at(STR_US_FM).sum() / scalers_map.at(STR_US_FM).size();
        }

        if (ds_ic_id > -1)
        {
            d_space = H5Dget_space(ds_ic_id);
            close_map.push({ d_space, H5O_DATASPACE });
            d_type = H5Dget_type(ds_ic_id);
            close_map.push({ d_type, H5O_DATATYPE });
            status = H5Dread(ds_ic_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)&ds_ic);
            if (status > -1)
            {
                override_values->DS_IC = (ds_ic);
            }
        }
        else if (scalers_map.count(STR_DS_IC) > 0)
        {
            override_values->DS_IC = scalers_map.at(STR_DS_IC).sum() / scalers_map.at(STR_DS_IC).size();
        }


        _close_h5_objects(close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        //std::time_t end_time = std::chrono::system_clock::to_time_t(end);

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool _load_calibration_curve_analyzed_h5(hid_t file_id, Fitting_Routines fitting_routine, data_struct::Detector<T_real>* detector)
    {
        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        hid_t cc_id;
        hid_t cc_labels_id;
        
        std::string str_cc_ds = "/MAPS/"+ STR_QUANTIFICATION + "/" + STR_CALIBRATION +  "/" + Fitting_Routine_To_Str.at(fitting_routine) + "/" + STR_CALIB_CURVE + STR_DS_IC;
        std::string str_cc_us = "/MAPS/" + STR_QUANTIFICATION + "/" + STR_CALIBRATION + "/" + Fitting_Routine_To_Str.at(fitting_routine) + "/" + STR_CALIB_CURVE + STR_US_IC;
        std::string str_cc_us_fm = "/MAPS/" + STR_QUANTIFICATION + "/" + STR_CALIBRATION + "/" + Fitting_Routine_To_Str.at(fitting_routine) + "/" + STR_CALIB_CURVE + STR_US_FM;
        std::string str_cc_sr = "/MAPS/" + STR_QUANTIFICATION + "/" + STR_CALIBRATION + "/" + Fitting_Routine_To_Str.at(fitting_routine) + "/" + STR_CALIB_CURVE + STR_SR_CURRENT;
        std::string str_cc_labels = "/MAPS/" + STR_QUANTIFICATION + "/" + STR_CALIBRATION + "/" + Fitting_Routine_To_Str.at(fitting_routine) + "/" + STR_CALIB_LABELS;
        
        if (false == _open_h5_object(cc_labels_id, H5O_DATASET, close_map, str_cc_labels, file_id))
        {
            logW << "Could not find labels " << str_cc_labels << ". Can not load the rest for " << Fitting_Routine_To_Str.at(fitting_routine) << "\n";
            return false;
        }

        std::map <std::string, std::string> ion_chambers = { {STR_DS_IC,str_cc_ds}, {STR_US_IC,str_cc_us}, {STR_US_FM,str_cc_us_fm}, {STR_SR_CURRENT,str_cc_sr} };

        // read ion chambers DS_IC, US_IC, and SR_Current
        for (auto& itr : ion_chambers)
        {
            //detector->fitting_quant_map.at(fitting_routine)[itr.first] = Quantification_Scaler_Struct<T_real>();
            if (true == _open_h5_object(cc_id, H5O_DATASET, close_map, itr.second.c_str(), file_id, false, false))
            {
                hsize_t offset[2] = { 0,0 };
                hsize_t count[2] = { 1,1 };
                hsize_t dims_in[2] = { 0,0 };
                hid_t d_space = H5Dget_space(cc_id);
                close_map.push({ d_space, H5O_DATASPACE });
                hid_t d_type = H5Dget_type(cc_id);
                close_map.push({ d_type, H5O_DATATYPE });

                int rank = H5Sget_simple_extent_ndims(d_space);
                if (rank != 2)
                {
                    _close_h5_objects(close_map);
                    logE << itr.second <<"  rank != 2. rank = " << rank << ". Can't load dataset. returning" << "\n";
                    return false;
                }

                int status_n = H5Sget_simple_extent_dims(d_space, &dims_in[0], nullptr);
                if (status_n < 0)
                {
                    _close_h5_objects(close_map);
                    logE << "getting dataset dims for "<< itr.second << "\n";
                    return false;
                }

                for (int i = 0; i < rank; i++)
                {

                    offset[i] = 0;
                    count[i] = dims_in[i];
                }

                //detector->fitting_quant_map.at(fitting_routine).at(itr.first).at(data_struct::Electron_Shell::K_SHELL)

                //detector->fitting_quant_map.at(fitting_routine).at(itr.first][data_struct::Electron_Shell::K_SHELL].resize(count[1]);
                //detector->fitting_quant_map.at(fitting_routine).at(itr.first][data_struct::Electron_Shell::L_SHELL].resize(count[1]);
                //detector->fitting_quant_map.at(fitting_routine)[itr.first][data_struct::Electron_Shell::M_SHELL].resize(count[1]);
                // vector < Element_Quant<T_real>
            }
            else
            {
                logW << "Could not find " << itr.second << "\n";
            }
        }

        _close_h5_objects(close_map);

        return true;
    }

    //-----------------------------------------------------------------------------

    /**
    * Loads whole quantification info
    */

    // TODO FINISH before using
    template<typename T_real>
    bool load_quantification_analyzed_h5(std::string path, data_struct::Detector<T_real>* detector)
    {

        std::lock_guard<std::mutex> lock(_mutex);

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        quantification::models::Quantification_Model<T_real> quantification_model;

        int num_standards = 0;

        hid_t file_id = -1;
        hid_t num_standards_id = -1;
        hid_t name_id = -1;
        hid_t ic_id = -1;

        herr_t status;

        char dset_name[2048] = { 0 };

        hsize_t count[1] = { 1 };

        std::vector< Fitting_Routines> fit_routes = { Fitting_Routines::NNLS, Fitting_Routines::GAUSS_MATRIX ,Fitting_Routines::ROI, Fitting_Routines::SVD };

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        if (detector == nullptr)
        {
            logW << "Can not load to nullptr detector\n";
            return false;
        }

        hid_t readwrite_space = H5Screate_simple(1, &count[0], &count[0]);
        close_map.push({ readwrite_space, H5O_DATASPACE });

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
        {
            logE << "Failed to open " << path << "\n";
            return false;
        }

        if (_open_h5_object(num_standards_id, H5O_DATASET, close_map, "/MAPS/Quantification/Number_Of_Standards", file_id))
        {
            hid_t d_space = H5Dget_space(num_standards_id);
            close_map.push({ d_space, H5O_DATASPACE });
            hid_t d_type = H5Dget_type(num_standards_id);
            close_map.push({ d_type, H5O_DATATYPE });
            status = H5Dread(num_standards_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)&num_standards);
        }
        else
        {
            logE << "Could not determine the number of standards! Returning\n";
            return false;
        }

        // Init fit quant map with all fitting routines
        for (auto& itr : fit_routes)
        {
            detector->fitting_quant_map[itr] = Fitting_Quantification_Struct<T_real>();
        }

        // iterate over number of standards and load 
        for (int i = 0; i < num_standards; i++)
        {
            // name
            std::string std_name = "";
            std::string str_std_name_loc = "/MAPS/Quantification/Standard" + std::to_string(0) + "/" + STR_STANDARD_NAME;
            if (_open_h5_object(name_id, H5O_DATASET, close_map, str_std_name_loc, file_id, false, false))
            {
                hid_t d_space = H5Dget_space(name_id);
                close_map.push({ d_space, H5O_DATASPACE });
                hid_t d_type = H5Dget_type(name_id);
                close_map.push({ d_type, H5O_DATATYPE });
                status = H5Dread(name_id, d_type, readwrite_space, d_space, H5P_DEFAULT, (void*)dset_name);
                if (status > -1)
                {
                    std_name = std::string(dset_name, 2047);
                    std_name.erase(std::remove_if(std_name.begin(), std_name.end(), ::isspace), std_name.end());
                }
            }
            
            if (std_name.length() == 0)
            {
                std_name = std::to_string(i);
            }

            detector->quantification_standards[std_name] = Quantification_Standard<T_real>();
            detector->quantification_standards[std_name].standard_filename = std_name;

            std::vector <std::string> ion_chambers = { STR_DS_IC, STR_US_IC, STR_US_FM, STR_SR_CURRENT };

            // read ion chambers DS_IC, US_IC, and SR_Current
            for (auto& itr : ion_chambers)
            {
                std::string str_loc = "/MAPS/Quantification/Standard" + std::to_string(0) + "/" + STR_SCALERS + "/" + itr;
                if (_open_h5_object(ic_id, H5O_DATASET, close_map, str_loc, file_id, false, false))
                {

                    hid_t d_space = H5Dget_space(ic_id);
                    close_map.push({ d_space, H5O_DATASPACE });
                    hid_t d_type = H5Dget_type(ic_id);
                    close_map.push({ d_type, H5O_DATATYPE });


                    T_real* val_addr = nullptr;

                    if (itr == STR_DS_IC)
                        val_addr = &(detector->quantification_standards[std_name].DS_IC);
                    else if (itr == STR_US_IC)
                        val_addr = &(detector->quantification_standards[std_name].US_IC);
                    else if (itr == STR_US_FM)
                        val_addr = &(detector->quantification_standards[std_name].US_FM);
                    else if (itr == STR_SR_CURRENT)
                        val_addr = &(detector->quantification_standards[std_name].sr_current);
                    else
                        val_addr = nullptr;
                    

                    if (val_addr != nullptr)
                    {
                        //read value
                        status = _read_h5d<T_real>(ic_id, readwrite_space, d_space, H5P_DEFAULT, (void*)val_addr);
                        if (status > -1)
                        {
                            Quantification_Standard<T_real>* quantification_standard = &(detector->quantification_standards[std_name]);
                            for (auto& fit_itr : fit_routes)
                            {
                                detector->update_element_quants(fit_itr, itr, quantification_standard, &quantification_model, *val_addr);
                            }
                        }
                        else
                        {
                            logE << "Could not read in " << str_loc << "\n";
                        }
                    }
                }
            }
        }

        for (auto& itr : fit_routes)
        {
            _load_calibration_curve_analyzed_h5(file_id, itr, detector);
        }

        _close_h5_objects(close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_scalers_analyzed_h5(std::string path, std::map<std::string, data_struct::ArrayXXr<T_real>> &scalers_map)
    {
        hid_t file_id = -1;
        hid_t maps_grp_id = -1;
        hid_t counts_dset_id, channels_dset_id, counts_dspace_id, channels_dspace_id;
        hid_t memoryspace_id, memoryspace_name_id, error;
        hsize_t offset[3] = { 0,0,0 };
        hsize_t count[3] = { 1,1,1 };
        hsize_t offset_name[1] = { 0 };
        hsize_t count_name[1] = { 1 };
        char tmp_name[255] = { 0 };
        hid_t   filetype, memtype, status;

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
        {
            logE << "Failed to open " << path << "\n";
            return false;
        }

        if (false == _open_h5_object(maps_grp_id, H5O_GROUP, close_map, "MAPS", file_id))
        {
            return false;
        }
        
        if (false == _open_h5_object(counts_dset_id, H5O_DATASET, close_map, "Scalers/Values", maps_grp_id, false, false))
        {
            if (false == _open_h5_object(counts_dset_id, H5O_DATASET, close_map, "scalers", maps_grp_id, false, false))
            {
                return false;
            }
        }
        counts_dspace_id = H5Dget_space(counts_dset_id);
        close_map.push({ counts_dspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(channels_dset_id, H5O_DATASET, close_map, "Scalers/Names", maps_grp_id, false, false))
        {
            if (false == _open_h5_object(channels_dset_id, H5O_DATASET, close_map, "scaler_names", maps_grp_id))
            {
                return false;
            }
        }
        
        channels_dspace_id = H5Dget_space(channels_dset_id);
        close_map.push({ channels_dspace_id, H5O_DATASPACE });

        int rank = H5Sget_simple_extent_ndims(counts_dspace_id);
        if (rank != 3)
        {
            logE << "Error getting rank for /MAPS/Scalers/Values\n";
        }
        hsize_t* dims_out = new hsize_t[rank];
        H5Sget_simple_extent_dims(counts_dspace_id, &dims_out[0], nullptr);

        filetype = H5Tcopy(H5T_C_S1);
        H5Tset_size(filetype, 256);
        memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

        for (int i = 0; i < 3; i++)
        {
            offset[i] = 0;
            count[i] = dims_out[i];
        }

        count[0] = 1;

        memoryspace_id = H5Screate_simple(3, count, nullptr);
        close_map.push({ memoryspace_id, H5O_DATASPACE });
        memoryspace_name_id = H5Screate_simple(1, count_name, nullptr);
        close_map.push({ memoryspace_name_id, H5O_DATASPACE });
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        H5Sselect_hyperslab(memoryspace_name_id, H5S_SELECT_SET, offset_name, nullptr, count_name, nullptr);

        for (hsize_t el_idx = 0; el_idx < dims_out[0]; el_idx++)
        {
            offset[0] = el_idx;
            offset_name[0] = el_idx;
            memset(&tmp_name[0], 0, 254);
            H5Sselect_hyperslab(channels_dspace_id, H5S_SELECT_SET, offset_name, nullptr, count_name, nullptr);
            error = H5Dread(channels_dset_id, memtype, memoryspace_name_id, channels_dspace_id, H5P_DEFAULT, (void*)&tmp_name[0]);
            std::string el_name = std::string(tmp_name);
            scalers_map.emplace(std::pair<std::string, data_struct::ArrayXXr<T_real>>(el_name, data_struct::ArrayXXr<T_real>()));
            scalers_map.at(el_name).resize(count[1], count[2]);

            H5Sselect_hyperslab(counts_dspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            status = _read_h5d<T_real>(counts_dset_id, memoryspace_id, counts_dspace_id, H5P_DEFAULT, (void*)(scalers_map.at(el_name).data()));
        }

        delete[]dims_out;

        _close_h5_objects(close_map);

        //nan to 0.f
        for (auto& itr : scalers_map)
        {
            itr.second = itr.second.unaryExpr([](T_real v) { return std::isfinite(v) ? v : 0.0f; });
        }

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_scan_info_analyzed_h5(std::string path, data_struct::Detector<T_real>* detector, data_struct::Scan_Info<T_real> &scan_info)
    {

        hid_t file_id = -1;
        hid_t x_axis_id = -1;
        hid_t y_axis_id = -1;
        hid_t amp_id = -1;
        hsize_t xdims[1] = { 0 };
        hsize_t ydims[1] = { 0 };
        herr_t status;
        hid_t amp_space = -1;
        char tmp_char[256] = { 0 };
        hid_t type;
        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
        {
            logE << "Failed to open " << path << "\n";
            return false;
        }

        if (false == _open_h5_object(x_axis_id, H5O_DATASET, close_map, "/MAPS/Scan/x_axis", file_id))
            return false;

        if (false == _open_h5_object(y_axis_id, H5O_DATASET, close_map, "/MAPS/Scan/y_axis", file_id))
            return false;

        hid_t x_space = H5Dget_space(x_axis_id);
        close_map.push({ x_space, H5O_DATASPACE });
        hid_t y_space = H5Dget_space(y_axis_id);
        close_map.push({ y_space, H5O_DATASPACE });
        status = H5Sget_simple_extent_dims(x_space, &xdims[0], nullptr);
        if (status < 0)
        {
            _close_h5_objects(close_map);
            logE << "getting dataset X dims" << "\n";
            return false;
        }

        status = H5Sget_simple_extent_dims(y_space, &ydims[0], nullptr);
        if (status < 0)
        {
            _close_h5_objects(close_map);
            logE << "getting dataset Y dims" << "\n";
            return false;
        }

        scan_info.meta_info.x_axis.resize(xdims[0]);
        scan_info.meta_info.y_axis.resize(ydims[0]);

        status = _read_h5d<T_real>(x_axis_id, x_space, x_space, H5P_DEFAULT, scan_info.meta_info.x_axis.data());
        if (status < 0)
        {
            logW << "Failed reading x_axis\n";
        }

        status = _read_h5d<T_real>(y_axis_id, y_space, y_space, H5P_DEFAULT, scan_info.meta_info.y_axis.data());
        if (status < 0)
        {
            logW << "Failed reading x_axis\n";
        }

        // US_AMPS_NUM
        if (false == _open_h5_object(amp_id, H5O_DATASET, close_map, "/MAPS/Scalers/us_amp_num", file_id))
            return false;

        amp_space = H5Dget_space(amp_id);
        close_map.push({ amp_space, H5O_DATASPACE });

        status = _read_h5d<T_real>(amp_id, amp_space, amp_space, H5P_DEFAULT, (void*)&(detector->fit_params_override_dict.us_amp_sens_num));
        if (status < 0)
        {
            logW << "Failed reading us_amp_sens_num\n";
        }

        // US_AMPS_UNIT
        if (false == _open_h5_object(amp_id, H5O_DATASET, close_map, "/MAPS/Scalers/us_amp_unit", file_id))
            return false;

        amp_space = H5Dget_space(amp_id);
        close_map.push({ amp_space, H5O_DATASPACE });

        type = H5Tget_native_type(H5Dget_type(amp_id), H5T_DIR_ASCEND);
        status = H5Dread(amp_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)tmp_char);
        if (status < 0)
        {
            logW << "Failed reading us_amp_sens_unit\n";
        }
        detector->fit_params_override_dict.us_amp_sens_unit = std::string(tmp_char, 255);
        detector->fit_params_override_dict.us_amp_sens_unit.erase(std::remove_if(detector->fit_params_override_dict.us_amp_sens_unit.begin(), detector->fit_params_override_dict.us_amp_sens_unit.end(), ::isspace), detector->fit_params_override_dict.us_amp_sens_unit.end());


        // DS_AMPS_NUM
        if (false == _open_h5_object(amp_id, H5O_DATASET, close_map, "/MAPS/Scalers/ds_amp_num", file_id))
            return false;

        amp_space = H5Dget_space(amp_id);
        close_map.push({ amp_space, H5O_DATASPACE });

        status = _read_h5d<T_real>(amp_id, amp_space, amp_space, H5P_DEFAULT, (void*)&(detector->fit_params_override_dict.ds_amp_sens_num));
        if (status < 0)
        {
            logW << "Failed reading ds_amp_sens_num\n";
        }

        // DS_AMPS_UNIT
        if (false == _open_h5_object(amp_id, H5O_DATASET, close_map, "/MAPS/Scalers/ds_amp_unit", file_id))
            return false;

        amp_space = H5Dget_space(amp_id);
        close_map.push({ amp_space, H5O_DATASPACE });

        type = H5Tget_native_type(H5Dget_type(amp_id), H5T_DIR_ASCEND);
        status = H5Dread(amp_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)tmp_char);
        if (status < 0)
        {
            logW << "Failed reading us_amp_sens_unit\n";
        }
        detector->fit_params_override_dict.ds_amp_sens_unit = std::string(tmp_char,255);
        detector->fit_params_override_dict.ds_amp_sens_unit.erase(std::remove_if(detector->fit_params_override_dict.ds_amp_sens_unit.begin(), detector->fit_params_override_dict.ds_amp_sens_unit.end(), ::isspace), detector->fit_params_override_dict.ds_amp_sens_unit.end());



        _close_h5_objects(close_map);
        
        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_quantification_scalers_gsecars(std::string path, data_struct::Params_Override<T_real> *override_values)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        int status = 0;
        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        logI << path << "\n";

        T_real us_ic, ds_ic;
        hid_t file_id, ds_ic_id, us_ic_id;
        hsize_t offset3[3] = { 0,0,0 };
        hsize_t count3[3] = { 1,1,1 };
        hsize_t dims3[3] = { 1,1,1 };
        hsize_t offset[1] = { 1 };
        hsize_t count[1] = { 1 };
        hid_t readwrite_space = H5Screate_simple(1, &count[0], &count[0]);

        GSE_CARS_SAVE_VER version = GSE_CARS_SAVE_VER::UNKNOWN;

        if (false == _open_h5_object(file_id, H5O_FILE, close_map, path, -1))
            return false;

        if (false == _open_h5_object(ds_ic_id, H5O_DATASET, close_map, "/xrmmap/scalars/I1_raw", file_id, false, false))
        {
            if (false == _open_h5_object(ds_ic_id, H5O_DATASET, close_map, "/xrfmap/roimap/det_name", file_id))
            {
                return false;
            }
            version = GSE_CARS_SAVE_VER::XRFMAP;
        }
        else
        {
            version = GSE_CARS_SAVE_VER::XRMMAP;
        }

        if (false == _open_h5_object(us_ic_id, H5O_DATASET, close_map, "/xrmmap/scalars/I0_raw", file_id, false, false))
        {
            if (false == _open_h5_object(us_ic_id, H5O_DATASET, close_map, "/xrfmap/roimap/det_raw", file_id))
            {
                return false;
            }
        }

        //read in scaler
        hid_t d_space = H5Dget_space(us_ic_id);
        close_map.push({ d_space, H5O_DATASPACE });
        status = H5Sget_simple_extent_dims(d_space, &dims3[0], nullptr);
        if (status < 0)
        {
            _close_h5_objects(close_map);
            logE << "getting dataset dims" << "\n";
            return false;
        }

        if (version == GSE_CARS_SAVE_VER::XRMMAP)
        {
            for (offset3[0] = 0; offset3[0] < dims3[0]; offset3[0]++)
            {
                for (offset3[1] = 0; offset3[1] < dims3[1]; offset3[1]++)
                {
                    H5Sselect_hyperslab(d_space, H5S_SELECT_SET, offset3, nullptr, count3, nullptr);

                    status = _read_h5d<T_real>(us_ic_id, readwrite_space, d_space, H5P_DEFAULT, (void*)&us_ic);
                    if (status > -1)
                    {
                        override_values->US_IC += (us_ic);
                    }

                    status = _read_h5d<T_real>(ds_ic_id, readwrite_space, d_space, H5P_DEFAULT, (void*)&ds_ic);
                    if (status > -1)
                    {
                        override_values->DS_IC += (ds_ic);
                    }
                }
            }
        }
        else if (version == GSE_CARS_SAVE_VER::XRFMAP)
        {
            int usIDX = -1;
            int dsIDX = -1;

            hid_t dtype = H5Tcopy(H5T_C_S1);
            H5Tset_size(dtype, 255);

            hid_t name_space = H5Dget_space(ds_ic_id);
            status = H5Sget_simple_extent_dims(name_space, &count[0], nullptr);
            int name_amt = count[0];
            count[0] = 1;

            for (offset[0] = 0; offset[0] < name_amt; offset[0]++)
            {
                char tmp_char[256] = { 0 };
                //read the detector names and find I0 and I1 indicies
                H5Sselect_hyperslab(name_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);
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
                    if (usIDX > -1 && dsIDX > -1)
                    {
                        break;
                    }
                }
            }

            if (usIDX == -1 || dsIDX == -1)
            {
                logW << "Could not find up stream ion or down stream ion chamber indicies\n";
            }

            for (offset3[0] = 0; offset3[0] < dims3[0]; offset3[0]++)
            {
                for (offset3[1] = 0; offset3[1] < dims3[1]; offset3[1]++)
                {
                    offset3[2] = usIDX;
                    H5Sselect_hyperslab(d_space, H5S_SELECT_SET, offset3, nullptr, count3, nullptr);

                    status = _read_h5d<T_real>(us_ic_id, readwrite_space, d_space, H5P_DEFAULT, (void*)&us_ic);
                    if (status > -1)
                    {
                        override_values->US_IC += (us_ic);
                    }

                    offset3[2] = dsIDX;
                    H5Sselect_hyperslab(d_space, H5S_SELECT_SET, offset3, nullptr, count3, nullptr);
                    status = _read_h5d<T_real>(us_ic_id, readwrite_space, d_space, H5P_DEFAULT, (void*)&ds_ic);
                    if (status > -1)
                    {
                        override_values->DS_IC += (ds_ic);
                    }
                }
            }
        }

        _close_h5_objects(close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool load_quantification_scalers_BNL(std::string path, data_struct::Params_Override<T_real>* override_values)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        logI << path << "\n";

        hid_t status;
        hid_t    file_id, src_maps_grp_id;

        hsize_t* det_dims_in = nullptr;
        hsize_t* val_dims_in = nullptr;;
        hsize_t scaler_offset[3] = { 0,0,0 };
        hsize_t single_offset[1] = { 0 };
        hsize_t mem_offset[2] = { 0,0 };
        hsize_t mem_count[2] = { 1,1 };
        double* buffer = nullptr;

        if (override_values == nullptr)
        {
            return false;
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

        //hid_t name_type = H5Dget_type(scaler_name_id);
        //close_map.push({ name_type, H5O_DATATYPE });
        hid_t scaler_type = H5Dget_type(scaler_val_id);
        close_map.push({ scaler_type, H5O_DATATYPE });
        hid_t scaler_val_space = H5Dget_space(scaler_val_id);
        close_map.push({ scaler_val_space, H5O_DATASPACE });
        hid_t val_rank = H5Sget_simple_extent_ndims(scaler_val_space);
        val_dims_in = new hsize_t[val_rank];
        H5Sget_simple_extent_dims(scaler_val_space, &val_dims_in[0], NULL);

        mem_count[0] = val_dims_in[0];
        mem_count[1] = val_dims_in[1];
        //names

        buffer = new double[val_dims_in[0] * val_dims_in[1]];
        hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);
        close_map.push({ mem_space, H5O_DATASPACE });
        size_t scaler_cnt = val_dims_in[2];
        val_dims_in[2] = 1;

        float famt = val_dims_in[0] * val_dims_in[1];

        char* tmp_char_arr[255];
        hid_t dtype2 = H5Tcopy(H5T_C_S1);
        H5Tset_size(dtype2, H5T_VARIABLE);
        _global_close_map.push({ dtype2, H5O_DATATYPE });
        hid_t rstatus = H5Dread(scaler_name_id, dtype2, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)tmp_char_arr);

        for (hsize_t i = 0; i < scaler_cnt; i++)
        {
            single_offset[0] = i;
            scaler_offset[2] = i;

            H5Sselect_hyperslab(scaler_val_space, H5S_SELECT_SET, scaler_offset, NULL, val_dims_in, NULL);
            if (rstatus > -1)
            {
                std::string read_name = std::string(tmp_char_arr[i]);
                read_name.erase(std::remove_if(read_name.begin(), read_name.end(), ::isspace), read_name.end());
                read_name.erase(std::find(read_name.begin(), read_name.end(), '\0'), read_name.end());
                std::string out_label = "";
                std::string beamline = "";
                bool tmpb;
                delete tmp_char_arr[i];
                if (data_struct::Scaler_Lookup::inst()->search_pv(read_name, out_label, tmpb, beamline))
                {
                    if (out_label == STR_DS_IC)
                    {
                        status = H5Dread(scaler_val_id, scaler_type, mem_space, scaler_val_space, H5P_DEFAULT, (void*)buffer);
                        override_values->DS_IC = 0.0;
                        if (status > -1)
                        {
                            for (hsize_t x = 0; x < val_dims_in[0] * val_dims_in[1]; x++)
                            {
                                override_values->DS_IC += (T_real)buffer[x];
                            }
                            override_values->DS_IC /= famt;
                        }
                        //break;
                    }
                    else if (out_label == STR_US_IC)
                    {
                        status = H5Dread(scaler_val_id, scaler_type, mem_space, scaler_val_space, H5P_DEFAULT, (void*)buffer);
                        override_values->US_IC = 0.0;
                        if (status > -1)
                        {
                            for (hsize_t x = 0; x < val_dims_in[0] * val_dims_in[1]; x++)
                            {
                                override_values->US_IC += (T_real)buffer[x];
                            }
                            override_values->US_IC /= famt;
                        }
                        //break;
                    }
                    else if (out_label == STR_SR_CURRENT)
                    {
                        status = H5Dread(scaler_val_id, scaler_type, mem_space, scaler_val_space, H5P_DEFAULT, (void*)buffer);
                        override_values->sr_current = 0.0;
                        if (status > -1)
                        {
                            for (hsize_t x = 0; x < val_dims_in[0] * val_dims_in[1]; x++)
                            {
                                override_values->sr_current += (T_real)buffer[x];
                            }
                            override_values->sr_current /= famt;
                        }
                        break;
                    }
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

    bool generate_avg(std::string avg_filename, std::vector<std::string> files_to_avg);

    bool generate_stream_dataset(std::string dataset_directory,
                                 std::string dataset_name,
                                 int detector_num,
                                 size_t height,
                                 size_t width);

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool get_scalers_and_metadata_emd(std::string path, data_struct::Scan_Info<T_real>* scan_info)
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
            if (nobj > 0)
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
                            close_map.push({ scalers_type, H5O_DATATYPE });

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

    template<typename T_real>
    bool get_scalers_and_metadata_confocal(std::string path, data_struct::Scan_Info<T_real>* scan_info)
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
            T_real* buffer = nullptr;
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
                    close_map.push({ scalers_type, H5O_DATATYPE });

                    if (first_save)
                    {
                        H5Sget_simple_extent_dims(scaler_space, &scalers_count[0], NULL);
                        scan_info->meta_info.requested_cols = scalers_count[0];
                        scan_info->meta_info.requested_rows = scalers_count[1];
                        first_save = false;
                    }
                    data_struct::Scaler_Map<T_real> sm;
                    sm.name = std::string(str_dset_name, len);
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
            close_map.push({ scalers_type, H5O_DATATYPE });
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
                data_struct::Scaler_Map<T_real> sm;
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
                    itr.name = std::string(detector_names[i], strlen(detector_names[i]));
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

    template<typename T_real>
    bool get_scalers_and_metadata_gsecars(std::string path, data_struct::Scan_Info<T_real>* scan_info)
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
                        close_map.push({ name_type, H5O_DATATYPE });
                        hid_t value_type = H5Dget_type(value_id);
                        close_map.push({ value_type, H5O_DATATYPE });
                        hid_t value_space = H5Dget_space(value_id);
                        close_map.push({ value_space, H5O_DATASPACE });
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
                data_struct::Scaler_Map sm;
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

    template<typename T_real>
    bool get_scalers_and_metadata_bnl(std::string path, data_struct::Scan_Info<T_real>* scan_info)
    {
        std::lock_guard<std::mutex> lock(_mutex);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        logI << path << "\n";
        hid_t    file_id, src_maps_grp_id;
        hid_t scaler_name_id, scaler_val_id, scaler_grp_id;
        hid_t tmp_id, status;
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
        close_map.push({ name_type, H5O_DATATYPE });
        hid_t scaler_type = H5Dget_type(scaler_val_id);
        close_map.push({ scaler_type, H5O_DATATYPE });
        hid_t scaler_val_space = H5Dget_space(scaler_val_id);
        close_map.push({ scaler_val_space, H5O_DATASPACE });
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

        hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);
        close_map.push({ mem_space, H5O_DATASPACE });
        size_t scaler_cnt = val_dims_in[2];
        val_dims_in[2] = 1;

        char* tmp_char_arr[255];
        hid_t dtype2 = H5Tcopy(H5T_C_S1);
        H5Tset_size(dtype2, H5T_VARIABLE);
        _global_close_map.push({ dtype2, H5O_DATATYPE });

        H5Sselect_hyperslab(name_space, H5S_SELECT_SET, single_offset, nullptr, single_count, nullptr);
        status = H5Dread(scaler_name_id, dtype2, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)tmp_char_arr);

        for (hsize_t i = 0; i < scaler_cnt; i++)
        {
            data_struct::Scaler_Map<T_real> scaler_map;
            scaler_map.values.resize(val_dims_in[1], val_dims_in[0]);
            scaler_map.unit = "cts";

            single_offset[0] = i;
            scaler_offset[2] = i;

            H5Sselect_hyperslab(scaler_val_space, H5S_SELECT_SET, scaler_offset, NULL, val_dims_in, NULL);
            if (status > -1)
            {
                scaler_map.name = std::string(tmp_char_arr[i]);
                scaler_map.name.erase(std::remove_if(scaler_map.name.begin(), scaler_map.name.end(), ::isspace), scaler_map.name.end());
                scaler_map.name.erase(std::find(scaler_map.name.begin(), scaler_map.name.end(), '\0'), scaler_map.name.end());
                delete tmp_char_arr[i];
            }
            status = _read_h5d<T_real>(scaler_val_id, mem_space, scaler_val_space, H5P_DEFAULT, scaler_map.values.data());

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

        if (_open_h5_object(tmp_id, H5O_GROUP, close_map, "scan_metadata", src_maps_grp_id))
        {

            // load attributes from this folder
            int na = H5Aget_num_attrs(tmp_id);

            double* ddata = new double[10];
            long* ldata = new long[10];

            for (int i = 0; i < na; i++)
            {
                char* adata[25];
                for (int a = 0; a < 25; a++)
                {
                    adata[a] = nullptr;
                }
                hid_t aid = H5Aopen_idx(tmp_id, (unsigned int)i);

                char buf[1000];
                ssize_t len = H5Aget_name(aid, 1000, buf);
                data_struct::Extra_PV e_pv;
                e_pv.name = std::string(buf, len);

                hid_t atype = H5Aget_type(aid);
                hid_t ntype = H5Tget_native_type(atype, H5T_DIR_ASCEND);
                //hid_t aspece = H5Aget_space(aid);
                H5T_class_t type_class = H5Tget_class(atype);
                //if (type_class == H5T_STRING)
                if (H5Tis_variable_str(atype) > 0)
                {
                    if (H5Aread(aid, ntype, &adata) > -1)
                    {
                        e_pv.value = std::string(adata[0]);
                    }
                    for (int q = 0; q < 25; q++)
                    {
                        if (adata[q] != nullptr)
                        {
                            delete adata[q];
                            adata[q] = nullptr;
                        }
                    }
                }
                else if (type_class == H5T_ARRAY)
                {
                    //
                }
                else if (type_class == H5T_FLOAT)
                {
                    if (H5Aread(aid, H5T_NATIVE_DOUBLE, (void*)ddata) > -1)
                    {
                        e_pv.value = std::to_string(ddata[0]);
                    }
                }
                else if (type_class == H5T_INTEGER)
                {
                    if (H5Aread(aid, H5T_NATIVE_LONG, (void*)ldata) > -1)
                    {
                        e_pv.value = std::to_string(ldata[0]);
                    }
                }


                scan_info->extra_pvs.push_back(e_pv);

                H5Tclose(atype);
                H5Aclose(aid);
            }
            delete[] ddata;
            delete[] ldata;
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
    /*
    template<typename T_real>
    bool save_stream_row(size_t d_hash, size_t detector_num, size_t row, std::vector< data_struct::Spectra<T_real>* >  *spectra_row)
    {
        return false;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_itegrade_spectra(data_struct::Spectra<T_real> * spectra)
    {
        return false;
    }
    
    //-----------------------------------------------------------------------------

    bool close_dataset(size_t d_hash);
    */
    bool start_save_seq(const std::string filename, bool force_new_file=false, bool open_file_only=false);

    bool start_save_seq(bool force_new_file=false){ return start_save_seq(_cur_filename, force_new_file, false);}

    void set_filename(std::string fname) {_cur_filename = fname;}

    const std::string& get_filename() { return _cur_filename; }

    bool polar_copy_raw(const std::string filename);

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_spectra_volume(const std::string path, data_struct::Spectra_Volume<T_real>* spectra_volume, const std::string &scan_type, size_t row_idx_start=0, int row_idx_end=-1, size_t col_idx_start=0, int col_idx_end=-1)
    {
        std::lock_guard<std::mutex> lock(_mutex);



        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        hid_t    dset_id, spec_grp_id, int_spec_grp_id, dataspace_id, memoryspace_id, memoryspace_time_id, maps_grp_id;
        hid_t dataspace_rt_id, dataspace_lt_id, dataspace_incr_id, dataspace_ocr_id;
        hid_t   dset_rt_id, dset_lt_id, incnt_dset_id, outcnt_dset_id;
        herr_t status = 0;

        hsize_t chunk_dims[3] = { 1,1,1 };
        hsize_t chunk_dims_times[2] = { 1,1 };
        hsize_t dims_out[3] = { 1,1,1 };
        hsize_t maxdims[3] = { H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED };
        hsize_t offset[3] = { 0,0,0 };
        hsize_t count[3] = { 1,1,1 };
        hsize_t dims_time_out[2] = { 0,0 };
        hsize_t offset_time[2] = { 0,0 };
        hsize_t count_time[2] = { 0,0 };
        hsize_t tmp_dims[3] = { 0,0,0 };

        if (row_idx_end < (int)row_idx_start || (size_t)row_idx_end > spectra_volume->rows() - 1)
        {
            row_idx_end = spectra_volume->rows();
        }
        if (col_idx_end < (int)col_idx_start || (size_t)col_idx_end > spectra_volume->cols() - 1)
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
        chunk_dims_times[0] = spectra_volume->rows();
        chunk_dims_times[1] = 1;

        offset_time[0] = 0;
        offset_time[1] = 0;
        count_time[0] = 1;
        count_time[1] = 1;

        _create_memory_space(3, count, memoryspace_id);
        _create_memory_space(2, count_time, memoryspace_time_id);

        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

        // open /MAPS
        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }

        // open /MAPS/Spectra
        if (false == _open_or_create_group(STR_SPECTRA, maps_grp_id, spec_grp_id))
        {
            return false;
        }

        // try to open mca dataset and expand before creating 
        if (false == _open_h5_dataset<T_real>(path, spec_grp_id, 3, dims_out, chunk_dims, dset_id, dataspace_id))
        {
            logE << "Error creating " << path << "\n";
            return false;
        }

        if (false == _open_h5_dataset<T_real>(STR_ELAPSED_REAL_TIME, spec_grp_id, 2, dims_time_out, chunk_dims_times, dset_rt_id, dataspace_rt_id))
        {
            logE << "Error creating " << path << "\n";
            return false;
        }
        if (false == _open_h5_dataset<T_real>(STR_ELAPSED_LIVE_TIME, spec_grp_id, 2, dims_time_out, chunk_dims_times, dset_lt_id, dataspace_lt_id))
        {
            logE << "Error creating " << path << "\n";
            return false;
        }
        if (false == _open_h5_dataset<T_real>(STR_INPUT_COUNTS, spec_grp_id, 2, dims_time_out, chunk_dims_times, incnt_dset_id, dataspace_incr_id))
        {
            logE << "Error creating " << path << "\n";
            return false;
        }
        if (false == _open_h5_dataset<T_real>(STR_OUTPUT_COUNTS, spec_grp_id, 2, dims_time_out, chunk_dims_times, outcnt_dset_id, dataspace_ocr_id))
        {
            logE << "Error creating " << path << "\n";
            return false;
        }

        H5Sselect_hyperslab(memoryspace_time_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

        T_real real_time;
        T_real life_time;
        T_real in_cnt;
        T_real out_cnt;
        for (size_t row = row_idx_start; row < (size_t)row_idx_end; row++)
        {
            offset[1] = row;
            offset_time[0] = row;
            for (size_t col = col_idx_start; col < (size_t)col_idx_end; col++)
            {
                const data_struct::Spectra<T_real>* spectra = &((*spectra_volume)[row][col]);
                offset[2] = col;
                offset_time[1] = col;
                H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);
                if (status < 0)
                {
                    logE << " H5Dwrite failed to write spectra\n";
                }

                H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_incr_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
                H5Sselect_hyperslab(dataspace_ocr_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

                real_time = spectra->elapsed_realtime();
                life_time = spectra->elapsed_livetime();
                in_cnt = spectra->input_counts();
                out_cnt = spectra->output_counts();
                status = _write_h5d<T_real>(dset_rt_id, memoryspace_time_id, dataspace_rt_id, H5P_DEFAULT, (void*)&real_time);
                if (status < 0)
                {
                    logE << " H5Dwrite failed to write " << STR_ELAPSED_REAL_TIME << "\n";
                }
                status = _write_h5d<T_real>(dset_lt_id, memoryspace_time_id, dataspace_lt_id, H5P_DEFAULT, (void*)&life_time);
                if (status < 0)
                {
                    logE << " H5Dwrite failed to write " << STR_ELAPSED_LIVE_TIME << "\n";
                }
                status = _write_h5d<T_real>(incnt_dset_id, memoryspace_time_id, dataspace_incr_id, H5P_DEFAULT, (void*)&in_cnt);
                if (status < 0)
                {
                    logE << " H5Dwrite failed to write " << STR_INPUT_COUNTS << "\n";
                }
                status = _write_h5d<T_real>(outcnt_dset_id, memoryspace_time_id, dataspace_ocr_id, H5P_DEFAULT, (void*)&out_cnt);
                if (status < 0)
                {
                    logE << " H5Dwrite failed to write " << STR_OUTPUT_COUNTS << "\n";
                }
            }
        }

        if (false == _open_or_create_group(STR_INT_SPEC, spec_grp_id, int_spec_grp_id))
        {
            return false;
        }

        //save integrated spectra
        data_struct::Spectra<T_real> spectra = spectra_volume->integrate();
        count[0] = spectra.size();
        _create_memory_space(1, count, memoryspace_id);
        if (false == _open_h5_dataset<T_real>(STR_SPECTRA, int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        offset[0] = 0;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&spectra[0]);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write " << STR_INT_SPEC << "/" << STR_SPECTRA << "\n";
        }

        if(scan_type == STR_SCAN_TYPE_POLAR_XANES)
        {
             //save POLAR left and right polarization integrated spectra
            data_struct::Spectra<T_real> lhcp_spectra;
            data_struct::Spectra<T_real> rhcp_spectra;
            if(spectra_volume->integrate_polar(lhcp_spectra, rhcp_spectra))
            {
                hid_t lhcp_dset_id, rhcp_dset_id;
                hid_t lhcp_space_id, rhcp_space_id;
                count[0] = lhcp_spectra.size();
                _create_memory_space(1, count, memoryspace_id);
                bool lhcp_bool = _open_h5_dataset<T_real>(STR_LHCP_SPECTRA, int_spec_grp_id, 1, count, count, lhcp_dset_id, lhcp_space_id);
                bool rhcp_bool = _open_h5_dataset<T_real>(STR_RHCP_SPECTRA, int_spec_grp_id, 1, count, count, rhcp_dset_id, rhcp_space_id);
                if (lhcp_bool && rhcp_bool)
                {
                    offset[0] = 0;
                    H5Sselect_hyperslab(lhcp_space_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                    H5Sselect_hyperslab(rhcp_space_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                    status = _write_h5d<T_real>(lhcp_dset_id, memoryspace_id, lhcp_space_id, H5P_DEFAULT, (void*)&lhcp_spectra[0]);
                    if (status < 0)
                    {
                        logE << " H5Dwrite failed to write left polarity int spectra" << STR_INT_SPEC << "/" << STR_LHCP_SPECTRA << "\n";
                    }
                    status = _write_h5d<T_real>(rhcp_dset_id, memoryspace_id, rhcp_space_id, H5P_DEFAULT, (void*)&rhcp_spectra[0]);
                    if (status < 0)
                    {
                        logE << " H5Dwrite failed to write right polarity int spectra " << STR_INT_SPEC << "/" << STR_RHCP_SPECTRA << "\n";
                    }
                }
            }
        }


        //save real_time
        count[0] = 1;
        T_real save_val = spectra.elapsed_realtime();
        _create_memory_space(1, count, memoryspace_id);
        if (false == _open_h5_dataset<T_real>(STR_ELAPSED_REAL_TIME, int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write " << STR_INT_SPEC << "/" << STR_ELAPSED_REAL_TIME << "\n";
        }

        //save life_time
        save_val = spectra.elapsed_livetime();
        if (false == _open_h5_dataset<T_real>(STR_ELAPSED_LIVE_TIME, int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write " << STR_INT_SPEC << "/" << STR_ELAPSED_LIVE_TIME << "\n";
        }

        //save input_counts
        save_val = spectra.input_counts();
        if (false == _open_h5_dataset<T_real>(STR_INPUT_COUNTS, int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write " << STR_INT_SPEC << "/" << STR_INPUT_COUNTS << "\n";
        }

        //save output_counts
        save_val = spectra.output_counts();
        if (false == _open_h5_dataset<T_real>(STR_OUTPUT_COUNTS, int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write " << STR_INT_SPEC << "/" << STR_OUTPUT_COUNTS << "\n";
        }

        //save file version
        save_val = HDF5_SAVE_VERSION;
        if (false == _open_h5_dataset<T_real>(STR_VERSION, maps_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write " << STR_MAPS << "/" << STR_VERSION << "\n";
        }
        _close_h5_objects(_global_close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        return true;

    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_energy_calib(int spectra_size, T_real energy_offset, T_real energy_slope, T_real energy_quad)
    {

        std::lock_guard<std::mutex> lock(_mutex);
        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        hid_t    dset_id, spec_grp_id, dataspace_id, memoryspace_id, maps_grp_id;
        herr_t status = 0;
        hsize_t offset[1] = { 0 };
        hsize_t count[1] = { 1 };

        data_struct::Range energy_range(0, spectra_size - 1);
        const data_struct::ArrayTr<T_real> ev = data_struct::generate_energy_array(energy_range, energy_offset, energy_slope, energy_quad);

        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }

        if (false == _open_or_create_group(STR_SPECTRA, maps_grp_id, spec_grp_id))
        {
            return false;
        }

        count[0] = ev.size();
        _create_memory_space(1, count, memoryspace_id);

        if (false == _open_h5_dataset<T_real>(STR_ENERGY, spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, ev.data());
        if (status < 0)
        {
            logE << " H5Dwrite failed to write " << STR_ENERGY << "\n";
        }

        // save energy calibration
        count[0] = 3;
        if (false == _open_h5_dataset<T_real>(STR_ENERGY_CALIB, spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        count[0] = 1;
        _create_memory_space(1, count, memoryspace_id);
        offset[0] = 0;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&energy_offset);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write energy offset\n";
        }
        offset[0] = 1;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&energy_slope);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write energy slope\n";
        }
        offset[0] = 2;
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&energy_quad);
        if (status < 0)
        {
            logE << " H5Dwrite failed to write energy quad\n";
        }


        _close_h5_objects(_global_close_map);

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_element_fits(const std::string path,
        const data_struct::Fit_Count_Dict<T_real>* const element_counts,
        [[maybe_unused]] size_t row_idx_start=0,
        [[maybe_unused]] int row_idx_end=-1,
        [[maybe_unused]] size_t col_idx_start=0,
        [[maybe_unused]] int col_idx_end=-1)
    {
        std::lock_guard<std::mutex> lock(_mutex);



        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        hid_t   dset_id, dset_ch_id, dset_un_id;
        hid_t   memoryspace, dataspace_id, dataspace_ch_id, dataspace_un_id, dataspace_ch_off_id;
        hid_t   filetype, memtype;
        herr_t  status;
        hid_t   xrf_grp_id, fit_grp_id, maps_grp_id;

        dset_id = -1;
        dset_ch_id = -1;
        hsize_t dims_out[3] = { 0, 0, 0 };
        hsize_t offset[1] = { 0 };
        hsize_t offset2[1] = { 0 };
        hsize_t offset_3d[3] = { 0, 0, 0 };
        hsize_t count[1] = { 1 };
        hsize_t count_3d[3] = { 1, 1, 1 };
        hsize_t chunk_dims[3];
        
        if (element_counts != nullptr)
        {
            if (element_counts->size() > 0)
            {
                auto iter = element_counts->begin();
                dims_out[1] = iter->second.rows();
                dims_out[2] = iter->second.cols();
            }
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

        _create_memory_space(1, count_3d, dataspace_ch_off_id);
        _create_memory_space(3, count_3d, memoryspace);

        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }

        if (false == _open_or_create_group(STR_XRF_ANALYZED, maps_grp_id, xrf_grp_id))
        {
            return false;
        }

        if (false == _open_or_create_group(path, xrf_grp_id, fit_grp_id))
        {
            return false;
        }

        if (false == _open_h5_dataset<T_real>(STR_COUNTS_PER_SEC, fit_grp_id, 3, dims_out, dims_out, dset_id, dataspace_id))
        {
            return false;
        }

        //filetype = H5Tcopy (H5T_FORTRAN_S1);
        filetype = H5Tcopy(H5T_C_S1);
        H5Tset_size(filetype, 256);
        memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

        if (false == _open_h5_dataset(STR_CHANNEL_NAMES, filetype, fit_grp_id, 1, dims_out, dims_out, dset_ch_id, dataspace_ch_id))
        {
            return false;
        }

        if (false == _open_h5_dataset(STR_CHANNEL_UNITS, filetype, fit_grp_id, 1, dims_out, dims_out, dset_un_id, dataspace_un_id))
        {
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
            element_lines.push_back(el_name + "_L");
        }
        for (std::string el_name : data_struct::Element_Symbols)
        {
            element_lines.push_back(el_name + "_M");
        }

        //add the rest 
        for (const auto& itr : *element_counts)
        {
            if (std::find(element_lines.begin(), element_lines.end(), itr.first) == element_lines.end())
            {
                element_lines.push_back(itr.first);
            }
        }

        //H5Sselect_hyperslab (memoryspace, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);

        int i = 0;
        //save by element Z order
        //for(const auto& iter : *element_counts)
        std::string units = "cts/s";
        for (std::string el_name : element_lines)
        {
            char tmp_char[256] = { 0 };
            if (element_counts->count(el_name) < 1)
            {
                continue;
            }
            offset[0] = i;
            offset_3d[0] = i;

            H5Sselect_hyperslab(dataspace_ch_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            H5Sselect_hyperslab(dataspace_un_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
            H5Sselect_hyperslab(dataspace_ch_off_id, H5S_SELECT_SET, offset2, nullptr, count, nullptr);

            el_name.copy(tmp_char, 254);

            status = H5Dwrite(dset_ch_id, memtype, dataspace_ch_off_id, dataspace_ch_id, H5P_DEFAULT, (void*)tmp_char);
            if (status < 0)
            {
                logE << " H5Dwrite failed to write " << STR_CHANNEL_NAMES << " at row " << i << "\n";
            }

            for (int z = 0; z < 256; z++)
            {
                tmp_char[z] = '\0';
            }
            if (el_name != STR_NUM_ITR && el_name != STR_RESIDUAL)
            {
                units.copy(tmp_char, 256);
            }
            status = H5Dwrite(dset_un_id, memtype, dataspace_ch_off_id, dataspace_un_id, H5P_DEFAULT, (void*)tmp_char);
            if (status < 0)
            {
                logE << " H5Dwrite failed to write " << STR_CHANNEL_UNITS << " at row " << i << "\n";
            }
            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_3d, nullptr, chunk_dims, nullptr);

            status = _write_h5d<T_real>(dset_id, memoryspace, dataspace_id, H5P_DEFAULT, (void*)element_counts->at(el_name).data());
            if (status < 0)
            {
                logE << " H5Dwrite failed to write " << STR_COUNTS_PER_SEC << " at row " << i << "\n";
            }

            i++;
        }

        _close_h5_objects(_global_close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";


        return true;

    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_fitted_int_spectra(const std::string path, const data_struct::Spectra<T_real>& spectra, const data_struct::Range& range, const data_struct::Spectra<T_real>& background, const size_t save_spectra_size)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        bool ret_val = true;
        hid_t   dset_id, dataspace_id, memoryspace_id;
        herr_t  status;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        std::string dset_name = "/" + STR_MAPS + "/" + STR_XRF_ANALYZED + "/" + path + "/" + STR_FIT_INT_SPEC;
        std::string background_name = "/" + STR_MAPS + "/" + STR_XRF_ANALYZED + "/" + path + "/" + STR_FIT_INT_BACKGROUND;

        hsize_t count[1] = { 1 };
        count[0] = save_spectra_size;

        // resize to the size of collected spectra
        data_struct::ArrayTr<T_real>   save_spectra;
        save_spectra.setZero(save_spectra_size);
        data_struct::ArrayTr<T_real>   save_background;
        save_background.setZero(save_spectra_size);

        int j = 0;
        for (size_t i = range.min; i <= range.max; i++)
        {
            if (std::isfinite(spectra[j]))
            {
                save_spectra[i] = spectra[j];
            }
            if (std::isfinite(background[j]))
            {
                save_background[i] = background[j];
            }
            j++;
        }

        _create_memory_space(1, count, memoryspace_id);

        // save spectra
        if (false == _open_h5_dataset<T_real>(dset_name, _cur_file_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)save_spectra.data());
        if (status < 0)
        {
            logW << "Failed to save " << dset_name << "\n";
            ret_val = false;
        }

        // save background
        if (false == _open_h5_dataset<T_real>(background_name, _cur_file_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)save_background.data());
        if (status < 0)
        {
            logW << "Failed to save " << background_name << "\n";
            ret_val = false;
        }

        _close_h5_objects(_global_close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";
        return ret_val;
    }
	
    //-----------------------------------------------------------------------------

    template<typename T_real>
	bool save_max_10_spectra(const data_struct::Range& range,
							const data_struct::Spectra<T_real>& max_spectra,
							const data_struct::Spectra<T_real>& max_10_spectra,
                            const data_struct::Spectra<T_real>& fit_int_background)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        bool ret_val = true;
        hid_t   dset_id, maps_grp_id, spec_grp_id, int_spec_grp_id, dataspace_id, memoryspace_id;
        herr_t  status;
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();


        hsize_t count[1] = { 1 };
        count[0] = max_spectra.size();

        _create_memory_space(1, count, memoryspace_id);

        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_SPECTRA, maps_grp_id, spec_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_INT_SPEC, spec_grp_id, int_spec_grp_id))
        {
            return false;
        }
        if (false == _open_h5_dataset<T_real>(STR_MAX_CHANNELS_INT_SPEC, int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)max_spectra.data());
        if (status < 0)
        {
            logW << "Failed to save " << STR_MAX_CHANNELS_INT_SPEC << "\n";
            ret_val = false;
        }

        if (false == _open_h5_dataset<T_real>(STR_MAX10_INT_SPEC, int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)max_10_spectra.data());
        if (status < 0)
        {
            logW << "Failed to save " << STR_MAX10_INT_SPEC << "\n";
            ret_val = false;
        }

        _close_h5_objects(_global_close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";
        return ret_val;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_quantification(data_struct::Detector<T_real>* detector)
    {

        std::lock_guard<std::mutex> lock(_mutex);

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        hid_t    dset_id, memoryspace_id, dataspace_id, filetype, dataspace_ch_id, dataspace_un_id, memtype, dset_ch_id, dset_un_id, q_int_spec_grp_id;
        hid_t   memtype_label, filetype_label, q_memoryspace_label_id, q_dataspace_label_id;
        hid_t   count_dataspace_id, count_dset_id, standard_grp_id, calib_grp_id;
        hid_t  memoryspace1_id, memoryspace2_id;

        hid_t q_dataspace_id, q_memoryspace_id, q_dset_id, q_grp_id, q_fit_grp_id, maps_grp_id, scalers_grp_id, xrf_fits_grp_id;
        hid_t dset_labels_id;
        hsize_t offset[3];
        hsize_t count[3];
        herr_t status;

        char unit_char[255] = "cts/s";

        hsize_t q_dims_out[2];

        if (_cur_file_id < 0)
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

        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }

        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
        memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);


        filetype_label = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype_label, 50);
        memtype_label = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype_label, 50);


        //--                        save calibration curve                  --

        if (detector != nullptr)
        {
            if (false == _open_or_create_group(STR_QUANTIFICATION, maps_grp_id, q_grp_id))
            {
                return false;
            }
            if (false == _open_or_create_group(STR_CALIBRATION, q_grp_id, calib_grp_id))
            {
                return false;
            }

            _create_memory_space(1, count, memoryspace1_id);
            //create dataset telling how many standards there are
            count[0] = 1;
            if (false == _open_h5_dataset(STR_NUMBER_OF_STANDARDS, H5T_NATIVE_INT, q_grp_id, 1, count, count, count_dset_id, count_dataspace_id))
            {
                return false;
            }
            size_t quant_size = detector->quantification_standards.size();
            status = H5Dwrite(count_dset_id, H5T_NATIVE_INT, memoryspace1_id, count_dataspace_id, H5P_DEFAULT, (void*)&quant_size);
            if (status < 0)
            {
                logE << "failed to write " << STR_NUMBER_OF_STANDARDS << "\n";
            }
            int standard_idx = 0;

            // ----------------------------------------- start per standard  ------------------------------------------------
            for (const auto& quant_itr : detector->quantification_standards)
            {
                //create group
                std::string standard_group_name = "Standard" + std::to_string(standard_idx);
                if (false == _open_or_create_group(standard_group_name, q_grp_id, standard_grp_id))
                {
                    return false;
                }
                if (false == _open_or_create_group(STR_SCALERS, standard_grp_id, scalers_grp_id))
                {
                    return false;
                }
                if (false == _open_or_create_group(STR_XRF_ANALYZED, standard_grp_id, xrf_fits_grp_id))
                {
                    return false;
                }

                //save quantification_standard element weights
                count[0] = quant_itr.second.element_standard_weights.size();
                _create_memory_space(1, count, memoryspace_id);

                if (false == _open_h5_dataset<T_real>(STR_ELEMENT_WEIGHTS, standard_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                if (false == _open_h5_dataset(STR_ELEMENT_WEIGHTS_NAMES, filetype, standard_grp_id, 1, count, count, dset_ch_id, dataspace_ch_id))
                {
                    return false;
                }

                offset[0] = 0;
                count[0] = 1;
                H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                int offset_idx = 0;
                for (auto itr : quant_itr.second.element_standard_weights)
                {
                    offset[0] = offset_idx;
                    H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                    H5Sselect_hyperslab(dataspace_ch_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

                    status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(itr.second));
                    if (status < 0)
                    {
                        logE << "failed to write " << STR_ELEMENT_WEIGHTS << "\n";
                    }
                    char tmp_char[256] = { 0 };
                    itr.first.copy(tmp_char, 254);
                    status = H5Dwrite(dset_ch_id, memtype, memoryspace_id, dataspace_ch_id, H5P_DEFAULT, (void*)tmp_char);
                    if (status < 0)
                    {
                        logE << "failed to write " << STR_ELEMENT_WEIGHTS_NAMES << "\n";
                    }
                    offset_idx++;
                }

                if (false == _open_or_create_group(STR_INT_SPEC, standard_grp_id, q_int_spec_grp_id))
                {
                    return false;
                }
                //save quantification_standard integrated spectra
                if (quant_itr.second.integrated_spectra.size() > 0)
                {
                    count[0] = quant_itr.second.integrated_spectra.size();
                    _create_memory_space(1, count, memoryspace2_id);

                    if (false == _open_h5_dataset<T_real>(STR_SPECTRA, q_int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
                    {
                        return false;
                    }
                    offset[0] = 0;
                    status = _write_h5d<T_real>(dset_id, memoryspace2_id, dataspace_id, H5P_DEFAULT, quant_itr.second.integrated_spectra.data());
                    if (status < 0)
                    {
                        logE << "failed to write " << STR_INT_SPEC << "\n";
                    }
                }

                //save standard name
                count[0] = 1;
                _create_memory_space(1, count, memoryspace_id);
                if (false == _open_h5_dataset(STR_STANDARD_NAME, filetype, standard_grp_id, 1, count, count, dset_ch_id, dataspace_ch_id))
                {
                    return false;
                }
                char tmp_char[255] = { 0 };
                quant_itr.second.standard_filename.copy(tmp_char, 254);
                status = H5Dwrite(dset_ch_id, memtype, memoryspace_id, dataspace_ch_id, H5P_DEFAULT, (void*)tmp_char);

                //save sr_current
                if (false == _open_h5_dataset<T_real>(STR_SR_CURRENT, scalers_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quant_itr.second.sr_current));
                if (status < 0)
                {
                    logE << "failed to write " << STR_SR_CURRENT << "\n";
                }

                //save us_ic
                if (false == _open_h5_dataset<T_real>(STR_US_IC, scalers_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quant_itr.second.US_IC));
                if (status < 0)
                {
                    logE << "failed to write " << STR_US_IC << "\n";
                }

                //save us_fm
                if (false == _open_h5_dataset<T_real>(STR_US_FM, scalers_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quant_itr.second.US_FM));
                if (status < 0)
                {
                    logE << "failed to write " << STR_US_FM << "\n";
                }

                //save ds_ic
                if (false == _open_h5_dataset<T_real>(STR_DS_IC, scalers_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(quant_itr.second.DS_IC));
                if (status < 0)
                {
                    logE << "failed to write " << STR_DS_IC << "\n";
                }

                //save real_time
                T_real save_val = quant_itr.second.integrated_spectra.elapsed_realtime();
                if (false == _open_h5_dataset<T_real>(STR_ELAPSED_REAL_TIME, q_int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
                if (status < 0)
                {
                    logE << "failed to write " << STR_ELAPSED_REAL_TIME << "\n";
                }

                //save life_time
                save_val = quant_itr.second.integrated_spectra.elapsed_livetime();
                if (false == _open_h5_dataset<T_real>(STR_ELAPSED_LIVE_TIME, q_int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
                if (status < 0)
                {
                    logE << "failed to write " << STR_ELAPSED_LIVE_TIME << "\n";
                }

                //save input counts
                save_val = quant_itr.second.integrated_spectra.input_counts();
                if (false == _open_h5_dataset<T_real>(STR_INPUT_COUNTS, q_int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
                if (status < 0)
                {
                    logE << "failed to write " << STR_INPUT_COUNTS << "\n";
                }

                //save output counts
                save_val = quant_itr.second.integrated_spectra.output_counts();
                if (false == _open_h5_dataset<T_real>(STR_OUTPUT_COUNTS, q_int_spec_grp_id, 1, count, count, dset_id, dataspace_id))
                {
                    return false;
                }
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
                if (status < 0)
                {
                    logE << "failed to write " << STR_OUTPUT_COUNTS << "\n";
                }

                for (const auto& fit_itr : detector->fitting_quant_map)
                {

                    if (false == _open_or_create_group(data_struct::Fitting_Routine_To_Str.at(fit_itr.first), xrf_fits_grp_id, q_fit_grp_id))
                    {
                        return false;
                    }

                    //save quantification_standard counts
                    std::unordered_map<std::string, T_real> element_counts;
                    if (quant_itr.second.element_counts.count(fit_itr.first) > 0)
                    {
                        element_counts = quant_itr.second.element_counts.at(fit_itr.first);
                    }
                    count[0] = element_counts.size();

                    _create_memory_space(1, count, memoryspace_id);

                    if (false == _open_h5_dataset<T_real>(STR_COUNTS_PER_SEC, q_fit_grp_id, 1, count, count, dset_id, dataspace_id))
                    {
                        return false;
                    }
                    if (false == _open_h5_dataset(STR_CHANNEL_NAMES, filetype, q_fit_grp_id, 1, count, count, dset_ch_id, dataspace_ch_id))
                    {
                        return false;
                    }
                    if (false == _open_h5_dataset(STR_CHANNEL_UNITS, filetype, q_fit_grp_id, 1, count, count, dset_un_id, dataspace_un_id))
                    {
                        return false;
                    }

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
                        T_real val = element_counts.at(el_name);

                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                        H5Sselect_hyperslab(dataspace_ch_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                        H5Sselect_hyperslab(dataspace_un_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                        std::memset(&tmp_char[0], 0, 255);
                        el_name.copy(tmp_char, 254);
                        status = H5Dwrite(dset_ch_id, memtype, memoryspace_id, dataspace_ch_id, H5P_DEFAULT, (void*)tmp_char);
                        if (status < 0)
                        {
                            logE << "failed to write " << STR_CHANNEL_NAMES << "\n";
                        }
                        status = H5Dwrite(dset_un_id, memtype, memoryspace_id, dataspace_un_id, H5P_DEFAULT, (void*)unit_char);
                        if (status < 0)
                        {
                            logE << "failed to write " << STR_CHANNEL_UNITS << "\n";
                        }
                        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&val);
                        if (status < 0)
                        {
                            logE << "failed to write " << STR_COUNTS_PER_SEC << "\n";
                        }

                        offset_idx++;
                    }
                }
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

            _create_memory_space(2, count, q_memoryspace_label_id);

            count[0] = 1;
            count[1] = q_dims_out[1];
            count[2] = 0;

            _create_memory_space(2, count, q_memoryspace_id);

            for (const auto& qitr : detector->fitting_quant_map)
            {
                if (false == _open_or_create_group(data_struct::Fitting_Routine_To_Str.at(qitr.first), calib_grp_id, q_fit_grp_id))
                {
                    return false;
                }
                if (false == _open_h5_dataset(STR_CALIB_LABELS, filetype_label, q_fit_grp_id, 2, q_dims_out, q_dims_out, dset_labels_id, q_dataspace_label_id))
                {
                    return false;
                }

                for (const auto& quant_scaler_itr : qitr.second.quant_scaler_map)
                {

                    std::string q_dset_name = "Calibration_Curve_" + quant_scaler_itr.first;

                    if (false == _open_h5_dataset<T_real>(q_dset_name, q_fit_grp_id, 2, q_dims_out, q_dims_out, q_dset_id, q_dataspace_id))
                    {
                        return false;
                    }

                    int j = 0;

                    //for(auto& shell_itr : quant_itr.second)
                    for (const auto& calib_itr : quant_scaler_itr.second.curve_quant_map)
                    {
                        char label[10] = { 0 };
                        //int element_offset = 0;
                        //create dataset for different shell curves
                        std::vector<T_real> calibration_curve;
                        std::vector<std::string> calibration_curve_labels;

                        if (calib_itr.first == data_struct::Electron_Shell::K_SHELL)
                        {
                            j = 0;
                        }
                        else if (calib_itr.first == data_struct::Electron_Shell::L_SHELL)
                        {
                            j = 1;
                        }
                        else if (calib_itr.first == data_struct::Electron_Shell::M_SHELL)
                        {
                            j = 2;
                        }

                        for (const auto& citr : calib_itr.second)
                        {
                            std::string name = citr.name;
                            if (calib_itr.first == data_struct::Electron_Shell::L_SHELL)
                            {
                                name += "_L";
                            }
                            if (calib_itr.first == data_struct::Electron_Shell::M_SHELL)
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

                        H5Sselect_hyperslab(q_dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                        status = _write_h5d<T_real>(q_dset_id, q_memoryspace_id, q_dataspace_id, H5P_DEFAULT, (void*)&calibration_curve[0]);
                        if (status < 0)
                        {
                            logE << "failed to write " << q_dset_name << "\n";
                        }

                        for (size_t k = 0; k < calibration_curve_labels.size(); k++)
                        {
                            memset(label, 0, 10);
                            calibration_curve_labels[k].copy(&label[0], 9);
                            offset[1] = k;
                            count[1] = 1;
                            status = H5Sselect_hyperslab(q_dataspace_label_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                            status = H5Dwrite(dset_labels_id, memtype_label, q_memoryspace_label_id, q_dataspace_label_id, H5P_DEFAULT, (void*)&label[0]);
                            if (status < 0)
                            {
                                logE << "failed to write " << STR_CALIB_LABELS << "\n";
                            }
                        }
                    }

                    bool index_ok = false;
                    bool name_ok = false;
                    bool prop_ok = false;
                    hsize_t p_offset[1] = { 0 };
                    hsize_t p_count[1] = { 0 };
                    hsize_t e_offset[2] = { 0 };
                    hsize_t e_count[2] = { 0 };
                    hid_t e_index_id = -1;
                    hid_t e_name_id = -1;
                    hid_t e_prop_id = -1;
                    hid_t e_index_dataspace_id = -1;
                    hid_t e_name_dataspace_id = -1;
                    hid_t e_prop_dataspace_id = -1;

                    e_count[0] = detector->all_element_quants[qitr.first][quant_scaler_itr.first].size();
                    e_count[1] = 10; // all properties of Element_Quant except for name
                    p_count[0] = 10; // all properties of Element_Quant except for name

                    // save quant generation info
                    std::string e_gen_idx = quant_scaler_itr.first + "_Element_Info_Index";
                    std::string e_gen_name = quant_scaler_itr.first + "_Element_Info_Names";
                    std::string e_gen_props = quant_scaler_itr.first + "_Element_Info_Values";
                    index_ok = _open_h5_dataset(e_gen_idx, filetype_label, q_fit_grp_id, 1, e_count, e_count, e_index_id, e_index_dataspace_id);
                    name_ok = _open_h5_dataset(e_gen_name, filetype_label, q_fit_grp_id, 1, p_count, p_count, e_name_id, e_name_dataspace_id);
                    prop_ok = _open_h5_dataset<T_real>(e_gen_props, q_fit_grp_id, 2, e_count, e_count, e_prop_id, e_prop_dataspace_id);
                    if (index_ok && name_ok && prop_ok)
                    {
                        int c = 0;
                        // proc_type          quant_scaler      element    quant_prop
                        for (auto& element_itr : detector->all_element_quants[qitr.first][quant_scaler_itr.first])
                        {

                            char label[50] = { 0 };
                            element_itr.first.copy(&label[0], 49);

                            e_offset[0] = c;
                            e_count[0] = 1;
                            e_count[1] = 1;
                            status = H5Sselect_hyperslab(e_index_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_index_id, memtype_label, q_memoryspace_label_id, e_index_dataspace_id, H5P_DEFAULT, (void*)&label);
                            if (status < 0)
                            {
                                logE << "failed to write Index name\n";
                            }

                            p_count[0] = 1;

                            // save properties 
                            p_offset[0] = 0;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"weight");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 0;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->weight));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }


                            // 
                            p_offset[0] = 1;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"absorption");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 1;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->absorption));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            // 
                            p_offset[0] = 2;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"transmission_Be");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 2;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->transmission_Be));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            // 
                            p_offset[0] = 3;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"transmission_Ge");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 3;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->transmission_Ge));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            // 
                            p_offset[0] = 4;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"yield");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 4;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->yield));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            // 
                            p_offset[0] = 5;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"transmission_through_Si_detector");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 5;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->transmission_through_Si_detector));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            // 
                            p_offset[0] = 6;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"transmission_through_air");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 6;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->transmission_through_air));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            // 
                            p_offset[0] = 7;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"Z");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 7;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_INT, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->Z));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            // 
                            p_offset[0] = 8;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"e_cal_ratio");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 8;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->e_cal_ratio));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            // 
                            p_offset[0] = 9;
                            status = H5Sselect_hyperslab(e_name_dataspace_id, H5S_SELECT_SET, p_offset, nullptr, p_count, nullptr);
                            status = H5Dwrite(e_name_id, memtype_label, q_memoryspace_label_id, e_name_dataspace_id, H5P_DEFAULT, (void*)"calib_curve_val");
                            if (status < 0)
                            {
                                logE << "failed to write prop name\n";
                            }

                            e_offset[1] = 9;
                            status = H5Sselect_hyperslab(e_prop_dataspace_id, H5S_SELECT_SET, e_offset, nullptr, e_count, nullptr);
                            status = H5Dwrite(e_prop_id, H5T_NATIVE_DOUBLE, q_memoryspace_label_id, e_prop_dataspace_id, H5P_DEFAULT, (void*)&(element_itr.second->calib_curve_val));
                            if (status < 0)
                            {
                                logE << "failed to write prop value\n";
                            }

                            c++;

                        }
                    }
                }
            }
        }

        _close_h5_objects(_global_close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_scan_scalers(data_struct::Scan_Info<T_real>* scan_info,
                           data_struct::Params_Override<T_real> * params_override,
          [[maybe_unused]] size_t row_idx_start=0,
          [[maybe_unused]] int row_idx_end=-1,
          [[maybe_unused]] size_t col_idx_start=0,
          [[maybe_unused]] int col_idx_end=-1)
    {

        std::lock_guard<std::mutex> lock(_mutex);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        hid_t scan_grp_id, maps_grp_id;

        if (scan_info == nullptr)
        {
            logW << "scalers_map == nullptr. Not returning from save_scan_scalers" << "\n";
            return false;
        }

        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        logI << "Saving scalers to hdf5" << "\n";

        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_SCAN, maps_grp_id, scan_grp_id))
        {
            return false;
        }
        
        _save_scan_meta_data(scan_grp_id, &(scan_info->meta_info));

        _save_extras(scan_grp_id, &(scan_info->extra_pvs));

        if (params_override != nullptr)
        {
            // if us_amps is not set in override file, search for it in extra_pvs
            if (params_override->us_amp_sens_num == -1 || params_override->ds_amp_sens_num == -1)
            {
                std::string label = "";
                std::string beamline = "";
                bool is_time_normalized = false;
                for (auto& itr : scan_info->extra_pvs)
                {
                    if (data_struct::Scaler_Lookup::inst()->search_pv(itr.name, label, is_time_normalized, beamline))
                    {
                        if (label == STR_US_AMP_NUM_UPPR)
                        {
                            params_override->us_amp_sens_num = parse_input_real<T_real>(itr.value);
                        }
                        else if (label == STR_US_AMP_UNIT_UPPR)
                        {
                            params_override->us_amp_sens_unit = itr.value;
                        }
                        if (label == STR_DS_AMP_NUM_UPPR)
                        {
                            params_override->ds_amp_sens_num = parse_input_real<T_real>(itr.value);
                        }
                        else if (label == STR_DS_AMP_UNIT_UPPR)
                        {
                            params_override->ds_amp_sens_unit = itr.value;
                        }
                    }
                }
            }

            _save_scalers(maps_grp_id, &(scan_info->scaler_maps), params_override->us_amp_sens_num, params_override->us_amp_sens_unit, params_override->ds_amp_sens_num, params_override->ds_amp_sens_unit);
        }
        else
        {
            _save_scalers(maps_grp_id, &(scan_info->scaler_maps), (T_real)0.0, "0.0", (T_real)0.0, "0.0");
        }
        _close_h5_objects(_global_close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_scan_scalers_confocal(std::string path,
                   [[maybe_unused]] size_t row_idx_start=0,
                   [[maybe_unused]] int row_idx_end=-1,
                   [[maybe_unused]] size_t col_idx_start=0,
                   [[maybe_unused]] int col_idx_end=-1)
    {

        std::lock_guard<std::mutex> lock(_mutex);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        hid_t scan_grp_id, maps_grp_id, scalers_grp_id, status, error;
        hid_t    file_id, src_maps_grp_id;
        hid_t    dataspace_detectors_id, dset_detectors_id, attr_detector_names_id;
        hid_t   xpos_dataspace_id, xpos_id, ypos_dataspace_id, ypos_id;
        hid_t dataspace_id, dataset_id;
        char* detector_names[256];
        int det_rank;
        hsize_t* det_dims_in;
        hid_t units_id;
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
        void* f_data;
        bool confocal_ver_2020 = false;

        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        logI << "Saving scalers to hdf5" << "\n";

        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_SCAN, maps_grp_id, scan_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_SCALERS, maps_grp_id, scalers_grp_id))
        {
            return false;
        }

        file_id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (file_id < 0)
            return false;
        _global_close_map.push({ file_id, H5O_FILE });

        if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, _global_close_map, "2D Scan", file_id))
        {
            return false;
        }

        if (false == _open_h5_object(xpos_id, H5O_DATASET, _global_close_map, "X Positions", src_maps_grp_id))
        {
            return false;
        }
        xpos_dataspace_id = H5Dget_space(xpos_id);
        _global_close_map.push({ xpos_dataspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(ypos_id, H5O_DATASET, _global_close_map, "Y Positions", src_maps_grp_id))
        {
            return false;
        }
        ypos_dataspace_id = H5Dget_space(ypos_id);
        _global_close_map.push({ ypos_dataspace_id, H5O_DATASPACE });

        hid_t xy_type = H5Dget_type(xpos_id);
        _global_close_map.push({ xy_type, H5O_DATATYPE });
        det_rank = H5Sget_simple_extent_ndims(xpos_dataspace_id);
        det_dims_in = new hsize_t[det_rank];
        H5Sget_simple_extent_dims(xpos_dataspace_id, &det_dims_in[0], NULL);

        if (false == _open_h5_dataset(STR_X_AXIS, xy_type, scan_grp_id, 1, &det_dims_in[1], &det_dims_in[1], dataset_id, dataspace_id))
        {
            logE << "Error creating " << STR_X_AXIS << "\n";
            return false;
        }
        x_count[1] = det_dims_in[1];
        f_data = malloc(sizeof(float) * det_dims_in[1]);
        H5Sselect_hyperslab(xpos_dataspace_id, H5S_SELECT_SET, x_offset, NULL, x_count, NULL);
        status = H5Dread(xpos_id, xy_type, dataspace_id, xpos_dataspace_id, H5P_DEFAULT, f_data);
        status = H5Dwrite(dataset_id, xy_type, H5S_ALL, dataspace_id, H5P_DEFAULT, f_data);
        free(f_data);

        if (false == _open_h5_dataset(STR_Y_AXIS, xy_type, scan_grp_id, 1, &det_dims_in[0], &det_dims_in[0], dataset_id, dataspace_id))
        {
            logE << "Error creating " << STR_Y_AXIS << "\n";
            return false;
        }
        y_count[1] = det_dims_in[0];
        f_data = malloc(sizeof(float) * det_dims_in[0]);
        H5Sselect_hyperslab(ypos_dataspace_id, H5S_SELECT_SET, y_offset, NULL, y_count, NULL);
        status = H5Dread(ypos_id, xy_type, dataspace_id, ypos_dataspace_id, H5P_DEFAULT, f_data);
        status = H5Dwrite(dataset_id, xy_type, H5S_ALL, dataspace_id, H5P_DEFAULT, f_data);
        free(f_data);
        delete[] det_dims_in;

        //Save scalers
        if (false == _open_h5_object(dset_detectors_id, H5O_DATASET, _global_close_map, "Detectors", src_maps_grp_id, false, false))
        {
            // try Detectors as group since it changed in 2020
            if (false == _open_h5_object(dset_detectors_id, H5O_GROUP, _global_close_map, "Detectors", src_maps_grp_id))
            {
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
            T_real* buffer = nullptr;
            hsize_t nobj = 0;
            H5Gget_num_objs(dset_detectors_id, &nobj);
            hid_t values_id;
            hid_t value_space;
            hid_t unit_space;
            hid_t mem_space;
            hid_t names_id;
            hid_t name_space;
            hid_t mem_name_space = H5Screate_simple(1, &names_cnt[0], &names_cnt[0]);
            _global_close_map.push({ mem_name_space, H5O_DATASPACE });
            hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
            H5Tset_size(filetype, 256);
            hid_t name_type = H5Tcopy(H5T_C_S1);
            status = H5Tset_size(name_type, 255);

            for (hsize_t i = 0; i < nobj; i++)
            {
                char str_dset_name[2048] = { 0 };
                hid_t dsid;
                hsize_t len = H5Gget_objname_by_idx(dset_detectors_id, i, str_dset_name, 2048);
                if (_open_h5_object(dsid, H5O_DATASET, _global_close_map, str_dset_name, dset_detectors_id))
                {
                    hid_t scaler_space = H5Dget_space(dsid);
                    _global_close_map.push({ scaler_space, H5O_DATASPACE });
                    hid_t scalers_type = H5Dget_type(dsid);
                    _global_close_map.push({ scalers_type, H5O_DATATYPE });
                    names_off[0] = i;
                    if (first_save)
                    {
                        H5Sget_simple_extent_dims(scaler_space, &scalers_count[0], NULL);
                        value_count[0] = nobj;
                        value_count[1] = scalers_count[0];
                        value_count[2] = scalers_count[1];
                        // create values
                        if (false == _open_h5_dataset(STR_VALUES, scalers_type, scalers_grp_id, 3, &value_count[0], &value_count[0], values_id, value_space))
                        {
                            logE << "Error creating " << STR_VALUES << "\n";
                            return false;
                        }
                        // create names 
                        value_count[0] = 1;
                        names_cnt[0] = nobj;
                        if (false == _open_h5_dataset(STR_NAMES, name_type, scalers_grp_id, 1, &names_cnt[0], &names_cnt[0], names_id, name_space))
                        {
                            logE << "Error creating " << STR_NAMES << "\n";
                            return false;
                        }

                        // create units but will be blank since we don't know what they are
                        if (false == _open_h5_dataset(STR_UNITS, name_type, scalers_grp_id, 1, &names_cnt[0], &names_cnt[0], units_id, unit_space))
                        {
                            logE << "Error creating " << STR_NAMES << "\n";
                            return false;
                        }
                        // reset name_cnt to 1 for writing 
                        names_cnt[0] = 1;
                        buffer = new T_real[scalers_count[0] * scalers_count[1]];
                        mem_count[0] = scalers_count[0];
                        mem_count[1] = scalers_count[1];
                        _create_memory_space(2, mem_count, mem_space);
                        first_save = false;
                    }
                    value_offset[0] = i;
                    H5Sselect_hyperslab(scaler_space, H5S_SELECT_SET, scalers_offset, nullptr, scalers_count, nullptr);
                    H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
                    H5Sselect_hyperslab(name_space, H5S_SELECT_SET, names_off, nullptr, names_cnt, nullptr);
                    H5Sselect_hyperslab(unit_space, H5S_SELECT_SET, names_off, nullptr, names_cnt, nullptr);
                    status = H5Dread(dsid, scalers_type, mem_space, scaler_space, H5P_DEFAULT, &buffer[0]);
                    status = H5Dwrite(values_id, scalers_type, mem_space, value_space, H5P_DEFAULT, &buffer[0]);
                    // Replace I0 with US_IC
                    if (strncmp("I0", str_dset_name, 2047) == 0)
                    {
                        STR_US_IC.copy(str_dset_name, 2047);
                    }
                    // Replace IT with  DS_IC
                    if (strncmp("IT", str_dset_name, 2047) == 0)
                    {
                        STR_DS_IC.copy(str_dset_name, 2047);
                    }
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
            hid_t values_id, value_space, dset_unit_id, dataspace_unit_id;
            hid_t scaler_space = H5Dget_space(dset_detectors_id);
            _global_close_map.push({ scaler_space, H5O_DATASPACE });
            H5Sget_simple_extent_dims(scaler_space, &scalers_count[0], NULL);
            value_count[0] = scalers_count[2];
            value_count[1] = scalers_count[0];
            value_count[2] = scalers_count[1];
            hid_t scalers_type = H5Dget_type(dset_detectors_id);
            _global_close_map.push({ scalers_type, H5O_DATATYPE });
            if (false == _open_h5_dataset(STR_VALUES, scalers_type, scalers_grp_id, 3, &value_count[0], &value_count[0], values_id, value_space))
            {
                logE << "Error creating " << STR_NAMES << "\n";
                return false;
            }
            T_real* buffer = new T_real[scalers_count[0] * scalers_count[1]];
            hsize_t scaler_amt = scalers_count[2];
            scalers_count[2] = 1;
            value_count[0] = 1;
            mem_count[0] = scalers_count[0];
            mem_count[1] = scalers_count[1];
            hid_t mem_space = H5Screate_simple(2, mem_count, mem_count);
            _global_close_map.push({ mem_space, H5O_DATASPACE });

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
            _global_close_map.push({ dataspace_detectors_id, H5O_DATASPACE });
            attr_detector_names_id = H5Aopen(dset_detectors_id, "Detector Names", H5P_DEFAULT);
            _global_close_map.push({ dataspace_detectors_id, H5O_ATTRIBUTE });

            det_rank = H5Sget_simple_extent_ndims(dataspace_detectors_id);
            det_dims_in = new hsize_t[det_rank];
            H5Sget_simple_extent_dims(dataspace_detectors_id, &det_dims_in[0], NULL);

            hid_t ftype = H5Aget_type(attr_detector_names_id);
            hid_t type = H5Tget_native_type(ftype, H5T_DIR_ASCEND);
            error = H5Aread(attr_detector_names_id, type, detector_names);

            if (error == 0)
            {
                if (false == _open_h5_dataset(STR_VALUES, type, scalers_grp_id, 1, &det_dims_in[2], &det_dims_in[2], dataset_id, dataspace_id))
                {
                    logE << "Error creating " << STR_NAMES << "\n";
                    return false;
                }
                if (false == _open_h5_dataset(STR_VALUES, type, scalers_grp_id, 1, &det_dims_in[2], &det_dims_in[2], dset_unit_id, dataspace_unit_id))
                {
                    logE << "Error creating " << STR_NAMES << "\n";
                    return false;
                }
                status = H5Dwrite(dataset_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, detector_names);

                //for (size_t z = 0; z < det_dims_in[2]; z++)
                //{
                    //TODO: look into why this is causing exception in windows

                //}
                //free(detector_names);
            }
            delete[] det_dims_in;
        }

        _close_h5_objects(_global_close_map);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        logI << "elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
	bool save_scan_scalers_gsecars(std::string path,	
               [[maybe_unused]] size_t row_idx_start = 0,
               [[maybe_unused]] int row_idx_end = -1,
               [[maybe_unused]] size_t col_idx_start = 0,
               [[maybe_unused]] int col_idx_end = -1)
    {

        std::lock_guard<std::mutex> lock(_mutex);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        hid_t scan_grp_id, maps_grp_id, scalers_grp_id, status;
        hid_t    file_id, src_maps_grp_id;
        hid_t   xypos_dataspace_id, xypos_id;
        hid_t x_dataspace_id, y_dataspace_id, x_dataset_id, y_dataset_id;
        int det_rank;
        hid_t names_id, name_space, units_id, unit_space, values_id, value_space, mem_single_space;
        hid_t scalers_ds_id, scalers_names_id;
        hid_t upstream_ic_id, downstream_ic_id, i2_id, tscaler_id;
        hid_t buffer_space;
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
        float* x_data;
        float* y_data;

        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        logI << "Saving scalers to hdf5" << "\n";

        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_SCAN, maps_grp_id, scan_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_SCALERS, maps_grp_id, scalers_grp_id))
        {
            return false;
        }

        file_id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (file_id < 0)
            return false;
        _global_close_map.push({ file_id, H5O_FILE });

        if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, _global_close_map, "xrmmap", file_id, true, false))
        {
            if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, _global_close_map, "xrfmap", file_id))
            {
                return false;
            }
            version = GSE_CARS_SAVE_VER::XRFMAP;
        }
        else
        {
            version = GSE_CARS_SAVE_VER::XRMMAP;
        }

        if (false == _open_h5_object(xypos_id, H5O_DATASET, _global_close_map, "positions/pos", src_maps_grp_id))
            return false;
        xypos_dataspace_id = H5Dget_space(xypos_id);
        _global_close_map.push({ xypos_dataspace_id, H5O_DATASPACE });

        hid_t xy_type = H5Dget_type(xypos_id);
        _global_close_map.push({ xy_type, H5O_DATATYPE });
        det_rank = H5Sget_simple_extent_ndims(xypos_dataspace_id);
        det_dims_in = new hsize_t[det_rank];
        H5Sget_simple_extent_dims(xypos_dataspace_id, &det_dims_in[0], NULL);

        if (false == _open_h5_dataset(STR_X_AXIS, xy_type, scan_grp_id, 1, &det_dims_in[1], &det_dims_in[1], x_dataset_id, x_dataspace_id))
        {
            logE << "Error creating " << STR_X_AXIS << "\n";
            return false;
        }
        if (false == _open_h5_dataset(STR_Y_AXIS, xy_type, scan_grp_id, 1, &det_dims_in[0], &det_dims_in[0], y_dataset_id, y_dataspace_id))
        {
            logE << "Error creating " << STR_Y_AXIS << "\n";
            return false;
        }

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
        _create_memory_space(1, &single_count[0], mem_single_space);
        single_count[0] = 4;
        single_count[0] = 1;
        hid_t dtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(dtype, 255);
        if (false == _open_h5_dataset(STR_NAMES, dtype, scan_grp_id, 1, &single_count[0], &single_count[0], names_id, name_space))
        {
            logE << "Error creating " << STR_NAMES << "\n";
            return false;
        }
        if (false == _open_h5_dataset(STR_UNITS, dtype, scan_grp_id, 1, &single_count[0], &single_count[0], units_id, unit_space))
        {
            logE << "Error creating " << STR_UNITS << "\n";
            return false;
        }

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

        if (false == _open_h5_dataset(STR_VALUES, dtype, scalers_grp_id, 3, &value_count[0], &value_count[0], values_id, value_space))
        {
            logE << "Error creating " << STR_VALUES << "\n";
            return false;
        }
        T_real* buffer = new T_real[det_dims_in[0] * det_dims_in[1]];
        value_count[0] = 1;


        if (version == GSE_CARS_SAVE_VER::XRMMAP)
        {
            if (false == _open_h5_object(upstream_ic_id, H5O_DATASET, _global_close_map, "scalars/I0_raw", src_maps_grp_id))
                return false;
            if (false == _open_h5_object(downstream_ic_id, H5O_DATASET, _global_close_map, "scalars/I1_raw", src_maps_grp_id))
                return false;
            if (false == _open_h5_object(i2_id, H5O_DATASET, _global_close_map, "scalars/I2_raw", src_maps_grp_id))
                return false;
            if (false == _open_h5_object(tscaler_id, H5O_DATASET, _global_close_map, "scalars/TSCALER_raw", src_maps_grp_id))
                return false;


            hid_t scaler_space = H5Dget_space(upstream_ic_id);

            status = _read_h5d<T_real>(upstream_ic_id, scaler_space, scaler_space, H5P_DEFAULT, &buffer[0]);
            if (status > -1)
            {
                value_offset[0] = 0;
                H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
                status = _write_h5d<T_real>(values_id, scaler_space, value_space, H5P_DEFAULT, &buffer[0]);
            }

            status = _read_h5d<T_real>(downstream_ic_id, scaler_space, scaler_space, H5P_DEFAULT, &buffer[0]);
            if (status > -1)
            {
                value_offset[0] = 1;
                H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
                status = _write_h5d<T_real>(values_id, scaler_space, value_space, H5P_DEFAULT, &buffer[0]);
            }

            status = _read_h5d<T_real>(i2_id, scaler_space, scaler_space, H5P_DEFAULT, &buffer[0]);
            if (status > -1)
            {
                value_offset[0] = 2;
                H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
                status = _write_h5d<T_real>(values_id, scaler_space, value_space, H5P_DEFAULT, &buffer[0]);
            }

            status = _read_h5d<T_real>(tscaler_id, scaler_space, scaler_space, H5P_DEFAULT, &buffer[0]);
            if (status > -1)
            {
                value_offset[0] = 3;
                H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
                status = _write_h5d<T_real>(values_id, scaler_space, value_space, H5P_DEFAULT, &buffer[0]);
            }

        }
        else if (version == GSE_CARS_SAVE_VER::XRFMAP)
        {
            std::map<std::string, int> save_name_idx = { {"I0", 0}, {"I1", 1}, {"I2", 2}, {"TSCALER", 3} };
            if (false == _open_h5_object(scalers_ds_id, H5O_DATASET, _global_close_map, "roimap/det_raw", src_maps_grp_id))
                return false;
            if (false == _open_h5_object(scalers_names_id, H5O_DATASET, _global_close_map, "roimap/det_name", src_maps_grp_id))
                return false;

            hid_t scaler_name_space = H5Dget_space(scalers_names_id);
            _global_close_map.push({ buffer_space, H5O_DATASPACE });
            hid_t scaler_space = H5Dget_space(scalers_ds_id);
            _global_close_map.push({ buffer_space, H5O_DATASPACE });

            _create_memory_space(2, &mem_count[0], buffer_space);


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
                        status = _read_h5d<T_real>(scalers_ds_id, buffer_space, scaler_space, H5P_DEFAULT, &buffer[0]);
                        if (status > -1)
                        {
                            value_offset[0] = save_name_idx[sname];
                            H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, nullptr, value_count, nullptr);
                            status = _write_h5d<T_real>(values_id, buffer_space, value_space, H5P_DEFAULT, &buffer[0]);
                        }
                    }
                }
            }
        }

        hid_t ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
        _global_close_map.push({ ocpypl_id, H5O_PROPERTY });
        status = H5Ocopy(src_maps_grp_id, "config", scan_grp_id, "config", ocpypl_id, H5P_DEFAULT);

        _close_h5_objects(_global_close_map);

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

    template<typename T_real>
    bool save_scan_scalers_bnl(std::string path,
        size_t detector_num,
        [[maybe_unused]] size_t row_idx_start = 0,
        [[maybe_unused]] int row_idx_end = -1,
        [[maybe_unused]] size_t col_idx_start = 0,
        [[maybe_unused]] int col_idx_end = -1)
    {

        ////std::lock_guard<std::mutex> lock(_mutex);
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        hid_t scan_grp_id, maps_grp_id, scalers_grp_id, status;
        hid_t    file_id, src_maps_grp_id;
        hid_t   xypos_dataspace_id, xypos_id;
        hid_t x_dataspace_id, y_dataspace_id, x_dataset_id, y_dataset_id;
        int det_rank;
        hsize_t* det_dims_in = nullptr;
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
        hid_t metadata_id;
        
        double* x_data = nullptr;
        double* y_data = nullptr;
        double* buffer = nullptr;

        if (_cur_file_id < 0)
        {
            logE << "hdf5 file was never initialized. Call start_save_seq() before this function." << "\n";
            return false;
        }

        logI << "Saving scalers to hdf5" << "\n";


        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_SCAN, maps_grp_id, scan_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_SCALERS, maps_grp_id, scalers_grp_id))
        {
            return false;
        }

        file_id = H5Fopen(path.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (file_id < 0)
            return false;
        _global_close_map.push({ file_id, H5O_FILE });

        if (false == _open_h5_object(src_maps_grp_id, H5O_GROUP, _global_close_map, "xrfmap", file_id))
        {
            return false;
        }

        if (false == _open_h5_object(xypos_id, H5O_DATASET, _global_close_map, "positions/pos", src_maps_grp_id))
            return false;
        xypos_dataspace_id = H5Dget_space(xypos_id);
        _global_close_map.push({ xypos_dataspace_id, H5O_DATASPACE });

        hid_t xy_type = H5Dget_type(xypos_id);
        _global_close_map.push({ xy_type, H5O_DATATYPE });
        det_rank = H5Sget_simple_extent_ndims(xypos_dataspace_id);
        det_dims_in = new hsize_t[det_rank];
        H5Sget_simple_extent_dims(xypos_dataspace_id, &det_dims_in[0], NULL);

        if (false == _open_h5_dataset(STR_X_AXIS, xy_type, scan_grp_id, 1, &det_dims_in[1], &det_dims_in[1], x_dataset_id, x_dataspace_id))
        {
            logE << "Error creating " << STR_X_AXIS << "\n";
            return false;
        }
        if (false == _open_h5_dataset(STR_Y_AXIS, xy_type, scan_grp_id, 1, &det_dims_in[2], &det_dims_in[2], y_dataset_id, y_dataspace_id))
        {
            logE << "Error creating " << STR_Y_AXIS << "\n";
            return false;
        }


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
        if (false == _open_h5_object(scaler_grp_id, H5O_GROUP, _global_close_map, "scalers", src_maps_grp_id))
        {
            return false;
        }
        if (false == _open_h5_object(scaler_name_id, H5O_DATASET, _global_close_map, "name", scaler_grp_id))
        {
            return false;
        }
        if (false == _open_h5_object(scaler_val_id, H5O_DATASET, _global_close_map, "val", scaler_grp_id))
        {
            return false;
        }

        hid_t scaler_name_space = H5Dget_space(scaler_name_id);
        _global_close_map.push({ scaler_name_space, H5O_DATASPACE });
        hid_t scaler_type = H5Dget_type(scaler_val_id);
        _global_close_map.push({ scaler_type, H5O_DATATYPE });
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
        hid_t mem_single_space;
        _create_memory_space(1, &single_count[0], mem_single_space);
        _global_close_map.push({ mem_single_space, H5O_DATASPACE });
        hid_t names_id, name_space, units_id, units_space, values_id, value_space;
        single_count[0] = val_dims_in[2];


        hid_t dtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(dtype, 255);
        _global_close_map.push({ dtype, H5O_DATATYPE });
        hid_t dtype2 = H5Tcopy(H5T_C_S1);
        H5Tset_size(dtype2, H5T_VARIABLE);
        _global_close_map.push({ dtype2, H5O_DATATYPE });

        if (false == _open_h5_dataset(STR_NAMES, dtype, scalers_grp_id, 1, &single_count[0], &single_count[0], names_id, name_space))
        {
            logW << "Error creating " << STR_NAMES << "\n";
        }
        if (false == _open_h5_dataset(STR_UNITS, dtype, scalers_grp_id, 1, &single_count[0], &single_count[0], units_id, units_space))
        {
            logW << "Error creating " << STR_UNITS << "\n";
        }
        if (false == _open_h5_dataset(STR_VALUES, scaler_type, scalers_grp_id, 3, &value_count[0], &value_count[0], values_id, value_space))
        {
            logW << "Error creating " << STR_VALUES << "\n";
        }

        char tmp_char[255] = { 0 };
        char *tmp_char_arr[255];
        char unit_char[255] = "cts";
        buffer = new double[val_dims_in[0] * val_dims_in[1]];
        hid_t mem_space;
        _create_memory_space(2, mem_count, mem_space);
        size_t scaler_cnt = val_dims_in[2];
        val_dims_in[2] = 1;
        value_count[0] = 1;
        single_count[0] = 1;

        hid_t rstatus = H5Dread(scaler_name_id, dtype2, H5S_ALL, H5S_ALL, H5P_DEFAULT, (void*)tmp_char_arr);

        for (hsize_t i = 0; i < scaler_cnt; i++)
        {
            single_offset[0] = i;
            value_offset[0] = i;
            scaler_offset[2] = i;

            H5Sselect_hyperslab(scaler_val_space, H5S_SELECT_SET, scaler_offset, NULL, val_dims_in, NULL);
            H5Sselect_hyperslab(value_space, H5S_SELECT_SET, value_offset, NULL, value_count, NULL);
            H5Sselect_hyperslab(name_space, H5S_SELECT_SET, single_offset, nullptr, single_count, nullptr);

            if (rstatus > -1)
            {
                std::string read_name = std::string(tmp_char_arr[i]);
                read_name.erase(std::remove_if(read_name.begin(), read_name.end(), ::isspace), read_name.end());
                read_name.erase(std::find(read_name.begin(), read_name.end(), '\0'), read_name.end());
                std::string out_label = "";
                std::string beamline = "";
                bool tmpb;
                for (int j = 0; j < 255; j++)
                {
                    tmp_char[j] = '\0';
                }
                if (data_struct::Scaler_Lookup::inst()->search_pv(read_name, out_label, tmpb, beamline))
                {
                    out_label.copy(tmp_char, 255);
                }
                else
                {
                    read_name.copy(tmp_char, 255);
                }
                H5Dwrite(names_id, dtype, mem_single_space, name_space, H5P_DEFAULT, (void*)tmp_char);
                delete tmp_char_arr[i];
            }

            H5Dwrite(units_id, dtype, mem_single_space, name_space, H5P_DEFAULT, (void*)unit_char);
            status = H5Dread(scaler_val_id, scaler_type, mem_space, scaler_val_space, H5P_DEFAULT, (void*)buffer);
            if (status > -1)
            {
                H5Dwrite(values_id, scaler_type, mem_space, value_space, H5P_DEFAULT, (void*)buffer);
            }
        }

        data_struct::Scan_Info<T_real> scan_info;
        if (_open_or_create_group(STR_SCAN_METADATA, src_maps_grp_id, metadata_id))
        {
            // load attributes from this folder
            int na = H5Aget_num_attrs(metadata_id);

            double* ddata = new double[10];
            long *ldata = new long[10];

            for (int i = 0; i < na; i++)
            {
                char* adata[25];
                for (int a = 0; a < 25; a++)
                {
                    adata[a] = nullptr;
                }
                hid_t aid = H5Aopen_idx(metadata_id, (unsigned int)i);

                char buf[1000];
                ssize_t len = H5Aget_name(aid, 1000, buf);
                data_struct::Extra_PV e_pv;
                e_pv.name = std::string(buf, len);

                hid_t atype = H5Aget_type(aid);
                hid_t ntype = H5Tget_native_type(atype, H5T_DIR_ASCEND);
                //hid_t aspece = H5Aget_space(aid);
                H5T_class_t type_class =  H5Tget_class(atype);
                //if (type_class == H5T_STRING)
                if (H5Tis_variable_str(atype) > 0)
                {
                    if (H5Aread(aid, ntype, &adata) > -1)
                    {
                        e_pv.value = std::string(adata[0]);
                    }
                    for (int q = 0; q < 25; q++)
                    {
                        if (adata[q] != nullptr)
                        {
                            delete adata[q];
                            adata[q] = nullptr;
                        }
                    }
                }
                else if (type_class == H5T_ARRAY)
                {
                   
                }
                else if(type_class == H5T_FLOAT)
                {
                    if (H5Aread(aid, H5T_NATIVE_DOUBLE, (void*)ddata) > -1)
                    {
                        e_pv.value = std::to_string(ddata[0]);
                    }
                }
                else if (type_class == H5T_INTEGER)
                {
                    if (H5Aread(aid, H5T_NATIVE_LONG, (void*)ldata) > -1)
                    {
                        e_pv.value = std::to_string(ldata[0]);
                    }
                } 


                scan_info.extra_pvs.push_back(e_pv);

                H5Tclose(atype);
                H5Aclose(aid);
            }
            delete [] ddata;
            delete [] ldata;
            save_scan_scalers<T_real>(&scan_info, nullptr);

        }

        

        _close_h5_objects(_global_close_map);

        if (x_data != nullptr)
            delete[] x_data;
        if (y_data != nullptr)
            delete[] y_data;
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

	// Add links to dataset and set version to 9 so legacy software can load it
    void add_v9_layout(std::string dataset_file);

	// Add exchange layout to be loadable by external software
    void add_exchange_layout(std::string dataset_file);

	// update theta value based on new pv name
    void update_theta(std::string dataset_file, std::string theta_pv_str);

	void update_amps(std::string dataset_file, std::string us_amp_str, std::string ds_amp_str);

	void update_quant_amps(std::string dataset_file, std::string us_amp_str, std::string ds_amp_str);

    //update scalers if maps_fit_parameters_override.txt has changes pv's and you don't want to refit
    //void update_scalers(std::string dataset_file, data_struct::Params_Override<T_real>* params_override);
    
    bool end_save_seq(bool loginfo = true);

    //-----------------------------------------------------------------------------

    //export integrated spec, fitted, background into csv
    template<typename T_real>
    void export_int_fitted_to_csv(std::string dataset_file)
    {
        std::lock_guard<std::mutex> lock(_mutex);

        logI << dataset_file << "\n";

        hid_t file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

        std::string exchange_fits[4] = { STR_FIT_GAUSS_MATRIX, STR_FIT_ROI, STR_FIT_NNLS };

        if (file_id > 0)
        {
            data_struct::ArrayTr<T_real> energy_array;
            data_struct::ArrayTr<T_real> int_spectra;
            data_struct::ArrayTr<T_real> model_spectra;
            data_struct::ArrayTr<T_real> background_array;
            std::string csv_path;


            std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;
            close_map.push( {file_id, H5O_FILE} );

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
                        _read_h5d<T_real>(dset_id, dataspace_id, dataspace_id, H5P_DEFAULT, int_spectra.data());
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
                    _read_h5d<T_real>(dset_id, dataspace_id, dataspace_id, H5P_DEFAULT, energy_array.data());
                }
            }

            for (const auto& fit_itr : exchange_fits)
            {

                size_t s_idx = dataset_file.find("img.dat");
                if (s_idx != std::string::npos)
                {
                    csv_path = dataset_file.substr(0, s_idx);
                    csv_path += "output";
                    csv_path += DIR_END_CHAR;
                    csv_path += dataset_file.substr(s_idx + 8);
                }
                else
                {
                    csv_path = dataset_file;
                }
                csv_path += "_";
                csv_path += fit_itr;
                csv_path += ".csv";

                std::string h5_path = "/MAPS/XRF_Analyzed/" + fit_itr + "/" + STR_FIT_INT_BACKGROUND;
                if (_open_h5_object(dset_id, H5O_DATASET, close_map, h5_path.c_str(), file_id))
                {
                    if (dims_in[0] > 0)
                    {
                        background_array.resize(dims_in[0]);
                        dataspace_id = H5Dget_space(dset_id);
                        close_map.push({ dataspace_id, H5O_DATASPACE });
                        _read_h5d<T_real>(dset_id, dataspace_id, dataspace_id, H5P_DEFAULT, background_array.data());
                    }
                }

                h5_path = "/MAPS/XRF_Analyzed/" + fit_itr + "/" + STR_FIT_INT_SPEC;
                if (_open_h5_object(dset_id, H5O_DATASET, close_map, h5_path.c_str(), file_id))
                {
                    if (dims_in[0] > 0)
                    {
                        model_spectra.resize(dims_in[0]);
                        dataspace_id = H5Dget_space(dset_id);
                        close_map.push({ dataspace_id, H5O_DATASPACE });
                        _read_h5d<T_real>(dset_id, dataspace_id, dataspace_id, H5P_DEFAULT, model_spectra.data());
                        csv::save_fit_and_int_spectra(csv_path, &energy_array, &int_spectra, &model_spectra, &background_array);
                    }
                }
            }
            _close_h5_objects(close_map);
        }
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool add_background(std::string directory, std::string filename, data_struct::Params_Override<T_real>& params)
    {

        std::lock_guard<std::mutex> lock(_mutex);

        std::string fullname = directory + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + filename;
        logI << fullname << "\n";
        hid_t file_id = H5Fopen(fullname.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
        if (file_id < 0)
        {
            logE << "Could not open file\n";
            return false;
        }
        hid_t mca_arr_id = H5Dopen(file_id, "/MAPS/mca_arr", H5P_DEFAULT);
        if (mca_arr_id < 0)
        {
            H5Fclose(file_id);
            logW << "Could not open /MAPS/mca_arr\n";
            return false;
        }

        hid_t back_arr_id = H5Dopen(file_id, "/MAPS/mca_background", H5P_DEFAULT);
        if (back_arr_id < 0)
        {
            hid_t ocpypl_id = H5Pcreate(H5P_OBJECT_COPY);
            hid_t status = H5Ocopy(file_id, "/MAPS/mca_arr", file_id, "/MAPS/mca_background", ocpypl_id, H5P_DEFAULT);
            if (status < 0)
            {
                H5Fclose(file_id);
                logW << "Could not copy /MAPS/mca_arr to /MAPS/mca_background\n";
                return false;
            }
            back_arr_id = H5Dopen(file_id, "/MAPS/mca_background", H5P_DEFAULT);
            if (back_arr_id < 0)
            {
                H5Fclose(file_id);
                logW << "Could not open mca_background after copy /MAPS/mca_arr to /MAPS/mca_background\n";
                return false;
            }
        }

        hid_t mca_arr_space = H5Dget_space(mca_arr_id);

        hsize_t dims_in[3];
        hsize_t offset[3] = { 0,0,0 };
        hsize_t count[3] = { 0,1,1 };
        int status_n = H5Sget_simple_extent_dims(mca_arr_space, &dims_in[0], NULL);
        if (status_n < 0)
        {
            H5Fclose(file_id);
            logW << "Can't get dims\n";
            return false;
        }
        count[0] = dims_in[0];
        hid_t memoryspace_id = H5Screate_simple(1, dims_in, nullptr);

        data_struct::ArrayTr<T_real>   buffer(count[0]);
        fitting::models::Range energy_range = data_struct::get_energy_range(dims_in[0], &(params.fit_params));

        logI << params.fit_params.value(STR_ENERGY_OFFSET) << " " << params.fit_params.value(STR_ENERGY_SLOPE) << " " << params.fit_params.value(STR_ENERGY_QUADRATIC) << " " << 0.0f << " " << params.fit_params.value(STR_SNIP_WIDTH) << " " << energy_range.min << " " << energy_range.max << "\n ";

        for (hsize_t x = 0; x < dims_in[1]; x++)
        {
            logI << fullname << " " << x << " " << dims_in[1] << "\n";
            for (hsize_t y = 0; y < dims_in[2]; y++)
            {
                offset[1] = x;
                offset[2] = y;
                H5Sselect_hyperslab(mca_arr_space, H5S_SELECT_SET, offset, nullptr, count, nullptr);
                hid_t error = _read_h5d<T_real>(mca_arr_id, memoryspace_id, mca_arr_space, H5P_DEFAULT, buffer.data());
                if (error > -1)
                {
                    data_struct::ArrayTr<T_real> background = data_struct::snip_background<T_real>((data_struct::Spectra<T_real>*) & buffer, params.fit_params.value(STR_ENERGY_OFFSET), params.fit_params.value(STR_ENERGY_SLOPE), params.fit_params.value(STR_ENERGY_QUADRATIC), params.fit_params.value(STR_SNIP_WIDTH), energy_range.min, energy_range.max);
                    error = _write_h5d<T_real>(back_arr_id, memoryspace_id, mca_arr_space, H5P_DEFAULT, background.data());
                    if (error < 0)
                    {
                        logE << x << " " << y << " bad write\n";
                    }
                }
            }
        }
        H5Dclose(mca_arr_id);
        H5Dclose(back_arr_id);
        H5Fclose(file_id);

        return true;
    }

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool save_params_override(data_struct::Params_Override<T_real>* params_override)
    {

        if (params_override == nullptr)
        {
            logE << " Params_override is null\n";
            return false;
        }
        hid_t   memtype_label, filetype_label, memoryspace_id, maps_grp_id, po_grp_id;
        hid_t br_grp_id, k_grp_id, l_grp_id, m_grp_id, l_fam_grp_id;
        hid_t fitp_names_id, fitp_values_id;
        hid_t fitp_names_space, fitp_values_space;

        filetype_label = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype_label, 50);
        memtype_label = H5Tcopy(H5T_C_S1);
        hid_t status = H5Tset_size(memtype_label, 50);

        hsize_t f_offset[1] = { 0 };
        hsize_t f_count[1] = { 0 };
        hsize_t fp_offset[2] = { 0 };
        hsize_t fp_count[2] = { 0 };
        f_count[0] = params_override->fit_params.size();
        fp_count[0] = f_count[0];
        fp_count[1] = 4; // value, min, max, step

        if (false == _open_or_create_group(STR_MAPS, _cur_file_id, maps_grp_id))
        {
            return false;
        }
        if (false == _open_or_create_group(STR_FIT_PARAMETERS_OVERRIDE, maps_grp_id, po_grp_id))
        {
            return false;
        }

        std::string str_name = STR_FIT_PARAMETERS + "_Names";
        std::string str_values = STR_FIT_PARAMETERS + "_Values";
        bool name_ok = _open_h5_dataset(str_name, filetype_label, po_grp_id, 1, f_count, f_count, fitp_names_id, fitp_names_space);
        bool prop_ok = _open_h5_dataset<T_real>(str_values, po_grp_id, 2, fp_count, fp_count, fitp_values_id, fitp_values_space);

        if (false == name_ok || false == prop_ok)
        {
            logE << " Failed to open or create datasetsl\n";
            return false;
        }

        f_count[0] = 1;
        _create_memory_space(1, f_count, memoryspace_id);
        fp_count[0] = 1;
        fp_count[1] = 1;

        int i = 0;
        for (typename std::unordered_map<std::string, Fit_Param<T_real>>::const_iterator itr = params_override->fit_params.begin(); itr != params_override->fit_params.end(); itr++)
        {
            f_offset[0] = i;
            fp_offset[0] = i;

            char label[50] = { 0 };
            itr->first.copy(&label[0], 49);

            status = H5Sselect_hyperslab(fitp_names_space, H5S_SELECT_SET, f_offset, nullptr, f_count, nullptr);

            status = H5Dwrite(fitp_names_id, memtype_label, memoryspace_id, fitp_names_space, H5P_DEFAULT, (void*)&label);
            if (status < 0)
            {
                logE << "failed to write name\n";
            }

            fp_offset[1] = 0;
            status = H5Sselect_hyperslab(fitp_values_space, H5S_SELECT_SET, fp_offset, nullptr, fp_count, nullptr);
            status = _write_h5d<T_real>(fitp_values_id, memoryspace_id, fitp_values_space, H5P_DEFAULT, (void*)&(itr->second.value));
            if (status < 0)
            {
                logE << "failed to write value\n";
            }

            fp_offset[1] = 1;
            status = H5Sselect_hyperslab(fitp_values_space, H5S_SELECT_SET, fp_offset, nullptr, fp_count, nullptr);
            status = _write_h5d<T_real>(fitp_values_id, memoryspace_id, fitp_values_space, H5P_DEFAULT, (void*)&(itr->second.min_val));
            if (status < 0)
            {
                logE << "failed to write min value\n";
            }

            fp_offset[1] = 2;
            status = H5Sselect_hyperslab(fitp_values_space, H5S_SELECT_SET, fp_offset, nullptr, fp_count, nullptr);
            status = _write_h5d<T_real>(fitp_values_id, memoryspace_id, fitp_values_space, H5P_DEFAULT, (void*)&(itr->second.max_val));
            if (status < 0)
            {
                logE << "failed to write max value\n";
            }

            fp_offset[1] = 3;
            status = H5Sselect_hyperslab(fitp_values_space, H5S_SELECT_SET, fp_offset, nullptr, fp_count, nullptr);
            status = _write_h5d<T_real>(fitp_values_id, memoryspace_id, fitp_values_space, H5P_DEFAULT, (void*)&(itr->second.step_size));
            if (status < 0)
            {
                logE << "failed to write step size\n";
            }

            i++;
        }
        
        if(params_override->branching_ratios.size() > 0
        || params_override->branching_ratio_K.size() > 0
        || params_override->branching_ratio_L.size() > 0
        || params_override->branching_ratio_M.size() > 0
        || params_override->branching_family_L.size() > 0)
        {
            // save branching ratio's 
            if (false == _open_or_create_group(STR_BRANCHING_RATIOS, _cur_file_id, br_grp_id))
            {
                return false;
            }
            if(params_override->branching_ratio_K.size() > 0)
            {
                // save K ratios 
                if (_open_or_create_group(STR_K_SHELL, br_grp_id, k_grp_id))
                {
                    
                }
            }
            if(params_override->branching_ratio_L.size() > 0)
            {
                // save l ratios 
                if (_open_or_create_group(STR_L_SHELL, br_grp_id, l_grp_id))
                {
                    
                }
            }
            if(params_override->branching_ratio_M.size() > 0)
            {
                // save K Shells 
                if (_open_or_create_group(STR_M_SHELL, br_grp_id, m_grp_id))
                {
                    
                }
            }
            if(params_override->branching_family_L.size() > 0)
            {
                // save K Shells 
                if (_open_or_create_group(STR_L_FAMILY, br_grp_id, l_fam_grp_id))
                {
                    
                }
            }
        }

        return true;
    }


    //-----------------------------------------------------------------------------

private:

    HDF5_IO();

    static HDF5_IO *_this_inst;

    //static std::mutex _mutex;
    std::mutex _mutex;

    //-----------------------------------------------------------------------------

    template<typename T_real>
	bool _load_integrated_spectra_analyzed_h5(hid_t file_id, data_struct::Spectra<T_real>* spectra)
    {

        hid_t    dset_id, dataspace_id, spec_grp_id, memoryspace_id, memoryspace_meta_id, dset_incnt_id, dset_outcnt_id, dset_rt_id, dset_lt_id;
        hid_t    dataspace_lt_id, dataspace_rt_id, dataspace_inct_id, dataspace_outct_id;
        herr_t   error;
        hsize_t dims_in[1] = { 0 };
        hsize_t offset[1] = { 0 };
        hsize_t count[1] = { 1 };
        hsize_t offset_time[1] = { 0 };
        hsize_t count_time[1] = { 1 };

        std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

        close_map.push({ file_id, H5O_FILE });

        if (false == _open_h5_object(spec_grp_id, H5O_GROUP, close_map, "/MAPS/Spectra/Integrated_Spectra", file_id))
            return false;

        if (false == _open_h5_object(dset_id, H5O_DATASET, close_map, "Spectra", spec_grp_id))
            return false;
        dataspace_id = H5Dget_space(dset_id);
        close_map.push({ dataspace_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_lt_id, H5O_DATASET, close_map, "Elapsed_Livetime", spec_grp_id))
            return false;
        dataspace_lt_id = H5Dget_space(dset_lt_id);
        close_map.push({ dataspace_lt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_rt_id, H5O_DATASET, close_map, "Elapsed_Realtime", spec_grp_id))
            return false;
        dataspace_rt_id = H5Dget_space(dset_rt_id);
        close_map.push({ dataspace_rt_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_incnt_id, H5O_DATASET, close_map, "Input_Counts", spec_grp_id))
            return false;
        dataspace_inct_id = H5Dget_space(dset_incnt_id);
        close_map.push({ dataspace_inct_id, H5O_DATASPACE });

        if (false == _open_h5_object(dset_outcnt_id, H5O_DATASET, close_map, "Output_Counts", spec_grp_id))
            return false;
        dataspace_outct_id = H5Dget_space(dset_outcnt_id);
        close_map.push({ dataspace_outct_id, H5O_DATASPACE });


        int rank = H5Sget_simple_extent_ndims(dataspace_id);
        if (rank != 1)
        {
            _close_h5_objects(close_map);
            logE << "Dataset /MAPS/Spectra/mca_arr  rank != 3. rank = " << rank << ". Can't load dataset. returning" << "\n";
            return false;
            //throw exception ("Dataset is not a volume");
        }

        int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], nullptr);
        if (status_n < 0)
        {
            _close_h5_objects(close_map);
            logE << "getting dataset dims for /MAPS/Spectra/mca_arr" << "\n";
            return false;
        }

        for (int i = 0; i < rank; i++)
        {

            offset[i] = 0;
            count[i] = dims_in[i];
        }

        spectra->resize(dims_in[0]);
        spectra->setZero(dims_in[0]);

        count[0] = dims_in[0];

        memoryspace_id = H5Screate_simple(1, count, nullptr);
        memoryspace_meta_id = H5Screate_simple(1, count_time, nullptr);
        H5Sselect_hyperslab(memoryspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);
        H5Sselect_hyperslab(memoryspace_meta_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

        T_real live_time = 1.0;
        T_real real_time = 1.0;
        T_real in_cnt = 1.0;
        T_real out_cnt = 1.0;

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

        error = _read_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(*spectra)[0]);

        H5Sselect_hyperslab(dataspace_lt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
        H5Sselect_hyperslab(dataspace_rt_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
        H5Sselect_hyperslab(dataspace_inct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);
        H5Sselect_hyperslab(dataspace_outct_id, H5S_SELECT_SET, offset_time, nullptr, count_time, nullptr);

        error = _read_h5d<T_real>(dset_rt_id, memoryspace_meta_id, dataspace_rt_id, H5P_DEFAULT, (void*)&real_time);
        error = _read_h5d<T_real>(dset_lt_id, memoryspace_meta_id, dataspace_lt_id, H5P_DEFAULT, (void*)&live_time);
        error = _read_h5d<T_real>(dset_incnt_id, memoryspace_meta_id, dataspace_inct_id, H5P_DEFAULT, (void*)&in_cnt);
        error = _read_h5d<T_real>(dset_outcnt_id, memoryspace_meta_id, dataspace_outct_id, H5P_DEFAULT, (void*)&out_cnt);

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

    template<typename T_real>
    bool _save_scan_meta_data(hid_t scan_grp_id, data_struct::Scan_Meta_Info<T_real>* meta_info)
    {

        hid_t dataspace_id = -1, memoryspace_id = -1;
        herr_t status = -1;
        hid_t filetype, memtype = -1;
        hid_t dset_id = -1;
        hsize_t count[1] = { 1 };

        try
        {
            filetype = H5Tcopy(H5T_FORTRAN_S1);
            H5Tset_size(filetype, 256);
            memtype = H5Tcopy(H5T_C_S1);
            status = H5Tset_size(memtype, 255);

            //save y axis
            count[0] = meta_info->y_axis.size();
            _create_memory_space(1, count, memoryspace_id);
            if (_open_h5_dataset<T_real>(STR_Y_AXIS, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)meta_info->y_axis.data());
                if (status < 0)
                {
                    logE << "failed to write " << STR_Y_AXIS << "\n";
                }
            }

            count[0] = meta_info->x_axis.size();
            _create_memory_space(1, count, memoryspace_id);
            if (_open_h5_dataset<T_real>(STR_X_AXIS, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)meta_info->x_axis.data());
                if (status < 0)
                {
                    logE << "failed to write " << STR_X_AXIS << "\n";
                }
            }

            //save requested rows
            count[0] = 1;
            _create_memory_space(1, count, memoryspace_id);
            if (_open_h5_dataset(STR_REQUESTED_ROWS, H5T_INTEL_I32, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {       
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(meta_info->requested_rows));
                if (status < 0)
                {
                    logE << "failed to write " << STR_REQUESTED_ROWS << "\n";
                }
            }

            //save requested cols
            if (_open_h5_dataset(STR_REQUESTED_COLS, H5T_INTEL_I32, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {
                status = H5Dwrite(dset_id, H5T_NATIVE_INT, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&(meta_info->requested_cols));
                if (status < 0)
                {
                    logE << "failed to write " << STR_REQUESTED_COLS << "\n";
                }
            }

            //Save theta
            if (_open_h5_dataset<T_real>(STR_THETA, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {   
                status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&meta_info->theta);
                if (status < 0)
                {
                    logE << "failed to write " << STR_THETA << "\n";
                }
            }

            //Save scan type  
            if ( _open_h5_dataset(STR_XRF_SCAN_TYPE, filetype, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {   
                char tmp_char[255] = { 0 };
                meta_info->scan_type.copy(tmp_char, 254);
                status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
                if (status < 0)
                {
                    logE << "failed to write " << STR_XRF_SCAN_TYPE << "\n";
                }
            }

            //Save polarity pattern
            if (_open_h5_dataset(STR_POLARITY_PATTERN, filetype, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {
                char tmp_char[255] = { 0 };
                meta_info->polarity_pattern.copy(tmp_char, 254);
                status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
                if (status < 0)
                {
                    logE << "failed to write " << STR_POLARITY_PATTERN << "\n";
                }
            }

            if (_open_h5_dataset(STR_SCAN_TIME_STAMP, filetype, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {
                char tmp_char[255] = { 0 };
                meta_info->scan_time_stamp.copy(tmp_char, 254);
                status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
                if (status < 0)
                {
                    logE << "failed to write " << STR_SCAN_TIME_STAMP << "\n";
                }
            }
            if (_open_h5_dataset(STR_NAME, filetype, scan_grp_id, 1, count, count, dset_id, dataspace_id))
            {   
                char tmp_char2[255] = { 0 };
                meta_info->name.copy(tmp_char2, 254);
                status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char2);
                if (status < 0)
                {
                    logE << "failed to write " << STR_NAME << "\n";
                }
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

	bool _save_extras(hid_t scan_grp_id, std::vector<data_struct::Extra_PV>* extra_pvs);

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool _save_scalers(hid_t maps_grp_id, std::vector<data_struct::Scaler_Map<T_real>>*scalers_map, T_real us_amps_val, std::string us_amps_unit, T_real ds_amps_val, std::string ds_amps_unit)
    {

        hid_t dataspace_values_id = -1, memoryspace_id = -1, dataspace_names_id = -1, memoryspace_str_id = -1;
        hid_t dataspace_units_id = -1;
        hid_t filetype, memtype = -1;
        hid_t dset_names_id = -1;
        hid_t dset_units_id = -1;
        hid_t dset_values_id = -1;
        hid_t scalers_grp_id = -1;
        hid_t dcpl_id = -1;
        herr_t status;

        hsize_t offset[1] = { 0 };
        hsize_t count[1] = { 1 };

        hsize_t count_2d[2] = { 1, 1 };

        hsize_t offset_3d[3] = { 0, 0, 0 };
        hsize_t count_3d[3] = { 1, 1, 1 };
        hsize_t chunk_3d[3] = { 1, 1, 1 };

        _create_memory_space(1, count, memoryspace_str_id);

        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
        memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);

        count[0] = 1;
        std::string units;

        if (false == _open_or_create_group(STR_SCALERS, maps_grp_id, scalers_grp_id))
        {
            return false;
        }

        _save_amps(scalers_grp_id, us_amps_val, us_amps_unit, ds_amps_val, ds_amps_unit);

        if (scalers_map != nullptr && scalers_map->size() > 0)
        {

            int cols = 0;
            int rows = 0;

            data_struct::Scaler_Map<T_real> map = scalers_map->front();
            rows = map.values.rows();
            cols = map.values.cols();
            if (rows > 0 && cols > 0)
            {
                // create calculated scalers
                data_struct::Scaler_Map<T_real> abs_ic_map, abs_cfg_map, H_dpc_cfg_map, V_dpc_cfg_map, dia1_dpc_cfg_map, dia2_dpc_cfg_map;
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
                data_struct::ArrayXXr<T_real>* us_ic_map = nullptr;
                data_struct::ArrayXXr<T_real>* ds_ic_map = nullptr;
                data_struct::ArrayXXr<T_real>* cfg_2_map = nullptr;
                data_struct::ArrayXXr<T_real>* cfg_3_map = nullptr;
                data_struct::ArrayXXr<T_real>* cfg_4_map = nullptr;
                data_struct::ArrayXXr<T_real>* cfg_5_map = nullptr;

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
                    if (upper_scaler_name == STR_CFG_2)
                    {
                        cfg_2_map = &(scaler.values);
                    }
                    if (upper_scaler_name == STR_CFG_3)
                    {
                        cfg_3_map = &(scaler.values);
                    }
                    if (upper_scaler_name == STR_CFG_4)
                    {
                        cfg_4_map = &(scaler.values);
                    }
                    if (upper_scaler_name == STR_CFG_5)
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
                    data_struct::ArrayXXr<T_real> t_abs_map;
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
                for (const auto& itr : *scalers_map)
                {
                    count_3d[1] = itr.values.rows();
                    count_2d[0] = count_3d[1];
                    count_3d[2] = itr.values.cols();
                    count_2d[1] = count_3d[2];
                    break;
                }

                chunk_3d[0] = 1;
                chunk_3d[1] = count_3d[1];
                chunk_3d[2] = count_3d[2];

                if (false == _open_h5_dataset<T_real>(STR_VALUES, scalers_grp_id, 3, count_3d, chunk_3d, dset_values_id, dataspace_values_id))
                {
                    return false;
                }

                count_3d[0] = 1;

                dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(dcpl_id, 3, count_3d);
                H5Pset_deflate(dcpl_id, 7);

                count_3d[0] = scalers_map->size();
                count[0] = count_3d[0];

                if (false == _open_h5_dataset(STR_NAMES, filetype, scalers_grp_id, 1, count, count, dset_names_id, dataspace_names_id))
                {
                    return false;
                }

                if (false == _open_h5_dataset(STR_UNITS, filetype, scalers_grp_id, 1, count, count, dset_units_id, dataspace_units_id))
                {
                    return false;
                }

                count_3d[0] = 1;
                count[0] = 1;

                _create_memory_space(2, count_2d, memoryspace_id);
                int idx = 0;
                for (auto& itr : *scalers_map)
                {
                    offset[0] = idx;
                    offset_3d[0] = idx;
                    idx++;
                    char tmp_char[255] = { 0 };
                    char tmp_char_units[255] = { 0 };
                    std::string out_label, out_beamline;
                    if (data_struct::Scaler_Lookup::inst()->search_pv(itr.name, out_label, itr.time_normalized, out_beamline))
                    {
                        out_label.copy(tmp_char, 254);
                    }
                    else
                    {
                        itr.name.copy(tmp_char, 254);
                    }
                    itr.unit.copy(tmp_char_units, 254);
                    H5Sselect_hyperslab(dataspace_names_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    H5Sselect_hyperslab(dataspace_units_id, H5S_SELECT_SET, offset, NULL, count, NULL);
                    status = H5Dwrite(dset_names_id, memtype, memoryspace_str_id, dataspace_names_id, H5P_DEFAULT, (void*)tmp_char);
                    if (status < 0)
                    {
                        logE << "failed to write " << STR_NAMES << "\n";
                    }
                    status = H5Dwrite(dset_units_id, memtype, memoryspace_str_id, dataspace_units_id, H5P_DEFAULT, (void*)tmp_char_units);
                    if (status < 0)
                    {
                        logE << "failed to write " << STR_UNITS << "\n";
                    }
                    H5Sselect_hyperslab(dataspace_values_id, H5S_SELECT_SET, offset_3d, NULL, count_3d, NULL);
                    itr.values = itr.values.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });
                    status = _write_h5d<T_real>(dset_values_id, memoryspace_id, dataspace_values_id, H5P_DEFAULT, (void*)itr.values.data());
                    if (status < 0)
                    {
                        logE << "failed to write " << STR_VALUES << "\n";
                    }
                }
            }
        }

        H5Tclose(filetype);
        H5Tclose(memtype);

        logI << "Done" << "\n";

        return true;
    }
    
    //-----------------------------------------------------------------------------

    template<typename T_real>
    void _save_amps(hid_t scalers_grp_id, T_real us_amp_sens_num_val, std::string us_amp_sens_unit_val, T_real ds_amp_sens_num_val, std::string ds_amp_sens_unit_val)
    {

        hid_t dataspace_us_id, dataspace_ds_id, memoryspace_id;
        hid_t dset_us_id, dset_ds_id;
        herr_t status;
        hid_t dset_id, dataspace_id;
        hid_t filetype, memtype;

        char tmp_char[255] = { 0 };
        hsize_t offset[1] = { 0 };
        hsize_t count[1] = { 3 };

        filetype = H5Tcopy(H5T_FORTRAN_S1);
        H5Tset_size(filetype, 256);
        memtype = H5Tcopy(H5T_C_S1);
        status = H5Tset_size(memtype, 255);
        std::string units;

        T_real trans_us_amp_sens_num_val = translate_back_sens_num<T_real>((int)us_amp_sens_num_val);
        T_real trans_us_amp_sens_unit_val = translate_back_sens_unit<T_real>(us_amp_sens_unit_val);;
        T_real trans_ds_amp_sens_num_val = translate_back_sens_num<T_real>((int)ds_amp_sens_num_val);;
        T_real trans_ds_amp_sens_unit_val = translate_back_sens_unit<T_real>(ds_amp_sens_unit_val);;

        if (false == _open_h5_dataset<T_real>(STR_US_AMP, scalers_grp_id, 1, count, count, dset_us_id, dataspace_us_id))
        {
            return;
        }
        if (false == _open_h5_dataset<T_real>(STR_DS_AMP, scalers_grp_id, 1, count, count, dset_ds_id, dataspace_ds_id))
        {
            return;
        }

        count[0] = 1;
        _create_memory_space(1, count, memoryspace_id);
        H5Sselect_hyperslab(dataspace_us_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = _write_h5d<T_real>(dset_us_id, memoryspace_id, dataspace_us_id, H5P_DEFAULT, (void*)&trans_us_amp_sens_num_val);

        offset[0] = 1;
        H5Sselect_hyperslab(dataspace_us_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = _write_h5d<T_real>(dset_us_id, memoryspace_id, dataspace_us_id, H5P_DEFAULT, (void*)&trans_us_amp_sens_unit_val);

        offset[0] = 2;
        H5Sselect_hyperslab(dataspace_us_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = _write_h5d<T_real>(dset_us_id, memoryspace_id, dataspace_us_id, H5P_DEFAULT, (void*)&us_amp_sens_num_val);


        offset[0] = 0;
        H5Sselect_hyperslab(dataspace_ds_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = _write_h5d<T_real>(dset_ds_id, memoryspace_id, dataspace_ds_id, H5P_DEFAULT, (void*)&trans_ds_amp_sens_num_val);

        offset[0] = 1;
        H5Sselect_hyperslab(dataspace_ds_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = _write_h5d<T_real>(dset_ds_id, memoryspace_id, dataspace_ds_id, H5P_DEFAULT, (void*)&trans_ds_amp_sens_unit_val);

        offset[0] = 2;
        H5Sselect_hyperslab(dataspace_ds_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = _write_h5d<T_real>(dset_ds_id, memoryspace_id, dataspace_ds_id, H5P_DEFAULT, (void*)&ds_amp_sens_num_val);

        offset[0] = 0;
        count[0] = 1;

        if (false == _open_h5_dataset<T_real>(STR_US_AMP_NUM, scalers_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return;
        }
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&us_amp_sens_num_val);

        if (false == _open_h5_dataset<T_real>(STR_DS_AMP_NUM, scalers_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return;
        }
        status = _write_h5d<T_real>(dset_id, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)&ds_amp_sens_num_val);


        us_amp_sens_unit_val.copy(tmp_char, 254);
        if (false == _open_h5_dataset(STR_US_AMP_UNIT, filetype, scalers_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return;
        }
        status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);

        ds_amp_sens_unit_val.copy(tmp_char, 254);
        if (false == _open_h5_dataset(STR_DS_AMP_UNIT, filetype, scalers_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return;
        }
        status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);


        H5Tclose(filetype);
        H5Tclose(memtype);

    }
    
    //-----------------------------------------------------------------------------

    void _gen_average(std::string full_hdf5_path, std::string dataset_name, hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids, bool avg=true);
    
    void _generate_avg_analysis(hid_t src_maps_grp_id, hid_t dst_maps_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids);
    
    void _generate_avg_integrated_spectra(hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids);

    void _add_v9_quant(hid_t file_id, hid_t chan_names, hid_t chan_space, int chan_amt, std::string quant_str, std::string new_loc);
    void _add_v9_scalers(hid_t file_id);
    void _add_extra_pvs(hid_t file_id, std::string group_name);

    bool _add_exchange_meta(hid_t file_id, std::string exchange_idx, std::string fits_link, std::string normalize_scaler);
	
    bool _open_h5_object(hid_t &id, H5_OBJECTS obj, std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map, std::string s1, hid_t id2, bool log_error=true, bool close_on_fail=true);
    bool _open_or_create_group(const std::string name, hid_t parent_id, hid_t& out_id, bool log_error = true, bool close_on_fail = true);
    bool _create_memory_space(int rank, const hsize_t* count, hid_t& out_id);
    bool _open_h5_dataset(const std::string& name, hid_t data_type, hid_t parent_id, int dims_size, const hsize_t* dims, const hsize_t* chunk_dims, hid_t& out_id, hid_t& out_dataspece);

    //-----------------------------------------------------------------------------

    template<typename T_real>
    bool _open_h5_dataset(const std::string& name, hid_t parent_id, int dims_size, const hsize_t* dims, const hsize_t* chunk_dims, hid_t& out_id, hid_t& out_dataspece)
    {
        if (std::is_same<T_real, float>::value)
        {
            return _open_h5_dataset(name, H5T_INTEL_F32, parent_id, dims_size, dims, chunk_dims, out_id, out_dataspece);
        }
        else if (std::is_same<T_real, double>::value)
        {
            return _open_h5_dataset(name, H5T_INTEL_F64, parent_id, dims_size, dims, chunk_dims, out_id, out_dataspece);
        }
        return false;
    }

    //-----------------------------------------------------------------------------

    void _close_h5_objects(std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map);

    //-----------------------------------------------------------------------------

    template<typename T_real>
    herr_t _read_h5d(hid_t dset_id, hid_t mem_space_id, hid_t file_space_id, hid_t plist_id, void* buf)
    {
        if (std::is_same<T_real, float>::value)
        {
            return H5Dread(dset_id, H5T_NATIVE_FLOAT, mem_space_id, file_space_id, plist_id, buf);
        }
        else if (std::is_same<T_real, double>::value)
        {
            return H5Dread(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, plist_id, buf);
        }
        return -1;
    }
    
    //-----------------------------------------------------------------------------

    template<typename T_real>
    herr_t _write_h5d(hid_t dset_id, hid_t mem_space_id, hid_t file_space_id, hid_t plist_id, const void* buf)
    {
        if (std::is_same<T_real, float>::value)
        {
            return H5Dwrite(dset_id, H5T_NATIVE_FLOAT, mem_space_id, file_space_id, plist_id, buf);
        }
        else if (std::is_same<T_real, double>::value)
        {
            return H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, file_space_id, plist_id, buf);
        }
        return -1;
    }

    //-----------------------------------------------------------------------------

    hid_t _cur_file_id;
    std::string _cur_filename;
    std::stack<std::pair<hid_t, H5_OBJECTS> > _global_close_map;

};


}// end namespace file
}// end namespace io

#endif // HDF5_IO_H
