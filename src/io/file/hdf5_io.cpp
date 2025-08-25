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
#include <filesystem>


const std::vector<std::string> xrf_analysis_save_names = { STR_FIT_ROI,
                                                          STR_FIT_GAUSS_MATRIX,
                                                          STR_FIT_SVD,
                                                          STR_FIT_NNLS
                                                         };





namespace io
{
namespace file
{

hsize_t max_dims_1d[1] = { H5S_UNLIMITED };
hsize_t max_dims_2d[2] = { H5S_UNLIMITED, H5S_UNLIMITED };
hsize_t max_dims_3d[3] = { H5S_UNLIMITED, H5S_UNLIMITED, H5S_UNLIMITED };

//std::mutex HDF5_IO::_mutex;

//-----------------------------------------------------------------------------

HDF5_IO* HDF5_IO::_this_inst(nullptr);


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
    //std::lock_guard<std::mutex> lock(_mutex);

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

//-----------------------------------------------------------------------------

bool HDF5_IO::_open_or_create_group(const std::string name, hid_t parent_id, hid_t& out_id, bool log_error, bool close_on_fail)
{
    out_id = H5Gopen(parent_id, name.c_str(), H5P_DEFAULT);
    if (out_id < 0)
        out_id = H5Gcreate(parent_id, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (out_id < 0)
    {
        if (log_error)
        {
            logE << "creating group " << name << "\n";
        }
        if (close_on_fail)
        {
            _close_h5_objects(_global_close_map);
        }
        return false;
    }
    _global_close_map.push({ out_id, H5O_GROUP });
    return true;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_create_memory_space(int rank, const hsize_t* count, hid_t& out_id)
{
    out_id = H5Screate_simple(rank, count, nullptr);
    if (out_id > -1)
    {
        _global_close_map.push({ out_id , H5O_DATASPACE });
        return true;
    }
    return false;
}

//-----------------------------------------------------------------------------

bool HDF5_IO::_open_h5_dataset(const std::string& name, hid_t data_type, hid_t parent_id, int dims_size, const hsize_t* dims, const hsize_t* chunk_dims, hid_t& out_id, hid_t& out_dataspece)
{
    out_id = H5Dopen(parent_id, name.c_str(), H5P_DEFAULT);
    if (out_id < 0)
    {
        // if doesn't exist, create new one
        hsize_t* max_dims;
        switch (dims_size)
        {
        case 1:
            max_dims = &max_dims_1d[0];
            break;
        case 2:
            max_dims = &max_dims_2d[0];
            break;
        case 3:
            max_dims = &max_dims_3d[0];
            break;
        default:
            return false;
        }
        out_dataspece = H5Screate_simple(dims_size, dims, max_dims);
        _global_close_map.push({ out_dataspece, H5O_DATASPACE });

        hid_t dcpl_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(dcpl_id, dims_size, chunk_dims);
        H5Pset_deflate(dcpl_id, 7);
        _global_close_map.push({ dcpl_id, H5O_PROPERTY });

        out_id = H5Dcreate(parent_id, name.c_str(), data_type, out_dataspece, H5P_DEFAULT, dcpl_id, H5P_DEFAULT);
        if (out_id > -1)
        {
            _global_close_map.push({ out_id, H5O_DATASET });
            return true;
        }
    }
    else
    {
        //if we are opening, check the dim size to see if we have to expand to fit new counts size
        _global_close_map.push({ out_id, H5O_DATASET });
        out_dataspece = H5Dget_space(out_id);
        _global_close_map.push({ out_dataspece, H5O_DATASPACE });
        hsize_t* tmp_dims = new hsize_t[dims_size];
        int status_n = H5Sget_simple_extent_dims(out_dataspece, tmp_dims, NULL);
        if (status_n > -1)
        {
            bool expand = false;
            for (size_t i = 0; i < dims_size; i++)
            {
                if (tmp_dims[i] != dims[i])
                {
                    expand = true;
                    break;
                }
            }
            if (expand)
            {
                herr_t err = H5Dset_extent(out_id, dims);
                if (err < 0)
                {
                    logE << "Failed to extend dim [" << name << "] from " << tmp_dims << " to " << dims << " \n";
                    delete[]tmp_dims;
                    return false;
                }
                out_dataspece = H5Dget_space(out_id);
                _global_close_map.push({ out_dataspece, H5O_DATASPACE });
            }
        }
        delete[]tmp_dims;
        return true;
    }

    _close_h5_objects(_global_close_map);
    return false;
}

//-----------------------------------------------------------------------------

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
                logW<<"Could not close h5 attribute id "<<obj.first<<"\n";
            }
        }
        else if (obj.second == H5O_DATATYPE)
        {
            err = H5Tclose(obj.first);
            if (err > 0)
            {
                logW << "Could not close h5 type id " << obj.first << "\n";
            }
        }
        else if (obj.second == H5O_PROPERTY)
        {
            err = H5Pclose(obj.first);
            if (err > 0)
            {
                logW << "Could not close h5 property id " << obj.first << "\n";
            }
        }
    }
}

//-----------------------------------------------------------------------------

bool HDF5_IO::start_save_seq(const std::string filename, bool force_new_file, bool open_file_only)
{

    if (_cur_file_id > -1)
    {
        logI<<" file already open, calling close() before opening new file. "<<"\n";
        end_save_seq();
    }

    if(false == force_new_file)
    {
        _cur_file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else
    {
        try 
        {
            if (std::filesystem::remove(filename))
            {
               logI << "file " << filename << " deleted.\n";
            }
        }
        catch(const std::filesystem::filesystem_error& err) 
        {
            logE << "filesystem error: " << err.what() << '\n';
        }
    }
    if (_cur_file_id < 1)
    {
        if (open_file_only == false)
        {
            logI << "Creating file " << filename << "\n";
            _cur_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        }
        else
        {
            return false;
        }
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
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_DATASET, (size_t)-1, objs ) > -1)
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
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_GROUP, (size_t)-1, objs ) > -1)
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
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_DATATYPE, (size_t)-1, objs ) > -1)
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
            if( H5Fget_obj_ids( _cur_file_id, H5F_OBJ_ATTR, (size_t)-1, objs ) > -1)
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

bool HDF5_IO::_save_extras(hid_t scan_grp_id, std::vector<data_struct::Extra_PV>* extra_pvs)
{
    hid_t memoryspace_id;
    herr_t status;
    hid_t dataspace_desc_id, dataspace_unit_id, dataspace_val_id, dataspace_id;
	hid_t filetype, memtype;
    hid_t dset_desc_id = -1, dset_unit_id = -1, dset_id = -1, dset_val_id = -1;
    hid_t extra_grp_id = -1;
	
	hsize_t offset[1] = { 0 };
	hsize_t count[1] = { 1 };

	if (extra_pvs == nullptr)
	{
		return false;
	}

    filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);
	memtype = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(memtype, 255);

    if (false == _open_or_create_group(STR_EXTRA_PVS, scan_grp_id, extra_grp_id))
    {
        return false;
    }
		
    if (extra_pvs->size() > 0)
    {
        //save extra pv's
        count[0] = extra_pvs->size();
        if (false == _open_h5_dataset(STR_NAMES, filetype, extra_grp_id, 1, count, count, dset_id, dataspace_id))
        {
            return false;
        }

        if (false == _open_h5_dataset(STR_VALUES, filetype, extra_grp_id, 1, count, count, dset_val_id, dataspace_val_id))
        {
            return false;
        }

        if (false == _open_h5_dataset(STR_DESCRIPTION, filetype, extra_grp_id, 1, count, count, dset_desc_id, dataspace_desc_id))
        {
            return false;
        }

        if (false == _open_h5_dataset(STR_UNIT, filetype, extra_grp_id, 1, count, count, dset_unit_id, dataspace_unit_id))
        {
            return false;
        }

        count[0] = 1;

        std::string str_val;

        count[0] = 1;
        _create_memory_space(1, count, memoryspace_id);

        for (int16_t i = 0; i < extra_pvs->size(); i++)
        {
            offset[0] = i;
            data_struct::Extra_PV pv = extra_pvs->at(i);

            H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(dataspace_val_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(dataspace_desc_id, H5S_SELECT_SET, offset, NULL, count, NULL);
            H5Sselect_hyperslab(dataspace_unit_id, H5S_SELECT_SET, offset, NULL, count, NULL);

            char tmp_char[255] = { 0 };
            pv.name.copy(tmp_char, 254);
            status = H5Dwrite(dset_id, memtype, memoryspace_id, dataspace_id, H5P_DEFAULT, (void*)tmp_char);
            if (status < 0)
            {
                logE << "failed to write " << STR_NAMES << "\n";
            }

            char tmp_char2[255] = { 0 };
            pv.value.copy(tmp_char2, 254);
            status = H5Dwrite(dset_val_id, memtype, memoryspace_id, dataspace_val_id, H5P_DEFAULT, (void*)tmp_char2);
            if (status < 0)
            {
                logE << "failed to write " << STR_VALUES << "\n";
            }

            char tmp_char3[255] = { 0 };
            pv.description.copy(tmp_char, 254);
            status = H5Dwrite(dset_desc_id, memtype, memoryspace_id, dataspace_desc_id, H5P_DEFAULT, (void*)tmp_char3);
            if (status < 0)
            {
                logE << "failed to write " << STR_DESCRIPTION << "\n";
            }

            char tmp_char4[255] = { 0 };
            pv.unit.copy(tmp_char, 254);
            status = H5Dwrite(dset_unit_id, memtype, memoryspace_id, dataspace_unit_id, H5P_DEFAULT, (void*)tmp_char4);
            if (status < 0)
            {
                logE << "failed to write " << STR_UNIT << "\n";
            }
        }
    }
    H5Tclose(filetype);
    H5Tclose(memtype);
    
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
                        dst_fit_grp_id = H5Gcreate(dst_calib_fit_grp_id, analysis_grp_name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                        std::string chan_name_loc = analysis_grp_name + "/" + STR_CALIB_LABELS;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        // copy element info index and names 
                        chan_name_loc = analysis_grp_name + "/" + STR_SR_CURRENT + STR_ELEMENT_INFO_INDEX;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        chan_name_loc = analysis_grp_name + "/" + STR_SR_CURRENT + STR_ELEMENT_INFO_NAMES;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        chan_name_loc = analysis_grp_name + "/" + STR_DS_IC + STR_ELEMENT_INFO_INDEX;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        chan_name_loc = analysis_grp_name + "/" + STR_DS_IC + STR_ELEMENT_INFO_NAMES;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        chan_name_loc = analysis_grp_name + "/" + STR_US_IC + STR_ELEMENT_INFO_INDEX;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        chan_name_loc = analysis_grp_name + "/" + STR_US_IC + STR_ELEMENT_INFO_NAMES;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        chan_name_loc = analysis_grp_name + "/" + STR_US_FM + STR_ELEMENT_INFO_INDEX;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        chan_name_loc = analysis_grp_name + "/" + STR_US_FM + STR_ELEMENT_INFO_NAMES;
                        status = H5Ocopy(cablib_grp_id, chan_name_loc.c_str(), dst_calib_fit_grp_id, chan_name_loc.c_str(), ocpypl_id, H5P_DEFAULT);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/" + STR_CALIB_CURVE + STR_SR_CURRENT, STR_CALIB_CURVE + STR_SR_CURRENT, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/" + STR_SR_CURRENT + STR_ELEMENT_INFO_VALUES, STR_SR_CURRENT + STR_ELEMENT_INFO_VALUES, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/" + STR_CALIB_CURVE + STR_DS_IC, STR_CALIB_CURVE + STR_DS_IC, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/" + STR_DS_IC + STR_ELEMENT_INFO_VALUES, STR_DS_IC + STR_ELEMENT_INFO_VALUES, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/" + STR_CALIB_CURVE + STR_US_IC, STR_CALIB_CURVE + STR_US_IC, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/" + STR_US_IC + STR_ELEMENT_INFO_VALUES, STR_US_IC + STR_ELEMENT_INFO_VALUES, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/" + STR_CALIB_CURVE + STR_US_FM, STR_CALIB_CURVE + STR_US_FM, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
                        _gen_average("MAPS/Quantification/Calibration/" + analysis_grp_name + "/" + STR_US_FM + STR_ELEMENT_INFO_VALUES, STR_US_FM + STR_ELEMENT_INFO_VALUES, src_fit_grp_id, dst_fit_grp_id, ocpypl_id, hdf5_file_ids);
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
                        std::string standard_group_name = "Standard" + std::to_string(i);
                        std::string whole_standard_name = "MAPS/Quantification/" + standard_group_name;
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

        H5Pclose(ocpypl_id);
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

void HDF5_IO::_gen_average(std::string full_hdf5_path, std::string dataset_name, hid_t src_fit_grp_id, hid_t dst_fit_grp_id, [[maybe_unused]] hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids, bool avg)
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
        if (dst_dset_id < 1)
        {
            //logI<<""<<full_hdf5_path<< " " <<dataset_name<<"\n";
            return;
        }
        int rank = H5Sget_simple_extent_ndims(dataspace_id);

        hsize_t* tmp_dims = new hsize_t[rank];
        hsize_t* dims_in = new hsize_t[rank];
        int status_n = H5Sget_simple_extent_dims(dataspace_id, &dims_in[0], NULL);
        if (status_n < 0)
        {
            logW << "could not get dataset dimensions for " << full_hdf5_path << " " << dataset_name << "\n";
            return;
        }

        //get all the other files dataset ids
        for (long unsigned int k = 1; k < hdf5_file_ids.size(); k++)
        {
            hid_t det_analysis_dset_id = H5Dopen2(hdf5_file_ids[k], full_hdf5_path.c_str(), H5P_DEFAULT);
            if (det_analysis_dset_id > -1)
            {
                // ran into a bug where detectors had different dims for counts per sec. need to check the min and use that to generate avg
                hid_t tdataspace_id = H5Dget_space(det_analysis_dset_id);
                status_n = H5Sget_simple_extent_dims(tdataspace_id, &tmp_dims[0], NULL);
                if (status_n > -1)
                {
                    for (int i = 0; i < rank; i++)
                    {
                        if (tmp_dims[i] < dims_in[i])
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
        for (int i = 0; i < rank; i++)
        {
            offset_rank[i] = 0;
            total *= dims_in[i];
        }

        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset_rank, nullptr, dims_in, nullptr);

        delete[]offset_rank;

        //long long avail_mem = get_total_mem() * .95;
        long long avail_mem = get_available_mem();
        //logI << "Available memory: " << avail_mem << "\n";
        
        if (H5Tequal(file_type, H5T_NATIVE_DOUBLE) || H5Tequal(file_type, H5T_INTEL_F64))
        {
            if ( ( (total * (long long)sizeof(double) * 2) > (avail_mem) ) && (rank == 3) ) //mca_arr
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
                data_struct::ArrayTr<double> buffer1(chunk_total); //don't need to zero because we are reading in full buffer
                data_struct::ArrayTr<double> buffer2(chunk_total); //don't need to zero because we are reading in full buffer

                hid_t memoryspace_id = H5Screate_simple(rank, chunk_dims, nullptr);
                for (int w = 0; w < dims_in[1]; w++)
                {
                    for (int h = 0; h < dims_in[2]; h++)
                    {
                        offset[1] = w;
                        offset[2] = h;
                        float divisor = 1.0;
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
                        error = H5Dread(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
                        if (error > -1)
                        {
                            buffer1 = buffer1.unaryExpr([](double v) { return std::isfinite(v) ? v : (double)0.0; });
                        }
                        for (long unsigned int k = 0; k < analysis_ids.size(); k++)
                        {
                            //hid_t dspace_id = H5Dget_space(analysis_ids[k]);
                            //H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
                            //error = H5Dread(analysis_ids[k], H5T_NATIVE_DOUBLE, memoryspace_id, dspace_id, H5P_DEFAULT, buffer2.data());
                            error = H5Dread(analysis_ids[k], H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer2.data());
                            if (error > -1)
                            {
                                buffer2 = buffer2.unaryExpr([](double v) { return std::isfinite(v) ? v : (double)0.0; });
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

                        error = H5Dwrite(dst_dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
                    }
                }
                delete[] offset;
                delete[] chunk_dims;
                H5Sclose(memoryspace_id);
            }
            else
            {
                //read in the whole dataset
                data_struct::ArrayTr<double> buffer1(total);
                data_struct::ArrayTr<double> buffer2(total);
                float divisor = 1.0;
                error = H5Dread(dset_id, H5T_NATIVE_DOUBLE, dataspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
                if (error > -1)
                {
                    buffer1 = buffer1.unaryExpr([](double v) { return std::isfinite(v) ? v : (double)0.0; });
                }
                for (long unsigned int k = 0; k < analysis_ids.size(); k++)
                {
                    error = H5Dread(analysis_ids[k], H5T_NATIVE_DOUBLE, dataspace_id, dataspace_id, H5P_DEFAULT, buffer2.data());
                    if (error > -1)
                    {
                        buffer2 = buffer2.unaryExpr([](double v) { return std::isfinite(v) ? v : (double)0.0; });
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
                error = H5Dwrite(dst_dset_id, H5T_NATIVE_DOUBLE, dataspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
            }
        }
        else  //else float
        {
            if ((total * (long long)sizeof(float) * 2) > (avail_mem) && rank == 3) //mca_arr
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
                data_struct::ArrayTr<float> buffer1(chunk_total); //don't need to zero because we are reading in full buffer
                data_struct::ArrayTr<float> buffer2(chunk_total); //don't need to zero because we are reading in full buffer

                hid_t memoryspace_id = H5Screate_simple(rank, chunk_dims, nullptr);
                for (int w = 0; w < dims_in[1]; w++)
                {
                    for (int h = 0; h < dims_in[2]; h++)
                    {
                        offset[1] = w;
                        offset[2] = h;
                        float divisor = 1.0;
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
                        error = H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
                        if (error > -1)
                        {
                            buffer1 = buffer1.unaryExpr([](float v) { return std::isfinite(v) ? v : (float)0.0; });
                        }
                        for (long unsigned int k = 0; k < analysis_ids.size(); k++)
                        {
                            //hid_t dspace_id = H5Dget_space(analysis_ids[k]);
                            //H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, nullptr, chunk_dims, nullptr);
                            //error = H5Dread(analysis_ids[k], H5T_NATIVE_REAL, memoryspace_id, dspace_id, H5P_DEFAULT, buffer2.data());
                            error = H5Dread(analysis_ids[k], H5T_NATIVE_FLOAT, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer2.data());
                            if (error > -1)
                            {
                                buffer2 = buffer2.unaryExpr([](float v) { return std::isfinite(v) ? v : (float)0.0; });
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

                        error = H5Dwrite(dst_dset_id, H5T_NATIVE_FLOAT, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
                    }
                }
                delete[] offset;
                delete[] chunk_dims;
                H5Sclose(memoryspace_id);
            }
            else
            {
                //read in the whole dataset
                data_struct::ArrayTr<float> buffer1(total);
                data_struct::ArrayTr<float> buffer2(total);
                float divisor = 1.0;
                error = H5Dread(dset_id, H5T_NATIVE_FLOAT, dataspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
                if (error > -1)
                {
                    buffer1 = buffer1.unaryExpr([](float v) { return std::isfinite(v) ? v : (float)0.0; });
                }
                for (long unsigned int k = 0; k < analysis_ids.size(); k++)
                {
                    error = H5Dread(analysis_ids[k], H5T_NATIVE_FLOAT, dataspace_id, dataspace_id, H5P_DEFAULT, buffer2.data());
                    if (error > -1)
                    {
                        buffer2 = buffer2.unaryExpr([](float v) { return std::isfinite(v) ? v : (float)0.0; });
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
                error = H5Dwrite(dst_dset_id, H5T_NATIVE_FLOAT, dataspace_id, dataspace_id, H5P_DEFAULT, buffer1.data());
            }
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
        H5Gclose(src_analyzed_grp_id);
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
		_gen_average(group_name + "/Integrated_Spectra/"+ STR_MAX_CHANNELS_INT_SPEC, STR_MAX_CHANNELS_INT_SPEC, src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids, false);
		_gen_average(group_name + "/Integrated_Spectra/"+ STR_MAX10_INT_SPEC, STR_MAX10_INT_SPEC, src_inner_grp_id, dst_inner_grp_id, ocpypl_id, hdf5_file_ids, false);

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
                     [[maybe_unused]] size_t height,
                     [[maybe_unused]] size_t width)
{

    std::string str_detector_num = std::to_string(detector_num);
    std::string full_save_path = dataset_directory+ DIR_END_CHAR+"img.dat"+ DIR_END_CHAR +dataset_name+".h5"+str_detector_num;
    //io::file::HDF5_IO::inst()->set_filename(full_save_path);


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
	float theta_value = 0;

	file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
		return;
	close_map.push({ file_id, H5O_FILE });
	
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
		value.erase(std::remove_if(value.begin(), value.end(), ::isspace), value.end());
		if (theta_pv_str == value)
		{
			for (int z = 0; z < 256; z++)
				tmp_char[z] = 0;
			rerror = H5Dread(extra_values, name_type, memoryspace_id, name_space, H5P_DEFAULT, (void*)tmp_char);
			theta_value = atof(tmp_char);
			rerror = H5Dwrite(theta_id, H5T_NATIVE_FLOAT, memoryspace_id, theta_space, H5P_DEFAULT, (void*)&theta_value);
			break;
		}

	}

	_close_h5_objects(close_map);

}

//-----------------------------------------------------------------------------

void HDF5_IO::update_amps(std::string dataset_file, std::string us_amp_str, std::string ds_amp_str)
{
	std::lock_guard<std::mutex> lock(_mutex);
	hid_t file_id, us_amp_id, us_amp_num_id, ds_amp_id, ds_amp_num_id;

	hsize_t offset_1d[1] = { 2 };
	hsize_t count_1d[1] = { 1 };
	hid_t rerror = 0;
	float us_amp_value = parse_input_real<float>(us_amp_str);
	float ds_amp_value = parse_input_real<float>(ds_amp_str);

	file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
		logW << "Could not open file " << dataset_file << "\n";
		return;
	}
    _global_close_map.push({ file_id, H5O_FILE });
    hid_t memoryspace_id, amp_space_id;
    _create_memory_space(1, count_1d, memoryspace_id);

    count_1d[0] = 3;
    if (false == _open_h5_dataset("/MAPS/Scalers/us_amp", H5T_NATIVE_FLOAT, file_id, 1, count_1d, count_1d, us_amp_id, amp_space_id))
    {
        logE << "Error creating " << "/MAPS/Scalers/us_amp" << "\n";
    }

	// try v9 layout 
	if (us_amp_id < 0)
	{
		_open_h5_object(us_amp_id, H5O_DATASET, _global_close_map, "/MAPS/us_amp", file_id, false, false);
	}

	if (us_amp_id > -1)
	{
		hid_t amp_space = H5Dget_space(us_amp_id);
        _global_close_map.push({ amp_space, H5O_DATASPACE });
		offset_1d[0] = 2;
		count_1d[0] = 1;
		H5Sselect_hyperslab(amp_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
		rerror = H5Dwrite(us_amp_id, H5T_NATIVE_FLOAT, memoryspace_id, amp_space, H5P_DEFAULT, (void*)&us_amp_value);
	}

    if (false == _open_h5_dataset("/MAPS/Scalers/ds_amp", H5T_NATIVE_FLOAT, file_id, 1, count_1d, count_1d, ds_amp_id, amp_space_id))
    {
        logE << "Error creating " << "/MAPS/Scalers/ds_amp" << "\n";
    }

	// try v9 layout 
	if (ds_amp_id < 0)
	{
		_open_h5_object(ds_amp_id, H5O_DATASET, _global_close_map, "/MAPS/ds_amp", file_id, false, false);
	}

	if (ds_amp_id > -1)
	{
		hid_t amp_space = H5Dget_space(ds_amp_id);
        _global_close_map.push({ amp_space, H5O_DATASPACE });
		offset_1d[0] = 2;
		count_1d[0] = 1;
		H5Sselect_hyperslab(amp_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
		rerror = H5Dwrite(ds_amp_id, H5T_NATIVE_FLOAT, memoryspace_id, amp_space, H5P_DEFAULT, (void*)&ds_amp_value);
	}


    if (false == _open_h5_dataset("/MAPS/Scalers/us_amp_num", H5T_NATIVE_FLOAT, file_id, 1, count_1d, count_1d, us_amp_num_id, amp_space_id))
    {
        logE << "Error creating " << "/MAPS/Scalers/ds_amp" << "\n";
    }

	if (us_amp_num_id > -1)
	{
		hid_t amp_space = H5Dget_space(us_amp_num_id);
        _global_close_map.push({ amp_space, H5O_DATASPACE });
		offset_1d[0] = 0;
		count_1d[0] = 1;
		H5Sselect_hyperslab(amp_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
		rerror = H5Dwrite(us_amp_num_id, H5T_NATIVE_FLOAT, memoryspace_id, amp_space, H5P_DEFAULT, (void*)&us_amp_value);
	}

    if (false == _open_h5_dataset("/MAPS/Scalers/ds_amp_num", H5T_NATIVE_FLOAT, file_id, 1, count_1d, count_1d, ds_amp_num_id, amp_space_id))
    {
        logE << "Error creating " << "/MAPS/Scalers/ds_amp_num" << "\n";
    }

	if (ds_amp_num_id > -1)
	{
		hid_t amp_space = H5Dget_space(ds_amp_num_id);
        _global_close_map.push({ amp_space, H5O_DATASPACE });
		offset_1d[0] = 0;
		count_1d[0] = 1;
		H5Sselect_hyperslab(amp_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
		rerror = H5Dwrite(ds_amp_num_id, H5T_NATIVE_FLOAT, memoryspace_id, amp_space, H5P_DEFAULT, (void*)&ds_amp_value);
	}

	_close_h5_objects(_global_close_map);

}

//-----------------------------------------------------------------------------

void HDF5_IO::update_quant_amps(std::string dataset_file, std::string us_amp_str, std::string ds_amp_str)
{
	std::lock_guard<std::mutex> lock(_mutex);
	hid_t file_id, us_amp_id, ds_amp_id, num_stand_id;
	std::stack<std::pair<hid_t, H5_OBJECTS> > close_map;

	hsize_t offset_1d[1] = { 2 };
	hsize_t count_1d[1] = { 1 };
	hid_t rerror = 0;
	float us_amp_value = parse_input_real<float>(us_amp_str);
	float ds_amp_value = parse_input_real<float>(ds_amp_str);
	int num_stands;
    std::string q_loc_pre_str = "/MAPS/Quantification/Standard";
    std::string q_loc_post_str = "/Scalers/";
    std::string q_us_str = "us_amp";
    std::string q_ds_str = "ds_amp";


	file_id = H5Fopen(dataset_file.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0)
	{
		logW << "Could not open file " << dataset_file << "\n";
		return;
	}
	close_map.push({ file_id, H5O_FILE });
	hid_t memoryspace_id = H5Screate_simple(1, count_1d, nullptr);
	close_map.push({ memoryspace_id, H5O_DATASPACE });

	bool update_v9 = false;

	if (_open_h5_object(num_stand_id, H5O_DATASET, close_map, "/MAPS/Quantification/Number_Of_Standards", file_id, false, false))
	{
		hid_t stand_space = H5Dget_space(num_stand_id);
		close_map.push({ stand_space, H5O_DATASPACE });
		rerror = H5Dread(num_stand_id, H5T_INTEL_I32, memoryspace_id, stand_space, H5P_DEFAULT, (void*)&num_stands);
		if (rerror > -1)
		{
			if (num_stands == 0)
			{
				update_v9 = true;
			}
			for (int i = 0; i < num_stands; ++i)
			{
                std::string q_loc = q_loc_pre_str + std::to_string(i) + q_loc_post_str + q_us_str;
				if (false == _open_h5_object(us_amp_id, H5O_DATASET, close_map, q_loc.c_str(), file_id, false, false))
				{
					count_1d[0] = 3;
                    hid_t amp_space_id;
                    if (false == _open_h5_dataset(q_loc, H5T_NATIVE_FLOAT, file_id, 1, count_1d, count_1d, us_amp_id, amp_space_id))
                    {
                        logE << "Error creating " << q_loc << "\n";
                    }
				}

				if (us_amp_id > -1)
				{
					hid_t amp_space = H5Dget_space(us_amp_id);
					close_map.push({ amp_space, H5O_DATASPACE });
					offset_1d[0] = 2;
					count_1d[0] = 1;
					H5Sselect_hyperslab(amp_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
					rerror = H5Dwrite(us_amp_id, H5T_NATIVE_FLOAT, memoryspace_id, amp_space, H5P_DEFAULT, (void*)&us_amp_value);
				}

				q_loc = q_loc_pre_str + std::to_string(i) + q_loc_post_str + q_ds_str;
				if (false == _open_h5_object(ds_amp_id, H5O_DATASET, close_map, q_loc.c_str(), file_id, false, false))
				{
					count_1d[0] = 3;
                    hid_t amp_space_id;
                    if (false == _open_h5_dataset(q_loc, H5T_NATIVE_FLOAT, file_id, 1, count_1d, count_1d, ds_amp_id, amp_space_id))
                    {
                        logE << "Error creating " << q_loc << "\n";
                    }
				}

				if (ds_amp_id > -1)
				{
					hid_t amp_space = H5Dget_space(ds_amp_id);
					close_map.push({ amp_space, H5O_DATASPACE });
					offset_1d[0] = 2;
					count_1d[0] = 1;
					H5Sselect_hyperslab(amp_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
					rerror = H5Dwrite(ds_amp_id, H5T_NATIVE_FLOAT, memoryspace_id, amp_space, H5P_DEFAULT, (void*)&ds_amp_value);
				}
			}
		}
	}
	else
	{
		update_v9 = true;
	}

	if(update_v9)
	{
		// try v9 layout 
		if (_open_h5_object(us_amp_id, H5O_DATASET, close_map, "/MAPS/make_maps_conf/nbs1832/us_amp", file_id, false, false))
		{
			hid_t amp_space = H5Dget_space(us_amp_id);
			close_map.push({ amp_space, H5O_DATASPACE });
			offset_1d[0] = 2;
			count_1d[0] = 1;
			H5Sselect_hyperslab(amp_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
			rerror = H5Dwrite(us_amp_id, H5T_NATIVE_FLOAT, memoryspace_id, amp_space, H5P_DEFAULT, (void*)&us_amp_value);
		}
		if (_open_h5_object(ds_amp_id, H5O_DATASET, close_map, "/MAPS/make_maps_conf/nbs1832/ds_amp", file_id, false, false))
		{
			hid_t amp_space = H5Dget_space(ds_amp_id);
			close_map.push({ amp_space, H5O_DATASPACE });
			offset_1d[0] = 2;
			count_1d[0] = 1;
			H5Sselect_hyperslab(amp_space, H5S_SELECT_SET, offset_1d, nullptr, count_1d, nullptr);
			rerror = H5Dwrite(ds_amp_id, H5T_NATIVE_FLOAT, memoryspace_id, amp_space, H5P_DEFAULT, (void*)&ds_amp_value);
		}
	}

	_close_h5_objects(close_map);

}

//-----------------------------------------------------------------------------
/*
void HDF5_IO::update_scalers(std::string dataset_file, data_struct::Params_Override<T_real>* params_override)
{
    std::lock_guard<std::mutex> lock(_mutex);
    if (params_override == nullptr)
    {
        return;
    }
    
    hid_t maps_grp_id;
    T_real time_scaler_clock = 1;
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
*/
//-----------------------------------------------------------------------------

void HDF5_IO::_add_v9_quant(hid_t file_id, 
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
    std::string currnt_quant_str = "/MAPS/Quantification/Calibration/" + quant_str + "/" + STR_CALIB_CURVE + STR_SR_CURRENT;
    std::string us_quant_str = "/MAPS/Quantification/Calibration/" + quant_str + "/" + STR_CALIB_CURVE + STR_US_IC;
    std::string ds_quant_str = "/MAPS/Quantification/Calibration/" + quant_str + "/" + STR_CALIB_CURVE + STR_DS_IC;
	hid_t cc_current = H5Dopen(file_id, currnt_quant_str.c_str(), H5P_DEFAULT);
	hid_t cc_us_ic = H5Dopen(file_id, us_quant_str.c_str(), H5P_DEFAULT);
	hid_t cc_ds_ic = H5Dopen(file_id, ds_quant_str.c_str(), H5P_DEFAULT);
    hid_t quant_dset, quant_space;
    hsize_t quant_dims[3] = { 3,1,1 };
    quant_dims[2] = chan_amt;
    if (false == _open_h5_dataset(new_loc, H5T_NATIVE_FLOAT, file_id, 3, quant_dims, quant_dims, quant_dset, quant_space))
    {
        logE << "Error creating " << new_loc << "\n";
    }
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
        float real_val = 0.0;
        hid_t err;

        STR_US_IC.copy(us_ic_carr, 254);
        STR_DS_IC.copy(ds_ic_carr, 254);

        count_1d[0] = 3;
        hid_t quant_names_dset;
        hid_t quant_name_space;
        count_1d[0] = 1;
        //save quant_names to know what each index is
        new_loc += "_names";

        if (false == _open_h5_dataset(new_loc, filetype, file_id, 1, count_1d, count_1d, quant_names_dset, quant_name_space))
        {
            logE << "Error creating " << new_loc << "\n";
        }
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
                size_t underscore_idx = el_name_str.find("_");
                //can check if > 0 instead of -1 since it shouldn't start with an '_'
                if (underscore_idx != std::string::npos && underscore_idx < el_name_str.length())
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
                auto element = data_struct::Element_Info_Map<float>::inst()->get_element(el_name_str);
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
                    if (H5Dread(cc_current, H5T_NATIVE_FLOAT, memoryspace_id, cc_space, H5P_DEFAULT, (void*)&real_val) > -1)
                    {
                        err = H5Dwrite(quant_dset, H5T_NATIVE_FLOAT, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);
                    }
                    offset_3d[0] = 1;
                    real_val = 0.0;
                    H5Sselect_hyperslab(quant_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                    if (H5Dread(cc_us_ic, H5T_NATIVE_FLOAT, memoryspace_id, us_space, H5P_DEFAULT, (void*)&real_val) > -1)
                    {
                        err = H5Dwrite(quant_dset, H5T_NATIVE_FLOAT, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);
                    }
                    offset_3d[0] = 2;
                    real_val = 0.0;
                    H5Sselect_hyperslab(quant_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
                    if (H5Dread(cc_ds_ic, H5T_NATIVE_FLOAT, memoryspace_id, ds_space, H5P_DEFAULT, (void*)&real_val) > -1)
                    {
                        err = H5Dwrite(quant_dset, H5T_NATIVE_FLOAT, memoryspace_id, quant_space, H5P_DEFAULT, (void*)&real_val);
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

    if (extra_names > -1)
        _global_close_map.push({extra_names, H5O_DATASET});
    if (extra_units > -1)
        _global_close_map.push({extra_units, H5O_DATASET});
    if (extra_values > -1)
        _global_close_map.push({extra_values, H5O_DATASET});
    if (extra_desc > -1)
        _global_close_map.push({extra_desc, H5O_DATASET});

    if(extra_names > -1 && extra_units > -1 && extra_values > -1 && extra_desc > -1)
    {
        hid_t name_space = H5Dget_space(extra_names);
        _global_close_map.push({name_space, H5O_DATASPACE});
        int rank = H5Sget_simple_extent_ndims(name_space);
        hsize_t* dims_in = new hsize_t[rank];
        H5Sget_simple_extent_dims(name_space, &dims_in[0], nullptr);
        hsize_t extra_pv_dims[2];
        extra_pv_dims[0] = 4;
        extra_pv_dims[1] = dims_in[0];

        hid_t file_space = H5Screate_simple(2, &extra_pv_dims[0], &extra_pv_dims[0]);
        _global_close_map.push({file_space, H5O_DATASPACE});

        hid_t name_type = H5Dget_type(extra_names);
        _global_close_map.push({name_type, H5O_DATATYPE});
        std::string extra_pvs_str = group_name + "/extra_pvs";
        std::string extra_pvs_as_csv_str = group_name + "/extra_pvs_as_csv";
        hid_t extra_pvs = H5Dcreate1(file_id, extra_pvs_str.c_str(), name_type, file_space, H5P_DEFAULT);
        if(extra_pvs > -1)
        {
            _global_close_map.push({extra_pvs, H5O_DATASET});
        }
        hid_t extra_pvs_as_csv = H5Dcreate1(file_id, extra_pvs_as_csv_str.c_str(), name_type, name_space, H5P_DEFAULT);
        if(extra_pvs_as_csv > -1)
        {
            _global_close_map.push({extra_pvs_as_csv, H5O_DATASET});
        }
        hsize_t offset_1d[1] = {0};
        hsize_t count_1d[1] = {1};
        hsize_t offset_2d[2] = {0,0};
        hsize_t count_2d[2] = {1,1};

        hid_t memoryspace_id = H5Screate_simple(1, count_1d, nullptr);
        _global_close_map.push({memoryspace_id, H5O_DATATYPE});

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

void HDF5_IO::add_v9_layout(std::string dataset_file)
{
    std::lock_guard<std::mutex> lock(_mutex);
    double* dbuf = nullptr;
    float* fbuf = nullptr;
    logI  << dataset_file << "\n";

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
        logW<<"Failed to open file "<<dataset_file<<"\n";
		return;
	}
    _global_close_map.push({file_id, H5O_FILE });
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
    std::string fit_int_back_name = "/MAPS/XRF_Analyzed/Fitted/" + STR_FIT_INT_BACKGROUND;
    std::string nnls_int_name = "/MAPS/XRF_Analyzed/NNLS/" + STR_FIT_INT_SPEC;
    std::string nnls_int_back_name = "/MAPS/XRF_Analyzed/NNLS/" + STR_FIT_INT_BACKGROUND;
	std::string max_name = "/MAPS/Spectra/Integrated_Spectra/"+ STR_MAX_CHANNELS_INT_SPEC;
	std::string max10_name = "/MAPS/Spectra/Integrated_Spectra/"+ STR_MAX10_INT_SPEC;
	std::string v9_max_name = "/MAPS/max_chan_spec";
	hid_t fit_int_id = -1, max_id = -1, max_10_id = -1, nnls_id = -1, back_id = -1, max_space = -1, max_type, v9_max_id = -1, v9_space = -1;
	
	max_id = H5Dopen(file_id, max_name.c_str(), H5P_DEFAULT);
    if(max_id > -1)
    {
        _global_close_map.push({max_id, H5O_DATASET });
    }
    fit_int_id = H5Dopen(file_id, fit_int_name.c_str(), H5P_DEFAULT);
    if(fit_int_id > -1)
    {
        _global_close_map.push({fit_int_id, H5O_DATASET });
    }
    max_10_id = H5Dopen(file_id, max10_name.c_str(), H5P_DEFAULT);
    if(max_10_id > -1)
    {
        _global_close_map.push({max_10_id, H5O_DATASET });
    }
    nnls_id = H5Dopen(file_id, nnls_int_name.c_str(), H5P_DEFAULT);
    if(nnls_id > -1)
    {
        _global_close_map.push({nnls_id, H5O_DATASET });
    }

    if (fit_int_id > -1)
    {
        back_id = H5Dopen(file_id, fit_int_back_name.c_str(), H5P_DEFAULT);
    }
    else if (nnls_id > -1)
    {
        back_id = H5Dopen(file_id, nnls_int_back_name.c_str(), H5P_DEFAULT);
    }

    if(back_id > -1)
    {
        _global_close_map.push({back_id, H5O_DATASET });
    }

    if (max_id > -1)
    {
        max_space = H5Dget_space(max_id);
        _global_close_map.push({max_space, H5O_DATASPACE });
        max_type = H5Dget_type(max_id);
        _global_close_map.push({max_type, H5O_DATATYPE });
        H5Sget_simple_extent_dims(max_space, &count2d[0], nullptr);
    }
    else if (nnls_id > -1)
    {
        max_space = H5Dget_space(nnls_id);
        _global_close_map.push({max_space, H5O_DATASPACE });
        max_type = H5Dget_type(nnls_id);
        _global_close_map.push({max_type, H5O_DATATYPE });
        H5Sget_simple_extent_dims(max_space, &count2d[0], nullptr);
    }

    if(max_space > -1)
    {
		//make 5 x spectra_size matrix
		count2d[1] = count2d[0];
		count2d[0] = 5;
		v9_space = H5Screate_simple(2, &count2d[0], &count2d[0]);

		v9_max_id = H5Dopen(file_id, v9_max_name.c_str(), H5P_DEFAULT);
		if (v9_max_id < 0)
		{
			v9_max_id = H5Dcreate(file_id, v9_max_name.c_str(), max_type, v9_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            if(v9_max_id > -1)
            {
                _global_close_map.push({v9_max_id, H5O_DATASET });
            }
		}
		if (v9_max_id > -1)
		{
            void* buf = nullptr;
            if (H5Tequal(max_type, H5T_NATIVE_DOUBLE) || H5Tequal(max_type, H5T_INTEL_F64))
            {
                dbuf = new double[count2d[1]];
                buf = dbuf;
            }
            else
            {
                fbuf = new float[count2d[1]];
                buf = fbuf;
            }
			count2d[0] = 1;
			
            if (max_id > -1)
            {
                if (H5Dread(max_id, max_type, max_space, max_space, H5P_DEFAULT, (void*)buf) > -1)
                {
                    offset2d[0] = 0;
                    H5Sselect_hyperslab(v9_space, H5S_SELECT_SET, offset2d, nullptr, count2d, nullptr);
                    H5Dwrite(v9_max_id, max_type, max_space, v9_space, H5P_DEFAULT, (void*)buf);
                }
            }
			if (max_10_id > -1)
			{
				if (H5Dread(max_10_id, max_type, max_space, max_space, H5P_DEFAULT, (void*)buf) > -1)
				{
					offset2d[0] = 1;
					H5Sselect_hyperslab(v9_space, H5S_SELECT_SET, offset2d, nullptr, count2d, nullptr);
					H5Dwrite(v9_max_id, max_type, max_space, v9_space, H5P_DEFAULT, (void*)buf);
				}
			}
			if (fit_int_id > -1)
			{
				if (H5Dread(fit_int_id, max_type, max_space, max_space, H5P_DEFAULT, (void*)buf) > -1)
				{
					offset2d[0] = 2;
					H5Sselect_hyperslab(v9_space, H5S_SELECT_SET, offset2d, nullptr, count2d, nullptr);
					H5Dwrite(v9_max_id, max_type, max_space, v9_space, H5P_DEFAULT, (void*)buf);
				}
			}
            if (nnls_id > -1)
            {
                if (H5Dread(nnls_id, max_type, max_space, max_space, H5P_DEFAULT, (void*)buf) > -1)
                {
                    offset2d[0] = 3;
                    H5Sselect_hyperslab(v9_space, H5S_SELECT_SET, offset2d, nullptr, count2d, nullptr);
                    H5Dwrite(v9_max_id, max_type, max_space, v9_space, H5P_DEFAULT, (void*)buf);
                }
            }
            if (back_id > -1)
            {
                if (H5Dread(back_id, max_type, max_space, max_space, H5P_DEFAULT, (void*)buf) > -1)
                {
                    offset2d[0] = 4;
                    H5Sselect_hyperslab(v9_space, H5S_SELECT_SET, offset2d, nullptr, count2d, nullptr);
                    H5Dwrite(v9_max_id, max_type, max_space, v9_space, H5P_DEFAULT, (void*)buf);
                }
            }

            if (dbuf != nullptr)
            {
                delete[] dbuf;
            }
            if(fbuf != nullptr)
            { 
                delete [] fbuf;
            }
		}
	}
	
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
    if (chan_names > -1)
    {
         _global_close_map.push({chan_names, H5O_DATASET });
    
        
        hid_t chan_space = H5Dget_space(chan_names);
        _global_close_map.push({chan_space, H5O_DATASPACE });
        int num_chan;
        if(chan_names > -1)
        {
            hsize_t chan_size = 1;
            H5Sget_simple_extent_dims(chan_space, &chan_size, nullptr);
            num_chan = chan_size; //num channel names
        }
        //hid_t quant_space = H5Screate_simple(3, &quant_dims[0], &quant_dims[0]);

        //Channel Units are a 4 x channels so we can't do a hardlink
        // the 4 are SR_current, US_IC, DS_IC, and cts/s
        
        
        hsize_t unit_dims[2];
        hsize_t offset_dims[2] = { 0,0 };
        unit_dims[0] = 4;
        unit_dims[1] =num_chan;
        hid_t ch_unit_id, units_space;
            
        if (false == _open_h5_dataset("/MAPS/channel_units", filetype, file_id, 2, &unit_dims[0], &unit_dims[0], ch_unit_id, units_space))
        {
            logE << "Error creating " << "/MAPS/channel_units" << "\n";
        }
        if (ch_unit_id > -1)
        {
            hsize_t mem_dims[1] = { 1 };
            hsize_t count_2d[2] = { 1,1 };
            hid_t mem_space = H5Screate_simple(1, &mem_dims[0], &mem_dims[0]);
            _global_close_map.push({mem_space, H5O_DATASPACE });
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
        }
        else
        {
            logW << "Couldn't create /MAPS/channel_units" << "\n";
        }
        

        if (H5Gget_objinfo(file_id, "/MAPS/XRF_roi", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/ROI/Counts_Per_Sec", 0, NULL) >= 0)
        {
            H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/ROI/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_roi", H5P_DEFAULT, H5P_DEFAULT);
        }
        if (H5Gget_objinfo(file_id, "/MAPS/XRF_roi", 0, NULL) >= 0)
        {
            _add_v9_quant(file_id, chan_names, chan_space, num_chan, "ROI", "/MAPS/XRF_roi_quant");
        }

        if (H5Gget_objinfo(file_id, "/MAPS/XRF_roi_plus", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec", 0, NULL) >= 0)
        {
            H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/NNLS/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_roi_plus", H5P_DEFAULT, H5P_DEFAULT);
        }
        if (H5Gget_objinfo(file_id, "/MAPS/XRF_roi_plus", 0, NULL) >= 0)
        {
            _add_v9_quant(file_id, chan_names, chan_space, num_chan, "NNLS", "/MAPS/XRF_roi_plus_quant");
        }

        if (H5Gget_objinfo(file_id, "/MAPS/XRF_fits", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec", 0, NULL) >= 0)
        {
            H5Lcreate_hard(file_id, "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec", H5L_SAME_LOC, "/MAPS/XRF_fits", H5P_DEFAULT, H5P_DEFAULT);
        }
        if (H5Gget_objinfo(file_id, "/MAPS/XRF_fits", 0, NULL) >= 0)
        {
            _add_v9_quant(file_id, chan_names, chan_space, num_chan, STR_FIT_GAUSS_MATRIX, "/MAPS/XRF_fits_quant");
        }

        // create links /MAPS/make_maps_conf/nbs1832/us_amp and ds_amp
        hid_t conf_id = H5Gcreate(file_id, "/MAPS/make_maps_conf", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (conf_id < 0)
        {
            conf_id = H5Gopen(file_id, "/MAPS/make_maps_conf", H5P_DEFAULT);
        }
        if (conf_id >-1)
        {
            _global_close_map.push({conf_id, H5O_GROUP });
            hid_t e_id = H5Gcreate(conf_id, "nbs1832", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);            
            if (e_id > -1)
            {
                _global_close_map.push({e_id, H5O_GROUP });
            }
        }

        if (H5Gget_objinfo(file_id, "/MAPS/make_maps_conf/nbs1832/us_amp", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/Quantification/Standard0/Scalers/us_amp", 0, NULL) >= 0)
        {
            H5Lcreate_hard(file_id, "/MAPS/Quantification/Standard0/Scalers/us_amp", H5L_SAME_LOC, "/MAPS/make_maps_conf/nbs1832/us_amp", H5P_DEFAULT, H5P_DEFAULT);
        }
        if (H5Gget_objinfo(file_id, "/MAPS/make_maps_conf/nbs1832/ds_amp", 0, NULL) < 0 && H5Gget_objinfo(file_id, "/MAPS/Quantification/Standard0/Scalers/ds_amp", 0, NULL) >= 0)
        {
            H5Lcreate_hard(file_id, "/MAPS/Quantification/Standard0/Scalers/ds_amp", H5L_SAME_LOC, "/MAPS/make_maps_conf/nbs1832/ds_amp", H5P_DEFAULT, H5P_DEFAULT);
        }
        _add_extra_pvs(file_id, "/MAPS");
    }

    //change version to 9
    hid_t version_id = H5Dopen(file_id, "/MAPS/version", H5P_DEFAULT);
    if(version_id > -1)
    {
        _global_close_map.push({version_id, H5O_DATASET });
    }
    hid_t ver_space = H5Dget_space(version_id);
    if(ver_space > -1)
    {
        _global_close_map.push({ver_space, H5O_DATASPACE });
    }
    hid_t ver_type = H5Dget_type(version_id);
    if(ver_type > -1)
    {
        _global_close_map.push({ver_type, H5O_DATATYPE });
    }
    if (H5Tget_class(ver_type) == H5T_FLOAT) 
    {
        float version = 9.0;
        H5Dwrite(version_id, ver_type, ver_space, ver_space, H5P_DEFAULT, (void*)&version);
    }
    else
    {
        double version = 9.0;
        H5Dwrite(version_id, ver_type, ver_space, ver_space, H5P_DEFAULT, (void*)&version);
    }
    if (H5Gget_objinfo(file_id, "/version", 0, NULL) < 0)
    {
        H5Lcreate_hard(file_id, "/MAPS/version", H5L_SAME_LOC, "/version", H5P_DEFAULT, H5P_DEFAULT);
    }

    _close_h5_objects(_global_close_map);
}

//-----------------------------------------------------------------------------

void HDF5_IO::_add_v9_scalers(hid_t file_id)
{

    std::map<std::string, int> scaler_map;
    hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
    H5Tset_size(filetype, 256);

    hid_t names_id, values_id, units_id;
    if (false == _open_h5_object(names_id, H5O_DATASET, _global_close_map, "/MAPS/Scalers/Names", file_id, true, false))
    {
        return;
    }
    if (false == _open_h5_object(values_id, H5O_DATASET, _global_close_map, "/MAPS/Scalers/Values", file_id, true, false))
    {
        return;
    }
    if (false == _open_h5_object(units_id, H5O_DATASET, _global_close_map, "/MAPS/Scalers/Units", file_id, true, false))
    {
        return;
    }
     
    hid_t name_space = H5Dget_space(names_id);
    _global_close_map.push({ name_space , H5O_DATASPACE });
    hid_t value_space = H5Dget_space(values_id);
    _global_close_map.push({ value_space , H5O_DATASPACE });
    hid_t unit_space = H5Dget_space(units_id);
    _global_close_map.push({ unit_space , H5O_DATASPACE });


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
            scaler_name_str.erase(std::remove_if(scaler_name_str.begin(), scaler_name_str.end(), ::isspace), scaler_name_str.end());
            size_t c_idx = scaler_name_str.find(':');

            if (c_idx == std::string::npos && scaler_name_str.length() > 0)
            {
                scaler_map[scaler_name_str] = i;
            }
        }
    }

    count_1d[0] = { scaler_map.size() };
    hid_t new_name_space, new_unit_space;

    count_3d[0] = { scaler_map.size() };
    hid_t new_value_space;// = H5Screate_simple(3, &count_3d[0], &new_max_3d[0]);
    hid_t new_names_id;
    hid_t new_units_id;
    hid_t new_values_id;
    if (false == _open_h5_dataset("/MAPS/scaler_names", filetype, file_id, 1, &count_1d[0], &count_1d[0], new_names_id, new_name_space))
    {
        logW << "Error creating /MAPS/scaler_names\n";
    }
    if (false == _open_h5_dataset("/MAPS/scaler_units", filetype, file_id, 1, &count_1d[0], &count_1d[0], new_units_id, new_unit_space))
    {
        logW << "Error creating /MAPS/scaler_units\n";
    }
    if (false == _open_h5_dataset("/MAPS/scalers", H5T_NATIVE_FLOAT, file_id, 3, &count_3d[0], &count_3d[0], new_values_id, new_value_space))
    {
        logW << "Error creating /MAPS/scalers\n";
    }


    if (new_names_id < 0 || new_units_id < 0 || new_values_id < 0)
    {
        logW << "Could not open /MAPS/scalers _names, or _units dataset to add v9 layout\n";
        return;
    }

    count_1d[0] = 1;
    count_3d[0] = 1;

    hid_t value_mem_space;
    _create_memory_space(3, &count_3d[0], value_mem_space);
    
    data_struct::ArrayXXr<float> tmp_values;
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
        H5Sselect_hyperslab(new_unit_space, H5S_SELECT_SET, new_offset_1d, nullptr, count_1d, nullptr);
        

        char tmp_char_name[256] = { 0 };
        itr.first.copy(tmp_char_name, 254);
        H5Dwrite(new_names_id, filetype, mem_space_1d, new_name_space, H5P_DEFAULT, (void*)&tmp_char_name[0]);

        char tmp_char[256] = { 0 };
        if (H5Dread(units_id, filetype, mem_space_1d, name_space, H5P_DEFAULT, (void*)tmp_char) > -1)
        {
            H5Dwrite(new_units_id, filetype, mem_space_1d, new_unit_space, H5P_DEFAULT, (void*)tmp_char);
        }

        H5Sselect_hyperslab(value_space, H5S_SELECT_SET, offset_3d, nullptr, count_3d, nullptr);
        if (H5Dread(values_id, H5T_NATIVE_FLOAT, value_mem_space, value_space, H5P_DEFAULT, (void*)tmp_values.data()) > -1)
        {
            H5Sselect_hyperslab(new_value_space, H5S_SELECT_SET, new_offset_3d, nullptr, count_3d, nullptr);
            H5Dwrite(new_values_id, H5T_NATIVE_FLOAT, value_mem_space, new_value_space, H5P_DEFAULT, (void*)tmp_values.data());
        }

        i++;
    }
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
    if(fits_grp > -1)
    {
        _global_close_map.push({fits_grp, H5O_GROUP});
    }



    hid_t filetype = H5Tcopy(H5T_FORTRAN_S1);
    _global_close_map.push({filetype, H5O_DATATYPE});

    H5Tset_size(filetype, 256);
    hsize_t count [1] = {1};
    hid_t dset_id;
    hid_t dataspace_id = H5Screate_simple (1, count, nullptr);
    _global_close_map.push({dataspace_id, H5O_DATASPACE});
    hid_t exchange_id;
    if (false == _open_or_create_group(exhange_str, file_id, exchange_id))
    {
        return false;
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
        dset_id = H5Dcreate(file_id, str_desc.c_str(), filetype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if(dset_id > -1)
        {
            _global_close_map.push({dset_id, H5O_DATASET});
            norm_desc.copy(desc, 256);
            H5Dwrite(dset_id, filetype, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&desc);
        }

        dset_id = H5Dopen(file_id, "/MAPS/XRF_Analyzed/Fitted/Counts_Per_Sec", H5P_DEFAULT);
        if(dset_id > -1)
        {
            _global_close_map.push({dset_id, H5O_DATASET});
        }
        //hid_t chan_units_id = H5Dopen(file_id, "/MAPS/XRF_Analyzed/Fitted/Channel_Units", H5P_DEFAULT);
        hid_t chan_names_id = H5Dopen(file_id, "/MAPS/XRF_Analyzed/Fitted/Channel_Names", H5P_DEFAULT);
        if(chan_names_id > -1)
        {
            _global_close_map.push({chan_names_id, H5O_DATASET});
        }
        hid_t scaler_dset_id = H5Dopen(file_id, "/MAPS/Scalers/Values", H5P_DEFAULT);
        if(scaler_dset_id > -1)
        {
            _global_close_map.push({scaler_dset_id, H5O_DATASET});
        }
        hid_t scaler_units_id = H5Dopen(file_id, "/MAPS/Scalers/Units", H5P_DEFAULT);
        if(scaler_units_id > -1)
        {
            _global_close_map.push({scaler_units_id, H5O_DATASET});
        }
        hid_t scaler_names_id = H5Dopen(file_id, "/MAPS/Scalers/Names", H5P_DEFAULT);
        if(scaler_names_id > -1)
        {
            _global_close_map.push({scaler_names_id, H5O_DATASET});
        }
        hid_t ds_ic_quant_id = H5Dopen(file_id, "/MAPS/Quantification/Calibration/Fitted/Calibration_Curve_DS_IC", H5P_DEFAULT);
        if(ds_ic_quant_id > -1)
        {
            _global_close_map.push({ds_ic_quant_id, H5O_DATASET});
        }

        hid_t quant_space = H5Dget_space(ds_ic_quant_id);
        if(quant_space > -1)
        {
            _global_close_map.push({quant_space, H5O_DATASPACE});
        }
        hid_t quant_type = H5Dget_type(ds_ic_quant_id);
        if(quant_type > -1)
        {
            _global_close_map.push({quant_type, H5O_DATATYPE});
        }
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
            hid_t image_names_dset_id;
            hid_t image_units_dset_id;
            hid_t unit_single_space;

            hid_t chan_type = H5Dget_type(dset_id);
            //hid_t scalername_type = H5Dget_type(scaler_names_id);

            hid_t chan_space = H5Dget_space(dset_id);
            hid_t chan_name_space = H5Dget_space(chan_names_id);
            hid_t scaler_space = H5Dget_space(scaler_dset_id);

            hid_t memtype = H5Tcopy(H5T_C_S1);
            H5Tset_size(memtype, 255);

            H5Sget_simple_extent_dims(chan_space, &chan_dims[0], nullptr);
            H5Sget_simple_extent_dims(scaler_space, &scaler_dims[0], nullptr);

            image_dims_single[0] = { 1 };
            hid_t readwrite_single_space = H5Screate_simple(1, &image_dims_single[0], &image_dims_single[0]);

            image_dims_single[0] = chan_dims[0] + scaler_dims[0];
            image_dims[0] = chan_dims[0] + scaler_dims[0];
            image_dims[1] = chan_dims[1];
            image_dims[2] = chan_dims[2];
            hid_t image_dset_id, image_space, image_single_space;

            if (false == _open_h5_dataset(exchange_images, chan_type, file_id, 3, &image_dims[0], &image_dims[0], image_dset_id, image_space))
            {
                logE << "Error creating " << exchange_images << "\n";
                return false;
            }
            if (false == _open_h5_dataset(exchange_image_names, filetype, file_id, 1, &image_dims_single[0], &image_dims_single[0], image_names_dset_id, image_single_space))
            {
                logE << "Error creating " << exchange_image_names << "\n";
                return false;
            }
            if (false == _open_h5_dataset(exchange_image_units, filetype, file_id, 1, &image_dims_single[0], &image_dims_single[0], image_units_dset_id, unit_single_space))
            {
                logE << "Error creating " << exchange_image_units << "\n";
                return false;
            }
            

            double *data = new double[chan_dims[1] * chan_dims[2]];
            double*ds_ic_data = new double[chan_dims[1] * chan_dims[2]];
            for (int z = 0; z < (chan_dims[1] * chan_dims[2]); z++)
            {
                data[z] = 0.;
                ds_ic_data[z] = 0.;
            }
            std::string scaler_name_str;
            char char_data[256]={0};
            char char_ug_data[256]="ug/cm2";
            int k =0;

            double quant_value = 1.0;

            std::transform(normalize_scaler.begin(), normalize_scaler.end(), normalize_scaler.begin(), [](unsigned char c) { return std::tolower(c); });

            image_dims[0] = 1;
            hid_t readwrite_space = H5Screate_simple(3, &image_dims[0], &image_dims[0]);
            image_dims_single[0] = 1;
            // save scalers first
            for(hsize_t i=0; i < scaler_dims[0]; i++)
            {
                offset[0] = i;
                offset_single[0] = i;
                offset_image[0] = k;
                k++;
                H5Sselect_hyperslab (image_space, H5S_SELECT_SET, offset, nullptr, image_dims, nullptr);
                //read write values
                hid_t status = H5Dread(scaler_dset_id, H5T_NATIVE_DOUBLE, readwrite_space, image_space, H5P_DEFAULT, (void*)&data[0]);
                if(status > -1)
                {
                    status = H5Dwrite(image_dset_id, H5T_NATIVE_DOUBLE, readwrite_space, image_space, H5P_DEFAULT, (void*)&data[0]);
                    if (status == -1)
                    {
                        logW << "Issue saving scaler index " << offset[0] << "\n";
                    }
                }

                //read write names
                H5Sselect_hyperslab (image_single_space, H5S_SELECT_SET, offset_single, nullptr, image_dims_single, nullptr);
                status = H5Dread(scaler_names_id, memtype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                if(status > -1)
                {
                    H5Dwrite(image_names_dset_id, memtype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                }


                scaler_name_str = std::string(char_data, 256);
                scaler_name_str.erase(std::find(scaler_name_str.begin(), scaler_name_str.end(), '\0'), scaler_name_str.end());
                scaler_name_str.erase(std::remove_if(scaler_name_str.begin(), scaler_name_str.end(), ::isspace), scaler_name_str.end());
                //to lower
                std::transform(scaler_name_str.begin(), scaler_name_str.end(), scaler_name_str.begin(), [](unsigned char c) { return std::tolower(c); });
                
                if(normalize_scaler.compare(scaler_name_str) == 0)
                {
                    for(hsize_t z=0; z < (chan_dims[1] * chan_dims[2]); z++)
                    {
                        ds_ic_data[z] = data[z];
                    }
                }

                //read write units
                status = H5Dread(scaler_units_id, memtype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                if(status > -1)
                {
                    H5Dwrite(image_units_dset_id, memtype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
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
                hid_t status = H5Dread(chan_names_id, memtype, readwrite_single_space, chan_name_space, H5P_DEFAULT, (void*)&char_data[0]);
                if(status > -1)
                {
                    H5Dwrite(image_names_dset_id, memtype, readwrite_single_space, image_single_space, H5P_DEFAULT, (void*)&char_data[0]);
                }

                // get quantification for ds_ic and store in quant_value
                if(ds_ic_quant_id > -1)
                {
                    std::string chan_name_str = std::string(char_data, 256);
                    chan_name_str.erase(std::find(chan_name_str.begin(), chan_name_str.end(), '\0'), chan_name_str.end());
                    chan_name_str.erase(std::remove_if(chan_name_str.begin(), chan_name_str.end(), ::isspace), chan_name_str.end());
                    data_struct::Element_Info<double>* element = data_struct::Element_Info_Map<double>::inst()->get_element(chan_name_str);
                    if(element != nullptr)
                    {
                        offset_quant[1] = element->number - 1;

                        H5Sselect_hyperslab (quant_space, H5S_SELECT_SET, offset_quant, nullptr, count_quant, nullptr);
                        status = H5Dread(ds_ic_quant_id, H5T_NATIVE_DOUBLE, readwrite_single_space, quant_space, H5P_DEFAULT, (void*)&quant_value);
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
                status = H5Dread(dset_id, H5T_NATIVE_DOUBLE, readwrite_space, chan_space, H5P_DEFAULT, (void*)&data[0]);
                if(status > -1)
                {
                    for(hsize_t z=0; z < (chan_dims[1] * chan_dims[2]); z++)
                    {
                        data[z] = data[z] / quant_value / ds_ic_data[z];
                    }
                    status = H5Dwrite(image_dset_id, H5T_NATIVE_DOUBLE, readwrite_space, image_space, H5P_DEFAULT, (void*)&data[0]);
                    if (status == -1)
                    {
                        logW << "Issue saving data index " << offset[0] << "\n";
                    }
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
            if(chan_names > -1)
            {
                _global_close_map.push({chan_names, H5O_DATASET});
                hid_t chan_space = H5Dget_space(chan_names);
                _global_close_map.push({chan_space, H5O_DATASPACE});
                hsize_t chan_size = 1;
                H5Sget_simple_extent_dims(chan_space, &chan_size, nullptr);
                quant_dims[2] = chan_size; //num channel names
            
                hid_t quant_space;
                _create_memory_space(3, &quant_dims[0], quant_space);
                _add_v9_quant(file_id, chan_names, chan_space, quant_dims[2], fits_link, exhange_str+"/quant");
            }
        }

        //Save description
        if (false == _open_h5_dataset(str_desc, filetype, file_id, 1, count, count, dset_id, dataspace_id))
        {
            logE << "Error creating " << str_desc << "\n";
            return false;
        }
        fits_link.copy(desc, 256);
        H5Dwrite(dset_id, filetype, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&desc);
        
    }




    //Add version dataset
    float save_val = HDF5_EXCHANGE_VERSION;
    
    if (false == _open_h5_dataset(str_version, H5T_NATIVE_DOUBLE, file_id, 1, count, count, dset_id, dataspace_id))
    {
        logE << "Error creating " << str_version << "\n";
        return false;
    }
    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, dataspace_id, dataspace_id, H5P_DEFAULT, (void*)&save_val);
    

    //Add theta
    if( H5Lcreate_hard(file_id, str_scan_theta.c_str(), H5L_SAME_LOC, str_theta.c_str(), H5P_DEFAULT, H5P_DEFAULT) < 0)
    {
        logW  << "Couldn't create soft link for"<< str_theta <<  "\n";
    }

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

    _close_h5_objects(_global_close_map);

    _cur_file_id = file_id;
    end_save_seq();
    logI<<"closing file"<<"\n";

    _cur_file_id = saved_file_id;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

herr_t h5_ext_file_info(hid_t loc_id, const char *name, [[maybe_unused]] const H5L_info2_t *linfo, void *opdata)
{
    hid_t group;


    logI<< "external link to "<< name << "\n";
    /*
     * Open the group using its name.
     */
    //hdf5_io_det_spec_line<T_real>* det_spec = (hdf5_io_det_spec_line<T_real>*)opdata;
    std::map<std::string, hid_t>* ext_links = (std::map<std::string, hid_t>*)opdata;

    if(ext_links !=  nullptr)
    {
        group = H5Gopen2(loc_id, name, H5P_DEFAULT);
        if(group >= 0)
        {
            std::string str_name = std::string(name);
            ext_links->insert({str_name, group});
        }
    }
    return 0;
}

//-----------------------------------------------------------------------------

} //end namespace file
}// end namespace io
