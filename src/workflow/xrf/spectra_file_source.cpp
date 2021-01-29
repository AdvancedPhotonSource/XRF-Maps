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

/// Initial Author <2017>: Arthur Glowacki



#include "spectra_file_source.h"
#include "io/file/hl_file_io.h"
#include "core/mem_info.h"

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

Spectra_File_Source::Spectra_File_Source() : Source<data_struct::Stream_Block*>()
{
    _analysis_job = nullptr;
    _current_dataset_directory = nullptr;
    _current_dataset_name = nullptr;
	_max_num_stream_blocks = -1;
    _cb_function = std::bind(&Spectra_File_Source::cb_load_spectra_data, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7);
}

//-----------------------------------------------------------------------------

Spectra_File_Source::Spectra_File_Source(data_struct::Analysis_Job* analysis_job) : Source<data_struct::Stream_Block*>()
{
    _analysis_job = analysis_job;
    _current_dataset_directory = nullptr;
    _current_dataset_name = nullptr;
    _init_fitting_routines = true;
	_max_num_stream_blocks = -1;
    _cb_function = std::bind(&Spectra_File_Source::cb_load_spectra_data, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7);
}

//-----------------------------------------------------------------------------

Spectra_File_Source::~Spectra_File_Source()
{

}

// ----------------------------------------------------------------------------

bool Spectra_File_Source::load_netcdf_line(std::string dirpath,
										   std::string filename,
                                           const std::vector<size_t>& detector_num_arr,
                                           size_t row,
                                           size_t row_size,
                                           size_t col_size)
{
	bool retVal;
	_current_dataset_directory = new std::string(dirpath);
	_current_dataset_name = new std::string(filename);
	retVal = io::file::NetCDF_IO::inst()->load_spectra_line_with_callback(dirpath+filename, detector_num_arr, row, row_size, col_size, _cb_function, nullptr);
	delete _current_dataset_directory;
	delete _current_dataset_name;
	return retVal;
}

// ----------------------------------------------------------------------------

data_struct::Stream_Block* Spectra_File_Source::_alloc_stream_block(int detector, size_t row, size_t col, size_t height, size_t width, size_t spectra_size)
{
	if (_max_num_stream_blocks == -1)
	{
		_max_num_stream_blocks = _analysis_job->mem_limit / (spectra_size * sizeof(real_t));
	}
	return new data_struct::Stream_Block(detector, row, col, height, width);
}

// ----------------------------------------------------------------------------

void Spectra_File_Source::cb_load_spectra_data(size_t row, size_t col, size_t height, size_t width, size_t detector_num, data_struct::Spectra* spectra, void* user_data)
{

    if(_output_callback_func != nullptr)
    {
		data_struct::Stream_Block * stream_block = _alloc_stream_block(detector_num, row, col, height, width, spectra->size());

        if(_init_fitting_routines && _analysis_job != nullptr)
        {
            _analysis_job->init_fit_routines(spectra->size());

            struct data_struct::Detector* cp = _analysis_job->get_detector(detector_num);

            if(cp == nullptr)
            {
                cp = _analysis_job->get_first_detector();
            }
            if(cp != nullptr)
            {
                stream_block->init_fitting_blocks(&(cp->fit_routines), &(cp->fit_params_override_dict.elements_to_fit));
                stream_block->model = cp->model;
            }
        }
        if(_analysis_job != nullptr)
        {
            stream_block->optimize_fit_params_preset = _analysis_job->optimize_fit_params_preset;
            stream_block->theta = _analysis_job->theta;
        }

        stream_block->spectra = spectra;
        stream_block->dataset_directory = _current_dataset_directory;
        stream_block->dataset_name = _current_dataset_name;

        _output_callback_func(stream_block);
    }
    else
    {
        delete spectra;
    }

}

// ----------------------------------------------------------------------------

void Spectra_File_Source::run()
{
    if(_analysis_job == nullptr)
    {
        logE<<"Class was not constructed with Analysis job. Don't know what to run?\n";
        return;
    }

	// if no memory limit is set, then query the system memory and make a stream block queue
	long long total_mem = get_available_mem();
	if (_analysis_job->mem_limit == -1)
	{
		_analysis_job->mem_limit = total_mem;
	}
	else
	{
		_analysis_job->mem_limit = std::min(_analysis_job->mem_limit, total_mem);
	}

    _netcdf_files = io::find_all_dataset_files(_analysis_job->dataset_directory + "flyXRF"+ DIR_END_CHAR, "_0.nc");
    _bnp_netcdf_files = io::find_all_dataset_files(_analysis_job->dataset_directory + "flyXRF"+ DIR_END_CHAR, "_001.nc");
    _hdf_files = io::find_all_dataset_files(_analysis_job->dataset_directory + "flyXRF.h5"+ DIR_END_CHAR, "_0.h5");

    for(std::string dataset_file : _analysis_job->dataset_files)
    {
        //load xfm dataset
        if (false == _load_spectra_volume_with_callback(_analysis_job->dataset_directory, dataset_file, _analysis_job->detector_num_arr, _cb_function))
        {
            logW << "Skipping dataset_file " << dataset_file << "\n";
            continue;
        }

		//send end of file stream block
		data_struct::Stream_Block* end_block = new data_struct::Stream_Block(-1, -1, -1, -1, -1);
		end_block->dataset_directory = _current_dataset_directory;
		end_block->dataset_name = _current_dataset_name;
		end_block->del_str_ptr = true;
		_output_callback_func(end_block);

		_current_dataset_directory = nullptr;
		_current_dataset_name = nullptr;
    }
}

//-----------------------------------------------------------------------------

bool Spectra_File_Source::_load_spectra_volume_with_callback(std::string dataset_directory,
                                                                 std::string dataset_file,
                                                                 const std::vector<size_t>& detector_num_arr,
																 data_struct::IO_Callback_Func_Def callback_fun)
{
    //Dataset importer
    io::file::MDA_IO mda_io;
    //data_struct::Detector detector;
    std::string tmp_dataset_file = dataset_file;

    logI<<"Loading dataset "<<dataset_directory+"mda"+ DIR_END_CHAR +dataset_file<<" \n";

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    bool hasNetcdf = false;
    bool hasBnpNetcdf = false;
    bool hasHdf = false;
    std::string file_middle = ""; //_2xfm3_ or dxpM...
    std::string bnp_netcdf_base_name = "bnp_fly_";
    for(auto &itr : _netcdf_files)
    {
        if (itr.find(tmp_dataset_file) == 0)
        {
            size_t slen = (itr.length()-4) - tmp_dataset_file.length();
            file_middle = itr.substr(tmp_dataset_file.length(), slen);
            hasNetcdf = true;
            break;
        }
    }
    if (hasNetcdf == false)
    {
        int idx = static_cast<int>(tmp_dataset_file.find("bnp_fly"));
        if (idx == 0)
        {
            std::string footer = tmp_dataset_file.substr(7, tmp_dataset_file.length() - 7);
            int file_index = std::atoi(footer.c_str());
            file_middle = std::to_string(file_index);
            bnp_netcdf_base_name = "bnp_fly_"+ file_middle + "_";
            for(auto &itr : _bnp_netcdf_files)
            {
                if (itr.find(bnp_netcdf_base_name) == 0)
                {
                    hasBnpNetcdf = true;
                    break;
                }
            }
        }
    }
    if (hasNetcdf == false && hasBnpNetcdf == false)
    {
        for(auto &itr : _hdf_files)
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length()-4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasHdf = true;
                break;
            }
        }
    }

    //TODO: add confocal and emd streaming
    //// load emd dataset
    //if(false == io::file::HDF5_IO::inst()->load_spectra_volume_emd_with_callback(_analysis_job->dataset_directory + dataset_file, _analysis_job->detector_num_arr, _cb_function, nullptr))
    //{
    //    logW << "Skipping dataset_file " << dataset_file << "\n";
    //    continue;
    //}

    size_t row_size = 0;
    size_t col_size = 0;
    _current_dataset_directory = new std::string(dataset_directory);
    _current_dataset_name = new std::string(dataset_file);
    //load spectra
    if (false == mda_io.load_spectra_volume_with_callback(dataset_directory+"mda"+ DIR_END_CHAR +dataset_file,
                                                        detector_num_arr,
                                                        hasNetcdf | hasBnpNetcdf | hasHdf,
                                                        _analysis_job,
                                                        row_size,
                                                        col_size,
                                                        callback_fun,
                                                        nullptr) )
    {
        logE<<"load spectra "<<dataset_directory+"mda"+ DIR_END_CHAR +dataset_file<<"\n";
        delete _current_dataset_directory;
        delete _current_dataset_name;
        return false;
    }
    else
    {
        if(hasNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF"+ DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for(int i=0; i<row_size; i++)
                {
                    full_filename = dataset_directory + "flyXRF"+ DIR_END_CHAR + tmp_dataset_file + file_middle + std::to_string(i) + ".nc";
                    //todo: add verbose option
                    //logI<<"Loading file "<<full_filename<<"\n";
                    io::file::NetCDF_IO::inst()->load_spectra_line_with_callback(full_filename, detector_num_arr, i, row_size, col_size, callback_fun, nullptr);
                }
            }
            else
            {
                logE<<"Did not find netcdf files "<<dataset_directory + "flyXRF"+ DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc"<<"\n";
                //return false;
            }
        }
        else if(hasBnpNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF"+ DIR_END_CHAR + bnp_netcdf_base_name + "001.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for(int i=0; i<row_size; i++)
                {
                    std::string row_idx_str = std::to_string(i+1);
                    int num_prepended_zeros = 3 - static_cast<int>(row_idx_str.size()); // 3 chars for num of rows, prepened with zeros if less than 100
                    std::string row_idx_str_full = "";
                    for(int z=0; z<num_prepended_zeros; z++)
                    {
                        row_idx_str_full += "0";
                    }
                    row_idx_str_full += row_idx_str;
                    full_filename = dataset_directory + "flyXRF"+ DIR_END_CHAR + bnp_netcdf_base_name + row_idx_str_full + ".nc";
                    io::file::NetCDF_IO::inst()->load_spectra_line_with_callback(full_filename, detector_num_arr, i, row_size, col_size, callback_fun, nullptr);
                }
            }
            else
            {
                logE<<"Did not find netcdf files "<<dataset_directory + "flyXRF"+ DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc"<<"\n";
                //return false;
            }
        }
        else if (hasHdf)
        {
            io::file::HDF5_IO::inst()->load_spectra_volume_with_callback(dataset_directory + "flyXRF.h5"+ DIR_END_CHAR + tmp_dataset_file + file_middle + "0.h5", detector_num_arr, callback_fun, nullptr);
        }

    }

    //move to stream_block so saver can deal with it
    mda_io.unload();
    logI<<"Finished Loading dataset "<<dataset_directory+"mda/"+dataset_file<<"\n";
    return true;
}

} //namespace xrf
} //namespace workflow
