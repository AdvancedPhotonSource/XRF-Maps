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



#include "spectra_net_source.h"

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

Spectra_Net_Source::Spectra_Net_Source(data_struct::xrf::Analysis_Job* analysis_job) : Source<data_struct::xrf::Stream_Block*>()
{
    _running = false;
    _receiver = nullptr;
    _analysis_job = analysis_job;
//    _cb_function = std::bind(&Spectra_Net_Source::cb_load_spectra_data, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7);
    _receiver = new zmq::socket_t(_context, ZMQ_PULL);
    _receiver->connect("tcp://localhost:5557");
}

//-----------------------------------------------------------------------------

Spectra_Net_Source::~Spectra_Net_Source()
{
    if (_receiver != nullptr)
    {
        delete _receiver;
    }
    _receiver = nullptr;
}

// ----------------------------------------------------------------------------
/*
void Spectra_Net_Source::cb_load_spectra_data(size_t row, size_t col, size_t height, size_t width, size_t detector_num, data_struct::xrf::Spectra* spectra, void* user_data)
{

    if(_output_callback_func != nullptr)
    {
        struct data_struct::xrf::Analysis_Sub_Struct* cp = _analysis_job->get_sub_struct(detector_num);

        data_struct::xrf::Stream_Block * stream_block = new data_struct::xrf::Stream_Block(row, col, height, width);
        stream_block->init_fitting_blocks(&(cp->fit_routines), &(cp->fit_params_override_dict.elements_to_fit));
        stream_block->spectra = spectra;
        stream_block->model = cp->model;
        stream_block->dataset_directory = _current_dataset_directory;
        stream_block->dataset_name = _current_dataset_name;
        stream_block->detector_number = detector_num;

        _output_callback_func(stream_block);
    }

}
*/
// ----------------------------------------------------------------------------

void Spectra_Net_Source::run()
{
    _running = true;

    while (_running)
    {
        zmq::message_t message;
        //int workload;           //  Workload in msecs

        _receiver->recv(&message);
        std::string smessage(static_cast<char*>(message.data()), message.size());

        if(_output_callback_func != nullptr)
        {
            //parse these from message
            size_t row;
            size_t col;
            size_t height;
            size_t width;
            size_t detector_num;
            data_struct::xrf::Spectra* spectra;

            struct data_struct::xrf::Analysis_Sub_Struct* cp = _analysis_job->get_sub_struct(detector_num);

            data_struct::xrf::Stream_Block * stream_block = new data_struct::xrf::Stream_Block(row, col, height, width);
            stream_block->init_fitting_blocks(&(cp->fit_routines), &(cp->fit_params_override_dict.elements_to_fit));
            stream_block->spectra = spectra;
            stream_block->model = cp->model;
            stream_block->dataset_directory = _current_dataset_directory;
            stream_block->dataset_name = _current_dataset_name;
            stream_block->detector_number = detector_num;

            _output_callback_func(stream_block);
        }
    }
}

//-----------------------------------------------------------------------------
/*
bool Spectra_Net_Source::_load_spectra_volume_with_callback(std::string dataset_directory,
                                                                 std::string dataset_file,
                                                                 size_t detector_num_start,
                                                                 size_t detector_num_end,
                                                                 io::file::IO_Callback_Func_Def callback_fun)
{
    //Dataset importer
    io::file::MDA_IO mda_io;
    //data_struct::xrf::Detector detector;
    std::string tmp_dataset_file = dataset_file;

    logit<<"Loading dataset "<<dataset_directory+"mda/"+dataset_file<<" detectors "<<detector_num_start<<":"<<detector_num_end<<std::endl;

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    bool hasNetcdf = false;
    bool hasHdf = false;
    std::string file_middle = ""; //_2xfm3_ or dxpM...
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

    _current_dataset_directory = new std::string(dataset_directory);
    _current_dataset_name = new std::string(dataset_file);
    //load spectra
    if (false == mda_io.load_spectra_volume_with_callback(dataset_directory+"mda/"+dataset_file,
                                                        detector_num_start,
                                                        detector_num_end,
                                                        hasNetcdf | hasHdf,
                                                        _analysis_job,
                                                        callback_fun,
                                                        nullptr) )
    {
        logit<<"Error load spectra "<<dataset_directory+"mda/"+dataset_file<<std::endl;
        delete _current_dataset_directory;
        delete _current_dataset_name;
        return false;
    }
    else
    {
        if(hasNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                int row_size = mda_io.rows();
                int col_size = mda_io.cols();
                for(size_t i=0; i<row_size; i++)
                {
                    full_filename = dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + std::to_string(i) + ".nc";
                    //todo: add verbose option
                    //logit<<"Loading file "<<full_filename<<std::endl;
                    io::file::NetCDF_IO::inst()->load_spectra_line_with_callback(full_filename, detector_num_start, detector_num_end, i, row_size, col_size, callback_fun, nullptr);
                }
            }
            else
            {
                logit<<"Did not find netcdf files "<<dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc"<<std::endl;
                //return false;
            }
        }
        else if (hasHdf)
        {
            io::file::HDF5_IO::inst()->load_spectra_volume_with_callback(dataset_directory + "flyXRF.h5/" + tmp_dataset_file + file_middle + "0.h5", detector_num_start, detector_num_end, callback_fun, nullptr);
        }

    }

    //move to stream_block so saver can deal with it
    mda_io.unload();
    logit<<"Finished Loading dataset "<<dataset_directory+"mda/"+dataset_file<<" detectors "<<detector_num_start<<":"<<detector_num_end<<std::endl;
    return true;
}
*/
} //namespace xrf
} //namespace workflow
