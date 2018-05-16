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



#include "spectra_stream_saver.h"

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

Spectra_Stream_Saver::Spectra_Stream_Saver() : Sink<data_struct::Stream_Block*>()
{
    _callback_func = std::bind(&Spectra_Stream_Saver::save_stream, this, std::placeholders::_1);
}

//-----------------------------------------------------------------------------

Spectra_Stream_Saver::~Spectra_Stream_Saver()
{

}

// ----------------------------------------------------------------------------

void Spectra_Stream_Saver::save_stream(data_struct::Stream_Block* stream_block)
{

    size_t d_hash = stream_block->dataset_hash();
    size_t detector_num = stream_block->detector_number;


    // Is this a new dataset
    if(_dataset_map.count(d_hash) < 1)
    {
        // Close any open datasets because we should not get any more data from them
        for (auto itr : _dataset_map)
        {
            _finalize_dataset(itr.second);
        }
        _dataset_map.clear();

        //insert new dataset
        _new_dataset(d_hash, stream_block);
    }
    else
    {
        // Get dataset and check if we have detector for it
        Dataset_Save *dataset = _dataset_map.at(d_hash);
        if(dataset->detector_map.count(detector_num) < 1)
        {
           _new_detector(dataset, stream_block);
        }
        else
        {
            Detector_Save *detector = dataset->detector_map.at(detector_num);

            if(detector->last_row > -1 && stream_block->row() > (size_t)detector->last_row)
            {
                io::file::HDF5_IO::inst()->save_stream_row(d_hash, detector_num, detector->last_row, &detector->spectra_line);
                for(size_t i=0; i<detector->spectra_line.size(); i++)
                {
                    if(detector->spectra_line[i] != nullptr)
                    {
                        delete detector->spectra_line[i];
                        detector->spectra_line[i] = nullptr;
                    }
                }
            }

            detector->last_row = stream_block->row();
            detector->integrated_spectra.add(*stream_block->spectra);
            //TODO: add limit checks to spectra_line
            if (detector->spectra_line[stream_block->col()]  != nullptr)
            {
                delete detector->spectra_line[stream_block->col()];
            }
            detector->spectra_line[stream_block->col()] = stream_block->spectra;
            stream_block->spectra = nullptr;
        }

        //TODO: this won't work because netcdf lines may be corrupt and last row/col won't match height/width
        if(stream_block->is_end_of_detector())
        {
            //_finalize_dataset();
            Detector_Save *detector = dataset->detector_map.at(detector_num);

            //save and close hdf5 for this detector
            ///io::file::HDF5_IO::inst()->save_scan_scalers(detector_num, stream_block->mda_io, params_override, false);
            io::file::HDF5_IO::inst()->save_itegrade_spectra(&detector->integrated_spectra);
            io::file::HDF5_IO::inst()->close_dataset(d_hash);
            ///delete stream_block->mda_io;

            delete detector;
            dataset->detector_map.erase(detector_num);
        }

        if(dataset->detector_map.size() == 0)
        {
            //done with dataset
            delete dataset;
            _dataset_map.erase(d_hash);
        }

    }

}

// ----------------------------------------------------------------------------

void Spectra_Stream_Saver::_new_dataset(size_t d_hash, data_struct::Stream_Block* stream_block)
{
    Dataset_Save *dataset = new Dataset_Save();
    dataset->dataset_directory = stream_block->dataset_directory;
    dataset->dataset_name = stream_block->dataset_name;
    _dataset_map.insert( {d_hash, dataset} );
    _new_detector(dataset, stream_block);
}

// ----------------------------------------------------------------------------

void Spectra_Stream_Saver::_new_detector(Dataset_Save *dataset, data_struct::Stream_Block* stream_block)
{
    Detector_Save *detector = new Detector_Save(stream_block->width());
    dataset->detector_map.insert( { stream_block->detector_number, detector } );

    detector->integrated_spectra = *stream_block->spectra;
    detector->spectra_line[stream_block->col()] = stream_block->spectra;

    //release ownership
    stream_block->spectra = nullptr;

    io::file::HDF5_IO::inst()->generate_stream_dataset(*dataset->dataset_directory, *dataset->dataset_name, stream_block->detector_number, stream_block->height(), stream_block->width());
}

// ----------------------------------------------------------------------------

void Spectra_Stream_Saver::_finalize_dataset(Dataset_Save *dataset)
{
    for(auto itr : dataset->detector_map)
    {
        Detector_Save *detector = itr.second;

        //save and close hdf5 for this detector
        ///io::file::HDF5_IO::inst()->save_scan_scalers(detector_num, stream_block->mda_io, params_override, false);
        io::file::HDF5_IO::inst()->save_itegrade_spectra(&detector->integrated_spectra);
//        io::file::HDF5_IO::inst()->close_dataset(d_hash);
        ///delete stream_block->mda_io;

        delete detector;
    }

    dataset->detector_map.clear();

    delete dataset;
}

// ----------------------------------------------------------------------------

} //namespace xrf
} //namespace workflow
