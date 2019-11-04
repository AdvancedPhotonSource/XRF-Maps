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



#ifndef Spectra_File_Source_H
#define Spectra_File_Source_H

#include "core/defines.h"

#include "workflow/source.h"
#include "data_struct/stream_block.h"
#include "data_struct/analysis_job.h"
#include "io/file/netcdf_io.h"
#include "io/file/mda_io.h"
#include "io/file/hdf5_io.h"
#include <functional>
#include <iostream>
#include <fstream>

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

class DLL_EXPORT Spectra_File_Source : public Source<data_struct::Stream_Block*>
{

public:

    Spectra_File_Source();

    //used with run function to process job
    Spectra_File_Source(data_struct::Analysis_Job* analysis_job);

    virtual ~Spectra_File_Source();

    virtual void cb_load_spectra_data(size_t row, size_t col, size_t height, size_t width, size_t detector_num, data_struct::Spectra* spectra, void* user_data);

    virtual void run();

    bool load_netcdf_line(std::string dirpath,
						  std::string filename,
							const std::vector<size_t>& detector_num_arr,
                          size_t row,
                          size_t row_size,
                          size_t col_size);

    void set_init_fitting_routines(bool val) {_init_fitting_routines = val;}

protected:

    virtual bool _load_spectra_volume_with_callback(std::string dataset_directory,
                                                    std::string dataset_file,
													const std::vector<size_t>& detector_num_arr,
													data_struct::IO_Callback_Func_Def callback_fun);


	data_struct::Stream_Block* _alloc_stream_block(size_t row, size_t col, size_t height, size_t width, size_t spectra_size);

	int _max_num_stream_blocks;
	int _allocated_stream_blocks;

    std::string *_current_dataset_directory;
    std::string *_current_dataset_name;

    data_struct::Analysis_Job* _analysis_job;

    std::vector<std::string> _netcdf_files;

    std::vector<std::string> _bnp_netcdf_files;

    std::vector<std::string> _hdf_files;

    std::function <void (size_t, size_t, size_t, size_t, size_t, data_struct::Spectra*, void*)> _cb_function;

    bool _init_fitting_routines;

};

} //namespace xrf
} //namespace workflow

#endif // Spectra_File_Source_H
