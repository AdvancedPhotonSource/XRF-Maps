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



#ifndef Spectra_Stream_Saver_H
#define Spectra_Stream_Saver_H

#include "core/defines.h"

#include "workflow/sink.h"
#include "data_struct/stream_block.h"
#include "io/file/mda_io.h"
#include "io/file/hdf5_io.h"
#include <functional>

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

class DLL_EXPORT Spectra_Stream_Saver : public Sink<data_struct::Stream_Block* >
{

public:

    Spectra_Stream_Saver();

    virtual ~Spectra_Stream_Saver();

    void save_stream(data_struct::Stream_Block* stream_block);

    virtual void set_function(std::function<void (data_struct::Stream_Block*)> func) { }

protected:

    class Detector_Save
    {
    public:
        Detector_Save(size_t width)
        {
            last_row = -1;
            spectra_line.resize(width);
            for(size_t i=0;i<width; i++)
            {
                spectra_line[i] = nullptr;
            }
        }
        ~Detector_Save()
        {
            spectra_line.resize(1);
        }

        int last_row;
        data_struct::Spectra integrated_spectra;
        std::vector< data_struct::Spectra* > spectra_line;
    };

    class Dataset_Save
    {
    public:
        Dataset_Save(){}
        ~Dataset_Save()
        {
            delete dataset_directory;
            delete dataset_name;
            for(auto& itr : detector_map)
            {
                if (itr.second != nullptr)
                {
                    delete itr.second;
                }
            }
            detector_map.clear();
        }

        std::string *dataset_directory;
        std::string *dataset_name;
        //by detector_num
        std::map<int, Detector_Save*> detector_map;
    };

    void _new_dataset(size_t d_hash, data_struct::Stream_Block* stream_block);

    void _new_detector(Dataset_Save *dataset, data_struct::Stream_Block* stream_block);

    void _finalize_dataset(Dataset_Save *dataset);

    //by detector_dir + dataset hash
    std::map<size_t, Dataset_Save*> _dataset_map;

};

//-----------------------------------------------------------------------------

} //namespace xrf
} //namespace workflow

#endif // Spectra_Stream_Saver_H
