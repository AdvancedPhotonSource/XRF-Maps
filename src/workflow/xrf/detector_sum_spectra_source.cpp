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



#include "detector_sum_spectra_source.h"

namespace workflow
{
namespace xrf
{

//-----------------------------------------------------------------------------

template<typename T_real>
Detector_Sum_Spectra_Source<T_real>::Detector_Sum_Spectra_Source() : Spectra_File_Source(nullptr)
{
    _cb_function = std::bind(&Detector_Sum_Spectra_Source<T_real>::cb_load_spectra_data, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7);
    _spectra = new data_struct::Spectra<T_real>(2048, 0.0, 0.0, 0.0, 0.0);
	for (size_t det = 0; det < 4; det++)
	{
		_detector_num_arr.push_back(det);
	}

}

//-----------------------------------------------------------------------------

template<typename T_real>
Detector_Sum_Spectra_Source<T_real>::Detector_Sum_Spectra_Source(data_struct::Analysis_Job<T_real>* analysis_job) : Spectra_File_Source(analysis_job)
{
    _cb_function = std::bind(&Detector_Sum_Spectra_Source<T_real>::cb_load_spectra_data, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7);
    _spectra = new data_struct::Spectra<T_real>(2048, 0.0, 0.0, 0.0, 0.0);
    if(analysis_job != nullptr)
    {
		for (size_t det : analysis_job->detector_num_arr)
		{
			_detector_num_arr.push_back(det);
		}
	}
    else
    {
		for (size_t det = 0; det < 4; det++)
		{
			_detector_num_arr.push_back(det);
		}
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
Detector_Sum_Spectra_Source<T_real>::~Detector_Sum_Spectra_Source()
{

}

// ----------------------------------------------------------------------------

template<typename T_real>
void Detector_Sum_Spectra_Source<T_real>::cb_load_spectra_data(size_t row, size_t col, size_t height, size_t width, size_t detector_num, data_struct::Spectra<T_real>* spectra, void* user_data)
{

    if(_spectra->size() < spectra->size())
    {
        _spectra->resize( spectra->size() );
    }

    _spectra->add(*spectra);

    if(detector_num == _detector_num_arr[_detector_num_arr.size()-1] && _output_callback_func != nullptr)
    {
        data_struct::Stream_Block<T_real>* stream_block = new data_struct::Stream_Block<T_real>(-1, row, col, height, width);

        if(_analysis_job != nullptr)
        {
            if(_init_fitting_routines)
            {
                _analysis_job->init_fit_routines(spectra->size());
            }
            data_struct::Detector<T_real>* cp = _analysis_job->get_detector(detector_num);

            if(cp == nullptr)
            {
                cp = _analysis_job->get_first_detector();
            }
            if(cp!=nullptr)
            {
                stream_block->init_fitting_blocks(&(cp->fit_routines), &(cp->fit_params_override_dict.elements_to_fit));
                stream_block->model = cp->model;
            }

            stream_block->theta = _analysis_job->theta;
            stream_block->optimize_fit_params_preset = _analysis_job->optimize_fit_params_preset;
        }

        stream_block->spectra = _spectra;
        stream_block->dataset_directory = _current_dataset_directory;
        stream_block->dataset_name = _current_dataset_name;
        
        _output_callback_func(stream_block);

        _spectra = new data_struct::Spectra<T_real>(spectra->size(), 0.0, 0.0, 0.0, 0.0);
    }

    delete spectra;

}


//-----------------------------------------------------------------------------

} //namespace xrf
} //namespace workflow
