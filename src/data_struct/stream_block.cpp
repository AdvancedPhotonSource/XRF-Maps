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



#include "stream_block.h"

namespace data_struct
{

template<typename T_real>
Stream_Block<T_real>::Stream_Block()
{
	_row = 0;
	_col = 0;
    _height = 0;
    _width = 0;
    theta = 0;
	elements_to_fit = nullptr;
    // by default we don't want to delete the string pointers becaues they are shared by stream blocks
    del_str_ptr = false;
	spectra = nullptr;
    optimize_fit_params_preset = fitting::models::Fit_Params_Preset::BATCH_FIT_NO_TAILS;
}

//-----------------------------------------------------------------------------

template<typename T_real>
Stream_Block<T_real>::Stream_Block(int detector,
                           size_t row,
                           size_t col,
                           size_t height,
                           size_t width)
{
    _row = row;
    _col = col;
    _height = height;
    _width = width;
    _detector = detector;
    theta = 0;
	elements_to_fit = nullptr;
    // by default we don't want to delete the string pointers becaues they are shared by stream blocks
    del_str_ptr = false;
    spectra = nullptr;
}

//-----------------------------------------------------------------------------

template<typename T_real>
Stream_Block<T_real>::~Stream_Block()
{
    if(del_str_ptr)
    {
        if(dataset_name != nullptr)
        {
            delete dataset_name;
            dataset_name = nullptr;
        }
        if(dataset_directory != nullptr)
        {
            delete dataset_directory;
            dataset_directory = nullptr;
        }
    }

    if(spectra != nullptr)
    {
        delete spectra;
        spectra = nullptr;
    }

    elements_to_fit = nullptr;

    model = nullptr;
}

//-----------------------------------------------------------------------------

template<typename T_real>
Stream_Block<T_real>::Stream_Block(const Stream_Block& stream_block)
{
	this->_col = stream_block._col;
	this->_row = stream_block._row;
	this->_height = stream_block._height;
	this->_width = stream_block._width;
	this->dataset_directory = stream_block.dataset_directory;
	this->dataset_name = stream_block.dataset_name;
	this->fitting_blocks = stream_block.fitting_blocks;
	this->_detector = stream_block._detector;
	this->spectra = stream_block.spectra;
	this->elements_to_fit = stream_block.elements_to_fit;
	this->model = stream_block.model;
}

//-----------------------------------------------------------------------------

template<typename T_real>
Stream_Block<T_real> &Stream_Block<T_real>::operator=(const Stream_Block<T_real>& stream_block)
{
	this->_col = stream_block._col;
	this->_row = stream_block._row;
	this->_height = stream_block._height;
	this->_width = stream_block._width;
	this->dataset_directory = stream_block.dataset_directory;
	this->dataset_name = stream_block.dataset_name;
	this->fitting_blocks = stream_block.fitting_blocks;
	this->_detector = stream_block._detector;
	this->spectra = stream_block.spectra;
	this->elements_to_fit = stream_block.elements_to_fit;
	this->model = stream_block.model;
	return *this;
}

//-----------------------------------------------------------------------------

template<typename T_real>
void Stream_Block<T_real>::init_fitting_blocks(std::unordered_map<Fitting_Routines, fitting::routines::Base_Fit_Routine<T_real> *> *fit_routines,
                                       Fit_Element_Map_Dict<T_real> * elements_to_fit_)
{
    elements_to_fit = elements_to_fit_;

    if(elements_to_fit == nullptr)
    {
        //throw Exception;
    }

    for(const auto &itr : *fit_routines)
    {
        fitting_blocks[itr.first] = Stream_Fitting_Block();
        fitting_blocks[itr.first].fit_routine = itr.second;
        for(auto& e_itr : *elements_to_fit)
        {
            fitting_blocks[itr.first].fit_counts.emplace(std::pair<std::string, T_real> (e_itr.first, (T_real)0.0));
        }
        fitting_blocks[itr.first].fit_counts.emplace(std::pair<std::string, T_real> (STR_NUM_ITR, (T_real)0.0));
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
size_t Stream_Block<T_real>::dataset_hash()
{
    if (dataset_directory != nullptr && dataset_name != nullptr)
    {
        return std::hash<std::string> {} ((*dataset_directory) + (*dataset_name)) + _detector;
    }
    return -1;
}

//-----------------------------------------------------------------------------

} //namespace data_struct
