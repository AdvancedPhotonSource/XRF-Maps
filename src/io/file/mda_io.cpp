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


#include "mda_io.h"

#include <string>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


#ifdef XDR_HACK
    #include "xdr_hack.h"
#else
  #include <rpc/types.h>
  #include <rpc/xdr.h>
#endif


namespace io
{
namespace file
{

//-----------------------------------------------------------------------------

MDA_IO::MDA_IO() : Base_File_IO()
{

    _mda_file = nullptr;
    _mda_file_info = nullptr;

}

//-----------------------------------------------------------------------------

MDA_IO::~MDA_IO()
{

    unload();

}

//-----------------------------------------------------------------------------

void MDA_IO::unload()
{
    if(_mda_file != nullptr)
    {
        mda_unload(_mda_file);
        _mda_file = nullptr;
    }
    if(_mda_file_info != nullptr)
    {
        mda_info_unload(_mda_file_info);
        _mda_file_info = nullptr;
    }
}

//-----------------------------------------------------------------------------

void MDA_IO::lazy_load()
{

    std::FILE *fptr = std::fopen(_filename.c_str(), "r");

    _mda_file_info = mda_info_load(fptr);

    std::fclose(fptr);

	if (_mda_file_info == nullptr)
	{
        logit << "Error loading mda file:" << _filename<<std::endl;
		return;
	}

    logit<<"mda info ver:"<<_mda_file_info->version<<" data rank:"<<_mda_file_info->data_rank;
    for(int16_t i = 0; i < _mda_file_info->data_rank; i++)
    {
        logit_s<<" dims["<<i<<"]:"<<_mda_file_info->dimensions[i];
    }
    logit_s<<std::endl;

}

//-----------------------------------------------------------------------------

bool MDA_IO::load_dataset(std::string path, Base_Dataset *dset)
{

    std::FILE *fptr = std::fopen(path.c_str(), "r");

    if (fptr == nullptr)
    {
        return false;
    }
    _mda_file = mda_load(fptr);
    std::fclose(fptr);
    logit<<"mda info ver:"<<_mda_file->header->version<<" data rank:"<<_mda_file->header->data_rank;
    //long total = 1;

    for(int16_t i = 0; i < _mda_file->header->data_rank; i++)
    {
        logit_s<<" dims["<<i<<"]:"<<_mda_file->header->dimensions[i];
        //total *= _mda_file->header->dimensions[i];
    }

    //dset->alloc(_mda_file->header->data_rank, _mda_file->header->dimensions);
    for(unsigned int j = 0; j < _mda_file->header->data_rank; j++)
    {

    }
    logit_s<<std::endl;
    logit<<"d "<<(_mda_file->scan->sub_scans[0]->sub_scans[0]->detectors_data[0][0])<<std::endl;
    logit<<"d "<<(_mda_file->scan->sub_scans[0]->sub_scans[0]->detectors_data[0][1])<<std::endl;


    return true;

}

//-----------------------------------------------------------------------------

int MDA_IO::find_2d_detector_index(struct mda_file* mda_file, std::string det_name, int detector_num, real_t& val)
{

    for(int k=0; k<mda_file->scan->number_detectors; k++)
    {
        //logit<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->name << std::endl;
        //logit<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->description << std::endl;
        if(strcmp(mda_file->scan->detectors[k]->name, det_name.c_str())  == 0)
        {
            val = mda_file->scan->detectors_data[k][detector_num];
            return k;
        }
    }

    for(int k=0; k<mda_file->scan->sub_scans[0]->number_detectors; k++)
    {
        //logit<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->name << std::endl;
        //logit<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->description << std::endl;
        if(strcmp(mda_file->scan->sub_scans[0]->detectors[k]->name, det_name.c_str())  == 0)
        {
            val = mda_file->scan->sub_scans[0]->detectors_data[k][detector_num];
            return k;
        }
    }
    return -1;

}

//-----------------------------------------------------------------------------
/*
bool MDA_IO::load_spectra_volume(std::string path,
                                 size_t detector_num,
                                 data_struct::xrf::Spectra_Volume* vol,
                                 bool hasNetCDF,
                                 data_struct::xrf::Params_Override *override_values,
                                 data_struct::xrf::Quantification_Standard * quantification_standard)
{
    //index per row and col
    int elt_idx = -1;
    int ert_idx = -1;
    int incnt_idx = -1;
    int outcnt_idx = -1;

    bool single_row_scan = false;

    std::FILE *fptr = std::fopen(path.c_str(), "rb");

    size_t cols = 1;
    size_t rows = 1;
    size_t samples = 1;

    if (fptr == nullptr)
    {
        return false;
    }


    _mda_file = mda_load(fptr);
    std::fclose(fptr);
    if (_mda_file == nullptr)
    {
        return false;
    }
    logit<<"mda info ver:"<<_mda_file->header->version<<" data rank:"<<_mda_file->header->data_rank;

    if (_mda_file->header->data_rank == 2)
    {
        logit_s<<" requested cols "<< _mda_file->header->dimensions[0] << " requested rows " << _mda_file->header->dimensions[1] <<
                  " acquired cols "<< _mda_file->scan->last_point << " acquired rows " << _mda_file->scan->sub_scans[0]->last_point <<std::endl;

        if(hasNetCDF)
        {
            if(_mda_file->scan->last_point == 0)
                rows = 1;
            else
                rows = _mda_file->scan->last_point;
            if(_mda_file->scan->sub_scans[0]->last_point == 0)
                cols = 1;
            else
                cols = _mda_file->scan->sub_scans[0]->last_point;
            vol->resize(rows, cols, 2048);
            return true;
        }
        else
        {
            if(_mda_file->header->dimensions[1] == 2000)
            {
                if(_mda_file->scan->sub_scans[0]->number_detectors-1 < detector_num)
                {
                    logit<<"Error: max detectors saved = "<<_mda_file->scan->sub_scans[0]->number_detectors<< std::endl;
                    unload();
                    return false;
                }

                rows = 1;
                if(_mda_file->scan->last_point == 0)
                    cols = 1;
                else
                cols = _mda_file->scan->last_point;
                samples = _mda_file->header->dimensions[1];
                vol->resize(rows, cols, samples);
                single_row_scan = true;
            }
            else
            {
                //if not then we don't know what is dataset is.
                unload();
                return false;
            }
        }
    }
    else if (_mda_file->header->data_rank == 3)
    {

        if(_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors-1 < detector_num)
        {
            logit<<"Error: max detectors saved = "<<_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors<< std::endl;
            unload();
            return false;
        }

        if(_mda_file->scan->last_point == 0)
            rows = 1;
        else
            rows = _mda_file->scan->last_point;
        if(_mda_file->scan->sub_scans[0]->last_point == 0)
            cols = 1;
        else
            cols = _mda_file->scan->sub_scans[0]->last_point;
        samples = _mda_file->header->dimensions[2];
        vol->resize(rows, cols, samples);
    }
    else
    {
        logit<<" Error: no support for data rank "<< _mda_file->header->data_rank <<std::endl;
        unload();
        return false;
    }

    //_load_detector_meta_data(detector);

    if (override_values != nullptr)
    {
        real_t tmp_val;

        if (override_values->elt_pv.length() > 0)
        {
            elt_idx = find_2d_detector_index(_mda_file, override_values->elt_pv, detector_num, tmp_val );
        }
        if (override_values->ert_pv.length() > 0)
        {
            ert_idx = find_2d_detector_index(_mda_file, override_values->ert_pv, detector_num, tmp_val );
        }
        if (override_values->in_cnt_pv.length() > 0)
        {
            incnt_idx = find_2d_detector_index(_mda_file, override_values->in_cnt_pv, detector_num, tmp_val );
        }
        if (override_values->out_cnt_pv.length() > 0)
        {
            outcnt_idx = find_2d_detector_index(_mda_file, override_values->out_cnt_pv, detector_num, tmp_val );
        }
        if(quantification_standard != nullptr)
        {
            if (override_values->scaler_pvs.count("SRCURRENT") > 0)
            {
                find_2d_detector_index(_mda_file, override_values->scaler_pvs.at("SRCURRENT"), detector_num, tmp_val );
                quantification_standard->sr_current(tmp_val);
            }
            if (override_values->scaler_pvs.count("US_IC") > 0)
            {
                find_2d_detector_index(_mda_file, override_values->scaler_pvs.at("US_IC"), detector_num, tmp_val );
                quantification_standard->US_IC(tmp_val);
            }
            if (override_values->scaler_pvs.count("DS_IC") > 0)
            {
                find_2d_detector_index(_mda_file, override_values->scaler_pvs.at("DS_IC"), detector_num, tmp_val );
                quantification_standard->DS_IC(tmp_val);
            }
        }

    }
    logit<<" elt_idx "<< elt_idx << " ert_idx " << ert_idx << " in cnt idx " << incnt_idx << " out cnt idx "<< outcnt_idx<<std::endl;


    try
    {
        if( single_row_scan )
        {
            if(_mda_file->scan->last_point < _mda_file->scan->requested_points)
            {
                cols = _mda_file->scan->last_point;
            }
        }
        else
        {
            if(_mda_file->scan->last_point < _mda_file->scan->requested_points)
            {
                rows = _mda_file->scan->last_point;
                //TODO: set a flag to return to tell that this is a bad scan
            }
        }

        for(size_t i=0; i<rows; i++)
        {
            // update num rows if header is incorrect and not single row scan

            if(false == single_row_scan)
            {
                if(_mda_file->scan->sub_scans[i]->last_point < _mda_file->scan->sub_scans[i]->requested_points)
                {
                    cols = _mda_file->scan->sub_scans[i]->last_point;
                    //TODO: set a flag to return to tell that this is a bad scan
                }
            }
            for(size_t j=0; j<cols; j++)
            {
// TODO: we might need to do the same check for samples size
//                if(_mda_file->scan->sub_scans[i]->sub_scan[j]->last_point < _mda_file->scan->sub_scans[i]->sub_scan[j]->requested_points)
//                {
//                    samples = _mda_file->scan->sub_scans[i]->sub_scan[j]->last_point;
//                }
//


                if (single_row_scan)
                {
                    if(elt_idx > -1)
                    {
                        //logit<<"eltm ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[elt_idx][j]<< std::endl;
                        (*vol)[i][j].elapsed_lifetime(_mda_file->scan->detectors_data[elt_idx][j]);
                    }
                    if(ert_idx > -1)
                    {
                        //logit<<"elrm ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[ert_idx][j]<< std::endl;
                        (*vol)[i][j].elapsed_realtime(_mda_file->scan->detectors_data[ert_idx][j]);
                    }
                    if(incnt_idx > -1)
                    {
                        //logit<<"incnt ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[incnt_idx][j]<< std::endl;
                        (*vol)[i][j].input_counts(_mda_file->scan->detectors_data[incnt_idx][j]);
                    }
                    if(outcnt_idx > -1)
                    {
                        //logit<<"outcnt ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[outcnt_idx][j]<< std::endl;
                        (*vol)[i][j].output_counts(_mda_file->scan->detectors_data[outcnt_idx][j]);
                    }
                    if(ert_idx > -1 && incnt_idx > -1 && outcnt_idx > -1)
                    {
                        (*vol)[i][j].recalc_elapsed_lifetime();
                    }


                    for(size_t k=0; k<samples; k++)
                    {

                        (*vol)[i][j][k] = (_mda_file->scan->sub_scans[j]->detectors_data[detector_num][k]);
                    }
                }
                else
                {
                    if(elt_idx > -1)
                    {
                        //logit<<"eltm ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[elt_idx][j]<< std::endl;
                        (*vol)[i][j].elapsed_lifetime(_mda_file->scan->sub_scans[i]->detectors_data[elt_idx][j]);
                    }
                    if(ert_idx > -1)
                    {
                        //logit<<"elrm ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[ert_idx][j]<< std::endl;
                        (*vol)[i][j].elapsed_realtime(_mda_file->scan->sub_scans[i]->detectors_data[ert_idx][j]);
                    }
                    if(incnt_idx > -1)
                    {
                        //logit<<"incnt ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[incnt_idx][j]<< std::endl;
                        (*vol)[i][j].input_counts(_mda_file->scan->sub_scans[i]->detectors_data[incnt_idx][j]);
                    }
                    if(outcnt_idx > -1)
                    {
                        //logit<<"outcnt ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[outcnt_idx][j]<< std::endl;
                        (*vol)[i][j].output_counts(_mda_file->scan->sub_scans[i]->detectors_data[outcnt_idx][j]);
                    }
                    if(ert_idx > -1 && incnt_idx > -1 && outcnt_idx > -1)
                    {
                        (*vol)[i][j].recalc_elapsed_lifetime();
                    }


                    for(size_t k=0; k<samples; k++)
                    {
                        (*vol)[i][j][k] = (_mda_file->scan->sub_scans[i]->sub_scans[j]->detectors_data[detector_num][k]);
                    }
                }
            }
        }
    }
    catch(std::exception& e)
    {
        logit<<"!!! Error Caught exception loading mda file."<<std::endl;
        std::cerr << "Exception catched : " << e.what() << std::endl;
        return false;
    }

    return true;
}
*/

void MDA_IO::find_scaler_indexes(size_t detector_num_start,
                                 size_t detector_num_end,
                                 data_struct::xrf::Analysis_Job *analysis_job)
{

    _spectra_scalers_map.clear();

//TODO: See if mda is 0 based or 1 for detector nums
    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        //_spectra_scalers_map.insert( std::pair<size_t, struct spectra_scalers_indexes> (detector_num, struct spectra_scalers_indexes()) );

        struct data_struct::xrf::Analysis_Sub_Struct* detector_struct = analysis_job->get_sub_struct(detector_num);

        real_t tmp_val;

        if (detector_struct->fit_params_override_dict.elt_pv.length() > 0)
        {
            _spectra_scalers_map[detector_num].elt_idx = find_2d_detector_index(_mda_file, detector_struct->fit_params_override_dict.elt_pv, detector_num, tmp_val );
        }
        if (detector_struct->fit_params_override_dict.ert_pv.length() > 0)
        {
            _spectra_scalers_map[detector_num].ert_idx = find_2d_detector_index(_mda_file, detector_struct->fit_params_override_dict.ert_pv, detector_num, tmp_val );
        }
        if (detector_struct->fit_params_override_dict.in_cnt_pv.length() > 0)
        {
            _spectra_scalers_map[detector_num].incnt_idx = find_2d_detector_index(_mda_file, detector_struct->fit_params_override_dict.in_cnt_pv, detector_num, tmp_val );
        }
        if (detector_struct->fit_params_override_dict.out_cnt_pv.length() > 0)
        {
            _spectra_scalers_map[detector_num].outcnt_idx = find_2d_detector_index(_mda_file, detector_struct->fit_params_override_dict.out_cnt_pv, detector_num, tmp_val );
        }


        if (detector_struct->fit_params_override_dict.scaler_pvs.count("SRCURRENT") > 0)
        {
            find_2d_detector_index(_mda_file, detector_struct->fit_params_override_dict.scaler_pvs.at("SRCURRENT"), detector_num, tmp_val );
            detector_struct->quant_standard.sr_current(tmp_val);
        }
        if (detector_struct->fit_params_override_dict.scaler_pvs.count("US_IC") > 0)
        {
            find_2d_detector_index(_mda_file, detector_struct->fit_params_override_dict.scaler_pvs.at("US_IC"), detector_num, tmp_val );
            detector_struct->quant_standard.US_IC(tmp_val);
        }
        if (detector_struct->fit_params_override_dict.scaler_pvs.count("DS_IC") > 0)
        {
            find_2d_detector_index(_mda_file, detector_struct->fit_params_override_dict.scaler_pvs.at("DS_IC"), detector_num, tmp_val );
            detector_struct->quant_standard.DS_IC(tmp_val);
        }
    }
}

//-----------------------------------------------------------------------------

bool MDA_IO::load_spectra_volume_with_callback(std::string path,
                                                 size_t detector_num_start,
                                                 size_t detector_num_end,
                                                 bool hasNetCDF,
                                                 data_struct::xrf::Analysis_Job *analysis_job,
                                                 IO_Callback_Func_Def callback_func,
                                                 void *user_data)
{
    //index per row and col
    int elt_idx = -1;
    int ert_idx = -1;
    int incnt_idx = -1;
    int outcnt_idx = -1;

    bool single_row_scan = false;

    std::FILE *fptr = std::fopen(path.c_str(), "rb");

    size_t samples = 1;

    if (fptr == nullptr)
    {
        return false;
    }


    _mda_file = mda_load(fptr);
    std::fclose(fptr);
    if (_mda_file == nullptr)
    {
        return false;
    }
    logit<<"mda info ver:"<<_mda_file->header->version<<" data rank:"<<_mda_file->header->data_rank;

    if (_mda_file->header->data_rank == 2)
    {
        logit_s<<" requested cols "<< _mda_file->header->dimensions[0] << " requested rows " << _mda_file->header->dimensions[1] <<
                  " acquired cols "<< _mda_file->scan->last_point << " acquired rows " << _mda_file->scan->sub_scans[0]->last_point <<std::endl;

        if(hasNetCDF)
        {
            if(_mda_file->scan->last_point == 0)
                _rows = 1;
            else
                _rows = _mda_file->scan->last_point;
            if(_mda_file->scan->sub_scans[0]->last_point == 0)
                _cols = 1;
            else
                _cols = _mda_file->scan->sub_scans[0]->last_point;
            return true;
        }
        else
        {
            if(_mda_file->header->dimensions[1] == 2000)
            {
                if(_mda_file->scan->sub_scans[0]->number_detectors-1 < detector_num_start || _mda_file->scan->sub_scans[0]->number_detectors-1 > detector_num_end)
                {
                    logit<<"Error: max detectors saved = "<<_mda_file->scan->sub_scans[0]->number_detectors<< std::endl;
                    unload();
                    return false;
                }

                _rows = 1;
                if(_mda_file->scan->last_point == 0)
                    _cols = 1;
                else
                _cols = _mda_file->scan->last_point;
                samples = _mda_file->header->dimensions[1];
                single_row_scan = true;
            }
            else
            {
                //if not then we don't know what is dataset is.
                unload();
                return false;
            }
        }
    }
    else if (_mda_file->header->data_rank == 3)
    {

        if(_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors-1 < detector_num_start || _mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors-1 > detector_num_end)
        {
            logit<<"Error: max detectors saved = "<<_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors<< std::endl;
            unload();
            return false;
        }

        if(_mda_file->scan->last_point == 0)
            _rows = 1;
        else
            _rows = _mda_file->scan->last_point;
        if(_mda_file->scan->sub_scans[0]->last_point == 0)
            _cols = 1;
        else
            _cols = _mda_file->scan->sub_scans[0]->last_point;
        samples = _mda_file->header->dimensions[2];
    }
    else
    {
        logit<<" Error: no support for data rank "<< _mda_file->header->data_rank <<std::endl;
        unload();
        return false;
    }

    //_load_detector_meta_data(detector);

    if (analysis_job != nullptr)
    {
        find_scaler_indexes(detector_num_start, detector_num_end, analysis_job);
    }
    logit<<" elt_idx "<< elt_idx << " ert_idx " << ert_idx << " in cnt idx " << incnt_idx << " out cnt idx "<< outcnt_idx<<std::endl;


    try
    {
        if( single_row_scan )
        {
            if(_mda_file->scan->last_point < _mda_file->scan->requested_points)
            {
                _cols = _mda_file->scan->last_point;
            }
        }
        else
        {
            if(_mda_file->scan->last_point < _mda_file->scan->requested_points)
            {
                _rows = _mda_file->scan->last_point;
                //TODO: set a flag to return to tell that this is a bad scan
            }
        }

        for(size_t i=0; i<_rows; i++)
        {
            // update num rows if header is incorrect and not single row scan

            if(false == single_row_scan)
            {
                if(_mda_file->scan->sub_scans[i]->last_point < _mda_file->scan->sub_scans[i]->requested_points)
                {
                    _cols = _mda_file->scan->sub_scans[i]->last_point;
                    //TODO: set a flag to return to tell that this is a bad scan
                }
            }
            for(size_t j=0; j<_cols; j++)
            {
/* TODO: we might need to do the same check for samples size
                if(_mda_file->scan->sub_scans[i]->sub_scan[j]->last_point < _mda_file->scan->sub_scans[i]->sub_scan[j]->requested_points)
                {
                    samples = _mda_file->scan->sub_scans[i]->sub_scan[j]->last_point;
                }
*/

                for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
                {
                    data_struct::xrf::Spectra* spectra = new data_struct::xrf::Spectra(samples);
                    if (single_row_scan)
                    {
                        if(elt_idx > -1)
                        {
                            spectra->elapsed_lifetime(_mda_file->scan->detectors_data[elt_idx][j]);
                        }
                        if(ert_idx > -1)
                        {
                            spectra->elapsed_realtime(_mda_file->scan->detectors_data[ert_idx][j]);
                        }
                        if(incnt_idx > -1)
                        {
                            spectra->input_counts(_mda_file->scan->detectors_data[incnt_idx][j]);
                        }
                        if(outcnt_idx > -1)
                        {
                            spectra->output_counts(_mda_file->scan->detectors_data[outcnt_idx][j]);
                        }
                        if(ert_idx > -1 && incnt_idx > -1 && outcnt_idx > -1)
                        {
                            spectra->recalc_elapsed_lifetime();
                        }


                        for(size_t k=0; k<samples; k++)
                        {

                            (*spectra)[k] = (_mda_file->scan->sub_scans[j]->detectors_data[detector_num][k]);
                        }
                        callback_func(i, j, _rows, _cols, detector_num, spectra, user_data);
                    }
                    else
                    {
                        if(elt_idx > -1)
                        {
                             spectra->elapsed_lifetime(_mda_file->scan->sub_scans[i]->detectors_data[elt_idx][j]);
                        }
                        if(ert_idx > -1)
                        {
                            spectra->elapsed_realtime(_mda_file->scan->sub_scans[i]->detectors_data[ert_idx][j]);
                        }
                        if(incnt_idx > -1)
                        {
                            spectra->input_counts(_mda_file->scan->sub_scans[i]->detectors_data[incnt_idx][j]);
                        }
                        if(outcnt_idx > -1)
                        {
                            spectra->output_counts(_mda_file->scan->sub_scans[i]->detectors_data[outcnt_idx][j]);
                        }
                        if(ert_idx > -1 && incnt_idx > -1 && outcnt_idx > -1)
                        {
                            spectra->recalc_elapsed_lifetime();
                        }


                        for(size_t k=0; k<samples; k++)
                        {
                            (*spectra)[k] = (_mda_file->scan->sub_scans[i]->sub_scans[j]->detectors_data[detector_num][k]);
                        }
                        callback_func(i, j, _rows, _cols, detector_num, spectra, user_data);
                    }
                }
            }
        }
    }
    catch(std::exception& e)
    {
        logit<<"!!! Error Caught exception loading mda file."<<std::endl;
        logit << "Exception catched : " << e.what() << std::endl;
        return false;
    }

    return true;
}

//-----------------------------------------------------------------------------

int MDA_IO::get_multiplied_dims(std::string path)
{
    int f_size = -1;


    std::FILE *fptr = std::fopen(path.c_str(), "rb");
    struct mda_header *header = mda_header_load(fptr);

    std::fclose(fptr);

    if (header == nullptr)
    {
        logit<<"Unable to open mda file "<< path <<std::endl;
        return f_size;
    }
    else if(header->data_rank == 1)
    {
        f_size = header->dimensions[0];
    }
    else if(header->data_rank == 2 || header->data_rank == 3)
    {
        f_size = header->dimensions[0] * header->dimensions[1];
    }
    else
    {
        logit<<"Unsupported mda data rank "<<header->data_rank<<" . Skipping file "<< path <<std::endl;
    }

    mda_header_unload(header);

    return f_size;
}

//-----------------------------------------------------------------------------

int MDA_IO::get_rank_and_dims(std::string path, int* dims)
{

    std::FILE *fptr = std::fopen(path.c_str(), "rb");
    struct mda_header *header = mda_header_load(fptr);
    int rank = -1;
    std::fclose(fptr);

    if (header == nullptr)
    {
        logit<<"Unable to open mda file "<< path <<std::endl;
        return -1;
    }
    else if(header->data_rank == 1)
    {
        dims[0] = header->dimensions[0];
    }
    else if(header->data_rank == 2)
    {
        dims[0] = header->dimensions[0];
        dims[1] = header->dimensions[1];
    }
    else if(header->data_rank == 3)
    {
        dims[0] = header->dimensions[0];
        dims[1] = header->dimensions[1];
        dims[2] = header->dimensions[2];
    }
    else
    {
        logit<<"Unsupported mda data rank "<<header->data_rank<<" . Skipping file "<< path <<std::endl;
    }
    rank = header->data_rank;

    mda_header_unload(header);

    return rank;
}

//-----------------------------------------------------------------------------

bool load_henke_from_xdr(std::string filename, data_struct::xrf::Element_Info_Map *element_map)
{
    std::ifstream fileStream(filename);

    if (false == fileStream.good())
    {
        logit<<"Error opening file "<<filename<<std::endl;
        return false;
    }

    std::FILE* xdr_file = fopen(filename.c_str(), "rb");

    XDR *xdrstream;
#ifndef XDR_HACK
    XDR xdrs;
    xdrstream = &xdrs;
    xdrstdio_create(xdrstream, xdr_file, XDR_DECODE);
#else
    xdrstream = xdr_file;
#endif

    int num_elements;
    int num_energies;
    int num_extra_energies;
    if( xdr_int32_t( xdrstream, &num_elements) == 0)
    {
        std::fclose(xdr_file);
        return false;
    }
    if( xdr_int32_t( xdrstream, &num_energies) == 0)
    {
        std::fclose(xdr_file);
        return false;
    }

    element_map->_energies.resize(num_energies);
    //float *energy_arr = new float[num_energies];
    if( xdr_vector( xdrstream, (char *) &(element_map->_energies)[0], num_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
    {
        std::fclose(xdr_file);
        return false;
    }

    //element_map->set_energies(energy_arr, num_energies);

    //delete [] energy_arr;

    for (int i=0; i<num_elements; i++)
    {
        data_struct::xrf::Element_Info* element = element_map->get_element(i+1);
        if (element == nullptr)
        {
            element = new data_struct::xrf::Element_Info();
            element->number = i+1;
            element->name = data_struct::xrf::Element_Symbols[i];
            element_map->add_element(element);
        }
        //element->init_f_energies(num_energies);
        element->f1_atomic_scattering_real.resize(num_energies);
        element->f2_atomic_scattering_imaginary.resize(num_energies);

        //element_information.
        if( xdr_vector( xdrstream, (char *) &(element->f1_atomic_scattering_real)[0], num_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
        if( xdr_vector( xdrstream, (char *) &(element->f2_atomic_scattering_imaginary)[0], num_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
    }

    if( xdr_int32_t( xdrstream, &num_extra_energies) == 0)
    {
        std::fclose(xdr_file);
        return false;
    }

    for (int i=0; i<num_elements; i++)
    {
        data_struct::xrf::Element_Info* element = element_map->get_element(i+1);
        element->init_extra_energies(num_extra_energies);

        int element_n;
        if( xdr_int32_t( xdrstream, &element_n) == 0)
        {
            std::fclose(xdr_file);
            return false;
        }

        if( xdr_vector( xdrstream, (char *) &(element->extra_energies)[0], num_extra_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
        if( xdr_vector( xdrstream, (char *) &(element->extra_f1)[0], num_extra_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
        if( xdr_vector( xdrstream, (char *) &(element->extra_f2)[0], num_extra_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
    }

    std::fclose(xdr_file);

    return true;
}

//-----------------------------------------------------------------------------

} //end namespace file
}// end namespace io
