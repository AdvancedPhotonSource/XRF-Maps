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



namespace io
{
namespace file
{

//-----------------------------------------------------------------------------

template<typename T_real>
MDA_IO<T_real>::MDA_IO()
{
    _mda_file = nullptr;
    _mda_file_info = nullptr;
    _hasNetcdf = false;
}

//-----------------------------------------------------------------------------

template<typename T_real>
MDA_IO<T_real>::~MDA_IO()
{

    unload();

}

//-----------------------------------------------------------------------------

template<typename T_real>
bool MDA_IO<T_real>::load_scalers(std::string path)
{
    if (_mda_file != nullptr)
    {
        unload();
    }

    std::FILE* fptr = fopen(path.c_str(), "rb");

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

    if (_mda_file->header->data_rank == 1)
    {
        logE << "Cannot load mda file data rank == 1" << "\n";
        return false;
    }

    _load_scalers(true, false); //not sure how to check if we are fly scan 
    _load_meta_info(false);
    _load_extra_pvs_vector();

    return true;
}

//-----------------------------------------------------------------------------

template<typename T_real>
void MDA_IO<T_real>::unload()
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

template<typename T_real>
bool MDA_IO<T_real>::load_quantification_scalers(std::string path, data_struct::Params_Override<T_real>* override_values, bool hasNetCDF)
{
    if (override_values == nullptr)
    {
        return false;
    }
    if (_scan_info.scaler_maps.size() == 0)
    {
        if (_mda_file == nullptr)
        {
            std::FILE* fptr = std::fopen(path.c_str(), "rb");

            if (fptr == nullptr)
            {
                return false;
            }

            _mda_file = mda_load(fptr);
            std::fclose(fptr);
        }
        if (_mda_file == nullptr)
        {
            return false;
        }

        if (_mda_file->header->data_rank == 1)
        {
            logE << "Cannot load mda file data rank == 1" << "\n";
            return false;
        }

        _load_scalers(false, hasNetCDF); //not sure how to check if we are fly scan
    }

    //const data_struct::ArrayXXr<T_real>* arr = nullptr;
    const data_struct::ArrayXXr<T_real>* arr_curr = _scan_info.scaler_values(STR_SR_CURRENT);
    const data_struct::ArrayXXr<T_real>* arr_us = _scan_info.scaler_values(STR_US_IC);
    const data_struct::ArrayXXr<T_real>* arr_us_fm = _scan_info.scaler_values(STR_US_FM);
    const data_struct::ArrayXXr<T_real>* arr_ds = _scan_info.scaler_values(STR_DS_IC);
    T_real cnt_curr = 0.;
    T_real sum_curr = 0.;
    T_real cnt_us = 0.;
    T_real sum_us = 0.;
    T_real cnt_us_fm = 0.;
    T_real sum_us_fm = 0.;
    T_real cnt_ds = 0.;
    T_real sum_ds = 0.;

    for (int i = 0; i < arr_us->rows(); i++)
    {
        for (int j = 0; j < arr_us->cols(); j++)
        {
            if (arr_curr && std::isfinite((* arr_curr)(i, j)) && (*arr_curr)(i, j) > 0.)
            {
                cnt_curr += 1.0;
                sum_curr += (*arr_curr)(i, j);
            }
            if (arr_us && std::isfinite((*arr_us)(i, j)) && (*arr_us)(i, j) > 0.)
            {
                cnt_us += 1.0;
                sum_us += (*arr_us)(i, j);
            }
            if (arr_us_fm && std::isfinite((*arr_us_fm)(i, j)) && (*arr_us_fm)(i, j) > 0.)
            {
                cnt_us_fm += 1.0;
                sum_us_fm += (*arr_us_fm)(i, j);
            }
            if (arr_ds && std::isfinite((*arr_ds)(i, j)) && (*arr_ds)(i, j) > 0.)
            {
                cnt_ds += 1.0;
                sum_ds += (*arr_ds)(i, j);
            }
        }
    }



    if (arr_curr != nullptr && cnt_curr > 0.)
    {
        override_values->sr_current = sum_curr / cnt_curr;
    }
    if (arr_us != nullptr && cnt_us > 0.)
    {
        override_values->US_IC = sum_us / cnt_us;
    }
    if (arr_us_fm != nullptr && cnt_us_fm > 0.)
    {
        override_values->US_FM = sum_us_fm / cnt_us_fm;
    }
    if (arr_ds != nullptr && cnt_ds > 0.)
    {
        override_values->DS_IC = sum_ds / cnt_ds;
    }

    return true;
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool MDA_IO<T_real>::load_spectra_volume(std::string path,
                                 size_t detector_num,
                                 data_struct::Spectra_Volume<T_real>* vol,
                                 bool hasNetCDF)
{
    bool is_single_row = false;
    const data_struct::ArrayXXr<T_real>* elt_arr = nullptr;
    const data_struct::ArrayXXr<T_real>* ert_arr = nullptr;
    const data_struct::ArrayXXr<T_real>* icr_arr = nullptr;
    const data_struct::ArrayXXr<T_real>* ocr_arr = nullptr;
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
    if (_mda_file == nullptr || vol == nullptr)
    {
        return false;
    }
    logI<<"mda info ver:"<<_mda_file->header->version<<" data rank:"<<_mda_file->header->data_rank<<"\n";

    if (_mda_file->header->data_rank == 1)
    {
        logE << "Cannot load mda file data rank == 1" << "\n";
        return false;
    }

    _load_scalers(false,hasNetCDF);
    _load_meta_info(hasNetCDF);
    _load_extra_pvs_vector();

    if (_mda_file->header->data_rank == 2)
    {
        logI<<" requested rows "<< _mda_file->header->dimensions[0] << " requested cols " << _mda_file->header->dimensions[1] <<
                  " acquired rows "<< _mda_file->scan->last_point << " acquired cols " << _mda_file->scan->sub_scans[0]->last_point <<"\n";

        if(hasNetCDF)
        {
            if(_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
            {
                rows = 1;
            }
            else
            {
                rows = _mda_file->scan->last_point;
            }
            if(_mda_file->scan->sub_scans[0]->requested_points == 0 || _mda_file->scan->sub_scans[0]->last_point == 0)
            {
                cols = 1;
            }
            else
            {
                cols = _mda_file->scan->sub_scans[0]->last_point - 2; // subtract 2 because hardwre trigger goofs up last 2 cols
            }

            if ( cols < 1)
            {
                cols = 1;
            }

            vol->resize_and_zero(rows, cols, 2048);
            return true;
        }
        else
        {
            // 2000 = APS step scan, 2048 = APS xanes scan
            if(_mda_file->header->dimensions[1] == 2000 || _mda_file->header->dimensions[1] == 2048)
            {
                if((size_t)_mda_file->scan->sub_scans[0]->number_detectors-1 < detector_num)
                {
                    logE<<"Max detectors saved = "<<_mda_file->scan->sub_scans[0]->number_detectors<< "\n";
                    unload();
                    return false;
                }

                rows = 1;
                if(_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
                {
                    cols = 1;
                }
                else
                {
                    cols = _mda_file->scan->last_point;
                }
                if (cols < 1)
                {
                    cols = 1;
                }
                samples = _mda_file->header->dimensions[1];
                vol->resize_and_zero(rows, cols, 2048); //default to 2048 since it is only 2000 saved
                is_single_row = true;
            }
            else
            {
                //if not then we don't know what is dataset is.
                logE << "Don't understand this dataset layout. Can not load it.\n";
                unload();
                return false;
            }
        }
    }
    else if (_mda_file->header->data_rank == 3)
    {
        // if num_detectors == 0 them it is new fly scan with tetram
        if(_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors == 0)
        {
            logI<<" requested rows "<< _mda_file->header->dimensions[0] << " requested cols " << _mda_file->header->dimensions[1] <<
            " acquired rows "<< _mda_file->scan->last_point << " acquired cols " << _mda_file->scan->sub_scans[0]->last_point <<"\n";

            if(hasNetCDF)
            {
                if(_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
                {
                    logW<<"mda_file->scan->requested_points == 0, Settign rows = 1\n";
                    rows = 1;
                }
                else
                {
                    rows = _mda_file->scan->last_point;
                }
                
                if(_mda_file->scan->sub_scans[0]->requested_points == 0 || _mda_file->scan->sub_scans[0]->last_point == 0)
                {
                    logW<<"mda_file->scan->>sub_scans[0]->requested_points == 0, Settign cols = 1\n";
                    cols = 1;
                }
                else
                {
                    cols = _mda_file->scan->sub_scans[0]->last_point - 2; // subtract 2 because hardwre trigger goofs up last 2 cols
                }
                if(cols < 1)
                {
                    cols = 1;
                }
                vol->resize_and_zero(rows, cols, 2048);
                return true;
            }
            // TODO: might need to check if is XANES like above 
        }
        if(_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors-1 < (int)detector_num)
        {
            logE<<"Max detectors saved = "<<_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors<< "\n";
            unload();
            return false;
        }

        if(_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
        {
            rows = 1;
        }
        else
        {
            rows = _mda_file->scan->last_point;
        }
        if(_mda_file->scan->sub_scans[0]->requested_points == 0 || _mda_file->scan->sub_scans[0]->last_point == 0)
        {
            cols = 1;
        }
        else
        {
            cols = _mda_file->scan->sub_scans[0]->last_point;
        }
        if(cols < 1)
        {
            cols = 1;
        }
        samples = _mda_file->header->dimensions[2];
        if(_mda_file->header->dimensions[2] == 2000)
        {
            vol->resize_and_zero(rows, cols, 2048); //default to 2048 since it is only 2000 saved
        }
        else if(_mda_file->header->dimensions[2] > 4096) // there can be a bug in mda files that the header has incorrect dimensions
        {
            samples = _mda_file->scan->sub_scans[0]->sub_scans[0]->requested_points;
            vol->resize_and_zero(rows, cols, samples);
        }
        else
        {
            vol->resize_and_zero(rows, cols, samples);
        }
    }
    else
    {
        logE<<" No support for data rank "<< _mda_file->header->data_rank <<"\n";
        unload();
        return false;
    }

    elt_arr = _scan_info.scaler_values(STR_ELT + std::to_string(detector_num + 1));
    ert_arr = _scan_info.scaler_values(STR_ERT + std::to_string(detector_num + 1));
    icr_arr = _scan_info.scaler_values(STR_ICR + std::to_string(detector_num + 1));
    ocr_arr = _scan_info.scaler_values(STR_OCR + std::to_string(detector_num + 1));

    logI<<" Found elt="<< (elt_arr != nullptr) << " ert=" << (ert_arr != nullptr) << " in cnt=" << (icr_arr != nullptr) << " out cnt="<< (ocr_arr != nullptr) <<"\n";

    try
    {
        for(size_t i=0; i<rows; i++)
        {
            for(size_t j=0; j<cols; j++)
            {
                if (is_single_row)
                {
                    if(elt_arr)
                    {
                        (*vol)[i][j].elapsed_livetime((*elt_arr)(i, j));
                    }
                    if(ert_arr)
                    {
                        (*vol)[i][j].elapsed_realtime((*ert_arr)(i, j));
                    }
                    if(icr_arr)
                    {
                        (*vol)[i][j].input_counts((*icr_arr)(i, j));
                    }
                    if(ocr_arr)
                    {
                        (*vol)[i][j].output_counts((*ocr_arr)(i, j));
                    }
                    if(ert_arr && icr_arr && ocr_arr)
                    {
                        (*vol)[i][j].recalc_elapsed_livetime();
                    }


                    for(size_t k=0; k<samples; k++)
                    {

                        (*vol)[i][j][k] = (_mda_file->scan->sub_scans[j]->detectors_data[detector_num][k]);
                    }
                }
                else
                {
                    if(elt_arr)
                    {
                        (*vol)[i][j].elapsed_livetime((*elt_arr)(i, j));
                    }
                    if(ert_arr)
                    {
                        (*vol)[i][j].elapsed_realtime((*ert_arr)(i, j));
                    }
                    if(icr_arr)
                    {
                        (*vol)[i][j].input_counts((*icr_arr)(i, j));
                    }
                    if(ocr_arr)
                    {
                        (*vol)[i][j].output_counts((*ocr_arr)(i, j));
                    }
                    if(ert_arr && icr_arr && ocr_arr)
                    {
                        (*vol)[i][j].recalc_elapsed_livetime();
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
        logE<<"Caught exception loading mda file."<<"\n";
        std::cerr << "Exception catched : " << e.what() << "\n";
        return false;
    }

    return true;
}


//-----------------------------------------------------------------------------

template<typename T_real>
bool MDA_IO<T_real>::load_spectra_volume_with_callback(std::string path,
												const std::vector<size_t>& detector_num_arr,
                                                 bool hasNetCDF,
                                                 data_struct::Analysis_Job<T_real>* analysis_job,
                                                 size_t &out_rows,
                                                 size_t &out_cols,
												 data_struct::IO_Callback_Func_Def<T_real> callback_func,
                                                 void *user_data)
{
	// detector , index
	std::map<size_t, const data_struct::ArrayXXr<T_real>*> elt_arr_map;
	std::map<size_t, const data_struct::ArrayXXr<T_real>*> ert_arr_map;
	std::map<size_t, const data_struct::ArrayXXr<T_real>*> incnt_arr_map;
	std::map<size_t, const data_struct::ArrayXXr<T_real>*> outcnt_arr_map;
    size_t max_detecotr_num = 0;
    bool is_single_row = false;

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
    logI<<"mda info ver:"<<_mda_file->header->version<<" data rank:"<<_mda_file->header->data_rank<<"\n";

    _load_scalers(false, hasNetCDF);

    for (size_t det : detector_num_arr)
	{
		max_detecotr_num = std::max(det, max_detecotr_num);
	}

    if (analysis_job != nullptr)
    {
        if (analysis_job->theta_pv.length() > 0)
        {
            _find_theta(analysis_job->theta_pv, &analysis_job->theta);
        }
    }

    if (_mda_file->header->data_rank == 2)
    {
        logit_s<<" requested rows "<< _mda_file->header->dimensions[0] << " requested cols " << _mda_file->header->dimensions[1] <<
                  " acquired rows "<< _mda_file->scan->last_point << " acquired cols " << _mda_file->scan->sub_scans[0]->last_point <<"\n";

        if(hasNetCDF)
        {
            if(_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
            {
                out_rows = 1;
            }
            else
            {
                out_rows = _mda_file->scan->last_point;
            }
            if(_mda_file->scan->sub_scans[0]->requested_points == 0 || _mda_file->scan->sub_scans[0]->last_point == 0)
            {
                out_cols = 1;
            }
            else
            {
                out_cols = _mda_file->scan->sub_scans[0]->last_point - 2; // subtract 2 because hardwre trigger goofs up last 2 cols
            }
            return true;
        }
        else
        {
            if(_mda_file->header->dimensions[1] == 2000 || _mda_file->header->dimensions[1] == 2048)
            {
                if((size_t)_mda_file->scan->sub_scans[0]->number_detectors-1 < max_detecotr_num)
                {
                    logE<<"Max detectors saved = "<<_mda_file->scan->sub_scans[0]->number_detectors<< "\n";
                    unload();
                    return false;
                }

                out_rows = 1;
                if (_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
                {
                    out_cols = 1;
                }
                else
                {
                    out_cols = _mda_file->scan->last_point;
                }
                samples = _mda_file->header->dimensions[1];
                is_single_row = true;
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

        if((size_t)_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors-1 < max_detecotr_num)
        {
            logE<<"Max detectors saved = "<<_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors<< "\n";
            unload();
            return false;
        }

        if(_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
        {
            out_rows = 1;
        }
        else
        {
            out_rows = _mda_file->scan->last_point;
        }
        if(_mda_file->scan->sub_scans[0]->requested_points == 0 || _mda_file->scan->sub_scans[0]->last_point == 0)
        {
            out_cols = 1;
        }
        else
        {
            out_cols = _mda_file->scan->sub_scans[0]->last_point;
        }
        samples = _mda_file->header->dimensions[2];
    }
    else
    {
        logE<<"No support for data rank "<< _mda_file->header->data_rank <<"\n";
        unload();
        return false;
    }

	for (size_t detector_num : detector_num_arr)
	{
        elt_arr_map[detector_num] = _scan_info.scaler_values(STR_ELT + std::to_string(detector_num + 1));
        ert_arr_map[detector_num] = _scan_info.scaler_values(STR_ERT + std::to_string(detector_num + 1));
        incnt_arr_map[detector_num] = _scan_info.scaler_values(STR_ICR + std::to_string(detector_num + 1));
        outcnt_arr_map[detector_num] = _scan_info.scaler_values(STR_OCR + std::to_string(detector_num + 1));

/*
		if (detector_struct->fit_params_override_dict.scaler_pvs.count("SRCURRENT") > 0)
		{
			find_scaler_index(_mda_file, detector_struct->fit_params_override_dict.scaler_pvs.at("SRCURRENT"), tmp_val, units);
			detector_struct->quant_standard.sr_current = (tmp_val);
		}
		if (detector_struct->fit_params_override_dict.scaler_pvs.count(STR_US_IC) > 0)
		{
			find_scaler_index(_mda_file, detector_struct->fit_params_override_dict.scaler_pvs.at(STR_US_IC), tmp_val, units);
			detector_struct->quant_standard.US_IC = (tmp_val);
		}
		if (detector_struct->fit_params_override_dict.scaler_pvs.count(STR_DS_IC) > 0)
		{
			find_scaler_index(_mda_file, detector_struct->fit_params_override_dict.scaler_pvs.at(STR_DS_IC), tmp_val, units);
			detector_struct->quant_standard.DS_IC = (tmp_val);
		}
*/
	}


    try
    {
        for(int i=0; i< out_rows; i++)
        {
            for(int j=0; j< out_cols; j++)
            {
                for(size_t detector_num : detector_num_arr)
                {
                    data_struct::Spectra<T_real>* spectra = new data_struct::Spectra<T_real>(samples);

                    if (is_single_row)
                    {
                        if(elt_arr_map.at(detector_num) != nullptr)
                        {
                            spectra->elapsed_livetime((*elt_arr_map.at(detector_num))(i,j));
                        }
                        if(ert_arr_map.at(detector_num) != nullptr)
                        {
                            spectra->elapsed_realtime((*ert_arr_map.at(detector_num))(i, j));
                        }
                        if(incnt_arr_map.at(detector_num) != nullptr)
                        {
                            spectra->input_counts((*incnt_arr_map.at(detector_num))(i, j));
                        }
                        if(outcnt_arr_map.at(detector_num) != nullptr)
                        {
                            spectra->output_counts((*outcnt_arr_map.at(detector_num))(i, j));
                        }
                        spectra->recalc_elapsed_livetime();


                        for(size_t k=0; k<samples; k++)
                        {

                            (*spectra)[k] = (_mda_file->scan->sub_scans[j]->detectors_data[detector_num][k]);
                        }
                        callback_func(i, j, out_rows, out_cols, detector_num, spectra, user_data);
                    }
                    else
                    {
                        if(elt_arr_map.at(detector_num) != nullptr)
                        {
                             spectra->elapsed_livetime((*elt_arr_map.at(detector_num))(i, j));
                        }
                        if(ert_arr_map.at(detector_num) != nullptr)
                        {
                            spectra->elapsed_realtime((*ert_arr_map.at(detector_num))(i, j));
                        }
                        if(incnt_arr_map.at(detector_num) != nullptr)
                        {
                            spectra->input_counts((*incnt_arr_map.at(detector_num))(i, j));
                        }
                        if(outcnt_arr_map.at(detector_num) != nullptr)
                        {
                            spectra->output_counts((*outcnt_arr_map.at(detector_num))(i, j));
                        }
                        spectra->recalc_elapsed_livetime();
                        

                        for(size_t k=0; k<samples; k++)
                        {
                            (*spectra)[k] = (_mda_file->scan->sub_scans[i]->sub_scans[j]->detectors_data[detector_num][k]);
                        }
                        callback_func(i, j, out_rows, out_cols, detector_num, spectra, user_data);
                    }
                }
            }
        }
    }
    catch(std::exception& e)
    {
        logE << "Exception catched : " << e.what() << "\n";
        return false;
    }

    return true;
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool MDA_IO<T_real>::load_integrated_spectra(std::string path,
		size_t detector_num,
		data_struct::Spectra<T_real>* out_integrated_spectra,
		bool hasNetCDF)
{
	//index per row and col
    const data_struct::ArrayXXr<T_real>* elt_arr = nullptr;
    const data_struct::ArrayXXr<T_real>* ert_arr = nullptr;
    const data_struct::ArrayXXr<T_real>* icr_arr = nullptr;
    const data_struct::ArrayXXr<T_real>* ocr_arr = nullptr;
	bool is_single_row = false;

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
	if (_mda_file == nullptr || out_integrated_spectra == nullptr)
	{
		logE << "_mda_file or out_integrated_spectra == nullptr\n";
		return false;
	}
	logI << "mda info ver:" << _mda_file->header->version << " data rank:" << _mda_file->header->data_rank << "\n";

	if (_mda_file->header->data_rank == 1)
	{
		logE << "Cannot load mda file data rank == 1\n";
		return false;
	}

	_load_scalers(false, hasNetCDF);
	_load_meta_info(hasNetCDF);
	_load_extra_pvs_vector();

	if (_mda_file->header->data_rank == 2)
	{
		logI << " requested rows " << _mda_file->header->dimensions[0] << " requested cols " << _mda_file->header->dimensions[1] <<
			" acquired rows " << _mda_file->scan->last_point << " acquired cols " << _mda_file->scan->sub_scans[0]->last_point << "\n";

		if (hasNetCDF)
		{
			if (_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
            {
				rows = 1;
            }
			else
            {
				rows = _mda_file->scan->last_point;
            }
			if (_mda_file->scan->sub_scans[0]->requested_points == 0 || _mda_file->scan->sub_scans[0]->last_point == 0)
            {
				cols = 1;
            }
			else
            {
				cols = _mda_file->scan->sub_scans[0]->last_point - 2; // subtract 2 because hardwre trigger goofs up last 2 cols
            }
			out_integrated_spectra->resize(2048);
			return true;
		}
		else
		{
			if (_mda_file->header->dimensions[1] == 2000 || _mda_file->header->dimensions[1] == 2048)
			{
				if ((size_t)_mda_file->scan->sub_scans[0]->number_detectors - 1 < detector_num)
				{
					logE << "Max detectors saved = " << _mda_file->scan->sub_scans[0]->number_detectors << "\n";
					unload();
					return false;
				}

				rows = 1;
				if (_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
                {
					cols = 1;
                }
				else
                {
					cols = _mda_file->scan->last_point;
                }
				samples = _mda_file->header->dimensions[1];
				out_integrated_spectra->resize(2048); //default to 2048 since it is only 2000 saved
				is_single_row = true;
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

		if ((size_t)_mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors - 1 < detector_num)
		{
			logE << "Max detectors saved = " << _mda_file->scan->sub_scans[0]->sub_scans[0]->number_detectors << "\n";
			unload();
			return false;
		}

		if (_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
        {
			rows = 1;
        }
		else
        {
			rows = _mda_file->scan->last_point;
        }
		if (_mda_file->scan->sub_scans[0]->requested_points == 0 ||  _mda_file->scan->sub_scans[0]->last_point == 0)
        {
			cols = 1;
        }
		else
        {
			cols = _mda_file->scan->sub_scans[0]->last_point;
        }
		samples = _mda_file->header->dimensions[2];
		if (_mda_file->header->dimensions[2] == 2000)
		{
			out_integrated_spectra->resize(2048); //default to 2048 since it is only 2000 saved
		}
		else if (_mda_file->header->dimensions[2] > 4096) // there can be a bug in mda files that the header has incorrect dimensions
		{
			samples = _mda_file->scan->sub_scans[0]->sub_scans[0]->requested_points;
			out_integrated_spectra->resize(samples);
		}
		else
		{
			out_integrated_spectra->resize(samples);
		}
	}
	else
	{
		logE << " No support for data rank " << _mda_file->header->data_rank << "\n";
		unload();
		return false;
	}

	out_integrated_spectra->setZero(samples);

    elt_arr = _scan_info.scaler_values(STR_ELT + std::to_string(detector_num + 1));
    ert_arr = _scan_info.scaler_values(STR_ERT + std::to_string(detector_num + 1));
    icr_arr = _scan_info.scaler_values(STR_ICR + std::to_string(detector_num + 1));
    ocr_arr = _scan_info.scaler_values(STR_OCR + std::to_string(detector_num + 1));

    logI << " Found elt=" << (elt_arr != nullptr) << " ert=" << (ert_arr != nullptr) << " in cnt=" << (icr_arr != nullptr) << " out cnt=" << (ocr_arr != nullptr) << "\n";

	try
	{
		for (size_t i = 0; i < rows; i++)
		{
			for (size_t j = 0; j < cols; j++)
			{
				if (is_single_row)
				{
					for (size_t k = 0; k < samples; k++)
					{
						(*out_integrated_spectra)[k] += (_mda_file->scan->sub_scans[j]->detectors_data[detector_num][k]);
					}
				}
				else
				{
					for (size_t k = 0; k < samples; k++)
					{
						(*out_integrated_spectra)[k] += (_mda_file->scan->sub_scans[i]->sub_scans[j]->detectors_data[detector_num][k]);
					}
				}
			}
		}
	}
	catch (std::exception& e)
	{
		logE << "Caught exception loading mda file." << "\n";
		std::cerr << "Exception catched : " << e.what() << "\n";
		return false;
	}
    
    if (elt_arr)
    {
        out_integrated_spectra->elapsed_livetime(elt_arr->sum());
    }
    if (ert_arr)
    {
        out_integrated_spectra->elapsed_realtime(ert_arr->sum());
    }
    if (icr_arr)
    {
        out_integrated_spectra->input_counts(icr_arr->sum());
    }
    if (ocr_arr)
    {
        out_integrated_spectra->output_counts(ocr_arr->sum());
    }
	out_integrated_spectra->recalc_elapsed_livetime();
	
	return true;
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool MDA_IO<T_real>::_find_theta(std::string pv_name, float* theta_out)
{
    float *tmpf;
    double *tmpd;
    for (int16_t i = 0; i < _mda_file->extra->number_pvs; i++)
    {
        struct mda_pv * pv = _mda_file->extra->pvs[i];
        if(pv == nullptr)
        {
            continue;
        }
        if(strcmp( pv->name, pv_name.c_str()) ==0 )
        {

            switch (pv->type)
            {
            case EXTRA_PV_FLOAT:
                tmpf = (float*)pv->values;
                *theta_out = *tmpf;
                break;
            case EXTRA_PV_DOUBLE:
                tmpd = (double*)pv->values;
                *theta_out = *tmpd;
                break;
            default:
                *theta_out = 0.0f;
                break;
            }
            return true;
        }
    }
    return false;
}


//-----------------------------------------------------------------------------

template<typename T_real>
void MDA_IO<T_real>::_load_scalers(bool load_int_spec, bool hasNetCDF)
{
    if (_mda_file == nullptr)
    {
        return;
    }
    size_t rows = 0;
    size_t cols = 0;
    bool single_row_scan = false;

    if (_mda_file->header->data_rank == 2)
    {
        if (_hasNetcdf == false && (_mda_file->header->dimensions[1] == 2000 || _mda_file->header->dimensions[1] == 2048))
        {
            single_row_scan = true;
        }
    }

    std::string beamline = "";

    //save scalers
    if (single_row_scan)
    {
        rows = 1;
        if (_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
        {
            rows = 1;
        }
        else
        {
            cols = _mda_file->scan->last_point;
            if(hasNetCDF)
            {
                cols = std::max(1, _mda_file->scan->last_point - 2);
            }
        }

        size_t iamt = std::min((size_t)_mda_file->scan->last_point, cols);
        for (size_t i = 0; i < iamt; i++)
        {
            for (int k = 0; k < _mda_file->scan->number_detectors; k++)
            {
                if (i == 0)
                {
                    data_struct::Scaler_Map<T_real> s_map;
                    s_map.values.resize(rows, cols);
                    s_map.values.setZero(rows, cols);
                    s_map.name = std::string(_mda_file->scan->detectors[k]->name);
                    std::string label = "";
                    bool is_time_normalized = false;
                    if (data_struct::Scaler_Lookup::inst()->search_pv(s_map.name, label, is_time_normalized, beamline))
                    {
                        s_map.name = label;
                    }
                    s_map.time_normalized = is_time_normalized;
                    s_map.unit = std::string(_mda_file->scan->detectors[k]->unit);
                    _scan_info.scaler_maps.push_back(s_map);
                }

                if (std::isfinite(_mda_file->scan->detectors_data[k][i]))
                {
                    _scan_info.scaler_maps[k].values(0, i) = _mda_file->scan->detectors_data[k][i];
                }
            }

            if (_mda_file->scan->sub_scans != nullptr && load_int_spec)
            {
                for (int32_t d = 0; d < _mda_file->scan->sub_scans[i]->number_detectors; d++)
                {
                    data_struct::ArrayTr<T_real>* int_spec;
                    if (_integrated_spectra_map.count(d) == 0)
                    {
                        // if this is the first one then zero it out
                        int_spec = &(_integrated_spectra_map[d]);
                        int_spec->resize(_mda_file->scan->sub_scans[i]->requested_points);
                        int_spec->setZero(_mda_file->scan->sub_scans[i]->requested_points);
                    }
                    else
                    {
                        int_spec = &(_integrated_spectra_map[d]);
                    }
                    for (int32_t m = 0; m < _mda_file->scan->sub_scans[i]->last_point; m++)
                    {
                        (*int_spec)[m] += (_mda_file->scan->sub_scans[i]->detectors_data[d][m]);
                    }
                }
            }
        }
        
    }
    else
    {
        if (_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
        {
            rows = 1;
        }
        else
        {
            rows = _mda_file->scan->last_point;
        }
        if (_mda_file->scan->sub_scans[0]->requested_points == 0 || _mda_file->scan->sub_scans[0]->last_point == 0)
        {
            cols = 1;
        }
        else
        {
            cols = _mda_file->scan->sub_scans[0]->last_point;
            if(hasNetCDF)
            {
                cols = std::max(1, _mda_file->scan->sub_scans[0]->last_point - 2);
            }
        }

        for (int32_t i = 0; i < _mda_file->scan->last_point; i++)
        {
            size_t jamt = std::min((size_t)_mda_file->scan->sub_scans[i]->last_point, cols);
            for (size_t j = 0; j < jamt; j++)
            {
                for (int k = 0; k < _mda_file->scan->sub_scans[i]->number_detectors; k++)
                {
                    if (i == 0 && j == 0)
                    {
                        data_struct::Scaler_Map<T_real> s_map;
                        s_map.values.resize(rows, cols);
                        s_map.values.setZero(rows, cols);
                        s_map.name = std::string(_mda_file->scan->sub_scans[0]->detectors[k]->name);
                        std::string label = "";
                        bool is_time_normalized = false;
                        if (data_struct::Scaler_Lookup::inst()->search_pv(s_map.name, label, is_time_normalized, beamline))
                        {
                            s_map.name = label;
                        }
                        s_map.time_normalized = is_time_normalized;
                        s_map.unit = std::string(_mda_file->scan->sub_scans[0]->detectors[k]->unit);
                        _scan_info.scaler_maps.push_back(s_map);
                    }

                    if (std::isfinite(_mda_file->scan->sub_scans[i]->detectors_data[k][j]))
                    {
                        _scan_info.scaler_maps[k].values(i, j) = _mda_file->scan->sub_scans[i]->detectors_data[k][j];
                    }
                }
                if (_mda_file->scan->sub_scans[i]->sub_scans != nullptr && load_int_spec)
                {
                    for (int32_t d = 0; d < _mda_file->scan->sub_scans[i]->sub_scans[j]->number_detectors; d++)
                    {
                        data_struct::ArrayTr<T_real>* int_spec;
                        if (_integrated_spectra_map.count(d) == 0)
                        {
                            // if this is the first one then zero it out
                            int_spec = &(_integrated_spectra_map[d]);
                            int_spec->resize(_mda_file->scan->sub_scans[i]->sub_scans[j]->requested_points);
                            int_spec->setZero(_mda_file->scan->sub_scans[i]->sub_scans[j]->requested_points);
                        }
                        else
                        {
                            int_spec = &(_integrated_spectra_map[d]);
                        }
                        for (int32_t m = 0; m < _mda_file->scan->sub_scans[i]->sub_scans[j]->last_point; m++)
                        {
                            (*int_spec)[m] += (_mda_file->scan->sub_scans[i]->sub_scans[j]->detectors_data[d][m]);
                        }
                    }
                }
            }
        }
    }

    std::vector<std::string> pv_names;
    for (const auto& itr : _scan_info.scaler_maps)
    {
        pv_names.push_back(itr.name);
    }
    std::string time_pv = "";
    double time_clock = 0.0;
    if (data_struct::Scaler_Lookup::inst()->search_for_timing_info(pv_names, time_pv, time_clock, beamline))
    {
        const data_struct::ArrayXXr<T_real>* time_array = _scan_info.scaler_values(time_pv);
        if (time_array != nullptr)
        {
            for (auto& itr : _scan_info.scaler_maps)
            {
                if (itr.time_normalized)
                {
                    itr.values = itr.values / (*time_array / time_clock);
                    itr.values = itr.values.unaryExpr([](T_real v) { return std::isfinite(v) ? v : (T_real)0.0; });
                }
            }
        }
    }

    auto summed_scalers = data_struct::Scaler_Lookup::inst()->get_summed_scaler_list(beamline);
    if (summed_scalers != nullptr)
    {
        for (const auto& itr : *summed_scalers)
        {
            data_struct::Scaler_Map<T_real> s_map;
            s_map.name = itr.scaler_name;
            s_map.values.resize(rows, cols);
            s_map.values.setZero(rows, cols);
            for (const auto& sitr : itr.scalers_to_sum)
            {
                const data_struct::ArrayXXr<T_real>* arr = _scan_info.scaler_values(sitr);
                if (arr != nullptr)
                {
                    s_map.values += (*arr);
                }
            }
            _scan_info.scaler_maps.push_back(s_map);
        }
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
void MDA_IO<T_real>::_load_extra_pvs_vector()
{
    std::map<std::string, std::string> name_hash;
    if (_mda_file == nullptr)
    {
        return;
    }

	if (_mda_file->extra != nullptr)
	{
		for (int16_t i = 0; i < _mda_file->extra->number_pvs; i++)
		{
			std::string str_val;
			short* s_val;
			int* i_val;
			float* f_val;
			double* d_val;

			struct mda_pv* pv = _mda_file->extra->pvs[i];
			if (pv == nullptr)
			{
				continue;
			}
			data_struct::Extra_PV e_pv;
			switch (pv->type)
			{

			case EXTRA_PV_STRING:
				e_pv.value = std::string(pv->values);
				break;
				//case EXTRA_PV_INT8:

				//    break;
			case EXTRA_PV_INT16:
				s_val = (short*)pv->values;
				e_pv.value = std::to_string(*s_val);
				break;
			case EXTRA_PV_INT32:
				i_val = (int*)pv->values;
				e_pv.value = std::to_string(*i_val);
				break;
			case EXTRA_PV_FLOAT:
				f_val = (float*)pv->values;
				e_pv.value = std::to_string(*f_val);
				break;
			case EXTRA_PV_DOUBLE:
				d_val = (double*)pv->values;
				e_pv.value = std::to_string(*d_val);
				break;

			}

			if (pv->name != nullptr)
			{
				e_pv.name = std::string(pv->name);
			}

            // rename extra pv's from scaler lookup table
            std::string beamline = "";
            std::string label = "";
            bool is_time_normalized = false;
            if (data_struct::Scaler_Lookup::inst()->search_pv(e_pv.name, label, is_time_normalized, beamline))
            {
                e_pv.name = label;
            }
            
            // check for tetramm column names and append to scaler_maps
            int t1idx = e_pv.name.find("tmm1:");
            int t2idx = e_pv.name.find("tmm2:");
            int ndidx = e_pv.name.find("NDArrayAddress");
            int nameidx = e_pv.name.find("Name");
            if((t1idx != std::string::npos || t2idx != std::string::npos) && nameidx != std::string::npos)
            {
                std::vector<std::string> tokens;
                std::stringstream ss(e_pv.name);
                std::string token;
                while (std::getline(ss, token, ':')) 
                {
                    tokens.push_back(token);
                }
                if(tokens.size() > 2)
                {
                    name_hash[tokens[1]] = e_pv.value;
                }
            }
            if((t1idx != std::string::npos || t2idx != std::string::npos) && ndidx != std::string::npos)
            {
                std::vector<std::string> tokens;
                std::stringstream ss(e_pv.name);
                std::string token;
                while (std::getline(ss, token, ':')) 
                {
                    tokens.push_back(token);
                }
                if(tokens.size() > 2 && _scan_info.scaler_maps.size() > 0)
                {
                    auto o_map = _scan_info.scaler_maps.at(0);
                    data_struct::Scaler_Map<T_real> s_map;
                    s_map.values.resize(o_map.values.rows(), o_map.values.cols());
                    s_map.values.setZero(o_map.values.rows(), o_map.values.cols());
                    if (data_struct::Scaler_Lookup::inst()->search_pv(tokens[0] + ":" + tokens[1], label, is_time_normalized, beamline))
                    {
                        s_map.name = label;
                    }
                    else
                    {
                        s_map.name = e_pv.name;
                    }
                    
                    if(t1idx != std::string::npos)
                    {
                        s_map.unit = "tetra1_" + e_pv.value;
                    }
                    else
                    {
                        s_map.unit = "tetra2_" + e_pv.value;
                    }
                    _scan_info.scaler_maps.push_back(s_map);  
                }
            }

			if (pv->description != nullptr)
			{
				e_pv.description = std::string(pv->description);
			}

			if (pv->unit != nullptr)
			{
				e_pv.unit = std::string(pv->unit);
			}
			_scan_info.extra_pvs.push_back(e_pv);
		}

        for(auto &itr : _scan_info.scaler_maps)
        {
            if(name_hash.count(itr.name) > 0)
            {
                itr.name = name_hash.at(itr.name);
            }
        }
	}
	else
	{
		data_struct::Extra_PV e_pv;
		e_pv.name = "NULL";
		e_pv.description = "NULL";
		e_pv.unit = "NULL";
		e_pv.value = "NULL";
		_scan_info.extra_pvs.push_back(e_pv);
	}

}

//-----------------------------------------------------------------------------

template<typename T_real>
void MDA_IO<T_real>::_load_meta_info(bool hasNetCDF)
{
    bool single_row_scan = false;

    if (_mda_file == nullptr)
    {
        return;
    }

    try
    {

        if (_mda_file->scan->scan_rank > 1)
        {
            if (_mda_file->header->data_rank == 2)
            {
                if (_mda_file->header->dimensions[1] == 2000 || _mda_file->header->dimensions[1] == 2048)
                {
                    single_row_scan = true;
                }
            }

            if (single_row_scan)
            {
                _scan_info.meta_info.requested_rows = 1;
                if(_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
                {
                    _scan_info.meta_info.requested_cols = 1;
                }
                else
                {
                    _scan_info.meta_info.requested_cols = _mda_file->scan->last_point;
                }
                if(hasNetCDF)
                {
                    _scan_info.meta_info.requested_cols = _scan_info.meta_info.requested_cols - 2; // subtract 2 because hardwre trigger goofs up last 2 cols                    
                }
                if (_scan_info.meta_info.requested_cols < 1)
                {
                    _scan_info.meta_info.requested_cols = 1;
                }
                _scan_info.meta_info.y_axis.resize(1);
                _scan_info.meta_info.y_axis.setZero(1);
                _scan_info.meta_info.x_axis.resize( _scan_info.meta_info.requested_cols);
                _scan_info.meta_info.x_axis.setZero( _scan_info.meta_info.requested_cols);
                size_t iamt = std::min((size_t)_mda_file->scan->last_point, (size_t)_scan_info.meta_info.requested_cols);
                for (size_t i = 0; i < iamt; i++)
                {
                    _scan_info.meta_info.x_axis(i) = _mda_file->scan->positioners_data[0][i];
                }
            }
            else
            {

                if(_mda_file->scan->requested_points == 0 || _mda_file->scan->last_point == 0)
                {
                    _scan_info.meta_info.requested_rows = 1;
                }
                else
                {
                    _scan_info.meta_info.requested_rows = _mda_file->scan->last_point;
                }
                if(_mda_file->scan->sub_scans[0]->requested_points == 0 || _mda_file->scan->sub_scans[0]->last_point == 0)
                {
                    _scan_info.meta_info.requested_cols = 1;
                }
                else
                {
                    _scan_info.meta_info.requested_cols = _mda_file->scan->sub_scans[0]->last_point;
                }
                
                if(hasNetCDF)
                {
                    _scan_info.meta_info.requested_cols = _scan_info.meta_info.requested_cols - 2; // subtract 2 because hardwre trigger goofs up last 2 cols
                    if (_scan_info.meta_info.requested_cols < 1)
                    {
                        _scan_info.meta_info.requested_cols = 1;
                    }
                }
                // resize and zero y axis
                _scan_info.meta_info.y_axis.resize(_scan_info.meta_info.requested_rows);
                _scan_info.meta_info.y_axis.setZero(_scan_info.meta_info.requested_rows);
                // if positioner exists, then read it
                if (_mda_file->scan->number_positioners > 0)
                {
                    // save y axis
                    for (int32_t i = 0; i < _mda_file->scan->last_point; i++)
                    {
                        _scan_info.meta_info.y_axis(i) = _mda_file->scan->positioners_data[0][i];
                    }
                }
                // resize and zero x axis
                _scan_info.meta_info.x_axis.resize(_scan_info.meta_info.requested_cols);
                _scan_info.meta_info.x_axis.setZero(_scan_info.meta_info.requested_cols);
                // if positioner exists, then read it
                if (_mda_file->scan->sub_scans[0]->number_positioners > 0)
                {
                    size_t iamt = std::min((size_t)_mda_file->scan->sub_scans[0]->last_point, (size_t)_scan_info.meta_info.requested_cols);
                    // save x axis
                    for (size_t i = 0; i < iamt; i++)
                    {
                        _scan_info.meta_info.x_axis(i) = _mda_file->scan->sub_scans[0]->positioners_data[0][i];
                    }
                }
            }
        }

        //set default theta to 0.0
        _scan_info.meta_info.theta = 0.0;

        if (_mda_file->extra != nullptr && _theta_pv_str.length() > 0)
        {    
            struct mda_pv* pv = nullptr;
            //find theta by param_override->theta_pv in extra names
            for (int16_t i = 0; i < _mda_file->extra->number_pvs; i++)
            {
                pv = _mda_file->extra->pvs[i];
                if (pv == nullptr)
                {
                    continue;
                }

                if (pv->name != nullptr && _theta_pv_str.compare(pv->name) == 0)
                {
                    break;
                }
            }
            if (pv != nullptr)
            {
                switch (pv->type)
                {
                //case EXTRA_PV_STRING:
                //    break;
                case EXTRA_PV_FLOAT:
                    _scan_info.meta_info.theta = *((float*)pv->values);
                    break;
                case EXTRA_PV_DOUBLE:
                    _scan_info.meta_info.theta = *((double*)pv->values);
                    break;
                default:
                    break;
                }
            }
        }

        if (_mda_file->scan->time != nullptr)
        {
            _scan_info.meta_info.scan_time_stamp = std::string(_mda_file->scan->time);
        }

        if (_mda_file->scan->name != nullptr)
        {
            _scan_info.meta_info.name = std::string(_mda_file->scan->name);
        }
    }
    catch (...)
    {
        logE << "loading meta data" << "\n";
    }
}

//-----------------------------------------------------------------------------

template<typename T_real>
data_struct::ArrayTr<T_real>* MDA_IO<T_real>::get_integrated_spectra(unsigned int detector)
{
    if (_integrated_spectra_map.count(detector) > 0)
    {
        return &(_integrated_spectra_map.at(detector));
    }
    return nullptr;
}

//-----------------------------------------------------------------------------

template<typename T_real>
bool MDA_IO<T_real>::load_henke_from_xdr(std::string filename)
{
    data_struct::Element_Info_Map<T_real>* element_map = data_struct::Element_Info_Map<T_real>::inst();

    std::ifstream fileStream(filename);

    if (false == fileStream.good())
    {
        logE << "Opening file " << filename << "\n";
        return false;
    }

    std::FILE* xdr_file = fopen(filename.c_str(), "rb");

    XDR* xdrstream;
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
    if (xdr_int32_t(xdrstream, &num_elements) == 0)
    {
        std::fclose(xdr_file);
        return false;
    }
    if (xdr_int32_t(xdrstream, &num_energies) == 0)
    {
        std::fclose(xdr_file);
        return false;
    }

    element_map->_energies.resize(num_energies);
    //float *energy_arr = new float[num_energies];
    if (xdr_vector(xdrstream, (char*)&(element_map->_energies)[0], num_energies, sizeof(float), (xdrproc_t)xdr_float) == false)
    {
        std::fclose(xdr_file);
        return false;
    }

    //element_map->set_energies(energy_arr, num_energies);

    //delete [] energy_arr;

    for (int i = 0; i < num_elements; i++)
    {
        data_struct::Element_Info<T_real>* element = element_map->get_element(i + 1);
        if (element == nullptr)
        {
            element = new data_struct::Element_Info<T_real>();
            element->number = i + 1;
            element->name = data_struct::Element_Symbols[i + 1];
            element_map->add_element(element);
        }
        //element->init_f_energies(num_energies);
        element->f1_atomic_scattering_real.resize(num_energies);
        element->f2_atomic_scattering_imaginary.resize(num_energies);

        //element_information.
        if (xdr_vector(xdrstream, (char*)&(element->f1_atomic_scattering_real)[0], num_energies, sizeof(float), (xdrproc_t)xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
        if (xdr_vector(xdrstream, (char*)&(element->f2_atomic_scattering_imaginary)[0], num_energies, sizeof(float), (xdrproc_t)xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
    }

    if (xdr_int32_t(xdrstream, &num_extra_energies) == 0)
    {
        std::fclose(xdr_file);
        return false;
    }

    for (int i = 0; i < num_elements; i++)
    {
        data_struct::Element_Info<T_real>* element = element_map->get_element(i + 1);
        element->init_extra_energies(num_extra_energies);

        int element_n;
        if (xdr_int32_t(xdrstream, &element_n) == 0)
        {
            std::fclose(xdr_file);
            return false;
        }

        if (xdr_vector(xdrstream, (char*)&(element->extra_energies)[0], num_extra_energies, sizeof(float), (xdrproc_t)xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
        if (xdr_vector(xdrstream, (char*)&(element->extra_f1)[0], num_extra_energies, sizeof(float), (xdrproc_t)xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
        if (xdr_vector(xdrstream, (char*)&(element->extra_f2)[0], num_extra_energies, sizeof(float), (xdrproc_t)xdr_float) == false)
        {
            std::fclose(xdr_file);
            return false;
        }
    }

    std::fclose(xdr_file);

    return true;
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int mda_get_multiplied_dims(std::string path)
{
    int f_size = -1;


    std::FILE* fptr = std::fopen(path.c_str(), "rb");
    if (fptr != nullptr)
    {
        struct mda_header* header = mda_header_load(fptr);

        std::fclose(fptr);

        if (header == nullptr)
        {
            logE << "Unable to open mda file " << path << "\n";
            return f_size;
        }
        else if (header->data_rank == 1)
        {
            f_size = header->dimensions[0];
        }
        else if (header->data_rank == 2 || header->data_rank == 3)
        {
            f_size = header->dimensions[0] * header->dimensions[1];
        }
        else
        {
            logW << "Unsupported mda data rank " << header->data_rank << " . Skipping file " << path << "\n";
        }

        mda_header_unload(header);
    }
    return f_size;
}

//-----------------------------------------------------------------------------

int mda_get_rank_and_dims(std::string path, size_t* dims)
{

    std::FILE* fptr = std::fopen(path.c_str(), "rb");
    struct mda_header* header = mda_header_load(fptr);
    int rank = -1;
    std::fclose(fptr);

    if (header == nullptr)
    {
        logE << "Unable to open mda file " << path << "\n";
        return -1;
    }
    else if (header->data_rank == 1)
    {
        dims[0] = header->dimensions[0];
    }
    else if (header->data_rank == 2)
    {
        dims[0] = header->dimensions[0];
        dims[1] = header->dimensions[1];
    }
    else if (header->data_rank == 3)
    {
        dims[0] = header->dimensions[0];
        dims[1] = header->dimensions[1];
        dims[2] = header->dimensions[2];
    }
    else
    {
        logW << "Unsupported mda data rank " << header->data_rank << " . Skipping file " << path << "\n";
    }
    rank = (int)header->data_rank;

    mda_header_unload(header);

    return rank;
}

//-----------------------------------------------------------------------------


TEMPLATE_CLASS_DLL_EXPORT MDA_IO<float>;
TEMPLATE_CLASS_DLL_EXPORT MDA_IO<double>;

} //end namespace file
}// end namespace io
