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

	if (_mda_file_info == nullptr)
	{
        std::cout << "Error loading mda file:" << _filename;
		return;
	}

    std::cout<<"mda info ver:"<<_mda_file_info->version<<" data rank:"<<_mda_file_info->data_rank;
    for(int16_t i = 0; i < _mda_file_info->data_rank; i++)
    {
        std::cout<<" dims["<<i<<"]:"<<_mda_file_info->dimensions[i];
    }
    std::cout<<std::endl;
    std::fclose(fptr);

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
    std::cout<<"mda info ver:"<<_mda_file->header->version<<" data rank:"<<_mda_file->header->data_rank;
    //long total = 1;

    for(int16_t i = 0; i < _mda_file->header->data_rank; i++)
    {
        std::cout<<" dims["<<i<<"]:"<<_mda_file->header->dimensions[i];
        //total *= _mda_file->header->dimensions[i];
    }

    //dset->alloc(_mda_file->header->data_rank, _mda_file->header->dimensions);
    for(unsigned int j = 0; j < _mda_file->header->data_rank; j++)
    {

    }
    std::cout<<std::endl;
    std::cout<<"d "<<(_mda_file->scan->sub_scans[0]->sub_scans[0]->detectors_data[0][0])<<std::endl;
    std::cout<<"d "<<(_mda_file->scan->sub_scans[0]->sub_scans[0]->detectors_data[0][1])<<std::endl;

    std::fclose(fptr);

    return true;

}

//-----------------------------------------------------------------------------

int MDA_IO::_find_2d_detector_index(std::string det_name, int detector_num, real_t& val)
{

    for(int k=0; k<_mda_file->scan->number_detectors; k++)
    {
        //std::cout<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->name << std::endl;
        //std::cout<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->description << std::endl;
        if(strcmp(_mda_file->scan->detectors[k]->name, det_name.c_str())  == 0)
        {
            val = _mda_file->scan->detectors_data[k][detector_num];
            return k;
        }
    }

    for(int k=0; k<_mda_file->scan->sub_scans[0]->number_detectors; k++)
    {
        //std::cout<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->name << std::endl;
        //std::cout<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->description << std::endl;
        if(strcmp(_mda_file->scan->sub_scans[0]->detectors[k]->name, det_name.c_str())  == 0)
        {
            val = _mda_file->scan->sub_scans[0]->detectors_data[k][detector_num];
            return k;
        }
    }
    return -1;

}

//-----------------------------------------------------------------------------

void MDA_IO::_load_detector_meta_data(data_struct::xrf::Detector * detector)
{
    detector->meta_array.resize(_mda_file->scan->number_detectors);

    for(int k=0; k<_mda_file->scan->number_detectors; k++)
    {
        detector->meta_array[k].name = std::string(_mda_file->scan->detectors[k]->name);
        detector->meta_array[k].description = std::string(_mda_file->scan->detectors[k]->description);
        detector->meta_array[k].unit = std::string(_mda_file->scan->detectors[k]->unit);
        detector->meta_array[k].number = _mda_file->scan->detectors[k]->number;
        detector->meta_array[k].value = (real_t)(_mda_file->scan->detectors_data[k][0]);
    }

}

//-----------------------------------------------------------------------------

bool MDA_IO::load_spectra_volume(std::string path,
                                 size_t detector_num,
                                 data_struct::xrf::Spectra_Volume* vol,
                                 bool hasNetCDF,
                                 std::unordered_map< std::string, std::string > *extra_override_values,
                                 data_struct::xrf::Quantification_Standard * quantification_standard)
{
    //index per row and col
    int elt_idx = -1;
    int ert_idx = -1;
    int incnt_idx = -1;
    int outcnt_idx = -1;


    //int sr_current_idx = -1;
    //int us_ic_idx = -1;
    //int ds_ic_idx = -1;
    int dpc1_ic_idx = -1;
    int dpc2_ic_idx = -1;

    int cfg1_idx = -1;
    int cfg2_idx = -1;
    int cfg3_idx = -1;
    int cfg4_idx = -1;
    int cfg5_idx = -1;
    int cfg6_idx = -1;
    int cfg7_idx = -1;
    int cfg8_idx = -1;
    int cfg9_idx = -1;


    bool single_row_scan = false;

    std::FILE *fptr = std::fopen(path.c_str(), "rb");

    size_t cols = 1;
    size_t rows = 1;
    size_t spectra = 1;

    if (fptr == nullptr)
    {
        return false;
    }

    struct mda_header *header = mda_header_load(fptr);

    if (header == nullptr)
    {
        return false;
    }

    std::cout<<"mda info ver:"<<header->version<<" data rank:"<<header->data_rank;
    //long total = 1;

    std::cout<<" cols "<< header->dimensions[0] << " rows " << header->dimensions[1] <<std::endl;


    _mda_file = mda_load(fptr);
    if (_mda_file == nullptr)
    {
        return false;
    }


    if (header->data_rank == 2)
    {
        if(hasNetCDF)
        {
            vol->resize(header->dimensions[0], header->dimensions[1], 2048);
            return true;
        }
        else
        {
            if(header->dimensions[1] == 2000)
            {
                rows = 1;
                cols = header->dimensions[0];
                spectra = header->dimensions[1];
                vol->resize(rows, cols, spectra);
                single_row_scan = true;
            }
            else
            {
                //if not then we don't know what is dataset is.
                 return false;
            }
        }
    }
    else if (header->data_rank == 3)
    {
        rows = header->dimensions[0];
        cols = header->dimensions[1];
        spectra = header->dimensions[2];
        vol->resize(header->dimensions[0], header->dimensions[1], header->dimensions[2]);
    }
    else
    {
        std::cout<<" Error: no support for data rank "<< header->data_rank <<std::endl;
        return false;
    }

    //_load_detector_meta_data(detector);

    if (extra_override_values != nullptr)
    {
        real_t tmp_val;
        if (extra_override_values->count("ELT1") > 0)
        {
            elt_idx = _find_2d_detector_index( extra_override_values->at("ELT1"), detector_num, tmp_val );
        }
        if (extra_override_values->count("ERT1") > 0)
        {
            ert_idx = _find_2d_detector_index( extra_override_values->at("ERT1"), detector_num, tmp_val );
        }
        if (extra_override_values->count("ICR1") > 0)
        {
            incnt_idx = _find_2d_detector_index( extra_override_values->at("ICR1"), detector_num, tmp_val );
        }
        if (extra_override_values->count("OCR1") > 0)
        {
            outcnt_idx = _find_2d_detector_index( extra_override_values->at("OCR1"), detector_num, tmp_val );
        }
        if(quantification_standard != nullptr)
        {
            if (extra_override_values->count("SRCURRENT") > 0)
            {
                _find_2d_detector_index( extra_override_values->at("SRCURRENT"), detector_num, tmp_val );
                quantification_standard->sr_current(tmp_val);
            }
            if (extra_override_values->count("US_IC") > 0)
            {
                _find_2d_detector_index( extra_override_values->at("US_IC"), detector_num, tmp_val );
                quantification_standard->US_IC(tmp_val);
            }
            if (extra_override_values->count("DS_IC") > 0)
            {
                _find_2d_detector_index( extra_override_values->at("DS_IC"), detector_num, tmp_val );
                quantification_standard->DS_IC(tmp_val);
            }
        }
/*
        if (extra_override_values->count("DPC1_IC") > 0)
        {
            dpc1_ic_idx = _find_2d_detector_index( extra_override_values->at("DPC1_IC") );
        }
        if (extra_override_values->count("DPC2_IC") > 0)
        {
            dpc2_ic_idx = _find_2d_detector_index( extra_override_values->at("DPC2_IC") );
        }

        if (extra_override_values->count("CFG_1") > 0)
        {
            cfg1_idx = _find_2d_detector_index( extra_override_values->at("CFG_1") );
        }
        if (extra_override_values->count("CFG_2") > 0)
        {
            cfg2_idx = _find_2d_detector_index( extra_override_values->at("CFG_2") );
        }
        if (extra_override_values->count("CFG_3") > 0)
        {
            cfg3_idx = _find_2d_detector_index( extra_override_values->at("CFG_3") );
        }
        if (extra_override_values->count("CFG_4") > 0)
        {
            cfg4_idx = _find_2d_detector_index( extra_override_values->at("CFG_4") );
        }
        if (extra_override_values->count("CFG_5") > 0)
        {
            cfg5_idx = _find_2d_detector_index( extra_override_values->at("CFG_5") );
        }
        if (extra_override_values->count("CFG_6") > 0)
        {
            cfg6_idx = _find_2d_detector_index( extra_override_values->at("CFG_6") );
        }
        if (extra_override_values->count("CFG_7") > 0)
        {
            cfg7_idx = _find_2d_detector_index( extra_override_values->at("CFG_7") );
        }
        if (extra_override_values->count("CFG_8") > 0)
        {
            cfg8_idx = _find_2d_detector_index( extra_override_values->at("CFG_8") );
        }
        if (extra_override_values->count("CFG_9") > 0)
        {
            cfg9_idx = _find_2d_detector_index( extra_override_values->at("CFG_9") );
        }
*/
    }
    std::cout<<" elt_idx "<< elt_idx << " ert_idx " << ert_idx << " in cnt idx " << incnt_idx << " out cnt idx "<< outcnt_idx<<std::endl;


    try
    {
        if( false == single_row_scan )
        {
            if(_mda_file->scan->last_point < _mda_file->scan->requested_points)
            {
                rows = _mda_file->scan->last_point + 1;
                //TODO: set a flag to return to tell that this is a bad scan
            }
        }

        for(size_t i=0; i<rows; i++)
        {
            // update num rows if header is incorrect and not single row scan

            if(_mda_file->scan->sub_scans[i]->last_point < _mda_file->scan->sub_scans[i]->requested_points)
            {
                cols = _mda_file->scan->sub_scans[i]->last_point + 1;
                //TODO: set a flag to return to tell that this is a bad scan
            }

            for(size_t j=0; j<cols; j++)
            {
/* TODO: we might need to do the same check for spectra
                if(_mda_file->scan->sub_scans[i]->last_point < _mda_file->scan->sub_scans[i]->requested_points)
                {
                    cols = _mda_file->scan->sub_scans[i]->last_point;
                }
*/
                if(elt_idx > -1)
                {
                    //std::cout<<"eltm ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[elt_idx][j]<< std::endl;
                    (*vol)[i][j].elapsed_lifetime(_mda_file->scan->sub_scans[i]->detectors_data[elt_idx][j]);
                }
                if(ert_idx > -1)
                {
                    //std::cout<<"elrm ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[ert_idx][j]<< std::endl;
                    (*vol)[i][j].elapsed_realtime(_mda_file->scan->sub_scans[i]->detectors_data[ert_idx][j]);
                }
                if(incnt_idx > -1)
                {
                    //std::cout<<"incnt ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[incnt_idx][j]<< std::endl;
                    (*vol)[i][j].input_counts(_mda_file->scan->sub_scans[i]->detectors_data[incnt_idx][j]);
                }
                if(outcnt_idx > -1)
                {
                    //std::cout<<"outcnt ["<<i<<"]["<<j<<"] = "<<_mda_file->scan->sub_scans[i]->detectors_data[outcnt_idx][j]<< std::endl;
                    (*vol)[i][j].output_counts(_mda_file->scan->sub_scans[i]->detectors_data[outcnt_idx][j]);
                }
                if(ert_idx > -1 && incnt_idx > -1 && outcnt_idx > -1)
                {
                    (*vol)[i][j].recalc_elapsed_lifetime();
                }

                if (single_row_scan)
                {
                    for(size_t k=0; k<spectra; k++)
                    {

                        (*vol)[i][j][k] = (_mda_file->scan->sub_scans[j]->detectors_data[detector_num][k]);
                    }
                }
                else
                {
                    for(size_t k=0; k<spectra; k++)
                    {
                        (*vol)[i][j][k] = (_mda_file->scan->sub_scans[i]->sub_scans[j]->detectors_data[detector_num][k]);
                    }
                }
            }
        }
        std::fclose(fptr);
    }
    catch(std::exception& e)
    {
        std::cout<<"!!! Error Caught exception loading mda file."<<std::endl;
        std::cerr << "Exception catched : " << e.what() << std::endl;
        return false;
    }


    return true;
}

//-----------------------------------------------------------------------------

bool load_henke_from_xdr(std::string filename, data_struct::xrf::Element_Info_Map *element_map)
{
    std::ifstream fileStream(filename);

    if (false == fileStream.good())
    {
        std::cout<<"Error opening file "<<filename<<std::endl;
        return false;
    }

    FILE* xdr_file = fopen(filename.c_str(), "rb");

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
        return false;
    if( xdr_int32_t( xdrstream, &num_energies) == 0)
        return false;

    element_map->_energies.resize(num_energies);
    //float *energy_arr = new float[num_energies];
    if( xdr_vector( xdrstream, (char *) &(element_map->_energies)[0], num_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
        return false;

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
            return false;
        if( xdr_vector( xdrstream, (char *) &(element->f2_atomic_scattering_imaginary)[0], num_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
            return false;
    }

    if( xdr_int32_t( xdrstream, &num_extra_energies) == 0)
        return false;

    for (int i=0; i<num_elements; i++)
    {
        data_struct::xrf::Element_Info* element = element_map->get_element(i+1);
        element->init_extra_energies(num_extra_energies);

        int element_n;
        if( xdr_int32_t( xdrstream, &element_n) == 0)
            return false;

        if( xdr_vector( xdrstream, (char *) &(element->extra_energies)[0], num_extra_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
            return false;
        if( xdr_vector( xdrstream, (char *) &(element->extra_f1)[0], num_extra_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
            return false;
        if( xdr_vector( xdrstream, (char *) &(element->extra_f2)[0], num_extra_energies, sizeof(float), (xdrproc_t) xdr_float) == false)
            return false;
    }

    fclose(xdr_file);

    return true;
}

//-----------------------------------------------------------------------------

} //end namespace file
}// end namespace io
