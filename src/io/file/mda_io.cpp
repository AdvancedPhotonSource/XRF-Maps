/***

Copyright (c) 2016 Arthur Glowacki

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

***/

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

MDA_IO::MDA_IO() : Base_File_IO()
{

    _mda_file = nullptr;
    _mda_file_info = nullptr;

}

MDA_IO::~MDA_IO()
{

    unload();

}

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

int MDA_IO::_find_2d_detector_index(std::string det_name)
{

    for(int k=0; k<_mda_file->scan->sub_scans[0]->number_detectors; k++)
    {
        //std::cout<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->name << std::endl;
        //std::cout<<"det name "<<_mda_file->scan->sub_scans[0]->detectors[k]->description << std::endl;
        if(strcmp(_mda_file->scan->sub_scans[0]->detectors[k]->name, det_name.c_str())  == 0)
        {
            return k;
        }
    }
    return -1;

}

bool MDA_IO::load_spectra_volume(std::string path,
                                 size_t detector_num,
                                 data_struct::xrf::Detector* detector,
                                 data_struct::xrf::Spectra_Volume* vol,
                                 std::unordered_map< std::string, std::string > *extra_override_values)
{
    int elt_idx = -1;
    int ert_idx = -1;
    int incnt_idx = -1;
    int outcnt_idx = -1;

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


    if (header->data_rank == 2)
    {
        vol->resize(header->dimensions[0], header->dimensions[1], 2048);
        return true;
    }
    else if (header->data_rank == 3)
    {
        cols = header->dimensions[0];
        rows = header->dimensions[1];
        spectra = header->dimensions[2];
        vol->resize(header->dimensions[0], header->dimensions[1], header->dimensions[2]);
    }
    else
    {
        std::cout<<" Error: no support for data rank "<< header->data_rank <<std::endl;
        return false;
    }

    _mda_file = mda_load(fptr);
    if (_mda_file == nullptr)
    {
        return false;
    }

    if (extra_override_values != nullptr)
    {
        if (extra_override_values->count("ELT1") > 0)
        {
            elt_idx = _find_2d_detector_index( extra_override_values->at("ELT1") );
        }
        if (extra_override_values->count("ERT1") > 0)
        {
            ert_idx = _find_2d_detector_index( extra_override_values->at("ERT1") );
        }
        if (extra_override_values->count("ICR1") > 0)
        {
            incnt_idx = _find_2d_detector_index( extra_override_values->at("ICR1") );
        }
        if (extra_override_values->count("OCR1") > 0)
        {
            outcnt_idx = _find_2d_detector_index( extra_override_values->at("OCR1") );
        }

    }
    std::cout<<" elt_idx "<< elt_idx << " ert_idx " << ert_idx << " in cnt idx " << incnt_idx << " out cnt idx "<< outcnt_idx<<std::endl;

    for(size_t i=0; i<cols; i++)
    {
        for(size_t j=0; j<rows; j++)
        {

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

            for(size_t k=0; k<spectra; k++)
            {
                (*vol)[i][j][k] = (_mda_file->scan->sub_scans[i]->sub_scans[j]->detectors_data[detector_num][k]);
            }
        }
    }

    std::fclose(fptr);

    return true;
}

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

} //end namespace file
}// end namespace io
