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

#ifndef HL_FILE_IO_H
#define HL_FILE_IO_H

#if defined _WIN32
#include "support/direct/dirent.h"
#else
#include <dirent.h>
#endif

#include <algorithm>

#include "io/file/netcdf_io.h"
#include "io/file/mda_io.h"
#include "io/file/mca_io.h"
#include "io/file/hdf5_io.h"
#include "io/file/csv_io.h"
#include "io/file/file_scan.h"
#include "io/file/esrf/edf_io.h"

#include "data_struct/spectra_volume.h"

#include "fitting/models/gaussian_model.h"

#include "data_struct/element_info.h"

#include "io/file/aps/aps_fit_params_import.h"

#include "fitting/routines/base_fit_routine.h"

#include "data_struct/fit_element_map.h"
#include "data_struct/params_override.h"
#include "data_struct/scaler_lookup.h"

#include "data_struct/quantification_standard.h"

#include "data_struct/stream_block.h"

#include "fitting/routines/roi_fit_routine.h"
#include "fitting/routines/svd_fit_routine.h"
#include "fitting/routines/nnls_fit_routine.h"

#include "fitting/models/gaussian_model.h"

#ifdef _BUILD_WITH_QT
  #include "visual/grapher.h"
#endif


namespace io
{

namespace file
{

DLL_EXPORT void check_and_create_dirs(std::string dataset_directory);

//DLL_EXPORT std::vector<std::string> find_all_dataset_files(std::string dataset_directory, std::string search_str);

DLL_EXPORT void generate_h5_averages(std::string dataset_directory, std::string dataset_file, const std::vector<size_t>& detector_num_arr);

DLL_EXPORT bool load_scalers_lookup(const std::string filename);

DLL_EXPORT bool load_quantification_standardinfo(std::string dataset_directory, std::string quantification_info_file, std::vector<Quantification_Standard<double>>& standard_element_weights);

//DLL_EXPORT void populate_netcdf_hdf5_files(std::string dataset_dir);

DLL_EXPORT void save_optimized_fit_params(std::string dataset_dir, std::string dataset_filename, int detector_num, std::string result, data_struct::Fit_Parameters<double>* fit_params, const data_struct::Spectra<double>* const spectra, const data_struct::Fit_Element_Map_Dict<double>* const elements_to_fit);

////DLL_EXPORT void save_roi_fit_params(std::string dataset_dir, std::string dataset_filename, int detector_num, string result, data_struct::Fit_Parameters<double>* fit_params, data_struct::Spectra<double>* spectra, data_struct::Fit_Element_Map_Dict<double>* elements_to_fit);

//DLL_EXPORT void sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string>* dataset_files);

// ----------------------------------------------------------------------------

template<typename T_real>
void cb_load_spectra_data_helper(size_t row, size_t col, size_t height, size_t width, size_t detector_num, data_struct::Spectra<T_real>* spectra, void* user_data)
{
    data_struct::Spectra<T_real>* integrated_spectra = nullptr;

    if (user_data != nullptr)
    {
        integrated_spectra = static_cast<data_struct::Spectra<T_real>*>(user_data);
    }

    if (integrated_spectra != nullptr && spectra != nullptr)
    {
        if (integrated_spectra->size() != spectra->size())
        {
            *integrated_spectra = *spectra;
        }
        else
        {
            integrated_spectra->add(*spectra);
        }
    }

    if (spectra != nullptr)
    {
        delete spectra;
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT fitting::routines::Base_Fit_Routine<T_real>* generate_fit_routine(data_struct::Fitting_Routines proc_type, fitting::optimizers::Optimizer<T_real>* optimizer)
{
    //Fitting routines
    fitting::routines::Base_Fit_Routine<T_real>* fit_routine = nullptr;
    switch (proc_type)
    {
    case data_struct::Fitting_Routines::GAUSS_TAILS:
        fit_routine = new fitting::routines::Param_Optimized_Fit_Routine<T_real>();
        ((fitting::routines::Param_Optimized_Fit_Routine<T_real>*)fit_routine)->set_optimizer(optimizer);
        break;
    case data_struct::Fitting_Routines::GAUSS_MATRIX:
        fit_routine = new fitting::routines::Matrix_Optimized_Fit_Routine<T_real>();
        ((fitting::routines::Matrix_Optimized_Fit_Routine<T_real>*)fit_routine)->set_optimizer(optimizer);
        break;
    case data_struct::Fitting_Routines::ROI:
        fit_routine = new fitting::routines::ROI_Fit_Routine<T_real>();
        break;
    case data_struct::Fitting_Routines::SVD:
        fit_routine = new fitting::routines::SVD_Fit_Routine<T_real>();
        break;
    case data_struct::Fitting_Routines::NNLS:
        fit_routine = new fitting::routines::NNLS_Fit_Routine<T_real>();
        break;
    default:
        break;
    }
    return fit_routine;
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT bool load_override_params(std::string dataset_directory,
    int detector_num,
    data_struct::Params_Override<T_real>* params_override,
    bool append_file_name = true)
{
    std::string det_num = "";
    std::string filename = dataset_directory;
    if (detector_num > -1)
        det_num = std::to_string(detector_num);


    if (append_file_name)
    {
        filename += "maps_fit_parameters_override.txt" + det_num;
    }

    if (false == io::file::aps::load_parameters_override(filename, params_override))
    {
        logE << "Loading fit param override file: " << filename << "\n";
        return false;
    }
    else
    {
        data_struct::Element_Info<T_real>* detector_element;
        if (params_override->detector_element.length() > 0)
        {
            // Get the element info class                                   // detector element as string "Si" or "Ge" usually

            detector_element = data_struct::Element_Info_Map<T_real>::inst()->get_element(params_override->detector_element);
        }
        else
        {
            //log error or warning
            logE << "No detector material defined in maps_fit_parameters_override.txt . Defaulting to Si" << "\n";
            detector_element = data_struct::Element_Info_Map<T_real>::inst()->get_element("Si");
        }

        if (params_override->elements_to_fit.count(STR_COMPTON_AMPLITUDE) == 0)
        {
            params_override->elements_to_fit.insert(std::pair<std::string, data_struct::Fit_Element_Map<T_real>*>(STR_COMPTON_AMPLITUDE, new data_struct::Fit_Element_Map<T_real>(STR_COMPTON_AMPLITUDE, nullptr)));
        }
        if (params_override->elements_to_fit.count(STR_COHERENT_SCT_AMPLITUDE) == 0)
        {
            params_override->elements_to_fit.insert(std::pair<std::string, data_struct::Fit_Element_Map<T_real>*>(STR_COHERENT_SCT_AMPLITUDE, new data_struct::Fit_Element_Map<T_real>(STR_COHERENT_SCT_AMPLITUDE, nullptr)));
        }

        logI << "Elements to fit:  ";
        //Update element ratios by detector element
        for (auto& itr : params_override->elements_to_fit)
        {
            itr.second->init_energy_ratio_for_detector_element(detector_element);
            logit_s << itr.first << " ";
        }
        logit_s << "\n";
    }

    return true;
}


// ----------------------------------------------------------------------------

/**
 * @brief init_analysis_job_detectors : Read in maps_fit_parameters_override.txt[0-3] and initialize data structres
 *                                      Override file contains information about which element to fit, updated branching
 *                                      conditions and other custom properties of the dataset.
 * @param analysis_job : data structure that holds information about the analysis to be perfomred.
 * @return True if successful
 */
template<typename T_real>
DLL_EXPORT bool init_analysis_job_detectors(data_struct::Analysis_Job<T_real>* analysis_job)
{

    //    _default_sub_struct.
    //    if( false == io::load_override_params(_dataset_directory, -1, override_params) )
    //    {
    //        return false;
    //    }

    analysis_job->detectors_meta_data.clear();



    //initialize models and fit routines for all detectors
    for (size_t detector_num : analysis_job->detector_num_arr)
    {
        if (analysis_job->detectors_meta_data.count(detector_num) < 1)
        {
            analysis_job->detectors_meta_data[detector_num] = data_struct::Detector<T_real>(detector_num);
        }
        data_struct::Detector<T_real>* detector = &analysis_job->detectors_meta_data[detector_num];

        if (detector->model == nullptr)
        {
            detector->model = new fitting::models::Gaussian_Model<T_real>();
        }
        data_struct::Params_Override<T_real>* override_params = &(detector->fit_params_override_dict);

        override_params->dataset_directory = analysis_job->dataset_directory;
        override_params->detector_num = detector_num;

        if (false == io::file::load_override_params(analysis_job->dataset_directory, detector_num, override_params))
        {
            if (false == io::file::load_override_params(analysis_job->dataset_directory, -1, override_params))
            {
                //last case, check current directory for override. This will be used for streaming
                if (false == io::file::load_override_params("./", detector_num, override_params))
                {
                    if (false == io::file::load_override_params("./", -1, override_params))
                    {
                        return false;
                    }
                }
            }
        }

        detector->update_from_fit_paramseters();

        if (override_params->elements_to_fit.size() < 1)
        {
            logE << "No elements to fit. Check  maps_fit_parameters_override.txt0 - 3 exist" << "\n";
            return false;
        }

        for (auto proc_type : analysis_job->fitting_routines)
        {
            //Fitting models
            detector->fit_routines[proc_type] = generate_fit_routine(proc_type, analysis_job->optimizer());

            //reset model fit parameters to defaults
            detector->model->reset_to_default_fit_params();
            //Update fit parameters by override values
            detector->model->update_fit_params_values(&(override_params->fit_params));

            //Fit_Element_Map_Dict *elements_to_fit = &(detector->fit_params_override_dict.elements_to_fit);
            //Initialize model
            //fit_routine->initialize(detector->model, elements_to_fit, energy_range);

        }
    }

    return true;
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT bool load_element_info(const std::string element_henke_filename, const std::string element_csv_filename)
{
    io::file::MDA_IO<T_real> mda_io;
    if (mda_io.load_henke_from_xdr(element_henke_filename) == false)
    {
        logE << "Could not load " << element_henke_filename << "\n";
        return false;
    }

    if (io::file::csv::load_element_info<T_real>(element_csv_filename) == false)
    {
        logE << "Could not load " << element_csv_filename << "\n";
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT bool load_and_integrate_spectra_volume(std::string dataset_directory,
													std::string dataset_file,
													size_t detector_num,
													data_struct::Spectra<T_real>* integrated_spectra,
													data_struct::Params_Override<T_real>* params_override)
{
    //Dataset importer
    io::file::MDA_IO<T_real> mda_io;
    //data_struct::Detector detector;
    std::string tmp_dataset_file = dataset_file;
    bool ret_val = true;
    std::vector<size_t> detector_num_arr{ detector_num };
    size_t out_rows = 0;
    size_t out_cols = 0;
    data_struct::IO_Callback_Func_Def<T_real>  cb_function = std::bind(&cb_load_spectra_data_helper<T_real>, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5, std::placeholders::_6, std::placeholders::_7);

    if (dataset_directory.back() != DIR_END_CHAR)
    {
        dataset_directory += DIR_END_CHAR;
    }
    //replace / with \ for windows, won't do anything for linux
    std::replace(dataset_directory.begin(), dataset_directory.end(), '/', DIR_END_CHAR);

    data_struct::Spectra_Volume<T_real> spectra_volume;

    logI << "Loading dataset " << dataset_directory + "mda" + DIR_END_CHAR + dataset_file << "\n";

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size() - 4);
    bool hasNetcdf = false;
    bool hasBnpNetcdf = false;
    bool hasHdf = false;
    bool hasXspress = false;
    std::string file_middle = ""; //_2xfm3_ or dxpM...
    std::string bnp_netcdf_base_name = "bnp_fly_";
    for (auto& itr : io::file::File_Scan::inst()->netcdf_files())
    {
        if (itr.find(tmp_dataset_file) == 0)
        {
            size_t slen = (itr.length() - 4) - tmp_dataset_file.length();
            file_middle = itr.substr(tmp_dataset_file.length(), slen);
            hasNetcdf = true;
            break;
        }
    }
    if (hasNetcdf == false)
    {
        int idx = static_cast<int>(tmp_dataset_file.find("bnp_fly"));
        if (idx == 0)
        {
            std::string footer = tmp_dataset_file.substr(7, tmp_dataset_file.length() - 7);
            int file_index = std::atoi(footer.c_str());
            file_middle = std::to_string(file_index);
            bnp_netcdf_base_name = "bnp_fly_" + file_middle + "_";
            for (auto& itr : io::file::File_Scan::inst()->bnp_netcdf_files())
            {
                if (itr.find(bnp_netcdf_base_name) == 0)
                {
                    hasBnpNetcdf = true;
                    break;
                }
            }
        }
    }
    if (hasNetcdf == false && hasBnpNetcdf == false)
    {
        for (auto& itr : io::file::File_Scan::inst()->hdf_files())
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length() - 4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasHdf = true;
                break;
            }
        }
    }
    if (hasNetcdf == false && hasBnpNetcdf == false && hasHdf == false)
    {
        for (auto& itr : io::file::File_Scan::inst()->hdf_xspress_files())
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length() - 4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasXspress = true;
                break;
            }
        }
    }

    if (hasNetcdf == false && hasBnpNetcdf == false && hasHdf == false && hasXspress == false)
    {

        int idx = static_cast<int>(tmp_dataset_file.find("bnp_fly"));
        if (idx == 0)
        {
            std::string footer = tmp_dataset_file.substr(7, tmp_dataset_file.length() - 7);
            int file_index = std::atoi(footer.c_str());
            file_middle = std::to_string(file_index);
            bnp_netcdf_base_name = "bnp_fly_" + file_middle + "_";
            for (auto& itr : io::file::File_Scan::inst()->hdf_xspress_files())
            {
                if (itr.find(bnp_netcdf_base_name) == 0)
                {
                    hasXspress = true;
                    break;
                }
            }
        }
    }


    bool ends_in_mca = false;
    size_t dlen = dataset_file.length();
    if (dataset_file[dlen - 4] == '.' && dataset_file[dlen - 3] == 'm' && dataset_file[dlen - 2] == 'c' && dataset_file[dlen - 1] == 'a')
    {
        ends_in_mca = true;
    }
    else
    {
        for (int i = 0; i < 8; i++)
        {
            std::string s1 = std::to_string(i);
            if (dataset_file[dlen - 5] == '.' && dataset_file[dlen - 4] == 'm' && dataset_file[dlen - 3] == 'c' && dataset_file[dlen - 2] == 'a' && dataset_file[dlen - 1] == s1[0])
            {
                ends_in_mca = true;
                break;
            }
        }
    }

    if (ends_in_mca)
    {
        Spectra<T_real> spec;
        std::unordered_map<std::string, T_real> pv_map;
        if (true == io::file::mca::load_integrated_spectra(dataset_directory + "mda" + DIR_END_CHAR + dataset_file, integrated_spectra, pv_map))
        {
            if (pv_map.count(STR_SR_CURRENT) > 0)
            {
                params_override->sr_current = pv_map.at(STR_SR_CURRENT);
            }
            if (pv_map.count(STR_US_IC) > 0)
            {
                params_override->US_IC = pv_map.at(STR_US_IC);
            }
            if (pv_map.count(STR_DS_IC) > 0)
            {
                params_override->DS_IC = pv_map.at(STR_DS_IC);
            }    
            return true;
        }
    }

    //  try to load from a pre analyzed file because they should contain the integrated spectra
    std::string fullpath = dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file + ".h5";
    if (detector_num != -1)
    {
        fullpath += std::to_string(detector_num);
    }
    bool loaded0 = io::file::HDF5_IO::inst()->load_integrated_spectra_analyzed_h5(fullpath, integrated_spectra, nullptr, false);
    if (false == loaded0)
    {
        // if failed, try just the name in img.dat folder
        fullpath = dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file;
        loaded0 = io::file::HDF5_IO::inst()->load_integrated_spectra_analyzed_h5(fullpath, integrated_spectra, nullptr, false);
    }
    if (true == loaded0)
    {
        logI << "Loaded integradted spectra from h5.\n";
        if (params_override != nullptr)
        {
            if (false == io::file::HDF5_IO::inst()->load_quantification_scalers_analyzed_h5(fullpath, params_override))
            {
                mda_io.load_quantification_scalers(dataset_directory + "mda" + DIR_END_CHAR + dataset_file, params_override);
            }
        }
        return true;
    }

    //try loading emd dataset
    if (io::file::HDF5_IO::inst()->load_spectra_volume_emd_with_callback(dataset_directory + DIR_END_CHAR + dataset_file, detector_num_arr, cb_function, integrated_spectra))
        //if (true == io::file::HDF5_IO::inst()->load_spectra_volume_emd(dataset_directory + DIR_END_CHAR + dataset_file, detector_num, &spectra_volume, false))
    {
        logI << "Loaded spectra volume confocal from h5.\n";
        return true;
    }

    //try loading confocal dataset
    if (true == io::file::HDF5_IO::inst()->load_spectra_volume_confocal(dataset_directory + DIR_END_CHAR + dataset_file, detector_num, &spectra_volume, false))
    {
        logI << "Loaded spectra volume confocal from h5.\n";
        *integrated_spectra = spectra_volume.integrate();
        return true;
    }

    //try loading gse cars dataset
    if (true == io::file::HDF5_IO::inst()->load_spectra_volume_gsecars(dataset_directory + DIR_END_CHAR + dataset_file, detector_num, &spectra_volume, false))
    {
        fullpath = dataset_directory + DIR_END_CHAR + dataset_file;
        if (false == io::file::HDF5_IO::inst()->load_quantification_scalers_gsecars(fullpath, params_override))
        {
            logW << "Failed to load ION chamber scalers from h5.\n";
        }

        logI << "Loaded spectra volume gse cars from h5.\n";
        *integrated_spectra = spectra_volume.integrate();
        return true;
    }

    if (true == io::file::HDF5_IO::inst()->load_integrated_spectra_bnl(dataset_directory + DIR_END_CHAR + dataset_file, detector_num, integrated_spectra, false))
    {
        fullpath = dataset_directory + DIR_END_CHAR + dataset_file;
        if (false == io::file::HDF5_IO::inst()->load_quantification_scalers_BNL(fullpath, params_override))
        {
            logW << "Failed to load ION chamber scalers from h5.\n";
        }
        return true;
    }


    //load spectra
    // load_spectra_volume will alloc memory for the whole vol, we don't want that for integrated spec
    bool has_external_files = hasNetcdf | hasBnpNetcdf | hasHdf | hasXspress;
    //if(false == mda_io.load_spectra_volume_with_callback(dataset_directory + "mda" + DIR_END_CHAR + dataset_file, detector_num_arr, has_external_files, analysis_job, out_rows, out_cols, cb_function, integrated_spectra))
    if (false == mda_io.load_integrated_spectra(dataset_directory + "mda" + DIR_END_CHAR + dataset_file, detector_num, integrated_spectra, has_external_files))

    {
        logE << "Load spectra " << dataset_directory + "mda" + DIR_END_CHAR + dataset_file << "\n";
        return false;
    }
    else
    {
        mda_io.load_quantification_scalers(dataset_directory + "mda" + DIR_END_CHAR + dataset_file, params_override);

        if (false == hasNetcdf && false == hasBnpNetcdf && false == hasHdf)
        {
            mda_io.unload();
        }
        else
        {
            int rank;
            size_t dims[10];
            dims[0] = 0;
            rank = io::file::mda_get_rank_and_dims(dataset_directory + "mda" + DIR_END_CHAR + dataset_file, &dims[0]);
            if (rank == 3)
            {
                integrated_spectra->resize(dims[2]);
                integrated_spectra->setZero(dims[2]);
            }
            else
            {
                integrated_spectra->resize(2048);
                integrated_spectra->setZero(2048);
            }


            if (hasNetcdf)
            {
                std::ifstream file_io(dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc");
                if (file_io.is_open())
                {
                    file_io.close();
                    std::string full_filename;
                    for (size_t i = 0; i < dims[0]; i++)
                    {
                        full_filename = dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + std::to_string(i) + ".nc";
                        //logI<<"Loading file "<<full_filename<<"\n";
                        size_t spec_size = io::file::NetCDF_IO<T_real>::inst()->load_spectra_line_integrated(full_filename, detector_num, dims[1], integrated_spectra);
                        if (detector_num > 3 && spec_size == -1) // this netcdf file only has 4 element detectors
                        {
                            return false;
                        }
                    }
                }
                else
                {
                    logW << "Did not find netcdf files " << dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc" << "\n";
                    //return false;
                }
            }
            else if (hasBnpNetcdf)
            {
                std::ifstream file_io(dataset_directory + "flyXRF" + DIR_END_CHAR + bnp_netcdf_base_name + "001.nc");
                if (file_io.is_open())
                {
                    file_io.close();
                    std::string full_filename;
                    for (size_t i = 0; i < dims[0]; i++)
                    {
                        std::string row_idx_str = std::to_string(i + 1);
                        int num_prepended_zeros = 3 - static_cast<int>(row_idx_str.size()); // 3 chars for num of rows, prepened with zeros if less than 100
                        std::string row_idx_str_full = "";
                        for (int z = 0; z < num_prepended_zeros; z++)
                        {
                            row_idx_str_full += "0";
                        }
                        row_idx_str_full += row_idx_str;
                        full_filename = dataset_directory + "flyXRF" + DIR_END_CHAR + bnp_netcdf_base_name + row_idx_str_full + ".nc";
                        size_t spec_size = io::file::NetCDF_IO<T_real>::inst()->load_spectra_line_integrated(full_filename, detector_num, dims[1], integrated_spectra);
                        if (detector_num > 3 && spec_size == -1) // this netcdf file only has 4 element detectors
                        {
                            return false;
                        }
                    }
                }
                else
                {
                    logW << "Did not find netcdf files " << dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc" << "\n";
                    //return false;
                }
            }
            else if (hasHdf)
            {
                ret_val = io::file::HDF5_IO::inst()->load_and_integrate_spectra_volume(dataset_directory + "flyXRF.h5" + DIR_END_CHAR + tmp_dataset_file + file_middle + "0.h5", detector_num, integrated_spectra);
            }
            else if (hasXspress)
            {
                std::string full_filename;
                data_struct::Spectra_Line<T_real> spectra_line;
                spectra_line.resize_and_zero(dims[1], integrated_spectra->size());
                for (size_t i = 0; i < dims[0]; i++)
                {
                    full_filename = dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + std::to_string(i) + ".hdf5";
                    if (io::file::HDF5_IO::inst()->load_spectra_line_xspress3(full_filename, detector_num, &spectra_line))
                    {
                        for (size_t k = 0; k < spectra_line.size(); k++)
                        {
                            *integrated_spectra += spectra_line[k];
                        }
                    }
                }
            }
        }
    }

    logI << "Finished Loading dataset " << dataset_directory + "mda" + DIR_END_CHAR + dataset_file << " detector " << detector_num << "\n";
    return ret_val;
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT bool load_spectra_volume(std::string dataset_directory,
                         std::string dataset_file,
                         size_t detector_num,
                         data_struct::Spectra_Volume<T_real>* spectra_volume,
                         data_struct::Params_Override<T_real>* params_override,
                         bool *is_loaded_from_analyazed_h5,
                         bool save_scalers)
{

    //Dataset importer
    io::file::MDA_IO<T_real> mda_io;
    //data_struct::Detector detector;
    std::string tmp_dataset_file = dataset_file;
    if (detector_num == -1)
    {
        logI << "Loading dataset " << dataset_directory << dataset_file << "\n";
    }
    else
    {
        logI << "Loading dataset " << dataset_directory << "mda" << DIR_END_CHAR << dataset_file << " detector " << detector_num << "\n";
    }
    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size() - 4);
    bool hasNetcdf = false;
    bool hasBnpNetcdf = false;
    bool hasHdf = false;
    bool hasXspress = false;
    std::string file_middle = ""; //_2xfm3_, dxpM, or file index in case of bnp...
    std::string bnp_netcdf_base_name = "bnp_fly_";
    std::vector<int> bad_rows;
    for (auto& itr : io::file::File_Scan::inst()->netcdf_files())
    {
        if (itr.find(tmp_dataset_file) == 0)
        {
            size_t slen = (itr.length() - 4) - tmp_dataset_file.length();
            file_middle = itr.substr(tmp_dataset_file.length(), slen);
            hasNetcdf = true;
            break;
        }
    }
    if (hasNetcdf == false)
    {
        int idx = static_cast<int>(tmp_dataset_file.find("bnp_fly"));
        if (idx == 0)
        {
            std::string footer = tmp_dataset_file.substr(7, tmp_dataset_file.length() - 7);
            int file_index = std::atoi(footer.c_str());
            file_middle = std::to_string(file_index);
            bnp_netcdf_base_name = "bnp_fly_" + file_middle + "_";
            for (auto& itr : io::file::File_Scan::inst()->bnp_netcdf_files())
            {
                if (itr.find(bnp_netcdf_base_name) == 0)
                {
                    hasBnpNetcdf = true;
                    break;
                }
            }
        }
    }
    if (hasNetcdf == false && hasBnpNetcdf == false)
    {
        for (auto& itr : io::file::File_Scan::inst()->hdf_files())
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length() - 4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasHdf = true;
                break;
            }
        }
    }
    if (hasNetcdf == false && hasBnpNetcdf == false && hasHdf == false)
    {
        for (auto& itr : io::file::File_Scan::inst()->hdf_xspress_files())
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length() - 6) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasXspress = true;
                break;
            }
        }
    }

    if (hasNetcdf == false && hasBnpNetcdf == false && hasHdf == false && hasXspress == false)
    {

        int idx = static_cast<int>(tmp_dataset_file.find("bnp_fly"));
        if (idx == 0)
        {
            std::string footer = tmp_dataset_file.substr(7, tmp_dataset_file.length() - 7);
            int file_index = std::atoi(footer.c_str());
            file_middle = std::to_string(file_index);
            bnp_netcdf_base_name = "bnp_fly_" + file_middle + "_";
            for (auto& itr : io::file::File_Scan::inst()->hdf_xspress_files())
            {
                if (itr.find(bnp_netcdf_base_name) == 0)
                {
                    hasXspress = true;
                    break;
                }
            }
        }
    }

    bool ends_in_h5 = false;
    bool ends_in_mca = false;
    size_t dlen = dataset_file.length();
    if (dataset_file[dlen - 3] == '.' && dataset_file[dlen - 2] == 'h' && dataset_file[dlen - 1] == '5')
    {
        ends_in_h5 = true;
    }
    else if (dataset_file[dlen - 4] == '.' && dataset_file[dlen - 3] == 'm' && dataset_file[dlen - 2] == 'c' && dataset_file[dlen - 1] == 'a')
    {
        ends_in_mca = true;
    }
    else
    {
        for (int i = 0; i < 8; i++)
        {
            std::string s1 = std::to_string(i);
            if (dataset_file[dlen - 4] == '.' && dataset_file[dlen - 3] == 'h' && dataset_file[dlen - 2] == '5' && dataset_file[dlen - 1] == s1[0])
            {
                ends_in_h5 = true;
                break;
            }
            else if (dataset_file[dlen - 5] == '.' && dataset_file[dlen - 4] == 'm' && dataset_file[dlen - 3] == 'c' && dataset_file[dlen - 2] == 'a' && dataset_file[dlen - 1] == s1[0])
            {
                ends_in_mca = true;
                break;
            }
        }
    }

    if (ends_in_mca)
    {
        Spectra<T_real> spec;
        std::unordered_map<std::string, T_real> pv_map;
        if (true == io::file::mca::load_integrated_spectra(dataset_directory + "mda" + DIR_END_CHAR + dataset_file, &spec, pv_map))
        {
            data_struct::Scan_Info<T_real> scan_info;

            for (auto& itr : pv_map)
            {
                data_struct::Scaler_Map<T_real> sm;
                sm.name = itr.first;
                sm.unit = " ";
                sm.time_normalized = false;
                sm.values.resize(1,1);
                sm.values(0, 0) = itr.second;
                scan_info.scaler_maps.push_back(sm);
            }
            
            scan_info.meta_info.requested_rows = 1;
            scan_info.meta_info.requested_cols = 1;
            scan_info.meta_info.x_axis.resize(1);
            scan_info.meta_info.x_axis.setZero(1);
            scan_info.meta_info.y_axis.resize(1);
            scan_info.meta_info.y_axis.setZero(1);
            scan_info.meta_info.theta = 0.0;

            data_struct::Extra_PV ep;
            ep.name = "Name";
            ep.description = "Filename";
            ep.unit = "Str";
            ep.value = dataset_file;

            scan_info.extra_pvs.push_back(ep);

            spectra_volume->resize_and_zero(1, 1, spec.size());
            (*spectra_volume)[0][0] = spec;
            io::file::HDF5_IO::inst()->start_save_seq(true);

            // add ELT, ERT, INCNT, OUTCNT to scaler map
            spectra_volume->generate_scaler_maps(&(scan_info.scaler_maps));
            io::file::HDF5_IO::inst()->save_scan_scalers(detector_num, &scan_info, params_override);
            return true;
        }
    }

    std::string fullpath = dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file;
    if (false == ends_in_h5)
    {
        fullpath += ".h5";

        if (detector_num != -1)
        {
            fullpath += std::to_string(detector_num);
        }
    }
    /*
    std::string fullpath;
    size_t dlen = dataset_file.length();
    if (dataset_file[dlen -4] == '.' && dataset_file[dlen - 3] == 'm' && dataset_file[dlen - 2] == 'd' && dataset_file[dlen - 1] == 'a')
    {
        fullpath = dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file + ".h5" + std::to_string(detector_num);
    }
    else
    {
        fullpath = dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file;
    }
    */
    //  try to load from a pre analyzed file because they should contain the whole mca_arr spectra volume
    if (true == io::file::HDF5_IO::inst()->load_spectra_vol_analyzed_h5(fullpath, spectra_volume))
    {
        logI << "Loaded spectra volume from h5.\n";
        *is_loaded_from_analyazed_h5 = true;
        return io::file::HDF5_IO::inst()->start_save_seq(false);
    }
    else
    {
        *is_loaded_from_analyazed_h5 = false;
    }

    // try to load esrf hdf5
    if (ends_in_h5)
    {
        fullpath = dataset_directory + DIR_END_CHAR + "mda" + DIR_END_CHAR + dataset_file;
        std::string file_title;
        if(true == io::file::HDF5_IO::inst()->load_spectra_vol_esrf(fullpath, file_title, spectra_volume))
        {
            logI << "Loaded spectra volume esrf.\n";
            std::string str_det_num = "";
            if (detector_num < 10)
            {
                // prepend 0 
                str_det_num = "0" + std::to_string(detector_num);
            }
            else
            {
                str_det_num = std::to_string(detector_num);
            }

            //  COX_4_50x_400nm_05_xia00_0001_0000_0000.edf
            //    title                 detector        row
            //   COX_4_50x_400nm_05      xia00          0000
           
            for (int r = 0; r < spectra_volume->rows(); r++)
            {
                std::string str_r = std::to_string(r);
                std::string str_row = std::string(4 - str_r.length(), '0') + std::to_string(r);
                fullpath = dataset_directory + DIR_END_CHAR + "edf" + DIR_END_CHAR + file_title + "_xia" + str_det_num + "_0001_0000_" + str_row + ".edf";
                io::file::edf::load_spectra_line(fullpath, &(*spectra_volume)[r]);
            }
            io::file::HDF5_IO::inst()->start_save_seq(true);
            //// io::file::HDF5_IO::inst()->save_scan_scalers_esrf<T_real>(dataset_directory + DIR_END_CHAR + dataset_file, detector_num);
            return true;
        }
    }

    //try loading emd dataset if it ends in .emd
    if (dataset_file.rfind(".emd") == dataset_file.length() - 4)
    {
        if (true == io::file::HDF5_IO::inst()->load_spectra_volume_emd(dataset_directory + DIR_END_CHAR + dataset_file, detector_num, spectra_volume))
        {
            //*is_loaded_from_analyazed_h5 = true;//test to not save volume
            std::string str_detector_num = "";
            if (detector_num != -1)
            {
                str_detector_num = std::to_string(detector_num);
            }
            std::string full_save_path = dataset_directory + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + dataset_file + "_frame_" + str_detector_num + ".h5";
            io::file::HDF5_IO::inst()->start_save_seq(full_save_path, true);
            return true;
        }
    }

    //try loading confocal dataset
    if (true == io::file::HDF5_IO::inst()->load_spectra_volume_confocal(dataset_directory + DIR_END_CHAR + dataset_file, detector_num, spectra_volume, false))
    {
        if (save_scalers)
        {
            io::file::HDF5_IO::inst()->start_save_seq(true);
            io::file::HDF5_IO::inst()->save_scan_scalers_confocal<T_real>(dataset_directory + DIR_END_CHAR + dataset_file, detector_num);
        }
        return true;
    }

    //try loading gse cars dataset
    if (true == io::file::HDF5_IO::inst()->load_spectra_volume_gsecars<T_real>(dataset_directory + DIR_END_CHAR + dataset_file, detector_num, spectra_volume, false))
    {
        if (save_scalers)
        {
            io::file::HDF5_IO::inst()->start_save_seq(true);
            io::file::HDF5_IO::inst()->save_scan_scalers_gsecars<T_real>(dataset_directory + DIR_END_CHAR + dataset_file, detector_num);
        }
        return true;
    }

    if (true == io::file::HDF5_IO::inst()->load_spectra_volume_bnl<T_real>(dataset_directory + DIR_END_CHAR + dataset_file, detector_num, spectra_volume, false))
    {
        if (save_scalers)
        {
            io::file::HDF5_IO::inst()->start_save_seq(true);
            io::file::HDF5_IO::inst()->save_scan_scalers_bnl<T_real>(dataset_directory + DIR_END_CHAR + dataset_file, detector_num);
        }
        return true;
    }

    // try to load spectra from mda file
    if (false == mda_io.load_spectra_volume(dataset_directory + "mda" + DIR_END_CHAR + dataset_file, detector_num, spectra_volume, hasNetcdf | hasBnpNetcdf | hasHdf | hasXspress))
    {
        logE << "Load spectra " << dataset_directory + "mda" + DIR_END_CHAR + dataset_file << "\n";
        return false;
    }
    else
    {
        if (hasNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc");
            if (file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for (size_t i = 0; i < spectra_volume->rows(); i++)
                {
                    full_filename = dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + std::to_string(i) + ".nc";
                    //todo: add verbose option
                    //logI<<"Loading file "<<full_filename<<"\n";
                    size_t spec_size = io::file::NetCDF_IO<T_real>::inst()->load_spectra_line(full_filename, detector_num, &(*spectra_volume)[i]);
                    if (detector_num > 0 && spec_size == -1) // this netcdf file only has 1 element detectors
                    {
                        return false;
                    }
                }
            }
            else
            {
                logW << "Did not find netcdf files " << dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc" << "\n";
                //return false;
            }
        }
        else if (hasBnpNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF" + DIR_END_CHAR + bnp_netcdf_base_name + "001.nc");
            if (file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for (size_t i = 0; i < spectra_volume->rows(); i++)
                {
                    std::string row_idx_str = std::to_string(i + 1);
                    int num_prepended_zeros = 3 - static_cast<int>(row_idx_str.size()); // 3 chars for num of rows, prepened with zeros if less than 100
                    std::string row_idx_str_full = "";
                    for (int z = 0; z < num_prepended_zeros; z++)
                    {
                        row_idx_str_full += "0";
                    }
                    row_idx_str_full += row_idx_str;
                    full_filename = dataset_directory + "flyXRF" + DIR_END_CHAR + bnp_netcdf_base_name + row_idx_str_full + ".nc";
                    size_t prev_size = 0;
                    size_t spec_size = io::file::NetCDF_IO<T_real>::inst()->load_spectra_line(full_filename, detector_num, &(*spectra_volume)[i]);
                    //
                    if (detector_num > 3 && spec_size == -1) // this netcdf file only has 4 element detectors
                    {
                        return false;
                    }
                    //if we failed to load and it isn't the first row, copy the previous one
                    if (i > 0)
                    {
                        prev_size = (*spectra_volume)[i - 1].size();
                    }
                    if (spec_size == 0 || spec_size < prev_size)
                    {
                        if (i > 0)
                        {
                            logW << "Bad row for file " << full_filename << " row " << i << ", using previous line\n";
                            bad_rows.push_back(i);
                            (*spectra_volume)[i] = (*spectra_volume)[i - 1];
                        }
                    }
                }
            }
            else
            {
                logW << "Did not find netcdf files " << dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + "0.nc" << "\n";
                //return false;
            }
        }
        else if (hasHdf)
        {
            if (false == io::file::HDF5_IO::inst()->load_spectra_volume(dataset_directory + "flyXRF.h5" + DIR_END_CHAR + tmp_dataset_file + file_middle + "0.h5", detector_num, spectra_volume))
            {
                std::string full_filename;
                for (size_t i = 0; i < spectra_volume->rows(); i++)
                {
                    //everyone else
                    full_filename = dataset_directory + "flyXRF.h5" + DIR_END_CHAR + tmp_dataset_file + file_middle + std::to_string(i) + ".h5";
                    io::file::HDF5_IO::inst()->load_spectra_line_xspress3(full_filename, detector_num, &(*spectra_volume)[i]);
                }

            }
        }
        else if (hasXspress)
        {
            std::string full_filename;
            for (size_t i = 0; i < spectra_volume->rows(); i++)
            //for (size_t i = 1; i <= spectra_volume->rows(); i++) //BNP hack of starting at 1 instead of 0
            {
                //bnp format
                //full_filename = dataset_directory + "flyXRF" + DIR_END_CHAR + bnp_netcdf_base_name + std::to_string(i) + ".hdf5";   
                //io::file::HDF5_IO::inst()->load_spectra_line_xspress3(full_filename, detector_num, &(*spectra_volume)[i - 1]);

                //everyone else
                full_filename = dataset_directory + "flyXRF" + DIR_END_CHAR + tmp_dataset_file + file_middle + std::to_string(i) + ".hdf5";
                io::file::HDF5_IO::inst()->load_spectra_line_xspress3(full_filename, detector_num, &(*spectra_volume)[i]);
            }
        }

    }

    if (save_scalers)
    {
        io::file::HDF5_IO::inst()->start_save_seq(true);
        data_struct::Scan_Info<T_real>* scan_info = mda_io.get_scan_info();
        // add ELT, ERT, INCNT, OUTCNT to scaler map
        if (spectra_volume != nullptr && scan_info != nullptr)
        {
            spectra_volume->generate_scaler_maps(&(scan_info->scaler_maps));
        }

        for (const auto& line : bad_rows)
        {
            for (auto& map : scan_info->scaler_maps)
            {
                // copy prev row
                for (Eigen::Index col = 0; col < map.values.cols(); col++)
                {
                    map.values(line, col) = map.values(line - 1, col);
                }
            }
        }
        io::file::HDF5_IO::inst()->save_scan_scalers(detector_num, scan_info, params_override);
    }

    mda_io.unload();
    logI << "Finished Loading dataset " << dataset_directory + "mda" + DIR_END_CHAR + dataset_file << " detector " << detector_num << "\n";
    return true;
}

// ----------------------------------------------------------------------------

// This is for HDF5 files only
template<typename T_real>
DLL_EXPORT bool get_scalers_and_metadata_h5(std::string dataset_directory, std::string dataset_file, data_struct::Scan_Info<T_real>* scan_info)
{
    if (true == io::file::HDF5_IO::inst()->get_scalers_and_metadata_emd(dataset_directory + DIR_END_CHAR + dataset_file, scan_info))
    {
        return true;
    }

    if (true == io::file::HDF5_IO::inst()->get_scalers_and_metadata_confocal(dataset_directory + DIR_END_CHAR + dataset_file, scan_info))
    {
        return true;
    }

    if (true == io::file::HDF5_IO::inst()->get_scalers_and_metadata_gsecars(dataset_directory + DIR_END_CHAR + dataset_file, scan_info))
    {
        return true;
    }

    if (true == io::file::HDF5_IO::inst()->get_scalers_and_metadata_bnl(dataset_directory + DIR_END_CHAR + dataset_file, scan_info))
    {
        return true;
    }

    logE << "Failed to load any meta info!\n";
    return false;
}

template<typename T_real>
DLL_EXPORT void save_quantification_plots(std::string path, Detector<T_real>* detector)
{
    std::string str_path = path + "/output/";
    //if(analysis_job->get_detector(detector_num))
    file::csv::save_quantification(str_path, detector);
#ifdef _BUILD_WITH_QT
    visual::SavePlotQuantificationFromConsole(str_path, detector);
#endif
}

}// end namespace file
}// end namespace io

#endif // HL_FILE_IO_H
