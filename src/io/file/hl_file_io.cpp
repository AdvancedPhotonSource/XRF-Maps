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

#include "hl_file_io.h"

#include <iostream>
#include <sstream>
#include <fstream>

namespace io
{

std::vector<std::string> netcdf_files;
std::vector<std::string> hdf_files;
std::vector<std::string> hdf_xspress_files;
std::vector<std::string> hdf_confocal_files;

// ----------------------------------------------------------------------------

bool compare_file_size (const file_name_size& first, const file_name_size& second)
{
    return ( first.total_rank_size > second.total_rank_size );
}

// ----------------------------------------------------------------------------

void populate_netcdf_hdf5_files(std::string dataset_dir)
{
    //populate netcdf and hdf5 files for fly scans
    netcdf_files = find_all_dataset_files(dataset_dir + "flyXRF/", "_0.nc");
    hdf_files = find_all_dataset_files(dataset_dir + "flyXRF.h5/", "_0.h5");
    hdf_xspress_files = find_all_dataset_files(dataset_dir + "flyXspress/", "_0.h5");
    hdf_confocal_files = find_all_dataset_files(dataset_dir , ".hdf5");
}

//-----------------------------------------------------------------------------

fitting::routines::Base_Fit_Routine* generate_fit_routine(data_struct::Fitting_Routines proc_type, fitting::optimizers::Optimizer* optimizer)
{
    //Fitting routines
    fitting::routines::Base_Fit_Routine *fit_routine = nullptr;
    switch(proc_type)
    {
		case data_struct::Fitting_Routines::GAUSS_TAILS:
            fit_routine = new fitting::routines::Param_Optimized_Fit_Routine();
            ((fitting::routines::Param_Optimized_Fit_Routine*)fit_routine)->set_optimizer(optimizer);
            break;
        case data_struct::Fitting_Routines::GAUSS_MATRIX:
            fit_routine = new fitting::routines::Matrix_Optimized_Fit_Routine();
            ((fitting::routines::Matrix_Optimized_Fit_Routine*)fit_routine)->set_optimizer(optimizer);
            break;
        case data_struct::Fitting_Routines::ROI:
            fit_routine = new fitting::routines::ROI_Fit_Routine();
            break;
        case data_struct::Fitting_Routines::SVD:
            fit_routine = new fitting::routines::SVD_Fit_Routine();
            break;
        case data_struct::Fitting_Routines::NNLS:
            fit_routine = new fitting::routines::NNLS_Fit_Routine();
            break;
        default:
            break;
    }
    return fit_routine;
}

// ----------------------------------------------------------------------------

bool init_analysis_job_detectors(data_struct::Analysis_Job* analysis_job)
{

//    _default_sub_struct.
//    if( false == io::load_override_params(_dataset_directory, -1, override_params) )
//    {
//        return false;
//    }

    analysis_job->detectors_meta_data.clear();



    //initialize models and fit routines for all detectors
    for(size_t detector_num = analysis_job->detector_num_start; detector_num <= analysis_job->detector_num_end; detector_num++)
    {
        analysis_job->detectors_meta_data[detector_num] = data_struct::Analysis_Sub_Struct();
        data_struct::Analysis_Sub_Struct *sub_struct = &analysis_job->detectors_meta_data[detector_num];

        sub_struct->model = new fitting::models::Gaussian_Model();
        data_struct::Params_Override * override_params = &(sub_struct->fit_params_override_dict);

        override_params->dataset_directory = analysis_job->dataset_directory;
        override_params->detector_num = detector_num;

        if( false == io::load_override_params(analysis_job->dataset_directory, detector_num, override_params) )
        {
            if( false == io::load_override_params(analysis_job->dataset_directory, -1, override_params) )
            {
                return false;
            }
        }

        if (override_params->elements_to_fit.size() < 1)
        {
            logit<<"Error, no elements to fit. Check  maps_fit_parameters_override.txt0 - 3 exist"<<"\n";
            return false;
        }

        for(auto proc_type : analysis_job->fitting_routines)
        {
            //Fitting models
            sub_struct->fit_routines[proc_type] = generate_fit_routine(proc_type, analysis_job->optimizer());

            //reset model fit parameters to defaults
            sub_struct->model->reset_to_default_fit_params();
            //Update fit parameters by override values
            sub_struct->model->update_fit_params_values(&(override_params->fit_params));

            //Fit_Element_Map_Dict *elements_to_fit = &(sub_struct->fit_params_override_dict.elements_to_fit);
            //Initialize model
            //fit_routine->initialize(sub_struct->model, elements_to_fit, energy_range);

        }        
    }

    return true;
}

// ----------------------------------------------------------------------------

bool load_element_info(std::string element_henke_filename, std::string element_csv_filename, data_struct::Element_Info_Map *element_info_map)
{

	if (io::file::load_henke_from_xdr(element_henke_filename, element_info_map) == false)
	{
		logit << "error loading " << element_henke_filename << "\n";
		return false;
	}

	if (io::file::load_element_info_from_csv(element_csv_filename, element_info_map) == false)
	{
		logit << "error loading " << element_csv_filename << "\n";
		return false;
	}

    return true;
}

// ----------------------------------------------------------------------------

bool save_results(std::string save_loc,
                  const data_struct::Fit_Count_Dict * const element_counts,
                  std::queue<std::future<bool> >* job_queue,
                  std::chrono::time_point<std::chrono::system_clock> start)
{


    //wait for queue to finish processing
    while(!job_queue->empty())
    {
        auto ret = std::move(job_queue->front());
        job_queue->pop();
        ret.get();
    }

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    logit << "Fitting [ "<< save_loc <<" ] elapsed time: " << elapsed_seconds.count() << "s"<<"\n";


    io::file::HDF5_IO::inst()->save_element_fits(save_loc, element_counts);

    delete job_queue;
    delete element_counts;

    return true;
}

// ----------------------------------------------------------------------------

bool save_volume(data_struct::Quantification_Standard * quantification_standard,
                 data_struct::Spectra_Volume *spectra_volume,
                 real_t energy_offset,
                 real_t energy_slope,
                 real_t energy_quad)
{
    io::file::HDF5_IO::inst()->save_quantification(quantification_standard);
    io::file::HDF5_IO::inst()->save_spectra_volume("mca_arr", spectra_volume, energy_offset, energy_slope, energy_quad);
    io::file::HDF5_IO::inst()->end_save_seq();

    delete spectra_volume;

    return true;
}

// ----------------------------------------------------------------------------

void save_optimized_fit_params(struct file_name_fit_params file_and_fit_params)
{
    io::file::CSV_IO csv_io;
    std::string full_path = file_and_fit_params.dataset_dir+"/output/"+file_and_fit_params.dataset_filename+std::to_string(file_and_fit_params.detector_num)+".csv";
    logit<<full_path<<"\n";

    fitting::models::Gaussian_Model model;
    //Range of energy in spectra to fit
    fitting::models::Range energy_range = data_struct::get_energy_range(file_and_fit_params.spectra.size(), &(file_and_fit_params.fit_params));

	data_struct::Spectra model_spectra = model.model_spectrum(&file_and_fit_params.fit_params, &file_and_fit_params.elements_to_fit, energy_range);

    if (file_and_fit_params.fit_params.contains(STR_SNIP_WIDTH))
	{
        data_struct::Fit_Param fit_snip_width = file_and_fit_params.fit_params[STR_SNIP_WIDTH];
        if (fit_snip_width.bound_type > data_struct::E_Bound_Type::FIXED)
		{
			real_t spectral_binning = 0.0;
			data_struct::ArrayXr s_background = data_struct::snip_background(&file_and_fit_params.spectra,
																			file_and_fit_params.fit_params.value(STR_ENERGY_OFFSET),
																			file_and_fit_params.fit_params.value(STR_ENERGY_SLOPE),
																			file_and_fit_params.fit_params.value(STR_ENERGY_QUADRATIC),
																			spectral_binning,
																			file_and_fit_params.fit_params.value(STR_SNIP_WIDTH),
																			energy_range.min,
																			energy_range.max);
			data_struct::ArrayXr background = s_background.segment(energy_range.min, energy_range.count());
			model_spectra += background;
		}
	}
	
	ArrayXr spec = file_and_fit_params.spectra.segment(energy_range.min, energy_range.count());
	
	data_struct::ArrayXr energy = data_struct::ArrayXr::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
	data_struct::ArrayXr ev = (file_and_fit_params.fit_params.value(STR_ENERGY_OFFSET) + energy) * file_and_fit_params.fit_params.value(STR_ENERGY_SLOPE) + (pow(energy, (real_t)2.0) * file_and_fit_params.fit_params.value(STR_ENERGY_QUADRATIC));
	data_struct::ArrayXr ev2 = file_and_fit_params.fit_params.at(STR_ENERGY_OFFSET).value + energy * file_and_fit_params.fit_params.at(STR_ENERGY_SLOPE).value + pow(energy, (real_t)2.0) * file_and_fit_params.fit_params.at(STR_ENERGY_QUADRATIC).value;
	int diff = (ev2 - ev).sum();

#ifdef _BUILD_WITH_VTK
    std::string str_path = file_and_fit_params.dataset_dir+"/output/fit_"+file_and_fit_params.dataset_filename+"_det"+std::to_string(file_and_fit_params.detector_num)+".png";
    visual::SavePlotSpectras(str_path, energy, spec, model_spectra, true);
#endif

    // TODO: save the spectra, model, and background to csv
    csv_io.save_fit_parameters(full_path, file_and_fit_params.fit_params );

    for(auto &itr : file_and_fit_params.elements_to_fit)
    {
        delete itr.second;
    }

}

// ----------------------------------------------------------------------------

void save_averaged_fit_params(std::string dataset_dir, std::unordered_map<int, data_struct::Fit_Parameters> fit_params_avgs, int detector_num_start, int detector_num_end)
{
    io::file::aps::APS_Fit_Params_Import aps_io;
    int i =0;
    std::string full_path = dataset_dir+"/maps_fit_parameters_override.txt";
    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        aps_io.save(full_path, fit_params_avgs[i], detector_num );
        i++;
    }
}

// ----------------------------------------------------------------------------

bool load_quantification_standard(std::string dataset_directory,
                                  std::string quantification_info_file,
                                  std::string *standard_file_name,
                                  std::unordered_map<std::string, real_t> *element_standard_weights)
                                  //data_struct::Quantification_Standard * quantification_standard)
{
    std::string path = dataset_directory + quantification_info_file;
    std::ifstream paramFileStream(path);

    if (paramFileStream.is_open() )
    {
        paramFileStream.exceptions(std::ifstream::failbit);
        bool has_filename = false;
        bool has_elements = false;
        bool has_weights = false;
        //std::string line;
        std::string tag;

        std::vector<std::string> element_names;
        std::vector<real_t> element_weights;

        try
        {

            for (std::string line; std::getline(paramFileStream, line); )
            //while(std::getline(paramFileStream, line))
            {
                std::istringstream strstream(line);
                std::getline(strstream, tag, ':');
                //logit<<"tag : "<<tag<<"\n";
                if (tag == "FILENAME")
                {
                    std::string standard_filename;
                    logit << line << "\n";
                    std::getline(strstream, standard_filename, ':');
                    standard_filename.erase(std::remove_if(standard_filename.begin(), standard_filename.end(), ::isspace), standard_filename.end());
                    logit << "Standard file name = "<< standard_filename << "\n";
                    //quantification_standard->standard_filename(standard_filename);
                    (*standard_file_name) = standard_filename;
                    has_filename = true;
                }
                else if (tag == "ELEMENTS_IN_STANDARD")
                {
                    std::string element_symb;
                    while(std::getline(strstream, element_symb, ','))
                    {
                        element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());
                        logit<<"Element : "<<element_symb<<"\n";
                        element_names.push_back(element_symb);
                    }
                    has_elements = true;
                }
                else if (tag == "WEIGHT")
                {
                    std::string element_weight_str;
                    while(std::getline(strstream, element_weight_str, ','))
                    {
                        element_weight_str.erase(std::remove_if(element_weight_str.begin(), element_weight_str.end(), ::isspace), element_weight_str.end());
                        logit<<"Element weight: "<<element_weight_str<<"\n";
                        real_t weight = std::stof(element_weight_str);
                        element_weights.push_back(weight);
                    }
                    has_weights = true;
                }

            }
        }
        catch(std::exception e)
        {
            if (paramFileStream.eof() == 0 && (paramFileStream.bad() || paramFileStream.fail()) )
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << paramFileStream.fail()
                    << "\neofbit: " << paramFileStream.eof()
                    << "\nbadbit: " << paramFileStream.bad() << "\n";
            }
        }


        paramFileStream.close();
        if(has_filename && has_elements && has_weights)
        {
            if(element_names.size() == element_weights.size())
            {
                for(size_t i=0; i<element_names.size(); i++)
                {
                    (*element_standard_weights)[element_names[i]] = element_weights[i];
                    //quantification_standard->append_element(element_names[i], element_weights[i]);
                }
            }
            else
            {
                logit<<"Error: number of element names ["<<element_names.size()<<"] does not match number of element weights ["<<element_weights.size()<<"]!"<<"\n";
            }

            return true;
        }

    }
    else
    {
        logit<<"Failed to open file "<<path<<"\n";
    }
    return false;

}

// ----------------------------------------------------------------------------

bool load_override_params(std::string dataset_directory,
                          int detector_num,
                          data_struct::Params_Override *params_override)
{
    //Importer for APS datasets
    io::file::aps::APS_Fit_Params_Import fit_param_importer;

    std::string det_num = "";
    if(detector_num > -1)
        det_num = std::to_string(detector_num);

    if(false == fit_param_importer.load(dataset_directory+"maps_fit_parameters_override.txt"+det_num,
                                        data_struct::Element_Info_Map::inst(),
                                        params_override))
    {
        logit<<"Error loading fit param override file: "<<"\n";
        return false;
    }
    else
    {
        data_struct::Element_Info* detector_element;
        if(params_override->detector_element.length() > 0)
        {
            // Get the element info class                                   // detector element as string "Si" or "Ge" usually

            detector_element = data_struct::Element_Info_Map::inst()->get_element(params_override->detector_element);
        }
        else
        {
         //log error or warning
            logit<<"Error, no detector material defined in maps_fit_parameters_override.txt . Defaulting to Si"<<"\n";
            detector_element = data_struct::Element_Info_Map::inst()->get_element("Si");
        }


        //add compton and coherant amp
        if(params_override->elements_to_fit.count(STR_COMPTON_AMPLITUDE) == 0)
        {
            params_override->elements_to_fit.insert(std::pair<std::string, data_struct::Fit_Element_Map*>(STR_COMPTON_AMPLITUDE, new data_struct::Fit_Element_Map(STR_COMPTON_AMPLITUDE, nullptr)) );
        }
        if(params_override->elements_to_fit.count(STR_COHERENT_SCT_AMPLITUDE) == 0)
        {
            params_override->elements_to_fit.insert(std::pair<std::string, data_struct::Fit_Element_Map*>(STR_COHERENT_SCT_AMPLITUDE, new data_struct::Fit_Element_Map(STR_COHERENT_SCT_AMPLITUDE, nullptr)) );
        }

        logit<<"Elements to fit:  ";
        //Update element ratios by detector element
        for(auto& itr : params_override->elements_to_fit)
        {
            itr.second->init_energy_ratio_for_detector_element(detector_element);
            logit_s<<itr.first<<" ";
        }
        logit_s<<"\n";

    }

    return true;
}

// ----------------------------------------------------------------------------

bool load_spectra_volume(std::string dataset_directory,
                         std::string dataset_file,
                         data_struct::Spectra_Volume *spectra_volume,
                         size_t detector_num,
                         data_struct::Params_Override * params_override,
                         data_struct::Quantification_Standard * quantification_standard,
                         bool *is_loaded_from_analyazed_h5,
                         bool save_scalers)
{

    //Dataset importer
    io::file::MDA_IO mda_io;
    //data_struct::Detector detector;
    std::string tmp_dataset_file = dataset_file;

    logit<<"Loading dataset "<<dataset_directory+"mda/"+dataset_file<<" detector "<<detector_num<<"\n";

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    bool hasNetcdf = false;
    bool hasHdf = false;
    bool hasXspress = false;
    std::string file_middle = ""; //_2xfm3_ or dxpM...
    for(auto &itr : netcdf_files)
    {
        if (itr.find(tmp_dataset_file) == 0)
        {
            size_t slen = (itr.length()-4) - tmp_dataset_file.length();
            file_middle = itr.substr(tmp_dataset_file.length(), slen);
            hasNetcdf = true;
            break;
        }
    }
    if (hasNetcdf == false)
    {
        for(auto &itr : hdf_files)
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length()-4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasHdf = true;
                break;
            }
        }
    }
    if (hasNetcdf == false && hasHdf == false)
    {
        for(auto &itr : hdf_xspress_files)
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length()-4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasXspress = true;
                break;
            }
        }
    }

    //  try to load from a pre analyzed file because they should contain the whole mca_arr spectra volume
    if(true == io::file::HDF5_IO::inst()->load_spectra_vol_analyzed_h5(dataset_directory+"img.dat/"+dataset_file, detector_num, spectra_volume))
    {
        *is_loaded_from_analyazed_h5 = true;
        io::file::HDF5_IO::inst()->start_save_seq(false);
        return true;
    }
    else
    {
        *is_loaded_from_analyazed_h5 = false;
    }

    //try loading confocal dataset
    if(true == io::file::HDF5_IO::inst()->load_spectra_volume_confocal(dataset_directory+"/"+dataset_file, detector_num, spectra_volume))
    {
        *is_loaded_from_analyazed_h5 = false;
        if(save_scalers)
        {
            io::file::HDF5_IO::inst()->start_save_seq(true);
            io::file::HDF5_IO::inst()->save_scan_scalers_confocal(dataset_directory+"/"+dataset_file, detector_num);
        }
        return true;
    }

    //load spectra
    if (false == mda_io.load_spectra_volume(dataset_directory+"mda/"+dataset_file, detector_num, spectra_volume, hasNetcdf | hasHdf | hasXspress, params_override, quantification_standard) )
    {
        logit<<"Error load spectra "<<dataset_directory+"mda/"+dataset_file<<"\n";
        return false;
    }
    else
    {
        if(hasNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for(size_t i=0; i<spectra_volume->rows(); i++)
                {
                    full_filename = dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + std::to_string(i) + ".nc";
                    //todo: add verbose option
                    //logit<<"Loading file "<<full_filename<<"\n";
                    io::file::NetCDF_IO::inst()->load_spectra_line(full_filename, detector_num, &(*spectra_volume)[i]);
                }
            }
            else
            {
                logit<<"Did not find netcdf files "<<dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc"<<"\n";
                //return false;
            }
        }
        else if (hasHdf)
        {
            io::file::HDF5_IO::inst()->load_spectra_volume(dataset_directory + "flyXRF.h5/" + tmp_dataset_file + file_middle + "0.h5", detector_num, spectra_volume);
        }
        else if (hasXspress)
        {
            std::string full_filename;
            for(size_t i=0; i<spectra_volume->rows(); i++)
            {
                full_filename = dataset_directory + "flyXspress/" + tmp_dataset_file + file_middle + std::to_string(i) + ".h5";
                io::file::HDF5_IO::inst()->load_spectra_line_xspress3(full_filename, detector_num, &(*spectra_volume)[i]);
            }
        }

    }

    if(save_scalers)
    {
        io::file::HDF5_IO::inst()->start_save_seq(true);
        io::file::HDF5_IO::inst()->save_scan_scalers(detector_num, mda_io.get_scan_ptr(), params_override, hasNetcdf | hasHdf | hasXspress);
    }

    mda_io.unload();
    logit<<"Finished Loading dataset "<<dataset_directory+"mda/"+dataset_file<<" detector "<<detector_num<<"\n";
    return true;
}

// ----------------------------------------------------------------------------

bool load_and_integrate_spectra_volume(std::string dataset_directory,
                                       std::string dataset_file,
                                       data_struct::Spectra *integrated_spectra,
                                       size_t detector_num,
                                       data_struct::Params_Override * params_override)
{
    //Dataset importer
    io::file::MDA_IO mda_io;
    //data_struct::Detector detector;
    std::string tmp_dataset_file = dataset_file;
    bool ret_val = true;

    data_struct::Spectra_Volume spectra_volume;

    logit<<"Loading dataset "<<dataset_directory+"mda/"+dataset_file<<"\n";

    //TODO: check if any analyized mda.h5 files are around to load. They should have integrated spectra saved already.

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    bool hasNetcdf = false;
    bool hasHdf = false;
    bool hasXspress = false;
    std::string file_middle = ""; //_2xfm3_ or dxpM...
    for(auto &itr : netcdf_files)
    {
        if (itr.find(tmp_dataset_file) == 0)
        {
            size_t slen = (itr.length()-4) - tmp_dataset_file.length();
            file_middle = itr.substr(tmp_dataset_file.length(), slen);
            hasNetcdf = true;
            break;
        }
    }
    if (hasNetcdf == false)
    {
        for(auto &itr : hdf_files)
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length()-4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasHdf = true;
                break;
            }
        }
    }
    if (hasNetcdf == false && hasHdf == false)
    {
        for(auto &itr : hdf_xspress_files)
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length()-4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasXspress = true;
                break;
            }
        }
    }

    //  try to load from a pre analyzed file because they should contain the integrated spectra
    if(true == io::file::HDF5_IO::inst()->load_integrated_spectra_analyzed_h5(dataset_directory+"img.dat/"+dataset_file, detector_num, integrated_spectra))
    {
        return true;
    }


    //try loading confocal dataset
    if(true == io::file::HDF5_IO::inst()->load_spectra_volume_confocal(dataset_directory+"/"+dataset_file, detector_num, &spectra_volume))
    {
        *integrated_spectra = spectra_volume.integrate();
        return true;
    }

    //load spectra
    if (false == hasNetcdf && false == hasHdf)
    {
        ret_val = mda_io.load_spectra_volume(dataset_directory+"mda/"+dataset_file, detector_num, &spectra_volume, hasNetcdf | hasHdf | hasXspress, params_override, nullptr);
        if(ret_val)
        {
            *integrated_spectra = spectra_volume.integrate();
        }
        mda_io.unload();
    }
    else
    {
        int rank;
        int dims[10];
        dims[0] = 0;
        rank = mda_io.get_rank_and_dims(dataset_directory + "mda/" + dataset_file, &dims[0]);
        if(rank == 3)
        {
            integrated_spectra->resize(dims[2]);
        }
        else
        {
            integrated_spectra->resize(2048);
        }


        if(hasNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for(size_t i=0; i<dims[0]; i++)
                {
                    data_struct::Spectra_Line spectra_line;
                    spectra_line.resize(dims[1], integrated_spectra->size());
                    full_filename = dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + std::to_string(i) + ".nc";
                    //logit<<"Loading file "<<full_filename<<"\n";
                    if( io::file::NetCDF_IO::inst()->load_spectra_line(full_filename, detector_num, &spectra_line) )
                    {
                        for(int k=0; k<spectra_line.size(); k++)
                        {
                            *integrated_spectra += spectra_line[k];
                        }
                    }
                }
            }
            else
            {
                logit<<"Did not find netcdf files "<<dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc"<<"\n";
                //return false;
            }
        }
        else if (hasHdf)
        {
            ret_val = io::file::HDF5_IO::inst()->load_and_integrate_spectra_volume(dataset_directory + "flyXRF.h5/" + tmp_dataset_file + file_middle + "0.h5", detector_num, integrated_spectra);
        }
        else if (hasXspress)
        {
            std::string full_filename;
            data_struct::Spectra_Line spectra_line;
            spectra_line.resize(dims[1], integrated_spectra->size());
            for(size_t i=0; i<dims[0]; i++)
            {
                full_filename = dataset_directory + "flyXspress/" + tmp_dataset_file + file_middle + std::to_string(i) + ".h5";
                if(io::file::HDF5_IO::inst()->load_spectra_line_xspress3(full_filename, detector_num, &spectra_line))
                {
                    for(int k=0; k<spectra_line.size(); k++)
                    {
                        *integrated_spectra += spectra_line[k];
                    }
                }
            }
        }

    }

    logit<<"Finished Loading dataset "<<dataset_directory+"mda/"+dataset_file<<" detector "<<detector_num<<"\n";
    return ret_val;
}

// ----------------------------------------------------------------------------

void generate_h5_averages(std::string dataset_directory,
                          std::string dataset_file,
                          size_t detector_num_start,
                          size_t detector_num_end)
{
    logit<<"\n";

    std::vector<std::string> hdf5_filenames;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    if (detector_num_start == detector_num_end)
    {
        logit << "Warning: detector range "<<detector_num_start<<":"<<detector_num_end<<" is only 1 detector. Nothing to avg."<<"\n";
        return;
    }


    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        hdf5_filenames.push_back(dataset_directory+"img.dat/"+dataset_file+".h5"+std::to_string(detector_num));
    }

    io::file::HDF5_IO::inst()->generate_avg(dataset_directory+"img.dat/"+dataset_file+".h5", hdf5_filenames);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";


}

// ----------------------------------------------------------------------------

std::vector<std::string> find_all_dataset_files(std::string dataset_directory, std::string search_str)
{
    std::vector<std::string> dataset_files;
    logit<<dataset_directory<<" searching for "<<search_str<<"\n";
    DIR *dir;
    struct dirent *ent;
    size_t search_str_size = search_str.length();
    if ((dir = opendir (dataset_directory.c_str())) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL)
        {
            std::string fname(ent->d_name);
            // check if extension is .mda
            if (fname.size() > 4)
            {
                if (fname.rfind(search_str) == fname.size() -search_str_size)
                {
                    dataset_files.push_back(fname);
                }
            }
        }
        closedir (dir);
    }
    else
    {
        /* could not open directory */
        logit<<"Error: could not open directory "<<dataset_directory<<"\n";
    }

    logit<<"found "<<dataset_files.size()<<"\n";
    return dataset_files;
}

// ----------------------------------------------------------------------------

void check_and_create_dirs(std::string dataset_directory)
{

    bool found_img_dat = false;
    logit<<dataset_directory<<"\n";
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dataset_directory.c_str())) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL)
        {
            if( strcmp(ent->d_name , "img.dat") == 0)
            {
                found_img_dat = true;
                break;
            }
        }
        closedir (dir);
    }
    else
    {
        /* could not open directory */
        logit<<"Error: could not open directory "<<dataset_directory<<"\n";
    }

    if (false == found_img_dat)
    {
        std::string cmd = "mkdir "+dataset_directory+"img.dat";
        system(cmd.c_str());
    }
    logit<<"done"<<"\n";

}

// ----------------------------------------------------------------------------

void sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string> *dataset_files)
{

    io::file::MDA_IO mda_io;
    logit<<dataset_directory<<" "<<dataset_files->size()<<" files"<<"\n";
    std::list<file_name_size> f_list;

    for (auto &itr : *dataset_files)
    {
        std::string full_path = dataset_directory + "/mda/"+itr;
        int fsize = mda_io.get_multiplied_dims(full_path);
        f_list.push_back(file_name_size(itr, fsize));
    }

    f_list.sort(compare_file_size);

    dataset_files->clear();

    for(auto &itr : f_list)
    {
        dataset_files->push_back(itr.filename);
    }

    logit<<"done"<<"\n";
}

}// end namespace io
