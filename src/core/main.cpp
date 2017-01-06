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


//#include <QCoreApplication>

#include <iostream>
#include <queue>
#include <string>
#include <array>
#include <vector>
#include <unordered_map>
#include <chrono>
#include <ctime>
#include <limits>
#include <sstream>
#include <fstream>

//#include "pybind11.h"

#include <dirent.h>

#include "threadpool.h"

#include "netcdf_io.h"
#include "mda_io.h"
#include "hdf5_io.h"
#include "csv_io.h"

#include "spectra_volume.h"
#include "detector.h"

//#include "task_queue.h"
//#include "general_function_task.h"

//#include "threadpool_distributor.h"

#include "gaussian_model.h"

#include "element_info.h"

#include "aps_fit_params_import.h"

#include "roi_fit_routine.h"
#include "svd_fit_routine.h"
#include "nnls_fit_routine.h"

#include "lmfit_optimizer.h"
//#include "mpfit_optimizer.h"

#include "fit_element_map.h"

#include "quantification_standard.h"

//#include "vtk_graph.h"

#include "command_line_parser.h"


using namespace std::placeholders; //for _1, _2,

// ----------------------------------------------------------------------------

enum Processing_Type { ROI=1 , GAUSS_TAILS=2, GAUSS_MATRIX=4, SVD=8, NNLS=16 };

struct file_name_size
{
    file_name_size(std::string name, int size) { filename = name; total_rank_size = size;}
    std::string filename;
    int total_rank_size;
};

struct file_name_fit_params
{
    std::string dataset_dir;
    std::string dataset_filename;
    int detector_num;
    data_struct::xrf::Fit_Parameters fit_params;
    bool success;
};

const std::unordered_map<int, std::string> save_loc_map = {
    {ROI, "ROI"},
    {GAUSS_TAILS, "Params"},
    {GAUSS_MATRIX, "Fitted"},
    {SVD, "SVD"},
    {NNLS, "NNLS"}
};

//Optimizers for fitting models
fitting::optimizers::LMFit_Optimizer lmfit_optimizer;

// ----------------------------------------------------------------------------

bool load_spectra_volume(std::string dataset_directory,
                         std::string dataset_file,
                         data_struct::xrf::Spectra_Volume *spectra_volume,
                         size_t detector_num,
                         std::unordered_map< std::string, std::string > *extra_override_values,
                         data_struct::xrf::Quantification_Standard * quantification_standard,
                         bool save_scalers);

// ----------------------------------------------------------------------------

bool load_element_info(std::string element_henke_filename, std::string element_csv_filename, data_struct::xrf::Element_Info_Map *element_info_map)
{

    if (io::file::load_element_info_from_csv(element_csv_filename, element_info_map) == false)
    {
        std::cout<<"error loading "<< element_csv_filename<<std::endl;
        return false;
    }

    if (io::file::load_henke_from_xdr(element_henke_filename, element_info_map) == false)
    {
        std::cout<<"error loading "<< element_henke_filename<<std::endl;
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

bool save_results(std::string full_path,
                  std::string save_loc,
                  const data_struct::xrf::Fit_Count_Dict * const element_counts,
                  fitting::routines::Base_Fit_Routine* fit_routine,
                  std::queue<std::future<bool> >* job_queue)
{

    //wait for queue to finish processing
    while(!job_queue->empty())
    {
        auto ret = std::move(job_queue->front());
        job_queue->pop();
        ret.get();
    }

    io::file::HDF5_IO hdf5_io;
    hid_t f_id = hdf5_io.start_save_seq(full_path);
    hdf5_io.save_element_fits(f_id, save_loc, element_counts);
    hdf5_io.end_save_seq(f_id);


    delete fit_routine;
    delete job_queue;
    delete element_counts;

    return true;
}

// ----------------------------------------------------------------------------

bool save_volume(std::string full_path,
                 data_struct::xrf::Quantification_Standard * quantification_standard,
                 data_struct::xrf::Spectra_Volume *spectra_volume,
                 real_t energy_offset,
                 real_t energy_slope,
                 real_t energy_quad,
                 std::queue<std::future<bool> >* job_queue)
{
    //wait for queue to finish processing
    while(!job_queue->empty())
    {
        auto ret = std::move(job_queue->front());
        job_queue->pop();
        ret.get();
    }

    io::file::HDF5_IO hdf5_io;
    hid_t f_id = hdf5_io.start_save_seq(full_path);
    hdf5_io.save_quantification(f_id, quantification_standard);
    hdf5_io.save_spectra_volume(f_id, "mca_arr", spectra_volume, energy_offset, energy_slope, energy_quad);
    hdf5_io.end_save_seq(f_id);

    delete job_queue;
    delete spectra_volume;

    return true;
}

// ----------------------------------------------------------------------------

void save_optimized_fit_params(struct file_name_fit_params file_and_fit_params)
{
    io::file::CSV_IO csv_io;
    std::string full_path = file_and_fit_params.dataset_dir+"/output/"+file_and_fit_params.dataset_filename+std::to_string(file_and_fit_params.detector_num)+".csv";
    std::cout<<"save_optimized_fit_params(): "<<full_path<<std::endl;
    csv_io.save_fit_parameters(full_path, file_and_fit_params.fit_params );
}

// ----------------------------------------------------------------------------

void save_averaged_fit_params(std::string dataset_dir, std::vector<data_struct::xrf::Fit_Parameters> fit_params_avgs, int detector_num_start, int detector_num_end)
{
    io::file::CSV_IO csv_io;
    int i =0;
    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        std::string full_path = dataset_dir+"/avrg_maps_fit_override_parameters.txt" + std::to_string(detector_num);
        std::cout<<"save_averaged_fit_params(): "<<full_path<<std::endl;
        csv_io.save_fit_parameters(full_path, fit_params_avgs[i] );
        i++;
    }
}

// ----------------------------------------------------------------------------

fitting::routines::Base_Fit_Routine * generate_fit_routine(Processing_Type proc_type)
{
    //Fitting routines
    fitting::routines::Base_Fit_Routine *fit_routine = nullptr;
    switch(proc_type)
    {
        case Processing_Type::GAUSS_TAILS:
            fit_routine = new fitting::routines::Param_Optimized_Fit_Routine();
            ((fitting::routines::Param_Optimized_Fit_Routine*)fit_routine)->set_optimizer(&lmfit_optimizer);
            break;
        case Processing_Type::GAUSS_MATRIX:
            fit_routine = new fitting::routines::Matrix_Optimized_Fit_Routine();
            ((fitting::routines::Matrix_Optimized_Fit_Routine*)fit_routine)->set_optimizer(&lmfit_optimizer);
            break;
        case Processing_Type::ROI:
            fit_routine = new fitting::routines::ROI_Fit_Routine();
            break;
        case Processing_Type::SVD:
            fit_routine = new fitting::routines::SVD_Fit_Routine();
            break;
        case Processing_Type::NNLS:
            fit_routine = new fitting::routines::NNLS_Fit_Routine();
            break;
        default:
            break;
    }
    return fit_routine;
}

// ----------------------------------------------------------------------------

//data_struct::xrf::Fit_Count_Dict* generate_fit_count_dict(data_struct::xrf::Fit_Element_Map_Dict *elements_to_fit, size_t width, size_t height )
template<typename T>
data_struct::xrf::Fit_Count_Dict* generate_fit_count_dict(std::unordered_map<std::string, T> *elements_to_fit, size_t width, size_t height )
{
    data_struct::xrf::Fit_Count_Dict* element_fit_counts_dict = new data_struct::xrf::Fit_Count_Dict();
    for(auto& e_itr : *elements_to_fit)
    {
        element_fit_counts_dict->emplace(std::pair<std::string, data_struct::xrf::Fit_Counts_Array>(e_itr.first, data_struct::xrf::Fit_Counts_Array()) );
        element_fit_counts_dict->at(e_itr.first).resize(width, height);
    }
    return element_fit_counts_dict;
}

// ----------------------------------------------------------------------------

bool fit_single_spectra(fitting::routines::Base_Fit_Routine * fit_routine,
                        const fitting::models::Base_Model * const model,
                        const data_struct::xrf::Spectra * const spectra,
                        const data_struct::xrf::Fit_Element_Map_Dict * const elements_to_fit,
                        data_struct::xrf::Fit_Count_Dict * out_fit_counts,
                        size_t i,
                        size_t j)
{
    std::unordered_map<std::string, real_t> counts_dict = fit_routine->fit_spectra(model, spectra, elements_to_fit);
    //save count / sec
    for (auto& el_itr : *elements_to_fit)
    {
        (*out_fit_counts)[el_itr.first][i][j] = counts_dict[el_itr.first] / spectra->elapsed_lifetime();
    }
    (*out_fit_counts)[data_struct::xrf::STR_NUM_ITR][i][j] = counts_dict[data_struct::xrf::STR_NUM_ITR];
    return true;
}

// ----------------------------------------------------------------------------

bool load_quantification_standard(std::string dataset_directory,
                                  std::string quantification_info_file,
                                  std::string *standard_file_name,
                                  std::unordered_map<std::string, real_t> *element_standard_weights)
                                  //data_struct::xrf::Quantification_Standard * quantification_standard)
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
                //std::cout<<"tag : "<<tag<<std::endl;
                if (tag == "FILENAME")
                {
                    std::string standard_filename;
                    std::cout << line << std::endl;
                    std::getline(strstream, standard_filename, ':');
                    standard_filename.erase(std::remove_if(standard_filename.begin(), standard_filename.end(), ::isspace), standard_filename.end());
                    std::cout << "Standard file name = "<< standard_filename << std::endl;
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
                        std::cout<<"Element : "<<element_symb<<std::endl;
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
                        std::cout<<"Element weight: "<<element_weight_str<<std::endl;
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
                    << "\nbadbit: " << paramFileStream.bad() << std::endl;
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
                std::cout<<"Error: number of element names ["<<element_names.size()<<"] does not match number of element weights ["<<element_weights.size()<<"]!"<<std::endl;
            }

            return true;
        }

    }
    else
    {
        std::cout<<"Failed to open file "<<path<<std::endl;
    }
    return false;

}


// ----------------------------------------------------------------------------

bool load_override_params(std::string dataset_directory,
                          int detector_num,
                          data_struct::xrf::Fit_Parameters *fit_params,
                          data_struct::xrf::Fit_Element_Map_Dict *elements_to_fit,
                          std::unordered_map< std::string, std::string > *extra_override_values)
{
    //Importer for APS datasets
    io::file::aps::APS_Fit_Params_Import fit_param_importer;

    std::string det_num = "";
    if(detector_num > -1)
        det_num = std::to_string(detector_num);

    if(false == fit_param_importer.load(dataset_directory+"maps_fit_parameters_override.txt"+det_num,
                                        data_struct::xrf::Element_Info_Map::inst(),
                                        fit_params,
                                        elements_to_fit,
                                        extra_override_values))
    {
        std::cout<<"Error loading fit param override file: "<<std::endl;
        return false;
    }
    else
    {
        data_struct::xrf::Element_Info* detector_element;
        if (extra_override_values->count(data_struct::xrf::STR_DETECTOR_ELEMENT) > 0)
        {
            // Get the element info class                                                                                 // detector element as string "Si" or "Ge" usually
            detector_element = data_struct::xrf::Element_Info_Map::inst()->get_element(extra_override_values->at(data_struct::xrf::STR_DETECTOR_ELEMENT));
        }
        else
        {
         //log error or warning
            std::cout<<"Error, no detector material defined in maps_fit_parameters_override.txt . Defaulting to Si";
            detector_element = data_struct::xrf::Element_Info_Map::inst()->get_element("Si");
        }

        std::cout<<"Elements to fit: "<<std::endl;
        //Update element ratios by detector element
        for(auto& itr : *elements_to_fit)
        {
            itr.second->init_energy_ratio_for_detector_element(detector_element);
            std::cout<<itr.first<<" ";
        }
        std::cout<<std::endl;

    }

    return true;
}

// ----------------------------------------------------------------------------

bool load_spectra_volume(std::string dataset_directory,
                         std::string dataset_file,
                         data_struct::xrf::Spectra_Volume *spectra_volume,
                         size_t detector_num,
                         std::unordered_map< std::string, std::string > *extra_override_values,
                         data_struct::xrf::Quantification_Standard * quantification_standard,
                         bool save_scalers)
{

    //Dataset importer
    io::file::MDA_IO mda_io;
    io::file::HDF5_IO hdf5_io;
    io::file::NetCDF_IO netcdf_io;
    //data_struct::xrf::Detector detector;
    std::string tmp_dataset_file = dataset_file;

    std::cout<<"Loading dataset "<<dataset_directory+"mda/"+dataset_file<<std::endl;

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    std::ifstream file_io_0(dataset_directory+"flyXRF/"+tmp_dataset_file+"_2xfm3__0.nc");
    std::ifstream file_io_hdf(dataset_directory+"flyXRF.h5/"+tmp_dataset_file+"_2xfm3__0.h5");
    bool hasNetcdf = file_io_0.is_open();
    bool hasHdf = file_io_hdf.is_open();
    if(hasNetcdf)
        file_io_0.close();
    if(hasHdf)
        file_io_hdf.close();

    //load spectra
    if (false == mda_io.load_spectra_volume(dataset_directory+"mda/"+dataset_file, detector_num, spectra_volume, hasNetcdf | hasHdf, extra_override_values, quantification_standard) )
    {
        std::cout<<"Error load spectra "<<dataset_directory+"mda/"+dataset_file<<std::endl;
        return false;
    }
    else
    {
        if(hasNetcdf)
        {
            std::ifstream file_io(dataset_directory+"flyXRF/"+tmp_dataset_file+"_2xfm3__0.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for(size_t i=0; i<spectra_volume->rows(); i++)
                {
                    full_filename = dataset_directory + "flyXRF/" + tmp_dataset_file + "_2xfm3__" + std::to_string(i) + ".nc";
                    std::cout<<"Loading file "<<full_filename<<std::endl;
                    netcdf_io.load_spectra_line(full_filename, detector_num, &(*spectra_volume)[i]);
                }
            }
            else
            {
                std::cout<<"Did not find netcdf files "<<dataset_directory+"flyXRF/"+tmp_dataset_file+"_2xfm3__0.nc"<<std::endl;
                //return false;
            }
        }
        else if (hasHdf)
        {
            hdf5_io.load_spectra_volume(dataset_directory+"flyXRF.h5/"+tmp_dataset_file+"_2xfm3__0.h5", detector_num, spectra_volume);
        }

    }
    if(save_scalers)
    {
        std::string str_detector_num = std::to_string(detector_num);
        std::string full_save_path = dataset_directory+"/img.dat/"+dataset_file+".h5"+str_detector_num;
        hdf5_io.save_scalers(full_save_path, detector_num, mda_io.get_scan_ptr(), extra_override_values);
    }
    mda_io.unload();

    return true;
}

// ----------------------------------------------------------------------------

bool load_and_integrate_spectra_volume(std::string dataset_directory,
                                       std::string dataset_file,
                                       data_struct::xrf::Spectra *integrated_spectra,
                                       size_t detector_num,
                                       std::unordered_map< std::string, std::string > *extra_override_values)
{
    //Dataset importer
    io::file::MDA_IO mda_io;
    io::file::HDF5_IO hdf5_io;
    io::file::NetCDF_IO netcdf_io;
    //data_struct::xrf::Detector detector;
    std::string tmp_dataset_file = dataset_file;
    bool ret_val = true;

    std::cout<<"load_spectra_volume_and_integrate(): Loading dataset "<<dataset_directory+"mda/"+dataset_file<<std::endl;

    //TODO: check if any analyized mda.h5 files are around to load. They should have integrated spectra saved already.

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    std::ifstream file_io_0(dataset_directory+"flyXRF/"+tmp_dataset_file+"_2xfm3__0.nc");
    std::ifstream file_io_hdf(dataset_directory+"flyXRF.h5/"+tmp_dataset_file+"_2xfm3__0.h5");
    bool hasNetcdf = file_io_0.is_open();
    bool hasHdf = file_io_hdf.is_open();
    if(hasNetcdf)
        file_io_0.close();
    if(hasHdf)
        file_io_hdf.close();

    //load spectra
    if (false == hasNetcdf && false == hasHdf)
    {
        data_struct::xrf::Spectra_Volume spectra_volume;
        ret_val = mda_io.load_spectra_volume(dataset_directory+"mda/"+dataset_file, detector_num, &spectra_volume, hasNetcdf | hasHdf, extra_override_values, nullptr);
        *integrated_spectra = spectra_volume.integrate();
        mda_io.unload();
    }
    else
    {
        int rank;
        int dims[10];
        dims[0] = 0;
        rank = mda_io.get_rank_and_dims(dataset_directory+"mda/"+dataset_file, &dims[0]);
        if(rank == 3)
            integrated_spectra->resize(dims[2]);
        else
            integrated_spectra->resize(2048);


        if(hasNetcdf)
        {
            std::ifstream file_io(dataset_directory+"flyXRF/"+tmp_dataset_file+"_2xfm3__0.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for(size_t i=0; i<dims[0]; i++)
                {
                    data_struct::xrf::Spectra_Line spectra_line;
                    full_filename = dataset_directory + "flyXRF/" + tmp_dataset_file + "_2xfm3__" + std::to_string(i) + ".nc";
                    std::cout<<"Loading file "<<full_filename<<std::endl;
                    netcdf_io.load_spectra_line(full_filename, detector_num, &spectra_line);
                    for(int k=0; k<spectra_line.size(); k++)
                    {
                        *integrated_spectra += spectra_line[k];
                    }
                }
            }
            else
            {
                std::cout<<"Did not find netcdf files "<<dataset_directory+"flyXRF/"+tmp_dataset_file+"_2xfm3__0.nc"<<std::endl;
                //return false;
            }
        }
        else if (hasHdf)
        {
            ret_val = hdf5_io.load_and_integrate_spectra_volume(dataset_directory+"flyXRF.h5/"+tmp_dataset_file+"_2xfm3__0.h5", detector_num, integrated_spectra);
        }

    }


    return ret_val;
}

// ----------------------------------------------------------------------------

 struct file_name_fit_params optimize_integrated_fit_params(std::string dataset_directory,
                                                            std::string  dataset_filename,
                                                            size_t detector_num)
{

    data_struct::xrf::Spectra spectra;
    //Fit parameter updates from user
    data_struct::xrf::Fit_Parameters override_fit_params;
    //Output of fits for elements specified
    data_struct::xrf::Fit_Element_Map_Dict elements_to_fit;
    std::unordered_map< std::string, std::string > extra_override_values;

    //return structure
    struct file_name_fit_params ret_struct;

    ret_struct.dataset_dir = dataset_directory;
    ret_struct.dataset_filename = dataset_filename;
    ret_struct.detector_num = detector_num;

    fitting::models::Gaussian_Model model;

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    //load override parameters
    load_override_params(dataset_directory, -1, &override_fit_params, &elements_to_fit, &extra_override_values);

    //load the quantification standard dataset
    if(false == load_and_integrate_spectra_volume(dataset_directory, dataset_filename, &spectra, detector_num, &extra_override_values) )
    {
        std::cout<<"Error in optimize_integrated_dataset loading dataset"<<dataset_filename<<" for detector"<<detector_num<<std::endl;
        ret_struct.success = false;
        return ret_struct;
    }

    energy_range.max = spectra.size() -1;


    //Fitting routines
    fitting::routines::Param_Optimized_Fit_Routine fit_routine;
    fit_routine.set_optimizer(&lmfit_optimizer);

    //reset model fit parameters to defaults
    model.reset_to_default_fit_params();
    //Update fit parameters by override values
    model.update_fit_params_values(override_fit_params);
    //Initialize the fit routine
    fit_routine.initialize(&model, &elements_to_fit, energy_range);
    //Fit the spectra saving the element counts in element_fit_count_dict
    ret_struct.fit_params = fit_routine.fit_spectra_parameters(&model, &spectra, &elements_to_fit);
    ret_struct.success = true;

    return ret_struct;

}

// ----------------------------------------------------------------------------

void generate_optimal_params(std::string dataset_directory,
                             std::vector<std::string> dataset_files,
                             ThreadPool* tp,
                             size_t detector_num_start,
                             size_t detector_num_end)
{
    bool first = true;
    std::queue<std::future<struct file_name_fit_params> > job_queue;

    std::vector<data_struct::xrf::Fit_Parameters> fit_params_avgs;
    fit_params_avgs.resize(detector_num_end - detector_num_start + 1);

    for(auto &itr : dataset_files)
    {
        for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
        {
            //data_struct::xrf::Fit_Parameters out_fitp;
            //out_fitp = optimize_integrated_fit_params(dataset_directory, itr, detector_num);
            job_queue.emplace( tp->enqueue(optimize_integrated_fit_params, dataset_directory, itr, detector_num) );
        }
    }


    while(!job_queue.empty())
    {
        auto ret = std::move(job_queue.front());
        job_queue.pop();
        struct file_name_fit_params f_struct = ret.get();
        if(f_struct.success)
        {
            if(first)
            {
                fit_params_avgs[f_struct.detector_num] = f_struct.fit_params;
            }
            else
            {
                fit_params_avgs[f_struct.detector_num].moving_average_with(f_struct.fit_params);
            }
            save_optimized_fit_params(f_struct);
        }
    }
    save_averaged_fit_params(dataset_directory, fit_params_avgs, detector_num_start, detector_num_end);

}

// ----------------------------------------------------------------------------

void process_dataset_file(std::string dataset_directory,
                          std::string dataset_file,
                          std::vector<Processing_Type> proc_types,
                          ThreadPool* tp,
                          std::vector<data_struct::xrf::Quantification_Standard>* quant_stand_list,
                          size_t detector_num_start,
                          size_t detector_num_end)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;

    fitting::models::Gaussian_Model model;

    //Output of fits for elements specified
    data_struct::xrf::Fit_Element_Map_Dict elements_to_fit;

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {

        data_struct::xrf::Quantification_Standard * quantification_standard = &(*quant_stand_list)[detector_num];

        data_struct::xrf::Fit_Parameters override_fit_params;

        std::string str_detector_num = std::to_string(detector_num);
        std::string full_save_path = dataset_directory+"/img.dat/"+dataset_file+".h5"+str_detector_num;

        //Spectra volume data
        data_struct::xrf::Spectra_Volume* spectra_volume = new data_struct::xrf::Spectra_Volume();

        //File job queue
        std::queue<std::future<bool> >* file_job_queue = new std::queue<std::future<bool> >();

        std::unordered_map< std::string, std::string > extra_override_values;
        extra_override_values.clear();
        //load override parameters
        load_override_params(dataset_directory, detector_num, &override_fit_params, &elements_to_fit, &extra_override_values);

        //add compton and coherant amp
        elements_to_fit.emplace(std::pair<std::string, data_struct::xrf::Fit_Element_Map*>(data_struct::xrf::STR_COMPTON_AMPLITUDE, new data_struct::xrf::Fit_Element_Map(data_struct::xrf::STR_COMPTON_AMPLITUDE, nullptr)) );
        elements_to_fit.emplace(std::pair<std::string, data_struct::xrf::Fit_Element_Map*>(data_struct::xrf::STR_COHERENT_SCT_AMPLITUDE, new data_struct::xrf::Fit_Element_Map(data_struct::xrf::STR_COHERENT_SCT_AMPLITUDE, nullptr)) );

        //load spectra volume
        if (false == load_spectra_volume(dataset_directory, dataset_file, spectra_volume, detector_num, &extra_override_values, nullptr, true) )
        {
            std::cout<<"Skipping detector "<<detector_num<<std::endl;
            continue;
        }

        energy_range.max = spectra_volume->samples_size() -1;

        for(auto proc_type : proc_types)
        {
            //Fit job queue
            std::queue<std::future<bool> >* fit_job_queue = new std::queue<std::future<bool> >();

            if (elements_to_fit.size() < 1)
            {
                std::cout<<"Error, no elements to fit. Check  maps_fit_parameters_override.txt0 - 3 exist"<<std::endl;
                return;
            }

            //Allocate memeory to save fit counts
            data_struct::xrf::Fit_Count_Dict  *element_fit_count_dict = generate_fit_count_dict(&elements_to_fit, spectra_volume->rows(), spectra_volume->cols());

            //Fitting models
            fitting::routines::Base_Fit_Routine *fit_routine = generate_fit_routine(proc_type);

            //for now we default to true to save iter count, in the future if we change the hdf5 layout we can store it per analysis.
            bool alloc_iter_count = true;

            //reset model fit parameters to defaults
            model.reset_to_default_fit_params();
            //Update fit parameters by override values
            model.update_fit_params_values(override_fit_params);

            if (alloc_iter_count)
            {
                //Allocate memeory to save number of fit iterations
                element_fit_count_dict->emplace(std::pair<std::string, data_struct::xrf::Fit_Counts_Array>(data_struct::xrf::STR_NUM_ITR, data_struct::xrf::Fit_Counts_Array()) );
                element_fit_count_dict->at(data_struct::xrf::STR_NUM_ITR).resize(spectra_volume->rows(), spectra_volume->cols());
            }

            //Initialize model
            fit_routine->initialize(&model, &elements_to_fit, energy_range);

            for(size_t i=0; i<spectra_volume->rows(); i++)
            {
                for(size_t j=0; j<spectra_volume->cols(); j++)
                {
                    //std::cout<< i<<" "<<j<<std::endl;
                    fit_job_queue->emplace( tp->enqueue(fit_single_spectra, fit_routine, &model, &(*spectra_volume)[i][j], &elements_to_fit, element_fit_count_dict, i, j) );
                }
            }

            //file_job_queue->emplace( tp->enqueue( save_results, full_save_path, save_loc, element_fit_count_dict, fit_job_queue) );
            save_results( full_save_path, save_loc_map.at(proc_type), element_fit_count_dict, fit_routine, fit_job_queue );
        }

        real_t energy_offset = 0.0;
        real_t energy_slope = 0.0;
        real_t energy_quad = 0.0;
        data_struct::xrf::Fit_Parameters fit_params = model.fit_parameters();
        if(fit_params.contains(fitting::models::STR_ENERGY_OFFSET))
        {
            energy_offset = fit_params[fitting::models::STR_ENERGY_OFFSET].value;
        }
        if(fit_params.contains(fitting::models::STR_ENERGY_SLOPE))
        {
            energy_slope = fit_params[fitting::models::STR_ENERGY_SLOPE].value;
        }
        if(fit_params.contains(fitting::models::STR_ENERGY_QUADRATIC))
        {
            energy_quad = fit_params[fitting::models::STR_ENERGY_QUADRATIC].value;
        }

        //tp->enqueue( save_volume, full_save_path, spectra_volume, file_job_queue );
        save_volume( full_save_path, quantification_standard, spectra_volume, energy_offset, energy_slope, energy_quad, file_job_queue);

    }

}


// ----------------------------------------------------------------------------

bool perform_quantification(std::string dataset_directory,
                            std::string quantification_info_file,
                            std::vector<Processing_Type> proc_types,
                            std::vector<data_struct::xrf::Quantification_Standard>* quant_stand_list,
                            size_t detector_num_start,
                            size_t detector_num_end)
{

    bool air_path = false;
    real_t detector_chip_thickness = 0.0;
    real_t beryllium_window_thickness = 0.0;
    real_t germanium_dead_layer = 0.0;
    real_t incident_energy = 10.0;

    fitting::models::Gaussian_Model model;

    data_struct::xrf::Element_Info* detector_element = data_struct::xrf::Element_Info_Map::inst()->get_element("Si");

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    std::string standard_file_name;
    std::unordered_map<std::string, real_t> element_standard_weights;
    //data_struct::xrf::Quantification_Standard* quantification_standard = &(*quant_stand_list)[0];

    if( load_quantification_standard(dataset_directory, quantification_info_file, &standard_file_name, &element_standard_weights) )
    {
        for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
        {
            data_struct::xrf::Quantification_Standard* quantification_standard = &(*quant_stand_list)[detector_num];
            quantification_standard->standard_filename(standard_file_name);
            for(auto& itr : element_standard_weights)
            {
                quantification_standard->append_element(itr.first, itr.second);
            }

            //Output of fits for elements specified
            std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> elements_to_fit;
            //Parameters for calibration curve
            data_struct::xrf::Fit_Parameters override_fit_params;
            std::unordered_map< std::string, std::string > extra_override_values;
            //load override parameters
            load_override_params(dataset_directory, detector_num, &override_fit_params, &elements_to_fit, &extra_override_values);

            if (extra_override_values.count(data_struct::xrf::STR_DETECTOR_ELEMENT) > 0)
            {
                // Get the element info class                                                                                 // detector element as string "Si" or "Ge" usually
                detector_element = data_struct::xrf::Element_Info_Map::inst()->get_element(extra_override_values.at(data_struct::xrf::STR_DETECTOR_ELEMENT));
            }
            if (extra_override_values.count(data_struct::xrf::STR_BE_WINDOW_THICKNESS) > 0)
            {
                beryllium_window_thickness = std::stof(extra_override_values.at(data_struct::xrf::STR_BE_WINDOW_THICKNESS));
            }
            if (extra_override_values.count(data_struct::xrf::STR_GE_DEAD_LAYER) > 0)
            {
                germanium_dead_layer = std::stof(extra_override_values.at(data_struct::xrf::STR_GE_DEAD_LAYER));
            }
            if (extra_override_values.count(data_struct::xrf::STR_DET_CHIP_THICKNESS) > 0)
            {
                detector_chip_thickness = std::stof(extra_override_values.at(data_struct::xrf::STR_DET_CHIP_THICKNESS));
            }

            if(override_fit_params.contains(data_struct::xrf::STR_COHERENT_SCT_ENERGY))
            {
                incident_energy = override_fit_params.at(data_struct::xrf::STR_COHERENT_SCT_ENERGY).value;
            }

            elements_to_fit.clear();
            for(auto& itr : element_standard_weights)
            {
                data_struct::xrf::Element_Info* e_info = data_struct::xrf::Element_Info_Map::inst()->get_element(itr.first);
                elements_to_fit[itr.first] = new data_struct::xrf::Fit_Element_Map(itr.first, e_info);
                elements_to_fit[itr.first]->init_energy_ratio_for_detector_element( detector_element );
            }

            data_struct::xrf::Spectra_Volume spectra_volume;
            //load the quantification standard dataset
            if(false == load_spectra_volume(dataset_directory, quantification_standard->standard_filename(), &spectra_volume, detector_num, &extra_override_values, quantification_standard, false) )
            {
                std::cout<<"Error in perform_quantification loading dataset for detector"<<detector_num<<std::endl;
                return false;
            }

            //First we integrate the spectra and get the elemental counts
            data_struct::xrf::Spectra integrated_spectra = spectra_volume.integrate();
            energy_range.max = integrated_spectra.size() -1;


            for(auto proc_type : proc_types)
            {

                //Fitting routines
                fitting::routines::Base_Fit_Routine *fit_routine = generate_fit_routine(proc_type);

                //reset model fit parameters to defaults
                model.reset_to_default_fit_params();
                //Update fit parameters by override values
                model.update_fit_params_values(override_fit_params);
                //Initialize the fit routine
                fit_routine->initialize(&model, &elements_to_fit, energy_range);
                //Fit the spectra
                std::unordered_map<std::string, real_t>counts_dict = fit_routine->fit_spectra(&model,
                                                                                              &integrated_spectra,
                                                                                              &elements_to_fit);

                if (fit_routine != nullptr)
                {
                    delete fit_routine;
                    fit_routine = nullptr;
                }

                for (auto& itr : elements_to_fit)
                {
                    counts_dict[itr.first] /= integrated_spectra.elapsed_lifetime();
                }
                quantification_standard->integrated_spectra(integrated_spectra);

                //save for each proc
                quantification_standard->quantifiy(&lmfit_optimizer,
                                                   save_loc_map.at(proc_type),
                                                  &counts_dict,
                                                  incident_energy,
                                                  detector_element,
                                                  air_path,
                                                  detector_chip_thickness,
                                                  beryllium_window_thickness,
                                                  germanium_dead_layer);

            }
        }
    }
    else
    {
        std::cout<<"Error loading quantification standard "<<quantification_info_file<<std::endl;
        return false;
    }

    return true;

}

// ----------------------------------------------------------------------------

bool compare_file_size (const file_name_size& first, const file_name_size& second)
{
    return ( first.total_rank_size > second.total_rank_size );
}

// ----------------------------------------------------------------------------

void sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string> *dataset_files)
{

    io::file::MDA_IO mda_io;

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

}

// ----------------------------------------------------------------------------

std::vector<std::string> find_all_dataset_files(std::string dataset_directory)
{
    std::vector<std::string> dataset_files;

    dataset_directory += "mda/";

    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dataset_directory.c_str())) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL)
        {
            std::string fname(ent->d_name);
            // check if extension is .mda
            if (fname.size() > 4)
            {
                if (fname.rfind(".mda") == fname.size() -4)
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
        std::cout<<"Error: could not open directory "<<dataset_directory<<std::endl;
    }


    return dataset_files;
}

// ----------------------------------------------------------------------------

void help()
{
    std::cout<<"Help: "<<std::endl;
    std::cout<<"Usage: xrf_maps [Options] [Fitting Routines] --dir [dataset directory] \n"<<std::endl;
    std::cout<<"Options: "<<std::endl;
    std::cout<<"--nthreads : <int> number of threads to use (default is all system threads) "<<std::endl;
    std::cout<<"--quantify-with : <standard.txt> File to use as quantification standard "<<std::endl;
    std::cout<<"--detector-range : <int:int> Start and end detector range. Defaults to 0:3 for 4 detector "<<std::endl;
    std::cout<<"--optimize-fit-override-params : Integrate the 8 largest mda datasets and fit with multiple params \n"<<std::endl;
    std::cout<<"Fitting Routines: "<<std::endl;
    std::cout<<"--roi : ROI "<<std::endl;
    std::cout<<"--roi_plus : SVD method "<<std::endl;
    std::cout<<"--nnls : Non-Negative Least Squares"<<std::endl;
    std::cout<<"--tails : Fit with multiple parameters "<<std::endl;
    std::cout<<"--matrix : Fit with locked parameters \n"<<std::endl;
    std::cout<<"Dataset: "<<std::endl;
    std::cout<<"--dir : Dataset directory "<<std::endl;
    std::cout<<"--files : Dataset files: comma (',') separated if multiple \n"<<std::endl;
    /*
    std::cout<<"Legacy Macros: "<<std::endl;
    std::cout<<"-A: Run roi and roi_plus"<<std::endl;
    std::cout<<"-B: Run optimize-fit-override-params"<<std::endl;
    std::cout<<"-C: Run roi and matrix"<<std::endl;
    */
    std::cout<<"Examples: "<<std::endl;
    std::cout<<"   Perform roi and matrix analysis on the directory /data/dataset1 "<<std::endl;
    std::cout<<"xrf_maps --roi --matrix --dir /data/dataset1 "<<std::endl;
    std::cout<<"   Perform roi and matrix analysis on the directory /data/dataset1 but only process scan1 and scan2 "<<std::endl;
    std::cout<<"xrf_maps --roi --matrix --dir /data/dataset1 --files scan1.mda,scan2.mda"<<std::endl;
    std::cout<<"   Perform roi, matrix, and nnls  analysis on the directory /data/dataset1, use maps_standard.txt information for quantification "<<std::endl;
    std::cout<<"xrf_maps --roi --matrix --nnls --quantify-with maps_standard.txt --dir /data/dataset1 "<<std::endl;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

    std::string dataset_dir;
    std::vector<std::string> dataset_files;
    std::vector<Processing_Type> proc_types;
    std::string quant_standard_filename = "";
    bool optimize_fit_override_params = false;

    //Default is to process detectors 0 through 3
    size_t detector_num_start = 0;
    size_t detector_num_end = 3;
    //ThreadPool tp(1);

    //Performance measure
    std::chrono::time_point<std::chrono::system_clock> start, end;

    //////// HENKE and ELEMENT INFO /////////////
	std::string element_csv_filename = "../reference/xrf_library.csv";
	std::string element_henke_filename = "../reference/henke.xdr";

    std::vector<data_struct::xrf::Quantification_Standard> quant_stand_list;

    Command_Line_Parser clp(argc, argv);

    if( clp.option_exists("-h") )
    {
       help();
       return 0;
    }


    size_t num_threads = std::thread::hardware_concurrency();
    if ( clp.option_exists("--nthreads") )
    {
        num_threads = std::stoi(clp.get_option("--nthreads"));
    }

    ThreadPool tp(num_threads);


    if ( clp.option_exists("--tails") )
    {
        proc_types.push_back(Processing_Type::GAUSS_TAILS);
    }
    if ( clp.option_exists("--matrix") )
    {
        proc_types.push_back(Processing_Type::GAUSS_MATRIX);
    }
    if ( clp.option_exists("--roi") )
    {
        proc_types.push_back(Processing_Type::ROI);
    }
    if ( clp.option_exists("--roi_plus") )
    {
        proc_types.push_back(Processing_Type::SVD);
    }
    if ( clp.option_exists("--nnls") )
    {
        proc_types.push_back(Processing_Type::NNLS);
    }

    if ( clp.option_exists("--quantify-with") )
    {
        quant_standard_filename = clp.get_option("--quantify-with");
    }

    if( clp.option_exists("--optimize-fit-override-params") )
    {
        optimize_fit_override_params = true;
    }

    //TODO: add --quantify-only option if you already did the fits and just want to add quantification

    if ( clp.option_exists("--detector-range") )
    {
        std::string detector_range_str = clp.get_option("--detector-range");
        if (detector_range_str.find(':') != std::string::npos )
        {
            // if we found a comma, split the string to get list of dataset files
            std::stringstream ss;
            ss.str(detector_range_str);
            std::string item;
            std::getline(ss, item, ':');
            detector_num_start = std::stoi(item);
            std::getline(ss, item, ':');
            detector_num_end = std::stoi(item);
        }
        else
        {
            detector_num_start = detector_num_end = std::stoi(detector_range_str);
        }
    }

    //quant_stand_list.resize( (detector_num_end - detector_num_start) + 1);
    //TODO: don't assume we can only have 4 detectors. The code also assumes quant_stand_list[0] is always detector 0. Maybe resize to detecotr_num_end?
    quant_stand_list.resize(4);

    dataset_dir = clp.get_option("--dir");
    if (dataset_dir.length() < 1)
    {
        help();
        return -1;
    }
    if (dataset_dir.back() != '/' && dataset_dir.back() != '\\')
    {
        dataset_dir += "/";
    }


    if (proc_types.size() == 0 && optimize_fit_override_params == false)
    {
        help();
        return -1;
    }

    std::string dset_file = clp.get_option("--files");
    if (dset_file.length() < 1)
    {
        // find all files in the dataset
        dataset_files = find_all_dataset_files(dataset_dir);
        if (dataset_files.size() == 0)
        {
            std::cout<<"Error: No mda files found in dataset directory "<<dataset_dir<<std::endl;
            return -1;
        }
    }
    else if (dset_file.find(',') != std::string::npos )
    {
        // if we found a comma, split the string to get list of dataset files
        std::stringstream ss;
        ss.str(dset_file);
        std::string item;
        while (std::getline(ss, item, ','))
        {
            dataset_files.push_back(item);
        }
    }
    else
    {
        dataset_files.push_back(dset_file);
    }

    std::cout << "Processing detectors " << detector_num_start << " - "<< detector_num_end <<" \n";

    start = std::chrono::system_clock::now();

    //load element information
    if(false == load_element_info(element_henke_filename, element_csv_filename, data_struct::xrf::Element_Info_Map::inst()))
    {
        std::cout<<"Error loading element information: "<<std::endl;
        return -1;
    }

    if(optimize_fit_override_params)
    {

        if (dataset_files.size() == 0)
        {
            std::cout<<"Error: No mda files found in dataset directory "<<dataset_dir<<std::endl;
            return -1;
        }
        sort_dataset_files_by_size(dataset_dir, &dataset_files);

        generate_optimal_params(dataset_dir, dataset_files, &tp, detector_num_start, detector_num_end);
    }

    if(proc_types.size() > 0)
    {
        if (quant_standard_filename.length() > 0)
        {
            perform_quantification(dataset_dir, quant_standard_filename, proc_types, &quant_stand_list, detector_num_start, detector_num_end);
        }

        for(std::string dataset_file : dataset_files)
        {
            process_dataset_file(dataset_dir, dataset_file, proc_types, &tp, &quant_stand_list, detector_num_start, detector_num_end);
        }
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "\n\n ---- Total elapsed time: " << elapsed_seconds.count() << "s -----\n\n";

    return 0;
    //return a.exec();
}
