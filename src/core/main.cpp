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

#include "element_info.h"

#include "aps_fit_params_import.h"

#include "roi_model.h"
#include "svd_model.h"
#include "nnls_model.h"

#include "lmfit_optimizer.h"
#include "mpfit_optimizer.h"

#include "fit_element_map.h"

#include "quantification_standard.h"

//#include "vtk_graph.h"

#include "command_line_parser.h"

using namespace std::placeholders; //for _1, _2,

// ----------------------------------------------------------------------------

enum Processing_Types { ROI=1 , GAUSS_TAILS=2, GAUSS_MATRIX=4, SVD=8, NNLS=16 };

// ----------------------------------------------------------------------------

bool load_spectra_volume(std::string dataset_directory,
                         std::string dataset_file,
                         data_struct::xrf::Spectra_Volume *spectra_volume,
                         size_t detector_num,
                         std::unordered_map< std::string, std::string > *extra_override_values);

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
    hdf5_io.save_element_fits(full_path, save_loc, element_counts);

    delete job_queue;
    delete element_counts;

    return true;
}

// ----------------------------------------------------------------------------

bool save_volume(std::string full_path,
                 data_struct::xrf::Spectra_Volume *spectra_volume,
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
    hdf5_io.save_spectra_volume(full_path, "mca_arr", spectra_volume);

    delete job_queue;
    delete spectra_volume;

    return true;
}

// ----------------------------------------------------------------------------
/*
void fit_and_plot_volume(fitting::models::Base_Model *model,
                         data_struct::xrf::Fit_Parameters fit_params,
                         data_struct::xrf::Spectra_Volume *spectra_volume,
                         data_struct::xrf::Quantification_Standard* calibration,
                         std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> *elements_to_fit)
{
    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = spectra_volume->samples_size() - 1;

    for(size_t i=0; i<spectra_volume->rows(); i++)
    {
        for(size_t j=0; j<spectra_volume->cols(); j++)
        {
            std::cout<< i<<" "<<j<<std::endl;
            data_struct::xrf::Spectra *spectra = &(*spectra_volume)[i][j];

            data_struct::xrf::Fit_Parameters result_fits = model->fit_spectra(fit_params, spectra, calibration, elements_to_fit, i, j);

            //fitp.print_non_fixed();
            data_struct::xrf::Spectra spectra_model = model->model_spectrum(&result_fits, spectra, calibration, elements_to_fit, energy_range);
//            spectra_model.log10();
//            spectra_model.set_min_zero();
//            spectra->log10();
            //visual::PlotSpectras(*spectra, spectra_model);

        }
    }

}
*/
// ----------------------------------------------------------------------------
/*
bool fit_volume(fitting::models::Base_Model *model,
                data_struct::xrf::Fit_Parameters fit_params,
                data_struct::xrf::Spectra_Volume *spectra_volume,
                data_struct::xrf::Quantification_Standard* calibration,
                std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> *elements_to_fit,
                ThreadPool *tp)
{

    std::queue<std::future<data_struct::xrf::Fit_Parameters> > res_q;
    for(size_t i=0; i<spectra_volume->rows(); i++)
    //for(size_t i=0; i<1; i++) //debug
    {
        for(size_t j=0; j<spectra_volume->cols(); j++)
        //for(size_t j=0; j<1; j++) //debug
        {
            std::cout<< i<<" "<<j<<std::endl;
            data_struct::xrf::Spectra *spectra = &(*spectra_volume)[i][j];

            res_q.emplace( tp->enqueue([](fitting::models::Base_Model* model, data_struct::xrf::Fit_Parameters fit_params, data_struct::xrf::Spectra* spectra, data_struct::xrf::Quantification_Standard* calibration, std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> *elements, size_t i, size_t j)
            {return model->fit_spectra(fit_params, spectra, calibration, elements, i, j);}, model, fit_params, spectra, calibration, elements_to_fit, i, j) );
            //std::function<void()> bound_f = std::bind(&fitting::models::Base_Model::fit_spectra, &model, fit_params, spectra, calibration, elements_to_fit, i, j);
            //res_q.emplace( tp->enqueue(  bound_f ) );
        }
    }

    return true;
}
*/
// ----------------------------------------------------------------------------

data_struct::xrf::Fit_Count_Dict* generate_fit_count_dict(data_struct::xrf::Fit_Element_Map_Dict *elements_to_fit, size_t width, size_t height )
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

bool fit_single_spectra(fitting::models::Base_Model *model,
                        data_struct::xrf::Fit_Parameters fit_params,
                        const data_struct::xrf::Spectra * const spectra,
                        const data_struct::xrf::Detector * const detector,
                        const data_struct::xrf::Fit_Element_Map_Dict * const elements_to_fit,
                        data_struct::xrf::Fit_Count_Dict * out_fit_counts,
                        size_t i,
                        size_t j)
{
    model->fit_spectra(fit_params, spectra, detector, elements_to_fit, out_fit_counts, i, j);
    return true;
}

// ----------------------------------------------------------------------------

bool load_quantification_standard(std::string dataset_directory,
                                  std::string quantification_info_file,
                                  std::vector<Processing_Types> proc_types,
                                  data_struct::xrf::Quantification_Standard * quantification_standard,
                                  size_t detector_num_start,
                                  size_t detector_num_end)
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
                    quantification_standard->standard_filename(standard_filename);
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
                    quantification_standard->append_element_weight(element_names[i], element_weights[i]);
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
                          data_struct::xrf::Detector *detector,
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
                                        detector,
                                        elements_to_fit,
                                        extra_override_values))
    {
        std::cout<<"Error loading fit param override file: "<<std::endl;
        return false;
    }
    else
    {

        if( fit_params->contains(data_struct::xrf::STR_E_OFFSET) )
            detector->energy_offset((*fit_params)[data_struct::xrf::STR_E_OFFSET].value);
        if( fit_params->contains(data_struct::xrf::STR_E_LINEAR) )
            detector->energy_slope((*fit_params)[data_struct::xrf::STR_E_LINEAR].value);
        if( fit_params->contains(data_struct::xrf::STR_E_QUADRATIC) )
            detector->energy_quadratic((*fit_params)[data_struct::xrf::STR_E_QUADRATIC].value);

        std::cout<<"Elements to fit: "<<std::endl;
        for (auto el_itr : *elements_to_fit)
        {
            std::cout<<el_itr.first<<" ";
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
                         std::unordered_map< std::string, std::string > *extra_override_values)
{

    //Dataset importer
    io::file::MDA_IO mda_io;
    io::file::HDF5_IO hdf5_io;
    io::file::NetCDF_IO netcdf_io;
    data_struct::xrf::Detector detector;
    std::string tmp_dataset_file = dataset_file;

    std::cout<<"Loading dataset "<<dataset_directory+"mda/"+dataset_file<<std::endl;

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    std::ifstream file_io_0(dataset_directory+"flyXRF/"+tmp_dataset_file+"_2xfm3__0.nc");
    bool hasNetcdf = file_io_0.is_open();
    if(hasNetcdf)
        file_io_0.close();

    //load spectra
    if (false == mda_io.load_spectra_volume(dataset_directory+"mda/"+dataset_file, detector_num, &detector, spectra_volume, hasNetcdf, extra_override_values) )
    {
        std::cout<<"Error load spectra "<<dataset_directory+"mda/"+dataset_file<<std::endl;
        return false;
    }
    else
    {
        //check to see if netcdf or hdf5 exist (fly scans)
        dataset_file = dataset_file.substr(0, dataset_file.size()-4);
        std::ifstream file_io(dataset_directory+"flyXRF/"+dataset_file+"_2xfm3__0.nc");
        if(file_io.is_open())
        {
            file_io.close();
            std::string full_filename;
            for(size_t i=0; i<spectra_volume->rows(); i++)
            {
                full_filename = dataset_directory + "flyXRF/" + dataset_file + "_2xfm3__" + std::to_string(i) + ".nc";
                std::cout<<"Loading file "<<full_filename<<std::endl;
                netcdf_io.load_spectra_line(full_filename, detector_num, &(*spectra_volume)[i]);
            }
        }
        else
        {
            std::cout<<"Did not find netcdf files "<<dataset_directory+"flyXRF/"+dataset_file+"_2xfm3__0.nc"<<std::endl;
            //return false;
        }

    }
    mda_io.unload(); //TODO optimize so we load straight to memeory

    return true;
}

// ----------------------------------------------------------------------------

bool process_integrated_dataset(std::string dataset_directory,
                                std::vector<Processing_Types> proc_types,
                                data_struct::xrf::Quantification_Standard * quantification_standard,
                                size_t detector_num_start,
                                size_t detector_num_end)
{

    //Spectra of quantification datasaet
    data_struct::xrf::Spectra_Volume spectra_volume;
    //Fit parameter updates from user
    data_struct::xrf::Fit_Parameters override_fit_params;
    //Element detector
    data_struct::xrf::Detector detector;
    //Output of fits for elements specified
    data_struct::xrf::Fit_Element_Map_Dict elements_to_fit;
    std::unordered_map< std::string, std::string > extra_override_values;
    //Optimizers for fitting models
    fitting::optimizers::LMFit_Optimizer lmfit_optimizer;

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {


        //load override parameters
        load_override_params(dataset_directory, detector_num, &override_fit_params, &detector, &elements_to_fit, &extra_override_values);

        //load the quantification standard dataset
        if(false == load_spectra_volume(dataset_directory, quantification_standard->standard_filename(), &spectra_volume, detector_num, &extra_override_values) )
        {
            std::cout<<"Error in process_integrated_dataset loading dataset for detector"<<detector_num<<std::endl;
            continue;
        }

        //First we integrate the spectra and get the elemental counts
        data_struct::xrf::Spectra integrated_spectra = spectra_volume.integrate();
        energy_range.max = integrated_spectra.size() -1;

        for(auto proc_type : proc_types)
        {
            //Fitting models
            fitting::models::Base_Model *model = nullptr;

            //for now we default to true to save iter count, in the future if we change the hdf5 layout we can store it per analysis.
            bool alloc_iter_count = true;
            switch(proc_type)
            {
            case Processing_Types::GAUSS_TAILS:
                model = new fitting::models::Gauss_Tails_Model();
                ((fitting::models::Gauss_Tails_Model*)model)->set_optimizer(&lmfit_optimizer);
                //save_loc = "XRF_tails_fits";
                alloc_iter_count = true;
                break;
            case Processing_Types::GAUSS_MATRIX:
                model = new fitting::models::Gauss_Matrix_Model();
                ((fitting::models::Gauss_Matrix_Model*)model)->set_optimizer(&lmfit_optimizer);
                //save_loc = "XRF_fits";
                alloc_iter_count = true;
                break;
            case Processing_Types::ROI:
                model = new fitting::models::ROI_Model();
                //save_loc = "XRF_roi";
                break;
            case Processing_Types::SVD:
                model = new fitting::models::SVD_Model();
                //save_loc = "XRF_roi_plus";
                break;
            case Processing_Types::NNLS:
                model = new fitting::models::NNLS_Model();
                //save_loc = "XRF_nnls";
                alloc_iter_count = true;
                break;
            default:
                continue;
            }

            data_struct::xrf::Fit_Count_Dict  *element_fit_count_dict = generate_fit_count_dict(&elements_to_fit, 1, 1);

            //Get required fit parameters from model
            data_struct::xrf::Fit_Parameters fit_params = model->get_fit_parameters();
            //Update fit parameters by override values
            fit_params.update_values(override_fit_params);

            model->initialize(&fit_params, &detector, &elements_to_fit, energy_range);

            data_struct::xrf::Fit_Parameters result_fits = model->fit_spectra(fit_params, &integrated_spectra, &detector, &elements_to_fit, element_fit_count_dict);

            //get new detector energy offset, slope, and quadratic

            //hdf5_io.save_element_fits(full_path, save_loc, element_counts);


            if (model != nullptr)
            {
                delete model;
                model = nullptr;
            }
        }
    }// end for
    return true;

}

// ----------------------------------------------------------------------------

void process_dataset_file(std::string dataset_directory,
                          std::string dataset_file,
                          std::vector<Processing_Types> proc_types,
                          ThreadPool* tp,
                          size_t detector_num_start,
                          size_t detector_num_end)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;

    //Fitting models
    fitting::models::Base_Model *model = nullptr;
    fitting::models::Gauss_Tails_Model gauss_tails_model;
    fitting::models::Gauss_Matrix_Model gauss_matrix_model;
    fitting::models::ROI_Model roi_model;
    fitting::models::SVD_Model svd_model;
    fitting::models::NNLS_Model nnls_model;

    //Optimizers for fitting models
    fitting::optimizers::LMFit_Optimizer lmfit_optimizer;

    //Fitting parameters
    data_struct::xrf::Fit_Parameters fit_params;

    //Element detector
    data_struct::xrf::Detector detector;

    //Output of fits for elements specified
    std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> elements_to_fit;

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    std::string save_loc;

    gauss_tails_model.set_optimizer(&lmfit_optimizer);
    gauss_matrix_model.set_optimizer(&lmfit_optimizer);

    //bool first_time = true;


    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {

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
        load_override_params(dataset_directory, detector_num, &override_fit_params, &detector, &elements_to_fit, &extra_override_values);

        //load spectra volume
        if (false == load_spectra_volume(dataset_directory, dataset_file, spectra_volume, detector_num, &extra_override_values) )
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

            //for now we default to true to save iter count, in the future if we change the hdf5 layout we can store it per analysis.
            bool alloc_iter_count = true;
            switch(proc_type)
            {
            case Processing_Types::GAUSS_TAILS:
                model = &gauss_tails_model;
                save_loc = "XRF_tails_fits";
                alloc_iter_count = true;
                break;
            case Processing_Types::GAUSS_MATRIX:
                model = &gauss_matrix_model;
                save_loc = "XRF_fits";
                alloc_iter_count = true;
                break;
            case Processing_Types::ROI:
                model = &roi_model;
                save_loc = "XRF_roi";
                break;
            case Processing_Types::SVD:
                model = &svd_model;
                save_loc = "XRF_roi_plus";
                break;
            case Processing_Types::NNLS:
                model = &nnls_model;
                save_loc = "XRF_nnls";
                alloc_iter_count = true;
                break;
            default:
                continue;
            }

            //Get required fit parameters from model
            fit_params = model->get_fit_parameters();
            //Update fit parameters by override values
            fit_params.update_values(override_fit_params);

            if (alloc_iter_count)
            {
                //Allocate memeory to save number of fit iterations
                element_fit_count_dict->emplace(std::pair<std::string, data_struct::xrf::Fit_Counts_Array>(data_struct::xrf::STR_NUM_ITR, data_struct::xrf::Fit_Counts_Array()) );
                element_fit_count_dict->at(data_struct::xrf::STR_NUM_ITR).resize(spectra_volume->rows(), spectra_volume->cols());
                //add num iters as a fit param
                fit_params.add_parameter(data_struct::xrf::STR_NUM_ITR, data_struct::xrf::Fit_Param(data_struct::xrf::STR_NUM_ITR, 0.0, std::numeric_limits<real_t>::max(), 0.0, 1.0, data_struct::xrf::FIXED));
            }

            //Initialize model
            model->initialize(&fit_params, &detector, &elements_to_fit, energy_range);

            for(size_t i=0; i<spectra_volume->rows(); i++)
            {
                for(size_t j=0; j<spectra_volume->cols(); j++)
                {
                    //std::cout<< i<<" "<<j<<std::endl;
                    fit_job_queue->emplace( tp->enqueue(fit_single_spectra, model, fit_params, &(*spectra_volume)[i][j], &detector, &elements_to_fit, element_fit_count_dict, i, j) );
                }
            }

            //file_job_queue->emplace( tp->enqueue( save_results, full_save_path, save_loc, element_fit_count_dict, fit_job_queue) );
            save_results( full_save_path, save_loc, element_fit_count_dict, fit_job_queue );

        }

        //tp->enqueue( save_volume, full_save_path, spectra_volume, file_job_queue );
        save_volume( full_save_path, spectra_volume, file_job_queue);

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
    std::cout<<"Usage: xrf_mapper [Options] [Fitting models] --dir [dataset directory] \n"<<std::endl;
    std::cout<<"Options: "<<std::endl;
    std::cout<<"--nthreads : <int> number of threads to use (default is all system threads) "<<std::endl;
    std::cout<<"--quantify-with : <standard.txt> File to use as quantification standard "<<std::endl;
    std::cout<<"--detector-range : <int:int> Start and end detector range. Defaults to 0:3 for 4 detector "<<std::endl;
    std::cout<<"Fitting models: "<<std::endl;
    std::cout<<"--roi : ROI "<<std::endl;
    std::cout<<"--roi_plus : SVD method "<<std::endl;
    std::cout<<"--nnls : Non-Negative Least Squares"<<std::endl;
    std::cout<<"--tails : Fit with multiple parameters "<<std::endl;
    std::cout<<"--matrix : Fit with locked parameters \n"<<std::endl;
    std::cout<<"Dataset: "<<std::endl;
    std::cout<<"--dir : Dataset directory "<<std::endl;
    std::cout<<"--files : Dataset files: comma (',') separated if multiple \n"<<std::endl;
    std::cout<<"Examples: "<<std::endl;
    std::cout<<"   Perform roi and matrix analysis on the directory /data/dataset1 "<<std::endl;
    std::cout<<"xrf_mapper --roi --matrix --dir /data/dataset1 "<<std::endl;
    std::cout<<"   Perform roi and matrix analysis on the directory /data/dataset1 but only process scan1 and scan2 "<<std::endl;
    std::cout<<"xrf_mapper --roi --matrix --dir /data/dataset1 --files scan1.mda,scan2.mda"<<std::endl;
    std::cout<<"   Perform roi, matrix, and nnls  analysis on the directory /data/dataset1, use maps_standard.txt information for quantification "<<std::endl;
    std::cout<<"xrf_mapper --roi --matrix --nnls --quantify-with maps_standard.txt --dir /data/dataset1 "<<std::endl;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

    std::string dataset_dir;
    std::vector<std::string> dataset_files;
    std::vector<Processing_Types> proc_types;
    std::string quant_standard_filename = "";

    //Default is to process detectors 0 through 3
    size_t detector_num_start = 0;
    size_t detector_num_end = 3;
    //ThreadPool tp(1);

    //Performance measure
    std::chrono::time_point<std::chrono::system_clock> start, end;

    //////// HENKE and ELEMENT INFO /////////////
	std::string element_csv_filename = "../reference/xrf_library.csv";
	std::string element_henke_filename = "../reference/henke.xdr";


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
        proc_types.push_back(Processing_Types::GAUSS_TAILS);
    }
    if ( clp.option_exists("--matrix") )
    {
        proc_types.push_back(Processing_Types::GAUSS_MATRIX);
    }
    if ( clp.option_exists("--roi") )
    {
        proc_types.push_back(Processing_Types::ROI);
    }
    if ( clp.option_exists("--roi_plus") )
    {
        proc_types.push_back(Processing_Types::SVD);
    }
    if ( clp.option_exists("--nnls") )
    {
        proc_types.push_back(Processing_Types::NNLS);
    }

    if ( clp.option_exists("--quantify-with") )
    {
        quant_standard_filename = clp.get_option("--quantify-with");
    }


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


    if (proc_types.size() == 0)
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

    if (quant_standard_filename.length() > 0)
    {
        data_struct::xrf::Quantification_Standard quantification_standard;
        if (load_quantification_standard(dataset_dir, quant_standard_filename, proc_types, &quantification_standard, detector_num_start, detector_num_end) )
        {
            process_integrated_dataset(dataset_dir, proc_types, &quantification_standard, detector_num_start, detector_num_end);
        }
    }

    for(std::string dataset_file : dataset_files)
    {
        process_dataset_file(dataset_dir, dataset_file, proc_types, &tp, detector_num_start, detector_num_end);
    }

    //process_integrated_dataset(dataset_dir, dataset_file, Processing_Types::GAUSS_TAILS);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "\n\n ---- Total elapsed time: " << elapsed_seconds.count() << "s -----\n\n";

    return 0;
    //return a.exec();
}
