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

void fit_and_plot_volume(fitting::models::Base_Model *model,
                         data_struct::xrf::Fit_Parameters fit_params,
                         data_struct::xrf::Spectra_Volume *spectra_volume,
                         data_struct::xrf::Calibration_Standard* calibration,
                         std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> *elements_to_fit)
{
    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = spectra_volume->samples_size() - 1;
/*
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
    */
}

// ----------------------------------------------------------------------------

bool fit_volume(fitting::models::Base_Model *model,
                data_struct::xrf::Fit_Parameters fit_params,
                data_struct::xrf::Spectra_Volume *spectra_volume,
                data_struct::xrf::Calibration_Standard* calibration,
                std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> *elements_to_fit,
                ThreadPool *tp)
{
/*
    std::queue<std::future<data_struct::xrf::Fit_Parameters> > res_q;
    for(size_t i=0; i<spectra_volume->rows(); i++)
    //for(size_t i=0; i<1; i++) //debug
    {
        for(size_t j=0; j<spectra_volume->cols(); j++)
        //for(size_t j=0; j<1; j++) //debug
        {
            std::cout<< i<<" "<<j<<std::endl;
            data_struct::xrf::Spectra *spectra = &(*spectra_volume)[i][j];

            res_q.emplace( tp->enqueue([](fitting::models::Base_Model* model, data_struct::xrf::Fit_Parameters fit_params, data_struct::xrf::Spectra* spectra, data_struct::xrf::Calibration_Standard* calibration, std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> *elements, size_t i, size_t j)
            {return model->fit_spectra(fit_params, spectra, calibration, elements, i, j);}, model, fit_params, spectra, calibration, elements_to_fit, i, j) );
            //std::function<void()> bound_f = std::bind(&fitting::models::Base_Model::fit_spectra, &model, fit_params, spectra, calibration, elements_to_fit, i, j);
            //res_q.emplace( tp->enqueue(  bound_f ) );
        }
    }
*/
    return true;
}

// ----------------------------------------------------------------------------

bool fit_single_spectra(fitting::models::Base_Model *model,
                        data_struct::xrf::Fit_Parameters fit_params,
                        const data_struct::xrf::Spectra * const spectra,
                        const data_struct::xrf::Calibration_Standard * const calibration,
                        const data_struct::xrf::Fit_Element_Map_Dict * const elements_to_fit,
                        data_struct::xrf::Fit_Count_Dict * out_fit_counts,
                        size_t i,
                        size_t j)
{
    model->fit_spectra(fit_params, spectra, calibration, elements_to_fit, out_fit_counts, i, j);
    return true;
}


// ----------------------------------------------------------------------------

bool load_override_params(std::string dataset_directory,
                          int detector_num,
                          data_struct::xrf::Fit_Parameters *fit_params,
                          data_struct::xrf::Calibration_Standard *calibration,
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

        if( fit_params->contains(data_struct::xrf::STR_E_OFFSET) )
            calibration->offset((*fit_params)[data_struct::xrf::STR_E_OFFSET].value);
        if( fit_params->contains(data_struct::xrf::STR_E_LINEAR) )
            calibration->slope((*fit_params)[data_struct::xrf::STR_E_LINEAR].value);
        if( fit_params->contains(data_struct::xrf::STR_E_QUADRATIC) )
            calibration->quad((*fit_params)[data_struct::xrf::STR_E_QUADRATIC].value);

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

void process_integrated_dataset(std::string dataset_directory, std::string dataset_file, Processing_Types proc_type)
{
    //for debugging
    bool plot = true;

    //Performance measure
    std::chrono::time_point<std::chrono::system_clock> start, end;

    //Fitting parameters
    data_struct::xrf::Fit_Parameters fit_params;
    data_struct::xrf::Fit_Parameters override_fit_params;

    //Standard calibration and importer for APS datasets
    data_struct::xrf::Calibration_Standard calibration;

    //Spectra volume data
    data_struct::xrf::Spectra_Volume spectra_volume;

    //Base fitting model
    fitting::models::Base_Model *model = nullptr;

    //Optimizers for fitting models
    fitting::optimizers::LMFit_Optimizer lmfit_optimizer;

    //Output of fits for elements specified
    std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> elements_to_fit;

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    std::unordered_map< std::string, std::string > extra_override_values;

    //load override parameters
    load_override_params(dataset_directory, 0, &override_fit_params, &calibration, &elements_to_fit, &extra_override_values);

    //load spectra volume
    if (false == load_spectra_volume(dataset_directory, dataset_file, &spectra_volume, 0, nullptr) ) // todo pass in extra override dict
    {
        return;
    }

    if ( proc_type == Processing_Types::GAUSS_TAILS )
    {
        model = new fitting::models::Gauss_Tails_Model();
        //Set the optimization technique to use
        ((fitting::models::Gauss_Tails_Model*) model)->set_optimizer(&lmfit_optimizer);
    }
    else if ( proc_type == Processing_Types::GAUSS_MATRIX )
    {
        model = new fitting::models::Gauss_Matrix_Model();
        //Set the optimization technique to use
        ((fitting::models::Gauss_Matrix_Model*) model)->set_optimizer(&lmfit_optimizer);
    }
    else if ( proc_type == Processing_Types::ROI )
    {
        model = new fitting::models::ROI_Model();
    }
    else if ( proc_type == Processing_Types::SVD )
    {
        model = new fitting::models::SVD_Model();
    }


    //Get required fit parameters from model
    fit_params = model->get_fit_parameters();
    //update fit parameters by override values
    fit_params.update_values(override_fit_params);

    energy_range.max = spectra_volume.samples_size() - 1;

    //Initialize model
    model->initialize(&fit_params, &calibration, &elements_to_fit, energy_range);

    //Allocate memeory to save fit counts
//    for(auto& e_itr : elements_to_fit)
//    {
//        e_itr.second.resize(1,1);
//    }

    data_struct::xrf::Spectra integrated_spectra = spectra_volume.integrate();

    fit_params.print();

    start = std::chrono::system_clock::now();

    data_struct::xrf::Fit_Parameters result_fits = model->fit_spectra(fit_params, &integrated_spectra, &calibration, &elements_to_fit, 0, 0);

    end = std::chrono::system_clock::now();

    //TODO: Save fit parameters

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

    result_fits.print();

    if(plot)
    {
        /*
        integrated_spectra = std::log10(integrated_spectra);
        integrated_spectra.set_min_zero();

        data_struct::xrf::Spectra spectra_model = model->model_spectrum(&result_fits, &integrated_spectra, &calibration, &elements_to_fit, energy_range);

        spectra_model = std::log10(spectra_model);
        spectra_model.set_min_zero();

        PlotSpectras(integrated_spectra, spectra_model);
        */
    }


    if (model != nullptr)
    {
        delete model;
        model = nullptr;
    }
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

    //Standard calibration and importer for APS datasets
    data_struct::xrf::Calibration_Standard calibration;

    //Output of fits for elements specified
    std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> elements_to_fit;

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    std::string save_loc;

    gauss_tails_model.set_optimizer(&lmfit_optimizer);
    gauss_matrix_model.set_optimizer(&lmfit_optimizer);

    bool first_time = true;


    for(size_t detector_num = detector_num_start; detector_num < detector_num_end; detector_num++)
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
        load_override_params(dataset_directory, detector_num, &override_fit_params, &calibration, &elements_to_fit, &extra_override_values);

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

            data_struct::xrf::Fit_Count_Dict  *element_fit_count_dict = new data_struct::xrf::Fit_Count_Dict();

            if (elements_to_fit.size() < 1)
            {
                std::cout<<"Error, no elements to fit. Check  maps_fit_parameters_override.txt0 - 3 exist"<<std::endl;
                return;
            }
            //Allocate memeory to save fit counts
            for(auto& e_itr : elements_to_fit)
            {

                element_fit_count_dict->emplace(std::pair<std::string, data_struct::xrf::Fit_Counts_Array>(e_itr.first, data_struct::xrf::Fit_Counts_Array()) );
                element_fit_count_dict->at(e_itr.first).resize(spectra_volume->rows(), spectra_volume->cols());
            }

            switch(proc_type)
            {
            case Processing_Types::GAUSS_TAILS:
                model = &gauss_tails_model;
                save_loc = "XRF_tails_fits";
                break;
            case Processing_Types::GAUSS_MATRIX:
                model = &gauss_matrix_model;
                save_loc = "XRF_fits";
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
                break;

            }

            //Get required fit parameters from model
            fit_params = model->get_fit_parameters();
            //Update fit parameters by override values
            fit_params.update_values(override_fit_params);

            if ( proc_type == Processing_Types::GAUSS_MATRIX )
            {
                //Set fit parameters to fixed so we only fit elements counts
                fit_params.set_all(data_struct::xrf::E_Bound_Type::FIXED);
            }

            //Initialize model
            model->initialize(&fit_params, &calibration, &elements_to_fit, energy_range);

            for(size_t i=0; i<spectra_volume->rows(); i++)
            {
                for(size_t j=0; j<spectra_volume->cols(); j++)
                {
                    //std::cout<< i<<" "<<j<<std::endl;
                    fit_job_queue->emplace( tp->enqueue(fit_single_spectra, model, fit_params, &(*spectra_volume)[i][j], &calibration, &elements_to_fit, element_fit_count_dict, i, j) );
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
    std::cout<<"xrf_mapper --roi --matrix --dir /data/dataset1 "<<std::endl;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

    std::string dataset_dir;
    std::vector<std::string> dataset_files;
    std::vector<Processing_Types> proc_types;

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


    start = std::chrono::system_clock::now();

    //load element information
    if(false == load_element_info(element_henke_filename, element_csv_filename, data_struct::xrf::Element_Info_Map::inst()))
    {
        std::cout<<"Error loading element information: "<<std::endl;
        return -1;
    }

    for(std::string dataset_file : dataset_files)
    {
        process_dataset_file(dataset_dir, dataset_file, proc_types, &tp, 0, 4);
    }

    //process_integrated_dataset(dataset_dir, dataset_file, Processing_Types::GAUSS_TAILS);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "\n\n ---- Total elapsed time: " << elapsed_seconds.count() << "s -----\n\n";

    return 0;
    //return a.exec();
}
