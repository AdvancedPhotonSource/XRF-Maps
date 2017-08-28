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


#include "core/process_streaming.h"
#include "core/process_whole.h"


// ----------------------------------------------------------------------------

void help()
{
    logit_s<<"Help: "<<std::endl;
    logit_s<<"Usage: xrf_maps [Options] [Fitting Routines] --dir [dataset directory] \n"<<std::endl;
    logit_s<<"Options: "<<std::endl;
    logit_s<<"--nthreads : <int> number of threads to use (default is all system threads) "<<std::endl;
    logit_s<<"--quantify-with : <standard.txt> File to use as quantification standard "<<std::endl;
    logit_s<<"--detector-range : <int:int> Start and end detector range. Defaults to 0:3 for 4 detector "<<std::endl;
    logit_s<<"--generate-avg-h5 : Generate .h5 file which is the average of all detectors .h50 - h.53 or range specified. "<<std::endl;
//    logit_s<<"--add-exchange : <us:ds:sr> Add exchange group into hdf5 file with normalized data.\n";
//    logit_s<<"    us = upstream ion chamber\n";
//    logit_s<<"    ds = downstream ion chamber\n";
//    logit_s<<"    sr = sr current. "<<std::endl;
    logit_s<<"--quick-and-dirty : Integrate the detector range into 1 spectra. "<<std::endl;
    logit_s<<"--optimize-fit-override-params : <int> Integrate the 8 largest mda datasets and fit with multiple params\n"<<
               "  1 = matrix batch fit\n  2 = batch fit without tails\n  3 = batch fit with tails\n  4 = batch fit with free E, everything else fixed"<<std::endl;
    logit_s<<"--optimizer <lmfit, mpfit> : Choose which optimizer to use for --optimize-fit-override-params or matrix fit routine \n"<<std::endl;
    logit_s<<"Fitting Routines: "<<std::endl;
    logit_s<<"--roi : ROI "<<std::endl;
    logit_s<<"--roi_plus : SVD method "<<std::endl;
    logit_s<<"--nnls : Non-Negative Least Squares"<<std::endl;
    logit_s<<"--tails : Fit with multiple parameters "<<std::endl;
    logit_s<<"--matrix : Fit with locked parameters \n"<<std::endl;
    logit_s<<"Dataset: "<<std::endl;
    logit_s<<"--dir : Dataset directory "<<std::endl;
    logit_s<<"--files : Dataset files: comma (',') separated if multiple \n"<<std::endl;
    logit_s<<"Examples: "<<std::endl;
    logit_s<<"   Perform roi and matrix analysis on the directory /data/dataset1 "<<std::endl;
    logit_s<<"xrf_maps --roi --matrix --dir /data/dataset1 "<<std::endl;
    logit_s<<"   Perform roi and matrix analysis on the directory /data/dataset1 but only process scan1 and scan2 "<<std::endl;
    logit_s<<"xrf_maps --roi --matrix --dir /data/dataset1 --files scan1.mda,scan2.mda"<<std::endl;
    logit_s<<"   Perform roi, matrix, and nnls  analysis on the directory /data/dataset1, use maps_standard.txt information for quantification "<<std::endl;
    logit_s<<"xrf_maps --roi --matrix --nnls --quantify-with maps_standard.txt --dir /data/dataset1 "<<std::endl;
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    //QCoreApplication a(argc, argv);

    std::string dataset_dir;
    std::vector<std::string> dataset_files;
    std::vector<data_struct::xrf::Fitting_Routines> proc_types;
    std::string quant_standard_filename = "";
    std::string whole_command_line = "";
    bool optimize_fit_override_params = false;
    bool quick_n_dirty = false;
    std::string opt = "";

    //Default is to process detectors 0 through 3
    size_t detector_num_start = 0;
    size_t detector_num_end = 3;

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


    size_t num_threads = std::thread::hardware_concurrency() - 2;
    if ( clp.option_exists("--nthreads") )
    {
        num_threads = std::stoi(clp.get_option("--nthreads"));
    }

    if ( clp.option_exists("--tails") )
    {
        proc_types.push_back(data_struct::xrf::GAUSS_TAILS);
    }
    if ( clp.option_exists("--matrix") )
    {
        proc_types.push_back(data_struct::xrf::GAUSS_MATRIX);
    }
    if ( clp.option_exists("--roi") )
    {
        proc_types.push_back(data_struct::xrf::ROI);
    }
    if ( clp.option_exists("--roi_plus") )
    {
        proc_types.push_back(data_struct::xrf::SVD);
    }
    if ( clp.option_exists("--nnls") )
    {
        proc_types.push_back(data_struct::xrf::NNLS);
    }

    if ( clp.option_exists("--quantify-with") )
    {
        quant_standard_filename = clp.get_option("--quantify-with");
    }

    if( clp.option_exists("--optimize-fit-override-params") )
    {
        optimize_fit_override_params = true;

        std::string opt = clp.get_option("--optimize-fit-override-params");
        if(opt == "1")
            optimize_fit_params_preset = fitting::models::MATRIX_BATCH_FIT;
        else if(opt == "2")
            optimize_fit_params_preset = fitting::models::BATCH_FIT_NO_TAILS;
        else if(opt == "3")
            optimize_fit_params_preset = fitting::models::BATCH_FIT_WITH_TAILS;
        else if(opt == "4")
            optimize_fit_params_preset = fitting::models::BATCH_FIT_WITH_FREE_ENERGY;
        else
            logit<<"Defaulting optimize_fit_params_preset to batch fit without tails"<<std::endl;
    }

    if( clp.option_exists("--optimizer"))
    {
        std::string opt = clp.get_option("--optimizer");
        /* TODO: connect to process_stream and process_whole
        if(opt == "mpfit")
        {
            optimizer = &mpfit_optimizer;
        }
        */
    }


    if( clp.option_exists("--add-exchange"))
    {
        //TODO:
    }
    if( clp.option_exists("--quick-and-dirty"))
    {
        quick_n_dirty = true;
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

    io::check_and_create_dirs(dataset_dir);

    if (proc_types.size() == 0 && optimize_fit_override_params == false && clp.option_exists("--generate-avg-h5") == false)
    {
        help();
        return -1;
    }

    std::string dset_file = clp.get_option("--files");
    if (dset_file.length() < 1)
    {
        // find all files in the dataset
        dataset_files = io::find_all_dataset_files(dataset_dir + "mda/", ".mda");
        if (dataset_files.size() == 0)
        {
            logit<<"Error: No mda files found in dataset directory "<<dataset_dir<<std::endl;
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

    //gen whole command line to save in hdf5 later
    for(int ic = 0; ic < argc; ic++)
    {
        whole_command_line += " " + std::string(argv[ic]);
    }

    logit<<"whole command line : "<<whole_command_line<<std::endl;

    logit << "Processing detectors " << detector_num_start << " - "<< detector_num_end <<std::endl;

    start = std::chrono::system_clock::now();

    //load element information
    if(false == io::load_element_info(element_henke_filename, element_csv_filename, data_struct::xrf::Element_Info_Map::inst()))
    {
        logit<<"Error loading element information: "<<std::endl;
        return -1;
    }

    io::populate_netcdf_hdf5_files(dataset_dir);

    if( clp.option_exists("--streaming"))
    {
        data_struct::xrf::Analysis_Job analysis_job(num_threads);
        analysis_job.set_optimizer(opt);

        if(optimize_fit_override_params)
        {
            std::vector<std::string> optim_dataset_files;
            if (dataset_files.size() == 0)
            {
                logit<<"Error: No mda files found in dataset directory "<<dataset_dir<<std::endl;
                return -1;
            }
            for (auto& itr : dataset_files)
            {
                optim_dataset_files.push_back(itr);
            }

            io::sort_dataset_files_by_size(dataset_dir, &optim_dataset_files);
            //if no files were specified only take the 8 largest datasets
            if (dset_file.length() < 1)
            {
                while (optim_dataset_files.size() > 9)
                {
                    optim_dataset_files.pop_back();
                }
            }

            if( analysis_job.load(dataset_dir, optim_dataset_files, proc_types, detector_num_start, detector_num_end) )
            {
                run_optimization_stream_pipeline(&analysis_job);
            }
            else
            {
                logit<<"Error initializing meta data. Exiting"<<std::endl;
                return -1;
            }
        }

        if(proc_types.size() > 0)
        {
            if (quant_standard_filename.length() > 0)
            {
                /*
                bool quant = false;
                quant = perform_quantification(dataset_dir, quant_standard_filename, proc_types, &quant_stand_list, &fit_params_override_dict, detector_num_start, detector_num_end);
                //if it is quick and dirty, average quants and save in first
                if( quant && quick_n_dirty)
                {
                    average_quantification(&quant_stand_list, detector_num_start, detector_num_end);
                }
                */
            }
            //relead to process all dataset files and in case optimizer changed fit parameters
            if(false == analysis_job.load(dataset_dir, dataset_files, proc_types, detector_num_start, detector_num_end) )
            {
                logit<<"Error initializing meta data. Exiting"<<std::endl;
                return -1;
            }

            if(quick_n_dirty)
            {
                run_quick_n_dirty_pipeline(&analysis_job);
            }
            else
            {
                run_stream_pipeline(&analysis_job);
            }
        }
        else
        {
            if(clp.option_exists("--generate-avg-h5"))
            {
                for(std::string dataset_file : dataset_files)
                {
                    io::generate_h5_averages(dataset_dir, dataset_file, detector_num_start, detector_num_end);
                }
            }
        }
    }
    else
    {
        //dict for override info and elements to fit.
        std::unordered_map<int, data_struct::xrf::Params_Override> fit_params_override_dict;
        ThreadPool tp(num_threads);
        std::vector<data_struct::xrf::Quantification_Standard> quant_stand_list(4);

        if(optimize_fit_override_params)
        {
            std::vector<std::string> optim_dataset_files;
            if (dataset_files.size() == 0)
            {
                logit<<"Error: No mda files found in dataset directory "<<dataset_dir<<std::endl;
                return -1;
            }
            for (auto& itr : dataset_files)
            {
                optim_dataset_files.push_back(itr);
            }

            io::sort_dataset_files_by_size(dataset_dir, &optim_dataset_files);
            //if no files were specified only take the 8 largest datasets
            if (dset_file.length() < 1)
            {
                while (optim_dataset_files.size() > 9)
                {
                    optim_dataset_files.pop_back();
                }
            }
            generate_optimal_params(dataset_dir, optim_dataset_files, &tp, detector_num_start, detector_num_end);
        }

        //try to load maps fit params override txt files for each detector. -1 is general one
        for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
        {
            data_struct::xrf::Params_Override params_override(dataset_dir, detector_num);
            if( io::load_override_params(dataset_dir, detector_num, &params_override) )
            {
                fit_params_override_dict[detector_num] = params_override;
            }
        }
        data_struct::xrf::Params_Override params(dataset_dir, -1);
        if( io::load_override_params(dataset_dir, -1, &params) )
        {
            fit_params_override_dict[-1] = params;
        }

        //check to make sure we have at least 1
        if(fit_params_override_dict.size() == 0)
        {
            logit<<"Error loading any maps_fit_params_override.txt "<<std::endl;
            return -1;
        }

        if(proc_types.size() > 0)
        {
            if (quant_standard_filename.length() > 0)
            {
                bool quant = false;
                quant = perform_quantification(dataset_dir, quant_standard_filename, proc_types, &quant_stand_list, &fit_params_override_dict, detector_num_start, detector_num_end);
                //if it is quick and dirty, average quants and save in first
                if( quant && quick_n_dirty)
                {
                    average_quantification(&quant_stand_list, detector_num_start, detector_num_end);
                }
            }

            for(std::string dataset_file : dataset_files)
            {
                if(quick_n_dirty)
                {
                    process_dataset_file_quick_n_dirty(dataset_dir, dataset_file, proc_types, &tp, &quant_stand_list, &fit_params_override_dict, detector_num_start, detector_num_end);
                }
                else
                {
                    process_dataset_file(dataset_dir, dataset_file, proc_types, &tp, &quant_stand_list, &fit_params_override_dict, detector_num_start, detector_num_end);
                    io::generate_h5_averages(dataset_dir, dataset_file, detector_num_start, detector_num_end);
                }

            }
        }
        else
        {
            if(clp.option_exists("--generate-avg-h5"))
            {
                for(std::string dataset_file : dataset_files)
                {
                    io::generate_h5_averages(dataset_dir, dataset_file, detector_num_start, detector_num_end);
                }
            }
        }
        //cleanup
        for (auto & itr : params.elements_to_fit)
        {
            delete itr.second;
        }
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    logit << "=-=-=-=-=-=- Total elapsed time: " << elapsed_seconds.count() << "s =-=-=-=-=-=-=-\n\n";

    data_struct::xrf::Element_Info_Map::inst()->clear();

    return 0;
    //return a.exec();
}
