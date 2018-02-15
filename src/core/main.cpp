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

#include "core/command_line_parser.h"
#include "core/process_streaming.h"
#include "core/process_whole.h"

// ----------------------------------------------------------------------------

void help()
{
    logit_s<<"Help: "<<"\n";
    logit_s<<"Usage: xrf_maps [Options] [Fitting Routines] --dir [dataset directory] \n"<<"\n";
    logit_s<<"Options: "<<"\n";
    logit_s<<"--nthreads : <int> number of threads to use (default is all system threads) "<<"\n";
    logit_s<<"--quantify-with : <standard.txt> File to use as quantification standard "<<"\n";
    logit_s<<"--detector-range : <int:int> Start and end detector range. Defaults to 0:3 for 4 detector "<<"\n";
    logit_s<<"--generate-avg-h5 : Generate .h5 file which is the average of all detectors .h50 - h.53 or range specified. "<<"\n";
//    logit_s<<"--add-exchange : <us:ds:sr> Add exchange group into hdf5 file with normalized data.\n";
//    logit_s<<"    us = upstream ion chamber\n";
//    logit_s<<"    ds = downstream ion chamber\n";
//    logit_s<<"    sr = sr current. "<<"\n";
    logit_s<<"--quick-and-dirty : Integrate the detector range into 1 spectra. "<<"\n";
    logit_s<<"--optimize-fit-override-params : <int> Integrate the 8 largest mda datasets and fit with multiple params\n"<<
               "  1 = matrix batch fit\n  2 = batch fit without tails\n  3 = batch fit with tails\n  4 = batch fit with free E, everything else fixed"<<"\n";
    logit_s<<"--optimizer <lmfit, mpfit> : Choose which optimizer to use for --optimize-fit-override-params or matrix fit routine \n"<<"\n";
    logit_s<<"Fitting Routines: "<<"\n";
    logit_s<<"--roi : ROI "<<"\n";
    logit_s<<"--roi_plus : SVD method "<<"\n";
    logit_s<<"--nnls : Non-Negative Least Squares"<<"\n";
    logit_s<<"--tails : Fit with multiple parameters "<<"\n";
    logit_s<<"--matrix : Fit with locked parameters \n"<<"\n";
    logit_s<<"Dataset: "<<"\n";
    logit_s<<"--dir : Dataset directory "<<"\n";
    logit_s<<"--files : Dataset files: comma (',') separated if multiple \n"<<"\n";
    logit_s<<"--confocal : load hdf confocal xrf datasets \n"<<"\n";
	logit_s << "--emd : load hdf electron microscopy FEI EMD xrf datasets \n" << "\n";
    logit_s<<"Examples: "<<"\n";
    logit_s<<"   Perform roi and matrix analysis on the directory /data/dataset1 "<<"\n";
    logit_s<<"xrf_maps --roi --matrix --dir /data/dataset1 "<<"\n";
    logit_s<<"   Perform roi and matrix analysis on the directory /data/dataset1 but only process scan1 and scan2 "<<"\n";
    logit_s<<"xrf_maps --roi --matrix --dir /data/dataset1 --files scan1.mda,scan2.mda"<<"\n";
    logit_s<<"   Perform roi, matrix, and nnls  analysis on the directory /data/dataset1, use maps_standard.txt information for quantification "<<"\n";
    logit_s<<"xrf_maps --roi --matrix --nnls --quantify-with maps_standard.txt --dir /data/dataset1 "<<"\n";
}

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{
    std::string dataset_dir;
    std::string whole_command_line = "";
    bool is_confocal = false;

    //Performance measure
    std::chrono::time_point<std::chrono::system_clock> start, end;

    //////// HENKE and ELEMENT INFO /////////////
    std::string element_csv_filename = "../reference/xrf_library.csv";
    std::string element_henke_filename = "../reference/henke.xdr";

    //main structure for analysis job information
    data_struct::Analysis_Job analysis_job;

    Command_Line_Parser clp(argc, argv);

    if( clp.option_exists("-h") )
    {
       help();
       return 0;
    }

    analysis_job.detector_num_start = 0;
    analysis_job.detector_num_end = 3;

    analysis_job.num_threads = std::thread::hardware_concurrency();
    if ( clp.option_exists("--nthreads") )
    {
        analysis_job.num_threads = std::stoi(clp.get_option("--nthreads"));
    }

    //Look for which analysis types we want to run
    if ( clp.option_exists("--tails") )
    {
        analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::GAUSS_TAILS);
    }
    if ( clp.option_exists("--matrix") )
    {
        analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::GAUSS_MATRIX);
    }
    if ( clp.option_exists("--roi") )
    {
        analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::ROI);
    }
    if ( clp.option_exists("--roi_plus") )
    {
        analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::SVD);
    }
    if ( clp.option_exists("--nnls") )
    {
        analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::NNLS);
    }

    //Check if we want to quantifiy with a standard
    if ( clp.option_exists("--quantify-with") )
    {
        analysis_job.quantificaiton_standard_filename = clp.get_option("--quantify-with");
    }

    //Do we want to optimize our fitting parameters
    if( clp.option_exists("--optimize-fit-override-params") )
    {
        analysis_job.optimize_fit_override_params = true;

        std::string opt = clp.get_option("--optimize-fit-override-params");
        if(opt == "1")
            analysis_job.optimize_fit_params_preset = fitting::models::MATRIX_BATCH_FIT;
        else if(opt == "2")
            analysis_job.optimize_fit_params_preset = fitting::models::BATCH_FIT_NO_TAILS;
        else if(opt == "3")
            analysis_job.optimize_fit_params_preset = fitting::models::BATCH_FIT_WITH_TAILS;
        else if(opt == "4")
            analysis_job.optimize_fit_params_preset = fitting::models::BATCH_FIT_WITH_FREE_ENERGY;
        else
            logit<<"Defaulting optimize_fit_params_preset to batch fit without tails"<<"\n";
    }

    //Which optimizer do we want to pick. Default is lmfit
    if( clp.option_exists("--optimizer"))
    {
        analysis_job.set_optimizer(clp.get_option("--optimizer"));
    }

    //Added exchange format to output file. Used as an interface to allow other analysis software to load out output file
    if( clp.option_exists("--add-exchange"))
    {
        //TODO:
    }

    //Should we sum up all the detectors and process it as one?
    if( clp.option_exists("--quick-and-dirty"))
    {
        analysis_job.quick_and_dirty = true;
        analysis_job.generate_average_h5 = false;
    }
    //Should create an averaged file of all detectors processed
    else if(clp.option_exists("--generate-avg-h5"))
    {
        analysis_job.generate_average_h5 = true;
    }


    if(clp.option_exists("--confocal"))
    {
        is_confocal = true;
    }
	if (clp.option_exists("--emd"))
	{
		analysis_job.is_emd = true;
        analysis_job.detector_num_start = 0;
        analysis_job.detector_num_end = 0;
	}


    //TODO: add --quantify-only option if you already did the fits and just want to add quantification

    //What detector range should we process. Usually there are 4 detectors.
    // 0:3 will do the first 4 detectors
    // 0:1 will do the first 2 detectors
    // 2:3 will do the last 2 detectors
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
            analysis_job.detector_num_start = std::stoi(item);
            std::getline(ss, item, ':');
            analysis_job.detector_num_end = std::stoi(item);
        }
        else
        {
            analysis_job.detector_num_start = analysis_job.detector_num_end = std::stoi(detector_range_str);
        }
    }

    //Get the dataset directory you want to process
    dataset_dir = clp.get_option("--dir");
    if (dataset_dir.length() < 1)
    {
        help();
        return -1;
    }
    //add a slash if missing at the end
    if (dataset_dir.back() != '/' && dataset_dir.back() != '\\')
    {
        dataset_dir += "/";
    }

    //We save our ouput file in $dataset_directory/img.dat  Make sure we create this directory if it doesn't exist
    io::check_and_create_dirs(dataset_dir);


    //Check to make sure we have something to do. If not then show the help screen
    if (analysis_job.fitting_routines.size() == 0 && analysis_job.optimize_fit_override_params == false && clp.option_exists("--generate-avg-h5") == false)
    {
        help();
        return -1;
    }


    //Look if files were specified
    std::string dset_file = clp.get_option("--files");
    //if they were not then look for them in $dataset_directory/mda
    if (dset_file.length() < 1)
    {
        if(is_confocal)
        {
            analysis_job.dataset_files = io::find_all_dataset_files(dataset_dir, ".hdf5");
        }
		else if (analysis_job.is_emd)
		{
			analysis_job.dataset_files = io::find_all_dataset_files(dataset_dir, ".emd");
		}
        else
        {
            // find all files in the dataset
            analysis_job.dataset_files = io::find_all_dataset_files(dataset_dir + "mda/", ".mda");
        }
        if (analysis_job.dataset_files.size() == 0)
        {
            logit<<"Error: No mda files found in dataset directory "<<dataset_dir<<"\n";
            return -1;
        }

        for (auto& itr : analysis_job.dataset_files)
        {
            analysis_job.optimize_dataset_files.push_back(itr);
        }

        if(!is_confocal && !analysis_job.is_emd)
            io::sort_dataset_files_by_size(dataset_dir, &analysis_job.optimize_dataset_files);

        //if no files were specified only take the 8 largest datasets
        while (analysis_job.optimize_dataset_files.size() > 9)
        {
            analysis_job.optimize_dataset_files.pop_back();
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
            analysis_job.dataset_files.push_back(item);
            analysis_job.optimize_dataset_files.push_back(item);
        }
    }
    else
    {
        analysis_job.dataset_files.push_back(dset_file);
        analysis_job.optimize_dataset_files.push_back(dset_file);
    }

    //gen whole command line to save in hdf5 later
    for(int ic = 0; ic < argc; ic++)
    {
        whole_command_line += " " + std::string(argv[ic]);
    }
    analysis_job.command_line = whole_command_line;
    logit<<"whole command line : "<<whole_command_line<<"\n";
    logit << "Processing detectors " << analysis_job.detector_num_start << " - "<< analysis_job.detector_num_end <<"\n";

    start = std::chrono::system_clock::now();

    //load element information
    if(false == io::load_element_info(element_henke_filename, element_csv_filename, data_struct::Element_Info_Map::inst()))
    {
        logit<<"Error loading element information: "<<"\n";
        return -1;
    }

    analysis_job.dataset_directory = dataset_dir;

    // init our job and run
    if(io::init_analysis_job_detectors(&analysis_job))
    {
        if(analysis_job.optimize_fit_override_params)
        {
            //run_optimization_stream_pipeline(&analysis_job);
            io::populate_netcdf_hdf5_files(dataset_dir);
            generate_optimal_params(&analysis_job);
        }

        if (analysis_job.quantificaiton_standard_filename.length() > 0)
        {
            perform_quantification(&analysis_job);
        }

        if( clp.option_exists("--stream") || analysis_job.is_emd)
        {
            //if we are streaming we use 1 thread for loading and 1 for saving
            analysis_job.num_threads = std::thread::hardware_concurrency() - 1;
            analysis_job.stream_over_network = true;
            //analysis_job.theta_pv = "2xfm:m53.VAL";
            analysis_job.theta_pv = clp.get_option("--theta_pv");
            run_stream_pipeline(&analysis_job);
        }
        else
        {
            io::populate_netcdf_hdf5_files(dataset_dir);
            process_dataset_files(&analysis_job);
            analysis_job.generate_average_h5 = true;
        }

        //average all detectors to one files
        if(analysis_job.generate_average_h5)
        {
            for(std::string dataset_file : analysis_job.dataset_files)
            {
                io::generate_h5_averages(analysis_job.dataset_directory, dataset_file, analysis_job.detector_num_start, analysis_job.detector_num_end);
            }
        }
    }
    else
    {
        logit<<"Error initializing analysis job"<<"\n";
    }

    data_struct::Element_Info_Map::inst()->clear();

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    logit << "=-=-=-=-=-=- Total elapsed time: " << elapsed_seconds.count() << "s =-=-=-=-=-=-=-\n\n";


    return 0;
}
