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
#include <cctype>

// ----------------------------------------------------------------------------

void help()
{
    logit_s<<"Help: \n";
    logit_s<<"Usage: xrf_maps [Options] --dir [dataset directory] \n\n";
    logit_s<<"Options: \n";
    logit_s<<"--nthreads : <int> number of threads to use (default is all system threads) \n";
    logit_s<<"--quantify-with : <standard.txt> File to use as quantification standard \n";
    logit_s<<"--detectors : <int,..> Detectors to process, Defaults to 0,1,2,3 for 4 detector \n";
    logit_s<<"--generate-avg-h5 : Generate .h5 file which is the average of all detectors .h50 - h.53 or range specified. \n";
    logit_s<<"--add-v9layout : Generate .h5 file which has v9 layout able to open in IDL MAPS software. \n";
    logit_s<<"--add-exchange : Add exchange group into hdf5 file with normalized data.\n";
    logit_s<< "--export-csv : Export Integrated spec, fitted, background to csv file.\n";
	logit_s<< "--update-theta : <theta_pv_string> Update the theta dataset value using theta_pv_string as new pv string ref.\n";
//    logit_s<< "--update-scalers : If scalers pv's have been changed in maps_fit_parameters_override.txt file, you can run this to just update scaler values without refitting.\n";
	logit_s << "--update-amps <us_amp>,<ds_amp>: Updates upstream and downstream amps if they changed inbetween scans.\n";
	logit_s << "--update-quant-amps <us_amp>,<ds_amp>: Updates upstream and downstream amps for quantification if they changed inbetween scans.\n";
    logit_s<<"--quick-and-dirty : Integrate the detector range into 1 spectra.\n";
//	logit_s<< "--mem-limit <limit> : Limit the memory usage. Append M for megabytes or G for gigabytes\n";
    logit_s<<"--optimize-fit-override-params : <int> Integrate the 8 largest mda datasets and fit with multiple params.\n"<<
               "  1 = matrix batch fit\n  2 = batch fit without tails\n  3 = batch fit with tails\n  4 = batch fit with free E, everything else fixed \n";
    logit_s<<"--optimize-fit-routine : <general,hybrid> General (default): passes elements amplitudes as fit parameters. Hybrid only passes fit parameters and fits element amplitudes using NNLS\n";
    logit_s<<"--optimizer <lmfit, mpfit> : Choose which optimizer to use for --optimize-fit-override-params or matrix fit routine \n";
    logit_s<<"Fitting Routines: \n";
	logit_s<< "--fit <routines,> comma seperated \n";
    logit_s<<"  roi : element energy region of interest \n";
    logit_s<<"  roi_plus : SVD method \n";
    logit_s<<"  nnls : Non-Negative Least Squares \n";
    logit_s<<"  tails : Fit with multiple parameters \n";
    logit_s<<"  matrix : Fit with locked parameters \n\n";
    logit_s<<"Dataset: "<<"\n";
    logit_s<<"--dir : Dataset directory \n";
    logit_s<<"--files : Dataset files: comma (',') separated if multiple \n";
#ifdef _BUILD_WITH_ZMQ
    logit_s<<"Network: \n";
    logit_s<<"--streamin [source ip] : Accept a ZMQ stream of spectra to process. Source ip defaults to localhost (must compile with -DBUILD_WITH_ZMQ option) \n";
    logit_s<<"--streamout : Streams the analysis counts over a ZMQ stream (must compile with -DBUILD_WITH_ZMQ option) \n\n";
#endif
    logit_s<<"Examples: \n";
    logit_s<<"   Perform roi and matrix analysis on the directory /data/dataset1 \n";
    logit_s<<"xrf_maps --fit roi,matrix --dir /data/dataset1 \n";
    logit_s<<"   Perform roi and matrix analysis on the directory /data/dataset1 but only process scan1 and scan2 \n";
    logit_s<<"xrf_maps --fit roi,matrix --dir /data/dataset1 --files scan1.mda,scan2.mda \n";
    logit_s<<"   Perform roi, matrix, and nnls  analysis on the directory /data/dataset1, use maps_standard.txt information for quantification \n";
    logit_s<<"xrf_maps --fit roi,matrix,nnls --quantify-with maps_standard.txt --dir /data/dataset1 \n";
}

// ----------------------------------------------------------------------------

template <typename T_real>
void set_optimizer(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    //Which optimizer do we want to pick. Default is lmfit
    if (clp.option_exists("--optimizer"))
    {
        analysis_job.set_optimizer(clp.get_option("--optimizer"));
    }
}

// ----------------------------------------------------------------------------

template <typename T_real>
void set_mem_limit(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    //check if we should set a ram memory limit
    if (clp.option_exists("--mem-limit"))
    {
        std::string memlimit = clp.get_option("--mem-limit");
        if (memlimit.rfind('M'))
        {

        }
        else if (memlimit.rfind('G'))
        {

        }
        else
        {
            logW << "Could not parse --mem-limit parameter. Make sure to use M for megabytes or G for gigabytes. ex 200M\n";
        }
        //analysis_job.mem_limit = ;
    }
}

// ----------------------------------------------------------------------------

template <typename T_real>
void set_num_threads(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    if (clp.option_exists("--nthreads"))
    {
        analysis_job.num_threads = std::stoi(clp.get_option("--nthreads"));
    }
}

// ----------------------------------------------------------------------------

template <typename T_real>
void set_detectors(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    //What detector range should we process. Usually there are 4 detectors.
    if (clp.option_exists("--detectors"))
    {
        std::string detector_range_str = clp.get_option("--detectors");
        if (detector_range_str.find(',') != std::string::npos)
        {
            // if we found a comma, split the string to get list of dataset files
            std::stringstream ss;
            ss.str(detector_range_str);
            std::string item;
            while (std::getline(ss, item, ','))
            {
                analysis_job.detector_num_arr.push_back(std::stoi(item));
            }
        }
        else
        {
            analysis_job.detector_num_arr.push_back(std::stoi(detector_range_str));
        }
    }
    else if (clp.option_exists("--detector-range"))
    {
        // 0:3 will do the first 4 detectors
        // 0:1 will do the first 2 detectors
        // 2:3 will do the last 2 detectors
        std::string detector_range_str = clp.get_option("--detector-range");
        if (detector_range_str.find(':') != std::string::npos)
        {
            // if we found a comma, split the string to get list of dataset files
            std::stringstream ss;
            ss.str(detector_range_str);
            std::string item;
            std::getline(ss, item, ':');
            size_t detector_num_start = std::stoi(item);
            std::getline(ss, item, ':');
            size_t detector_num_end = std::stoi(item);
            for (size_t det = detector_num_start; det <= detector_num_end; det++)
            {
                analysis_job.detector_num_arr.push_back(det);
            }
        }
        else
        {
            analysis_job.detector_num_arr.push_back(std::stoi(detector_range_str));
        }
    }
    else
    {
        for (size_t det = 0; det < 7; det++)
        {
            analysis_job.detector_num_arr.push_back(det);
        }
    }
}

// ----------------------------------------------------------------------------

template <typename T_real>
void set_fit_routines(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    //Look for which analysis types we want to run
    if (clp.option_exists("--fit"))
    {
        std::string fitting = clp.get_option("--fit");
        if (fitting.find(',') != std::string::npos)
        {
            // if we found a comma, split the string to get list of dataset files
            std::stringstream ss;
            ss.str(fitting);
            std::string item;
            while (std::getline(ss, item, ','))
            {
                for (std::string::size_type x = 0; x < item.length(); ++x)
                {
                    item[x] = std::toupper(item[x]);
                }

                if (item == STR_FIT_ROI)
                {
                    analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::ROI);
                }
                else if (item == STR_FIT_SVD)
                {
                    analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::SVD);
                }
                else if (item == STR_FIT_NNLS)
                {
                    analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::NNLS);
                }
                else if (item == STR_FIT_GAUSS_MATRIX || item == "MATRIX")
                {
                    analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::GAUSS_MATRIX);
                }
            }
        }
        else
        {
            for (std::string::size_type x = 0; x < fitting.length(); ++x)
            {
                fitting[x] = std::toupper(fitting[x]);
            }

            if (fitting == STR_FIT_ROI)
            {
                analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::ROI);
            }
            else if (fitting == STR_FIT_SVD)
            {
                analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::SVD);
            }
            else if (fitting == STR_FIT_NNLS)
            {
                analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::NNLS);
            }
            else if (fitting == STR_FIT_GAUSS_MATRIX || fitting == "MATRIX")
            {
                analysis_job.fitting_routines.push_back(data_struct::Fitting_Routines::GAUSS_MATRIX);
            }
        }
    }
}

// ----------------------------------------------------------------------------

template <typename T_real>
int set_dir_and_files(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    //Get the dataset directory you want to process
    std::string dataset_dir = clp.get_option("--dir");
    if (dataset_dir.length() < 1)
    {
        help();
        return -1;
    }
    //add a slash if missing at the end
    if (dataset_dir.back() != DIR_END_CHAR)
    {
        dataset_dir += DIR_END_CHAR;
    }
    //replace / with \ for windows, won't do anything for linux
    std::replace(dataset_dir.begin(), dataset_dir.end(), '/', DIR_END_CHAR);

    //We save our ouput file in $dataset_directory/img.dat  Make sure we create this directory if it doesn't exist
    io::file::check_and_create_dirs(dataset_dir);


    //Look if files were specified
    std::string dset_file = clp.get_option("--files");
    //if they were not then look for them in $dataset_directory/mda
    if (dset_file.length() < 1)
    {
        for (auto& itr : io::file::File_Scan::inst()->find_all_dataset_files(dataset_dir, ".hdf5"))
        {
            analysis_job.dataset_files.push_back(itr);
        }
        for (auto& itr : io::file::File_Scan::inst()->find_all_dataset_files(dataset_dir, ".h5"))
        {
            analysis_job.dataset_files.push_back(itr);
        }
        for (auto& itr : io::file::File_Scan::inst()->find_all_dataset_files(dataset_dir, ".emd"))
        {
            analysis_job.dataset_files.push_back(itr);
        }
        for (auto& itr : io::file::File_Scan::inst()->find_all_dataset_files(dataset_dir, ".mda"))
        {
            analysis_job.dataset_files.push_back(itr);
        }
        for (auto& itr : io::file::File_Scan::inst()->find_all_dataset_files(dataset_dir + "mda" + DIR_END_CHAR, ".mda"))
        {
            analysis_job.dataset_files.push_back(itr);
        }
        // don't want to open h5 avg files for optimize
        for (auto& itr : analysis_job.dataset_files)
        {
            analysis_job.optimize_dataset_files.push_back(itr);
        }

        for (auto& itr : io::file::File_Scan::inst()->find_all_dataset_files(dataset_dir + "img.dat" + DIR_END_CHAR, ".h5"))
        {
            analysis_job.dataset_files.push_back(itr);
        }

        if (analysis_job.dataset_files.size() == 0)
        {
            logE << "No mda files found in dataset directory " << dataset_dir << "\n";
            return -1;
        }

        io::file::File_Scan::inst()->sort_dataset_files_by_size(dataset_dir, &analysis_job.optimize_dataset_files);

        //if no files were specified only take the 8 largest datasets
        while (analysis_job.optimize_dataset_files.size() > 9)
        {
            analysis_job.optimize_dataset_files.pop_back();
        }

    }
    else if (dset_file.find(',') != std::string::npos)
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

    analysis_job.dataset_directory = dataset_dir;

    return 0;
}

// ----------------------------------------------------------------------------

template <typename T_real>
void set_whole_command(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    /*
    //gen whole command line to save in hdf5 later
    for (int ic = 0; ic < argc; ic++)
    {
        whole_command_line += " " + std::string(argv[ic]);
    }
    analysis_job.command_line = whole_command_line;
    logI << "whole command line : " << whole_command_line << "\n";
    */
}

// ----------------------------------------------------------------------------

template <typename T_real>
int set_general_options(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    set_num_threads(clp, analysis_job);
    set_detectors(clp, analysis_job);
    set_optimizer(clp, analysis_job);
    set_whole_command(clp, analysis_job);
    return set_dir_and_files(clp, analysis_job);
}

// ----------------------------------------------------------------------------

int run_optimization(Command_Line_Parser& clp)
{
    //main structure for analysis job information
    // Force optimization to double prec
    data_struct::Analysis_Job<double> analysis_job; 

    if (set_general_options(clp, analysis_job) == -1)
    {
        return -1;
    }

    //Do we want to optimize our fitting parameters
    if (clp.option_exists("--optimize-fit-override-params"))
    {
        std::string opt = clp.get_option("--optimize-fit-override-params");
        if (opt == "1")
        {
            analysis_job.optimize_fit_params_preset = fitting::models::Fit_Params_Preset::MATRIX_BATCH_FIT;
        }
        else if (opt == "2")
        {
            analysis_job.optimize_fit_params_preset = fitting::models::Fit_Params_Preset::BATCH_FIT_NO_TAILS;
        }
        else if (opt == "3")
        {
            analysis_job.optimize_fit_params_preset = fitting::models::Fit_Params_Preset::BATCH_FIT_WITH_TAILS;
        }
        else if (opt == "4")
        {
            analysis_job.optimize_fit_params_preset = fitting::models::Fit_Params_Preset::BATCH_FIT_WITH_FREE_ENERGY;
        }
        else
        {
            logI << "Defaulting optimize_fit_params_preset to batch fit without tails" << "\n";
        }
    }

    if (clp.option_exists("--optimize-fit-routine"))
    {
        std::string opt = clp.get_option("--optimize-fit-routine");
        if (opt == "hybrid")
        {
            analysis_job.optimize_fit_routine = OPTIMIZE_FIT_ROUTINE::HYBRID;
        }
    }

    if (io::file::init_analysis_job_detectors(&analysis_job))
    {
        io::file::File_Scan::inst()->populate_netcdf_hdf5_files(analysis_job.dataset_directory);
        
        generate_optimal_params(&analysis_job);
        // reload override files
        //io::init_analysis_job_detectors(&analysis_job);
        //iterate_datasets_and_update(analysis_job);
    }
    else
    {
        logE << "Error initalizing detectors!\n";
        return -1;
    }
    return 0;
}

// ----------------------------------------------------------------------------

int run_quantification(Command_Line_Parser& clp)
{
    //main structure for analysis job information
    // Force optimization to double prec
    data_struct::Analysis_Job<double> analysis_job;

    if (set_general_options(clp, analysis_job) == -1)
    {
        return -1;
    }
    set_fit_routines(clp, analysis_job);

    //Check if we want to quantifiy with a standard
    if (clp.option_exists("--quantify-with"))
    {
        analysis_job.quantification_standard_filename = clp.get_option("--quantify-with");
    }

    if (io::file::init_analysis_job_detectors(&analysis_job))
    {
        io::file::File_Scan::inst()->populate_netcdf_hdf5_files(analysis_job.dataset_directory);

        if (analysis_job.quantification_standard_filename.length() > 0)
        {
            perform_quantification(&analysis_job);
        }
        else
        {
            logE << "Please specify filename with quantification information (usually maps_standardinfo.txt)\n";
            return -1;
        }


        iterate_datasets_and_update(analysis_job);
    }
    else
    {
        logE << "Error initalizing detectors!\n";
        return -1;
    }
    return 0;
}

// ----------------------------------------------------------------------------

int run_streaming(Command_Line_Parser& clp)
{
    data_struct::Analysis_Job<float> analysis_job;

    if (set_general_options(clp, analysis_job) == -1)
    {
        return -1;
    }
    set_fit_routines(clp, analysis_job);

    // init our job and run
    if (io::file::init_analysis_job_detectors(&analysis_job))
    {
        io::file::File_Scan::inst()->populate_netcdf_hdf5_files(analysis_job.dataset_directory);
        
        // If we have fitting routines then stream the counts per sec
        if (analysis_job.fitting_routines.size() > 0)
        {
            //if we are streaming we use 1 thread for loading and 1 for saving
            //analysis_job.num_threads = std::thread::hardware_concurrency() - 1;
            //analysis_job.theta_pv = "2xfm:m53.VAL";
            analysis_job.theta_pv = clp.get_option("--theta_pv");
            run_stream_pipeline(&analysis_job);
        }
        // otherwise just stream the spectra
        else
        {
            stream_spectra(&analysis_job);
        }
    }
    else
    {
        logE << "Error initalizing detectors!\n";
        return -1;
    }
    return 0;
}

// ----------------------------------------------------------------------------

int run_fits(Command_Line_Parser& clp)
{
    //main structure for analysis job information
    data_struct::Analysis_Job<float> analysis_job;

    if (set_general_options(clp, analysis_job) == -1)
    {
        return -1;
    }
    set_fit_routines(clp, analysis_job);
    set_mem_limit(clp, analysis_job);


    //Should we sum up all the detectors and process it as one?
    if (clp.option_exists("--quick-and-dirty"))
    {
        analysis_job.quick_and_dirty = true;
        analysis_job.generate_average_h5 = false;
    }
        
    /*
    bool update_h5_without_fitting = analysis_job.generate_average_h5 ||
        analysis_job.add_v9_layout ||
        analysis_job.add_exchange_layout ||
        analysis_job.update_theta_str.length() > 0 ||
        //analysis_job.update_scalers || 
        analysis_job.export_int_fitted_to_csv ||
        analysis_job.update_us_amps_str.length() > 0 ||
        analysis_job.update_quant_us_amps_str.length() > 0 ||
        analysis_job.add_background;

    bool update_h5_fit = analysis_job.fitting_routines.size() > 0 ||
        optimize_fit_override_params ||
        analysis_job.stream_over_network ||
        update_h5_without_fitting;
        
    //Check to make sure we have something to do. If not then show the help screen
    if (update_h5_fit == false)
    {
        help();
        return -1;
    }
    */

    // init our job and run
    if (io::file::init_analysis_job_detectors(&analysis_job))
    {
        io::file::File_Scan::inst()->populate_netcdf_hdf5_files(analysis_job.dataset_directory);
        
        if (analysis_job.fitting_routines.size() > 0)
        {
            process_dataset_files(&analysis_job);
            analysis_job.generate_average_h5 = true;
        }
        else
        {
            logW << "No fitting routines picked! Please select from [--fit roi,nnls,matrix]\n";
        }

        iterate_datasets_and_update(analysis_job);
    }
    else
    {
        logE << "Error initalizing detectors!\n";
        return -1;
    }
    return 0;
}

// ----------------------------------------------------------------------------

int run_h5_file_updates(Command_Line_Parser& clp)
{
    data_struct::Analysis_Job<float> analysis_job;

    //Should create an averaged file of all detectors processed
    if (clp.option_exists("--generate-avg-h5"))
    {
        analysis_job.generate_average_h5 = true;
    }

    if (clp.option_exists("--add-v9layout"))
    {
        analysis_job.add_v9_layout = true;
    }

    if (clp.option_exists("--update-theta"))
    {
        analysis_job.update_theta_str = clp.get_option("--update-theta");
    }

    /*
    if (clp.option_exists("--update-scalers"))
    {
        analysis_job.update_scalers = true;
    }
    */

    if (clp.option_exists("--update-amps"))
    {
        string amps = clp.get_option("--update-amps");
        size_t idx = amps.find(',');
        if (idx != string::npos)
        {
            analysis_job.update_us_amps_str = amps.substr(0, idx);
            analysis_job.update_ds_amps_str = amps.substr(idx + 1);
        }
        else
        {
            logW << "Could not find ',' while parsing --update-amps\n";
        }
    }

    if (clp.option_exists("--update-quant-amps"))
    {
        string amps = clp.get_option("--update-quant-amps");
        size_t idx = amps.find(',');
        if (idx != string::npos)
        {
            analysis_job.update_quant_us_amps_str = amps.substr(0, idx);
            analysis_job.update_quant_ds_amps_str = amps.substr(idx + 1);
        }
        else
        {
            logW << "Could not find ',' while parsing --update-quant-amps\n";
        }
    }

    //Added exchange format to output file. Used as an interface to allow other analysis software to load out output file
    if (clp.option_exists("--add-exchange"))
    {
        analysis_job.add_exchange_layout = true;
    }

    if (clp.option_exists("--export-csv"))
    {
        analysis_job.export_int_fitted_to_csv = true;
    }

    if (clp.option_exists("--add_background"))
    {
        analysis_job.add_background = true;
    }

    iterate_datasets_and_update(analysis_job);

    return 0;
}

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    
    //std::string whole_command_line = "";

    //Performance measure
    std::chrono::time_point<std::chrono::system_clock> start, end;

    //////// HENKE and ELEMENT INFO /////////////
    const std::string element_csv_filename = "../reference/xrf_library.csv";
    const std::string element_henke_filename = "../reference/henke.xdr";
    const std::string scaler_lookup_yaml = "../reference/Scaler_to_PV_map.yaml";

    start = std::chrono::system_clock::now();


    if (false == io::file::load_scalers_lookup(scaler_lookup_yaml))
    {
        logE << " Could not load " << scaler_lookup_yaml << ". Won't be able to translate from PV to Label for scalers!\n";
    }

    //load element information
    if (false == io::file::load_element_info<float>(element_henke_filename, element_csv_filename))
    {
        logE << "loading element information: " << "\n";
        return -1;
    }
    if (false == io::file::load_element_info<double>(element_henke_filename, element_csv_filename))
    {
        logE << "loading element information: " << "\n";
        return -1;
    }

    Command_Line_Parser clp(argc, argv);

    if (clp.option_exists("-h"))
    {
        help();
        return 0;
    }

    if (clp.option_exists("--optimize-fit-override-params"))
    {
        run_optimization(clp);
    }

    if (clp.option_exists("--streamin") || clp.option_exists("--streamout") )
    {
        run_streaming(clp);
    }
    else if (clp.option_exists("--fit") )
    {
        run_fits(clp);
    }

    if (clp.option_exists("--quantify-with"))
    {
        run_quantification(clp);
    }

    run_h5_file_updates(clp);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "=-=-=-=-=-=- Total elapsed time: " << elapsed_seconds.count() << "s =-=-=-=-=-=-=-\n" << std::endl; //endl will flush the print.


    data_struct::Element_Info_Map<float>::inst()->clear();
    data_struct::Element_Info_Map<double>::inst()->clear();

}
