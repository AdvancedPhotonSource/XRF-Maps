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


#define MAX_DETECTORS 7

// ----------------------------------------------------------------------------

void help()
{
    logit_s<<"Help: \n";
    logit_s<<"Usage: xrf_maps [Options] --dir [dataset directory] \n\n";
    logit_s<<"Options: \n";
    logit_s<<"--nthreads : <int> number of threads to use (default is all system threads) \n";
    logit_s<<"--quantify-with : <standard.txt> File to use as quantification standard \n";
    logit_s<<"--quantify-fit <routines,>: If you want to perform quantification without having to re-fit all datasets. See --fit for routine options \n";
    logit_s<<"--detectors : <int,..> Detectors to process, Defaults to 0,1,2,3 for 4 detector \n";
    logit_s<<"--detector-range : <int:int> Detector range 0:3 will process 0,1,2,3 detectors \n";
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
               "  0 = use override file\n  1 = matrix batch fit\n  2 = batch fit without tails\n  3 = batch fit with tails\n  4 = batch fit with free E, everything else fixed \n  5 = batch fit without tails, and fit energy quadratic\n";
    logit_s<<"--optimize-fit-routine : <general,hybrid> General (default): passes elements amplitudes as fit parameters. Hybrid only passes fit parameters and fits element amplitudes using NNLS\n";
    //logit_s<<"--optimizer <lmfit, mpfit> : Choose which optimizer to use for --optimize-fit-override-params or matrix fit routine \n";
   // logit_s<<"--optimizer-fx-tols <tol_override_val> : F_TOL, X_TOL, Default is LM_FIT = " << DP_LM_USERTOL << " , MP_FIT = " << 1.192e-10 << "\n";
   // logit_s<<"--optimizer-fxg-tols <tol_override_val> : F_TOL, X_TOL, G_TOL, Default is LM_FIT = " << DP_LM_USERTOL << " , MP_FIT = " << 1.192e-10 << "\n";
    logit_s<<"--optimizer-use-weights : Calculate and use weights for residual error function.\n";
    logit_s<<"--optimize-rois : Looks in 'rois' directory and performs --optimize-fit-override-params on each roi separately. Needs to have --quantify-rois-with <maps_standardinfo.txt> and --quantify-fit <routines,>  \n";
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
    logit_s<<"--streamout [port]: Streams the analysis counts over a ZMQ stream (must compile with -DBUILD_WITH_ZMQ option) \n\n";
#endif
    logit_s<<"Examples: \n";
    logit_s<<"   Perform roi and matrix analysis on the directory /data/dataset1 \n";
    logit_s<<"xrf_maps --fit roi,matrix --dir /data/dataset1 \n";
    logit_s<<"   Perform roi and matrix analysis on the directory /data/dataset1 but only process scan1 and scan2 \n";
    logit_s<<"xrf_maps --fit roi,matrix --dir /data/dataset1 --files scan1.mda,scan2.mda \n";
    logit_s<<"   Perform roi, matrix, and nnls  analysis on the directory /data/dataset1, use maps_standard.txt information for quantification \n";
    logit_s<<"xrf_maps --fit roi,matrix,nnls --quantify-with maps_standard.txt --dir /data/dataset1 \n";
    logit_s<<"   Perform optimization of an integrated roi spectra \n";
    logit_s<<"xrf_maps  --optimize-rois --quantify-rois-with maps_standardinfo.txt --quantify-fit roi,matrix,nnls --dir /data/dataset1 \n";
}

// ----------------------------------------------------------------------------

template <typename T_real>
void set_optimizer(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    bool fx_exists = clp.option_exists("--optimizer-fx-tols");
    bool fxg_exists = clp.option_exists("--optimizer-fxg-tols");

    if(clp.option_exists("--optimizer-use-weights"))
    {
        std::string val = clp.get_option("--optimizer-use-weights");
        std::transform(val.begin(), val.end(), val.begin(), [](unsigned char c) { return std::tolower(c); });
        if(val == "on" || val == "talse")
        {
            analysis_job.use_weights = true;
        }
        if(val == "off" || val == "false")
        {
            analysis_job.use_weights = false;
        }
        logI<<" Use weights = "<<analysis_job.use_weights<<"\n";
    }
    //Which optimizer do we want to pick. Default is lmfit
    if (clp.option_exists("--optimizer"))
    {
        analysis_job.set_optimizer(clp.get_option("--optimizer"));
    }

    if (clp.option_exists("--optimize-fit-routine"))
    {
        std::string opt = clp.get_option("--optimize-fit-routine");
        if (opt == "hybrid")
        {
            logI << "Using Hybrid optimizer\n";
            analysis_job.optimize_fit_routine = OPTIMIZE_FIT_ROUTINE::HYBRID;
        }
    }

    if (fx_exists || fxg_exists)
    {
        T_real fxg_tol = 0.00000000000000001; 
        if (std::is_same<T_real, float>::value)
        {
            if (fxg_exists)
            {
                fxg_tol = std::stof(clp.get_option("--optimizer-fxg-tols"));
            }
            else if (fx_exists)
            {
                fxg_tol = std::stof(clp.get_option("--optimizer-fx-tols"));
            }
        }
        else if (std::is_same<T_real, double>::value)
        {
            if (fxg_exists)
            {
                fxg_tol = std::stod(clp.get_option("--optimizer-fxg-tols"));
            }
            else if (fx_exists)
            {
                fxg_tol = std::stod(clp.get_option("--optimizer-fx-tols"));
            }
        }
        
        std::unordered_map<std::string, T_real> opt_map;
        opt_map[STR_OPT_FTOL] = fxg_tol;
        opt_map[STR_OPT_XTOL] = fxg_tol;
        if (fxg_exists)
        {
            opt_map[STR_OPT_GTOL] = fxg_tol;
            logI << "Setting FTOL, XTOL, GTOL to " << fxg_tol << "\n";
        }
        else
        {
            logI << "Setting FTOL, XTOL to " << fxg_tol << "\n";
        }
        analysis_job.optimizer()->set_options(opt_map);
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
        for (size_t det = 0; det < MAX_DETECTORS; det++)
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
    if (clp.option_exists("--fit") || clp.option_exists("--quantify-fit"))
    {
        std::string fitting = clp.get_option("--fit");
        if (fitting.length() < 1)
        {
            fitting = clp.get_option("--quantify-fit");
        }
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
    std::vector<std::string> ignore_dir_list = { "img.dat", "rois", "mda", "output" };
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

    analysis_job.dataset_directory = dataset_dir;

    if(clp.option_exists("--output-dir"))
    {
        analysis_job.output_dir = clp.get_option("--output-dir");
        if (analysis_job.output_dir.length() < 1)
        {
            logE<<"Please enter a proper directory path with --output-dir option\n";
            return -1;
        }
        //add a slash if missing at the end
        if (analysis_job.output_dir.back() != DIR_END_CHAR)
        {
            analysis_job.output_dir += DIR_END_CHAR;
        }
        //replace / with \ for windows, won't do anything for linux
        std::replace(analysis_job.output_dir.begin(), analysis_job.output_dir.end(), '/', DIR_END_CHAR);
        logI<<"Setting output directory: "<<analysis_job.output_dir<<"\n";
    }
    else
    {
        analysis_job.output_dir = dataset_dir;
    }
    
    //We save our ouput file in $dataset_directory/img.dat  Make sure we create this directory if it doesn't exist
    io::file::check_and_create_dirs(analysis_job.output_dir);

    //Look if files were specified
    std::string dset_file = clp.get_option("--files");
    //if they were not then look for them in $dataset_directory/mda
    if (dset_file.length() < 1 || dset_file == "all")
    {
        std::vector<std::string> flist;

        std::vector<std::string> search_ext_list;
        std::vector<std::string> search_ext_mda_list;
        std::vector<std::string> search_ext_h5_list;

        search_ext_list.push_back(".hdf5");
        search_ext_list.push_back(".h5");
        search_ext_list.push_back(".emd");
        search_ext_list.push_back(".mda");

        flist = io::file::File_Scan::inst()->find_all_dataset_files_by_list(dataset_dir, search_ext_list);
        for (const auto& fitr : flist)
        {
            if (fitr == "flyXRF.h5") // ignore folder
            {
                continue;
            }
            analysis_job.dataset_files.push_back(fitr);
        }

        search_ext_mda_list.push_back(".mda");
        search_ext_mda_list.push_back(".mca");
        for (int i = 0; i < MAX_DETECTORS; i++)
        {
            std::string mca_str = ".mca" + std::to_string(i);
            search_ext_mda_list.push_back(mca_str);
        }

        flist = io::file::File_Scan::inst()->find_all_dataset_files_by_list(dataset_dir + "mda" + DIR_END_CHAR, search_ext_mda_list);
        for (const auto& fitr : flist)
        {
            analysis_job.dataset_files.push_back(fitr);
        }

        // search recursivly for esrf datasets
        std::vector<std::string> root_dir_list = io::file::File_Scan::inst()->find_all_dirs(dataset_dir + DIR_END_CHAR, ignore_dir_list, false);
        
        search_ext_h5_list.push_back(".h5"); // added h5 for esrf datasets

        for (const auto& itr : root_dir_list)
        {
            if (itr == "flyXRF.h5" || itr == "rois" || itr == "output" || itr == "vlm" || itr == "flyXRF") // ignore these 
            {
                continue;
            }
            std::vector<std::string> flist = io::file::File_Scan::inst()->find_all_dataset_files_by_list(dataset_dir + DIR_END_CHAR + itr + DIR_END_CHAR, search_ext_h5_list);
            for (const auto& fitr : flist)
            {
                analysis_job.dataset_files.push_back(itr + DIR_END_CHAR + fitr);
            }
        }

        // don't want to open h5 avg files for optimize because we optimize by detector , not avg
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

        if(dset_file != "all")
        {
            //if no files were specified only take the 8 largest datasets
            while (analysis_job.optimize_dataset_files.size() > 9)
            {
                analysis_job.optimize_dataset_files.pop_back();
            }
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
void set_streaming_options(Command_Line_Parser& clp, data_struct::Analysis_Job<T_real>& analysis_job)
{
    if (clp.option_exists("--streamin"))
    {
        analysis_job.is_network_source = true;
        std::string source_addr = clp.get_option("--streamin");
        if (source_addr.length() > 0)
        {
            size_t idx = source_addr.find(':');
            if (idx != std::string::npos)
            {
                analysis_job.network_source_ip = source_addr.substr(0, idx);
                analysis_job.network_source_port = source_addr.substr(idx+1);
            }
            else
            {
                analysis_job.network_source_ip = source_addr;
                analysis_job.network_source_port = "43434";
            }
        }
    }
    if (clp.option_exists("--streamout"))
    {
        analysis_job.stream_over_network = true;
        // set default streaming out port
        analysis_job.network_stream_port = "43434";
        // if port is specified, then update it
        std::string out_port = clp.get_option("--streamout");
        if (out_port.length() > 0)
        {
            analysis_job.network_stream_port = out_port;
        }
    }
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
        if (opt == "0")
        {
            analysis_job.optimize_fit_params_preset = fitting::models::Fit_Params_Preset::NOT_SET; // use tags from override file
        }
        else if (opt == "1")
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
        else if (opt == "5")
        {
            analysis_job.optimize_fit_params_preset = fitting::models::Fit_Params_Preset::BATCH_FIT_NO_TAILS_E_QUAD;
        }
        else
        {
            logI << "Defaulting optimize_fit_params_preset to batch fit without tails" << "\n";
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
    if (analysis_job.fitting_routines.size() == 0)
    {
        logE << "Please specify fit routines with --quantify-fit roi,nnls,matrix \n";
        return -1;
    }
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
            perform_quantification(&analysis_job, true);
        }
        else
        {
            logE << "Please specify filename with quantification information (usually maps_standardinfo.txt)\n";
            return -1;
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

int run_optimize_rois(Command_Line_Parser& clp)
{
    data_struct::Analysis_Job<double> analysis_job;

    // Force perform quantification until we have proper loading quant from hdf5 finished
    if (set_general_options(clp, analysis_job) == -1)
    {
        return -1;
    }
    set_fit_routines(clp, analysis_job);
    /*
    if (analysis_job.fitting_routines.size() == 0)
    {
        logE << "Please specify fit routines with --quantify-fit roi,nnls,matrix \n";
        return -1;
    }
    */
    //Check if we want to quantifiy with a standard
    if (clp.option_exists("--quantify-rois-with"))
    {
        analysis_job.quantification_standard_filename = clp.get_option("--quantify-rois-with");
    }

    if (io::file::init_analysis_job_detectors(&analysis_job))
    {
        io::file::File_Scan::inst()->populate_netcdf_hdf5_files(analysis_job.dataset_directory);

        if (analysis_job.quantification_standard_filename.length() > 0)
        {
            perform_quantification(&analysis_job, false);
        }
        optimize_rois(analysis_job);
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

    set_streaming_options(clp, analysis_job);
    set_general_options(clp, analysis_job);
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

int run_h5_file_updates(Command_Line_Parser& clp)
{
    data_struct::Analysis_Job<float> analysis_job;

    set_detectors(clp, analysis_job);
    set_whole_command(clp, analysis_job);
    if (set_dir_and_files(clp, analysis_job) != -1)
    {


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
            std::string amps = clp.get_option("--update-amps");
            size_t idx = amps.find(',');
            if (idx != std::string::npos)
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
            std::string amps = clp.get_option("--update-quant-amps");
            size_t idx = amps.find(',');
            if (idx != std::string::npos)
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
    }
    return 0;
}

// ----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    
    std::string whole_command_line = "";
    for (int i = 0; i < argc; i++)
    {
        whole_command_line += std::string(argv[i]) + " ";
    }
    logI << whole_command_line << "\n";

    //Performance measure
    std::chrono::time_point<std::chrono::system_clock> start, end;

    // get location of where we are running from and use it to find ref files
    std::string exe_loc = std::string(argv[0]);
    int prog_idx = exe_loc.find("xrf_maps");
    if (prog_idx > 0)
    {
        exe_loc = exe_loc.substr(0, prog_idx);
    }

    //////// HENKE and ELEMENT INFO /////////////
    const std::string element_csv_filename = exe_loc + "../reference/xrf_library.csv";
    const std::string element_henke_filename = exe_loc + "../reference/henke.xdr";
    const std::string scaler_lookup_yaml = exe_loc + "../reference/Scaler_to_PV_map.yaml";

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

    if (clp.option_exists("--optimize-rois"))
    {
        run_optimize_rois(clp);
    }

    run_h5_file_updates(clp);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    logI << "=-=-=-=-=-=- Total elapsed time: " << elapsed_seconds.count() << "s =-=-=-=-=-=-=-\n" << std::endl; //endl will flush the print.


    data_struct::Element_Info_Map<float>::inst()->clear();
    data_struct::Element_Info_Map<double>::inst()->clear();

}
