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

#include "core/process_whole.h"

using namespace std::placeholders; //for _1, _2,

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

bool optimize_integrated_fit_params(data_struct::Analysis_Job<double> * analysis_job,
                                    data_struct::Spectra<double>& int_spectra,
                                    size_t detector_num,
                                    data_struct::Params_Override<double>* params_override,
                                    std::string save_filename,
                                    data_struct::Fit_Parameters<double>& out_fitp,
                                    Callback_Func_Status_Def* status_callback)
{
    fitting::models::Gaussian_Model<double> model;
    bool ret_val = false;

    if (params_override != nullptr)
    {
        //Range of energy in spectra to fit
        fitting::models::Range energy_range = data_struct::get_energy_range<double>(int_spectra.size(), &(params_override->fit_params));

        //Fitting routines
        fitting::routines::Param_Optimized_Fit_Routine<double>* fit_routine;


        if (analysis_job->optimize_fit_routine == OPTIMIZE_FIT_ROUTINE::HYBRID)
        {
            fit_routine = new fitting::routines::Hybrid_Param_NNLS_Fit_Routine<double>();
        }
        else
        {
            fit_routine = new fitting::routines::Param_Optimized_Fit_Routine<double>();
        }

        fit_routine->set_optimizer(analysis_job->optimizer());
        fit_routine->set_update_coherent_amplitude_on_fit(false);

        //reset model fit parameters to defaults
        model.reset_to_default_fit_params();
        //Update fit parameters by override values
        model.update_fit_params_values(&(params_override->fit_params));
        if (analysis_job->optimize_fit_params_preset != fitting::models::Fit_Params_Preset::NOT_SET)
        {
            //set fixed/fit preset
            model.set_fit_params_preset(analysis_job->optimize_fit_params_preset);
        }
        //Initialize the fit routine
        fit_routine->initialize(&model, &params_override->elements_to_fit, energy_range);

        //Fit the spectra saving the element counts in element_fit_count_dict
        fitting::optimizers::OPTIMIZER_OUTCOME outcome = fit_routine->fit_spectra_parameters(&model, &int_spectra, &params_override->elements_to_fit, analysis_job->use_weights, out_fitp, status_callback);
        out_fitp.add_parameter(data_struct::Fit_Param<double>(STR_OUTCOME, double(outcome)));
        std::string result = optimizer_outcome_to_str(outcome);
        logI << "Outcome = " << result << "\n";
        switch (outcome)
        {
        case fitting::optimizers::OPTIMIZER_OUTCOME::CONVERGED:
        case fitting::optimizers::OPTIMIZER_OUTCOME::F_TOL_LT_TOL:
        case fitting::optimizers::OPTIMIZER_OUTCOME::X_TOL_LT_TOL:
        case fitting::optimizers::OPTIMIZER_OUTCOME::G_TOL_LT_TOL:
        case fitting::optimizers::OPTIMIZER_OUTCOME::EXHAUSTED:
            // if we have a good fit, update our fit parameters so we are closer for the next fit
            //params_override->fit_params.update_values(&out_fitp);
            ret_val = true;
            break;
        case fitting::optimizers::OPTIMIZER_OUTCOME::CRASHED:
        case fitting::optimizers::OPTIMIZER_OUTCOME::EXPLODED:
        case fitting::optimizers::OPTIMIZER_OUTCOME::FAILED:
        case fitting::optimizers::OPTIMIZER_OUTCOME::FOUND_NAN:
        case fitting::optimizers::OPTIMIZER_OUTCOME::FOUND_ZERO:
        case fitting::optimizers::OPTIMIZER_OUTCOME::STOPPED:
        case fitting::optimizers::OPTIMIZER_OUTCOME::TRAPPED:
            ret_val = false;
            break;
        }
        io::file::save_optimized_fit_params(analysis_job->dataset_directory, save_filename, detector_num, result, &out_fitp, &int_spectra, &(params_override->elements_to_fit));

        delete fit_routine;
    }
    
    return ret_val;

}

// ----------------------------------------------------------------------------

void generate_optimal_params(data_struct::Analysis_Job<double>* analysis_job)
{
    std::unordered_map<int, data_struct::Fit_Parameters<double>> fit_params_avgs;
    std::unordered_map<int, data_struct::Params_Override<double>*> params;
    std::unordered_map<int, float> detector_file_cnt;
    data_struct::Params_Override<double>* params_override = nullptr;
    data_struct::Spectra<double> int_spectra;

    std::string full_path = analysis_job->dataset_directory + DIR_END_CHAR + "maps_fit_parameters_override.txt";

    for (size_t detector_num : analysis_job->detector_num_arr)
    {
        detector_file_cnt[detector_num] = 0.0;
    }

    for (auto& itr : analysis_job->optimize_dataset_files)
    {
        for (size_t detector_num : analysis_job->detector_num_arr)
        {
            //reuse previous param override if it exists
            if (params.count(detector_num) > 0)
            {
                params_override = params[detector_num];
            }
            else
            {
                params_override = new data_struct::Params_Override<double>();
                //load override parameters
                if (false == io::file::load_override_params(analysis_job->dataset_directory, detector_num, params_override))
                {
                    if (false == io::file::load_override_params(analysis_job->dataset_directory, -1, params_override))
                    {
                        logE << "Loading maps_fit_parameters_override.txt\n";
                        delete params_override;
                        params_override = nullptr;
                        continue;
                    }
                }
                if (params_override != nullptr)
                {
                    params[detector_num] = params_override;
                }
				
            }

            //load the int spectra from the dataset.
            if (io::file::load_and_integrate_spectra_volume(analysis_job->dataset_directory, itr, detector_num, &int_spectra, params_override))
            {
                data_struct::Fit_Parameters<double> out_fitp;
                if (optimize_integrated_fit_params(analysis_job, int_spectra, detector_num, params_override, itr, out_fitp, nullptr))
                {
                    detector_file_cnt[detector_num] += 1.0;
                    if (fit_params_avgs.count(detector_num) > 0)
                    {
                        fit_params_avgs[detector_num].sum_values(out_fitp);
                    }
                    else
                    {
                        fit_params_avgs[detector_num] = out_fitp;
                    }
                }
            }
            else
            {
                logE << "In optimize_integrated_dataset loading dataset" << itr << " for detector" << detector_num << "\n";
            }
        }
    }

    for(size_t detector_num : analysis_job->detector_num_arr)
	{
        if (detector_file_cnt[detector_num] > 1.)
        {
            fit_params_avgs[detector_num].divide_fit_values_by(detector_file_cnt[detector_num]);
        }

		io::file::aps::create_detector_fit_params_from_avg(full_path, fit_params_avgs[detector_num], detector_num);

        if (params.count(detector_num) > 0)
        {
            params_override = params[detector_num];
            if (params_override != nullptr)
            {
                delete params_override;
            }
            params.erase(detector_num);
        }
    }
}

// ----------------------------------------------------------------------------

void load_and_fit_quatification_datasets(data_struct::Analysis_Job<double>* analysis_job, size_t detector_num)
{
    fitting::models::Gaussian_Model<double> model;
    quantification::models::Quantification_Model<double> quantification_model;
    // Z number : count
    std::unordered_map<int, float> element_amt_in_all_standards;

    data_struct::Detector<double>* detector = analysis_job->get_detector(detector_num);
    data_struct::Params_Override<double>* override_params = &(detector->fit_params_override_dict);

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;


    for (Quantification_Standard<double>& standard_itr : analysis_job->standard_element_weights)
    {
        // detecotr_struct descructor will delete this memory
        detector->quantification_standards[standard_itr.standard_filename] = Quantification_Standard<double>(standard_itr.standard_filename, standard_itr.element_standard_weights);
        Quantification_Standard<double>* quantification_standard = &(detector->quantification_standards[standard_itr.standard_filename]);

        //Output of fits for elements specified
        std::unordered_map<std::string, data_struct::Fit_Element_Map<double>*> elements_to_fit;
        for (auto& itr : standard_itr.element_standard_weights)
        {
            data_struct::Element_Info<double>* e_info = data_struct::Element_Info_Map<double>::inst()->get_element(itr.first);
            elements_to_fit[itr.first] = new data_struct::Fit_Element_Map<double>(itr.first, e_info);
            elements_to_fit[itr.first]->init_energy_ratio_for_detector_element(detector->detector_element, standard_itr.disable_Ka_for_quantification, standard_itr.disable_La_for_quantification);

            if (element_amt_in_all_standards.count(e_info->number) > 0)
            {
                element_amt_in_all_standards[e_info->number] += 1.0;
            }
            else
            {
                element_amt_in_all_standards[e_info->number] = 1.0;
            }
        }

        std::unordered_map<std::string, double> pv_map;
        //load the quantification standard dataset
        size_t fn_str_len = quantification_standard->standard_filename.length();
        if (fn_str_len > 5 &&
            quantification_standard->standard_filename[fn_str_len - 4] == '.' &&
            quantification_standard->standard_filename[fn_str_len - 3] == 'm' &&
            quantification_standard->standard_filename[fn_str_len - 2] == 'c' &&
            quantification_standard->standard_filename[fn_str_len - 1] == 'a')
        {
            //try with adding detector_num on the end for 2ide datasets
            std::string qfilepath = analysis_job->dataset_directory + quantification_standard->standard_filename;
            if (detector_num != -1)
            {
                qfilepath += std::to_string(detector_num);
            }
            if (false == io::file::mca::load_integrated_spectra(qfilepath, &quantification_standard->integrated_spectra, pv_map))
            {
                //try without detector number on end 2idd
                if (false == io::file::mca::load_integrated_spectra(analysis_job->dataset_directory + quantification_standard->standard_filename, &quantification_standard->integrated_spectra, pv_map))
                {

                    //legacy code would load mca files, check for mca and replace with mda
                    size_t std_str_len = standard_itr.standard_filename.length();
                    if (standard_itr.standard_filename[std_str_len - 4] == '.' && standard_itr.standard_filename[std_str_len - 3] == 'm' && standard_itr.standard_filename[std_str_len - 2] == 'c' && standard_itr.standard_filename[std_str_len - 1] == 'a')
                    {
                        standard_itr.standard_filename[std_str_len - 2] = 'd';
                        quantification_standard->standard_filename = standard_itr.standard_filename;
                        if (false == io::file::load_and_integrate_spectra_volume(analysis_job->dataset_directory, quantification_standard->standard_filename, detector_num, &quantification_standard->integrated_spectra, override_params))
                        {
                            logE << "Could not load file " << standard_itr.standard_filename << " for detector" << detector_num << "\n";
                            continue;
                        }
                        else
                        {
                            quantification_standard->sr_current = override_params->sr_current;
                            quantification_standard->US_IC = override_params->US_IC;
                            quantification_standard->DS_IC = override_params->DS_IC;
                        }
                    }
                    else
                    {
                        logE << "Could not load file " << standard_itr.standard_filename << " for detector" << detector_num << "\n";
                        continue;
                    }
                }
                else
                {
                    find_quantifier_scalers(pv_map, quantification_standard);
                }
            }
            else
            {
                find_quantifier_scalers(pv_map, quantification_standard);
            }
        }
        else
        {
            if (false == io::file::load_and_integrate_spectra_volume(analysis_job->dataset_directory, quantification_standard->standard_filename, detector_num, &quantification_standard->integrated_spectra, override_params))
            {
                logE << "Could not load file " << standard_itr.standard_filename << " for detector " << detector_num << "\n";
                continue;
            }
            else
            {
                quantification_standard->sr_current = override_params->sr_current;
                quantification_standard->US_IC = override_params->US_IC;
                quantification_standard->DS_IC = override_params->DS_IC;
            }
        }
        
        if (quantification_standard->integrated_spectra.size() == 0)
        {
            logE << "Spectra size == 0! Can't process it!\n";
            continue;
        }

        //This is what IDL MAPS DID
        if (quantification_standard->sr_current == 0.0 )
        {
            quantification_standard->sr_current = 100.0;
        }

        energy_range = get_energy_range(quantification_standard->integrated_spectra.size(), &(override_params->fit_params));
        //First we integrate the spectra and get the elemental counts
        for (auto& fit_itr : detector->fit_routines)
        {

            fitting::routines::Base_Fit_Routine<double>* fit_routine = fit_itr.second;
            for (auto& el_itr : standard_itr.element_standard_weights)
            {
                quantification_standard->element_counts[fit_itr.first][el_itr.first] = 0;

                detector->append_element(fit_itr.first, STR_SR_CURRENT, el_itr.first, el_itr.second);
                detector->append_element(fit_itr.first, STR_US_IC, el_itr.first, el_itr.second);
                detector->append_element(fit_itr.first, STR_DS_IC, el_itr.first, el_itr.second);
            }

            //reset model fit parameters to defaults
            model.reset_to_default_fit_params();
            //Update fit parameters by override values
            model.update_fit_params_values(&(override_params->fit_params));
            //Initialize the fit routine
            fit_routine->initialize(&model, &elements_to_fit, energy_range);
            //Fit the spectra
            fit_routine->fit_spectra(&model, &quantification_standard->integrated_spectra, &elements_to_fit, quantification_standard->element_counts[fit_itr.first]);

            quantification_standard->normalize_counts_by_time(fit_itr.first);

            //Save csv and png if matrix or nnls
            if (fit_itr.first == Fitting_Routines::GAUSS_MATRIX || fit_itr.first == Fitting_Routines::NNLS)
            {
                Fit_Parameters<double> fit_params = model.fit_parameters();
                
                //add elements to fit parameters if they don't exist
                for (auto& itr2 : elements_to_fit)
                {
                    if (false == override_params->fit_params.contains(itr2.first))
                    {
                        fit_params.add_parameter(Fit_Param<double>(itr2.first, quantification_standard->element_counts.at(fit_itr.first).at(itr2.first)));
                    }
                }
                fitting::routines::Matrix_Optimized_Fit_Routine<double>* f_routine = (fitting::routines::Matrix_Optimized_Fit_Routine<double>*)fit_routine;
                double energy_offset = fit_params.value(STR_ENERGY_OFFSET);
                double energy_slope = fit_params.value(STR_ENERGY_SLOPE);
                double energy_quad = fit_params.value(STR_ENERGY_QUADRATIC);

                data_struct::ArrayTr<double> energy = data_struct::ArrayTr<double>::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
                data_struct::ArrayTr<double> ev = energy_offset + (energy * energy_slope) + (Eigen::pow(energy, (double)2.0) * energy_quad);
                data_struct::ArrayTr<double> sub_spectra = quantification_standard->integrated_spectra.segment(energy_range.min, energy_range.count());

                std::string full_path = analysis_job->dataset_directory + DIR_END_CHAR + "output" + DIR_END_CHAR + "calib_" + fit_routine->get_name() + "_" + standard_itr.standard_filename;
                if (detector_num != -1)
                {
                    full_path += std::to_string(detector_num);
                }
                logI << full_path << "\n";
                #ifdef _BUILD_WITH_QT
                visual::SavePlotSpectrasFromConsole(full_path + ".png", &ev, &sub_spectra, (ArrayTr<double>*)(&f_routine->fitted_integrated_spectra()), (ArrayTr<double>*)(&f_routine->fitted_integrated_background()), true);
                #endif

                io::file::csv::save_fit_and_int_spectra(full_path + ".csv", &ev, &sub_spectra, (ArrayTr<double>*)(&f_routine->fitted_integrated_spectra()), (ArrayTr<double>*)(&f_routine->fitted_integrated_background()));
            }

            detector->update_element_quants(fit_itr.first, STR_SR_CURRENT, quantification_standard, &quantification_model, quantification_standard->sr_current);
            detector->update_element_quants(fit_itr.first, STR_US_IC, quantification_standard, &quantification_model, quantification_standard->US_IC);
            detector->update_element_quants(fit_itr.first, STR_DS_IC, quantification_standard, &quantification_model, quantification_standard->DS_IC);
        }

        //cleanup
        for (auto& itr3 : elements_to_fit)
        {
            delete itr3.second;
        }
        elements_to_fit.clear();
    }

    float divisor = (float)analysis_job->standard_element_weights.size();    
    
    if (divisor > 1.0)
    {
        for (auto& fit_itr : detector->fit_routines)
        {
            detector->avg_element_quants(fit_itr.first, STR_SR_CURRENT, element_amt_in_all_standards);
        }
    }
}

// ----------------------------------------------------------------------------

bool perform_quantification(data_struct::Analysis_Job<double>* analysis_job, bool save_when_done)
{
    quantification::models::Quantification_Model<double> quantification_model;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    logI << "Perform_quantification()"<<"\n";

    std::vector<std::string> quant_scaler_name_list = { STR_SR_CURRENT, STR_US_IC, STR_DS_IC };

    if( io::file::load_quantification_standardinfo(analysis_job->dataset_directory, analysis_job->quantification_standard_filename, analysis_job->standard_element_weights) )
    {
        for(size_t detector_num : analysis_job->detector_num_arr)
        {
            data_struct::Detector<double>* detector = analysis_job->get_detector(detector_num);
            data_struct::Params_Override<double>* override_params = &(detector->fit_params_override_dict);

            
            load_and_fit_quatification_datasets(analysis_job, detector_num);
            detector->generage_avg_quantification_scalers();

            for(auto &fit_itr : detector->fit_routines)
            {
               fitting::optimizers::Optimizer<double>* optimizer = analysis_job->optimizer();
               for (auto& quant_itr : detector->avg_quantification_scaler_map)
               {
                    fitting::routines::Base_Fit_Routine<double>* fit_routine = fit_itr.second;

                    logI << Fitting_Routine_To_Str.at(fit_itr.first) << " " << quant_itr.first << "\n";
                    Fit_Parameters<double> fit_params;
                    // min, and max values doen't matter because we are free fitting in mpfit and lmfit
                    fit_params.add_parameter(Fit_Param<double>("quantifier", 0.0, std::numeric_limits<double>::max(), 1.0, 0.0001, E_Bound_Type::FIT));
                    //initial guess: parinfo_value[0] = 100000.0 / factor
                    fit_params["quantifier"].value = (double)100000.0 / quant_itr.second;
                    optimizer->minimize_quantification(&fit_params, &detector->all_element_quants[fit_itr.first][quant_itr.first], &quantification_model);
                    double val = fit_params["quantifier"].value;

                    if(false == std::isfinite(val))
                    {
                        logW<<"Quantifier Value = Inf. setting it to 0.\n";
                        val = 0;
                    }
                    else
                    {
                        logI<<"Quantifier Value = "<<val<<"\n";
                    }

                    detector->update_calibration_curve(fit_itr.first, quant_itr.first, &quantification_model, val);
                }
            }

            if (save_when_done)
            {
                io::file::save_quantification_plots(analysis_job->dataset_directory, detector);
            }
        }
    }
    else
    {
        logE<<"Loading quantification standard "<<analysis_job->quantification_standard_filename<<"\n";
        return false;
    }

    if (save_when_done)
    {
        //Save quantification to each file
        for (auto& dataset_file : analysis_job->dataset_files)
        {
            for (size_t detector_num : analysis_job->detector_num_arr)
            {
                data_struct::Detector<double>* detector = analysis_job->get_detector(detector_num);
                std::string str_detector_num = "";
                if (detector_num != -1)
                {
                    str_detector_num = std::to_string(detector_num);
                }

                std::string full_save_path = analysis_job->dataset_directory + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + dataset_file;

                size_t fn_str_len = dataset_file.length();
                if (fn_str_len > 5 &&
                    dataset_file[fn_str_len - 4] == '.' &&
                    dataset_file[fn_str_len - 3] == 'm' &&
                    dataset_file[fn_str_len - 2] == 'd' &&
                    dataset_file[fn_str_len - 1] == 'a')
                {
                    full_save_path += ".h5" + str_detector_num;
                }
                if (io::file::HDF5_IO::inst()->start_save_seq(full_save_path, false, true))
                {
                    if (false == io::file::HDF5_IO::inst()->save_quantification(detector))
                    {
                        logE << "Error saving quantification to hdf5\n";
                    }
                    io::file::HDF5_IO::inst()->end_save_seq();
                }
            }
        }
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logI << "quantification elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;

}

// ----------------------------------------------------------------------------

void find_and_optimize_roi(data_struct::Analysis_Job<double>& analysis_job,
                            int detector_num,
                            std::map<std::string, std::vector<std::pair<int, int>>>& rois,
                            std::unordered_map<std::string, data_struct::Spectra<double> >& int_spectra_map,
                            std::string search_filename,
                            std::map<std::string, data_struct::Fit_Parameters<double>>& out_roi_fit_params,
                            Callback_Func_Status_Def* status_callback)
{
    int cnt = int_spectra_map.count(search_filename);
    std::vector<std::string> files = io::file::File_Scan::inst()->find_all_dataset_files(analysis_job.dataset_directory + "img.dat", search_filename);
    if (files.size() > 0 || cnt > 0)
    {
        std::string sfile_name;
        for (auto& roi_itr : rois)
        {
            data_struct::Spectra<double> int_spectra;

            std::string file_path = analysis_job.dataset_directory + "img.dat" + DIR_END_CHAR; 
            if (cnt > 0)
            {
                int_spectra = int_spectra_map.at(search_filename);
                sfile_name = search_filename;
                file_path += sfile_name;
            }
            else
            {
                if (files.size() == 1) // v9 will find just 1, 
                {
                    sfile_name = files[0];
                }
                else // v10 finds more so we just use search name
                {
                    sfile_name = search_filename;
                }
                file_path += sfile_name;
                if (false == io::file::HDF5_IO::inst()->load_integrated_spectra_analyzed_h5_roi(file_path, &int_spectra, roi_itr.second))
                {
                    logE << "Could not load int spectra for " << file_path << ".  skipping..\n";
                    continue;
                }
            }

            data_struct::Detector<double>* detector = analysis_job.get_detector(detector_num);
            if (detector == nullptr)
            {
                logE << "Detector " << detector_num << " is not initialized, skipping..\n";
                continue;
            }
                
            data_struct::Params_Override<double>* params_override = &(detector->fit_params_override_dict);
                
            /// If quant is not done also, Need to save more info from Quantification first, then we can finish implementing loading quant from hdf5 instead of rerunning each time roi's are done
            ///io::file::HDF5_IO::inst()->load_quantification_analyzed_h5(file_path, detector);
            // load scalers and scan info
            data_struct::Scan_Info<double> scan_info;
            io::file::HDF5_IO::inst()->load_scan_info_analyzed_h5(file_path, detector, scan_info);

            double sr_current = 1.0;
            double us_ic = 1.0;
            double ds_ic = 1.0;

            std::map<std::string, data_struct::ArrayXXr<double>> scalers_map;
            if (io::file::HDF5_IO::inst()->load_scalers_analyzed_h5(file_path, scalers_map))
            {
                for (auto& in_itr : roi_itr.second)
                {
                    hsize_t xoffset = in_itr.first;
                    hsize_t yoffset = in_itr.second;
                    if (scalers_map.count(STR_DS_IC) > 0)
                    {
                        ds_ic += scalers_map.at(STR_DS_IC)(yoffset, xoffset);
                    }
                    if (scalers_map.count(STR_US_IC) > 0)
                    {
                        us_ic += scalers_map.at(STR_US_IC)(yoffset, xoffset);
                    }
                    if (scalers_map.count(STR_SR_CURRENT) > 0)
                    {
                        sr_current += scalers_map.at(STR_SR_CURRENT)(yoffset, xoffset);
                    }
                }
            }
            else
            {
                sr_current = detector->fit_params_override_dict.sr_current;
                us_ic = detector->fit_params_override_dict.US_IC;
                ds_ic = detector->fit_params_override_dict.DS_IC;
            }


            data_struct::Fit_Parameters<double> out_fitp;
            std::string roi_name = roi_itr.first;
            if (false == optimize_integrated_fit_params(&analysis_job, int_spectra, detector_num, params_override, sfile_name + "_roi_" + roi_name + "_det_", out_fitp, status_callback))
            {
                logE << "Failed to optimize ROI "<< file_path<<" : "<< roi_name<<".\n";
            }

                
            double total_counts = 0.0;
            for (auto& itr : out_fitp)
            {
                if (data_struct::Element_Info_Map<double>::inst()->is_element(itr.first))
                {
                    total_counts += itr.second.value;
                }
            }
                

            double abs_err = std::abs(out_fitp.at(STR_RESIDUAL).value);
            double rel_err = abs_err / int_spectra.sum();;
            double roi_area = 0;
            if (scan_info.meta_info.x_axis.rows() > 0 && scan_info.meta_info.x_axis.cols() > 0
                && scan_info.meta_info.y_axis.rows() > 0 && scan_info.meta_info.y_axis.cols() > 0)
            {
                roi_area = roi_itr.second.size() * 1000.0 * 1000.0 * (scan_info.meta_info.x_axis.maxCoeff() - scan_info.meta_info.x_axis.minCoeff()) / (scan_info.meta_info.x_axis.size() - 1) * (scan_info.meta_info.y_axis.maxCoeff() - scan_info.meta_info.y_axis.minCoeff()) / (scan_info.meta_info.y_axis.size() - 1);
            }
            // add in other properties that will be saved to csv
            out_fitp.add_parameter(data_struct::Fit_Param<double>("real_time", int_spectra.elapsed_realtime()));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("live_time", int_spectra.elapsed_livetime()));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("SRcurrent", sr_current));
            out_fitp.add_parameter(data_struct::Fit_Param<double>(STR_US_IC, us_ic));
            out_fitp.add_parameter(data_struct::Fit_Param<double>(STR_DS_IC, ds_ic));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("total_counts", total_counts));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("status", out_fitp.at(STR_OUTCOME).value));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("niter", out_fitp.at(STR_NUM_ITR).value));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("total_perror", out_fitp.at(STR_RESIDUAL).value));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("abs_error", abs_err));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("relative_error", rel_err));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("roi_areas", roi_area));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("roi_pixels", roi_itr.second.size()));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("US_num", params_override->us_amp_sens_num));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("US_unit", io::file::translate_back_sens_unit<double>(params_override->us_amp_sens_unit)));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("US_sensfactor", params_override->us_amp_sens_num));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("DS_num", params_override->ds_amp_sens_num));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("DS_unit", io::file::translate_back_sens_unit<double>(params_override->ds_amp_sens_unit)));
            out_fitp.add_parameter(data_struct::Fit_Param<double>("DS_sensfactor", params_override->ds_amp_sens_num));
            out_roi_fit_params[sfile_name + "_roi_" + roi_name] = out_fitp;
        }
    }
    else
    {
        logE << "Could not find dataset file for " << search_filename << "\n";
    }
}

// ----------------------------------------------------------------------------

void optimize_single_roi(data_struct::Analysis_Job<double>& analysis_job,
    std::string roi_file_name,
    std::map<int, std::map<std::string,
    data_struct::Fit_Parameters<double>>> & out_roi_fit_params,
    Callback_Func_Status_Def* status_callback)
{
    //std::map<int, std::vector<std::pair<unsigned int, unsigned int>>> rois;
    std::map<std::string, std::vector<std::pair<int, int>>> rois;
    std::string search_filename;
    // int specs loaded from v10 version
    std::unordered_map<std::string, data_struct::Spectra<double> > int_specs;

    int slen = roi_file_name.size();
    if (slen < 6)
    {
        logE << "Roi file name too short " << roi_file_name << ". Skipping file.\n";
        return;
    }
    std::string base_file_name = roi_file_name.substr(0, slen - 4);

    if (roi_file_name[slen - 4] == '.' && roi_file_name[slen - 3] == 'r' && roi_file_name[slen - 2] == '0' && roi_file_name[slen - 1] == 'i')
    {
        if (io::file::aps::load_v10_rois(analysis_job.dataset_directory + "rois" + DIR_END_CHAR + roi_file_name, rois, int_specs))
        {
            for (size_t detector_num : analysis_job.detector_num_arr)
            {
                if (detector_num != -1)
                { 
                    search_filename = "";
                    std::string str_detector_num = std::to_string(detector_num);
                    // search for detector in list of int specs loaded
                    for (const auto& spec_itr : int_specs)
                    {
                        int slen = spec_itr.first.length();
                        if (slen > 0 && spec_itr.first[slen - 1] == str_detector_num[0])
                        {
                            search_filename = spec_itr.first;
                            break;
                        }
                    }
                    if (search_filename.length() == 0)
                    {
                        search_filename = base_file_name + ".mda.h5" + str_detector_num;
                    }
                    find_and_optimize_roi(analysis_job, detector_num, rois, int_specs, search_filename, out_roi_fit_params[detector_num], status_callback);
                }
            }
            //now for the avg h5
            search_filename = base_file_name + ".mda.h5";
            find_and_optimize_roi(analysis_job, -1, rois, int_specs, search_filename, out_roi_fit_params[-1], status_callback);
        }
        else
        {
            logE << "Error loading roi file " << roi_file_name << "\n";
        }
    }
    else
    {
        if (io::file::aps::load_v9_rois(analysis_job.dataset_directory + "rois" + DIR_END_CHAR + roi_file_name, rois))
        {
            logI << "Loaded roi file " << roi_file_name << "\n";
            //get dataset file number
            std::string dataset_num;
            std::istringstream strstream(roi_file_name);
            if (roi_file_name.find('_') != std::string::npos)
            {
                std::getline(strstream, dataset_num, '_');
            }
            else if (roi_file_name.find('-') != std::string::npos)
            {
                std::getline(strstream, dataset_num, '-');
            }
            if (dataset_num.length() > 0)
            {
                // iterate through all 
                for (size_t detector_num : analysis_job.detector_num_arr)
                {
                    //detector = analysis_job.get_detector(detector_num);
                    if (detector_num != -1)
                    {
                        std::string str_detector_num = std::to_string(detector_num);
                        search_filename = dataset_num + ".h5" + str_detector_num; // v9 save
                        find_and_optimize_roi(analysis_job, detector_num, rois, int_specs, search_filename, out_roi_fit_params[detector_num], status_callback);
                        search_filename = dataset_num + ".mda.h5" + str_detector_num;
                        find_and_optimize_roi(analysis_job, detector_num, rois, int_specs, search_filename, out_roi_fit_params[detector_num], status_callback);
                    }
                }
                //now for the avg h5
                //detector = analysis_job.get_detector(-1); 
                search_filename = dataset_num + ".h5"; // v9 save
                find_and_optimize_roi(analysis_job, -1, rois, int_specs, search_filename, out_roi_fit_params[-1], status_callback);
                search_filename = dataset_num + ".mda.h5";
                find_and_optimize_roi(analysis_job, -1, rois, int_specs, search_filename, out_roi_fit_params[-1], status_callback);
            }
            else
            {
                logE << "Error parsing dataset number for " << roi_file_name << "\n";
            }

        }
        else
        {
            logE << "Error loading roi file " << roi_file_name << "\n";
        }
    }
}

// ----------------------------------------------------------------------------

void optimize_rois(data_struct::Analysis_Job<double>& analysis_job)
{
     //       detector_num        file_name_roi           fit_params
    std::map<int, std::map<std::string, data_struct::Fit_Parameters<double>>> roi_fit_params;
    // old v9 format
    std::vector<std::string> files_to_proc = io::file::File_Scan::inst()->find_all_dataset_files(analysis_job.dataset_directory + "rois", ".roi");
    // new roi format
    std::vector<std::string> r0i_files = io::file::File_Scan::inst()->find_all_dataset_files(analysis_job.dataset_directory + "rois", ".r0i");

    for (auto itr : r0i_files)
    {
        files_to_proc.push_back(itr);
    }

    for (auto& itr : files_to_proc)
    {
        optimize_single_roi(analysis_job, itr, roi_fit_params);
    }

    // save all to csv
    for (auto detector_num : analysis_job.detector_num_arr)
    {
        std::string save_path = analysis_job.dataset_directory + "output/specfit_results" + std::to_string(detector_num) + ".csv";
        std::string quant_save_path = analysis_job.dataset_directory + "output/specfit_results" + std::to_string(detector_num) + "_quantified.csv";
        io::file::csv::save_v9_specfit(save_path, roi_fit_params.at(detector_num));
        if (analysis_job.get_detector(detector_num)->quantification_standards.size() > 0)
        {
            io::file::csv::save_v9_specfit_quantified(quant_save_path, analysis_job.get_detector(detector_num), analysis_job.fitting_routines, roi_fit_params.at(detector_num));
        }
    }
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
