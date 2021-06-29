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

template<typename T>
data_struct::Fit_Count_Dict* generate_fit_count_dict(std::unordered_map<std::string, T> *elements_to_fit, size_t height, size_t width, bool alloc_iter_count)
{
    data_struct::Fit_Count_Dict* element_fit_counts_dict = new data_struct::Fit_Count_Dict();
    for(auto& e_itr : *elements_to_fit)
    {
        element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr >(e_itr.first, data_struct::ArrayXXr()) );
        element_fit_counts_dict->at(e_itr.first).resize(height, width);
    }

    if (alloc_iter_count)
    {
        //Allocate memeory to save number of fit iterations
        element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr >(STR_NUM_ITR, data_struct::ArrayXXr() ));
        element_fit_counts_dict->at(STR_NUM_ITR).resize(height, width);
        element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr >(STR_RESIDUAL, data_struct::ArrayXXr()));
        element_fit_counts_dict->at(STR_RESIDUAL).resize(height, width);
    }

	//  TOTAL_FLUORESCENCE_YIELD
	element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr >(STR_TOTAL_FLUORESCENCE_YIELD, data_struct::ArrayXXr()));
	element_fit_counts_dict->at(STR_TOTAL_FLUORESCENCE_YIELD).resize(height, width);

	//SUM_ELASTIC_INELASTIC
	element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr >(STR_SUM_ELASTIC_INELASTIC_AMP , data_struct::ArrayXXr()));
	element_fit_counts_dict->at(STR_SUM_ELASTIC_INELASTIC_AMP).resize(height, width);


    return element_fit_counts_dict;
}

// ----------------------------------------------------------------------------

bool fit_single_spectra(fitting::routines::Base_Fit_Routine * fit_routine,
                        const fitting::models::Base_Model * const model,
                        const data_struct::Spectra * const spectra,
                        const data_struct::Fit_Element_Map_Dict * const elements_to_fit,
                        data_struct::Fit_Count_Dict * out_fit_counts,
                        size_t i,
                        size_t j)
{
    std::unordered_map<std::string, real_t> counts_dict;
    fit_routine->fit_spectra(model, spectra, elements_to_fit, counts_dict);
    //save count / sec
    for (auto& el_itr : *elements_to_fit)
    {
        (*out_fit_counts)[el_itr.first](i,j) = counts_dict[el_itr.first] / spectra->elapsed_livetime();
    }
    (*out_fit_counts)[STR_NUM_ITR](i,j) = counts_dict[STR_NUM_ITR];

    (*out_fit_counts)[STR_RESIDUAL](i, j) = counts_dict[STR_RESIDUAL];
	// add total fluorescense yield
	if (out_fit_counts->count(STR_TOTAL_FLUORESCENCE_YIELD))
	{
		(*out_fit_counts)[STR_TOTAL_FLUORESCENCE_YIELD](i, j) = spectra->sum();
	}
	// add sum coherent and compton
	if (out_fit_counts->count(STR_SUM_ELASTIC_INELASTIC_AMP) > 0  && counts_dict.count(STR_COHERENT_SCT_AMPLITUDE) > 0 && counts_dict.count(STR_COMPTON_AMPLITUDE) > 0)
	{
		(*out_fit_counts)[STR_SUM_ELASTIC_INELASTIC_AMP](i, j) = counts_dict[STR_COHERENT_SCT_AMPLITUDE] + counts_dict[STR_COMPTON_AMPLITUDE];
        // add total fluorescense yield
        if (out_fit_counts->count(STR_TOTAL_FLUORESCENCE_YIELD))
        {                   //                                      (sum - (elastic + inelastic)) / live time
            (*out_fit_counts)[STR_TOTAL_FLUORESCENCE_YIELD](i, j) = ( spectra->sum() - (*out_fit_counts)[STR_SUM_ELASTIC_INELASTIC_AMP](i, j) ) / spectra->elapsed_livetime();
        }
	}
    else
    {
        // add total fluorescense yield
        if (out_fit_counts->count(STR_TOTAL_FLUORESCENCE_YIELD))
        {
            (*out_fit_counts)[STR_TOTAL_FLUORESCENCE_YIELD](i, j) = spectra->sum() / spectra->elapsed_livetime();
        }
    }

    return true;
}

// ----------------------------------------------------------------------------

bool optimize_integrated_fit_params(std::string dataset_directory,
                                    std::string  dataset_filename,
                                    size_t detector_num,
                                    data_struct::Params_Override* params_override,
                                    fitting::models::Fit_Params_Preset optimize_fit_params_preset,
                                    fitting::optimizers::Optimizer* optimizer,
                                    data_struct::Fit_Parameters& out_fitp)
{
    fitting::models::Gaussian_Model model;
    bool ret_val = false;
    data_struct::Spectra int_spectra;

    if (params_override != nullptr)
    {
        //load the quantification standard dataset
        if (false == io::load_and_integrate_spectra_volume(dataset_directory, dataset_filename, detector_num, &int_spectra, params_override))
        {
            logE << "In optimize_integrated_dataset loading dataset" << dataset_filename << " for detector" << detector_num << "\n";
            return false;
        }
        //Range of energy in spectra to fit
        fitting::models::Range energy_range = data_struct::get_energy_range(int_spectra.size(), &(params_override->fit_params));

        //Fitting routines
        fitting::routines::Param_Optimized_Fit_Routine fit_routine;
        fit_routine.set_optimizer(optimizer);
		fit_routine.set_update_coherent_amplitude_on_fit(false);

        //reset model fit parameters to defaults
        model.reset_to_default_fit_params();
        //Update fit parameters by override values
        model.update_fit_params_values(&(params_override->fit_params));
        //set fixed/fit preset
        model.set_fit_params_preset(optimize_fit_params_preset);

        //Initialize the fit routine
        fit_routine.initialize(&model, &params_override->elements_to_fit, energy_range);

        data_struct::Fit_Parameters out_fitp;
        //Fit the spectra saving the element counts in element_fit_count_dict
        fitting::optimizers::OPTIMIZER_OUTCOME outcome = fit_routine.fit_spectra_parameters(&model, &int_spectra, &params_override->elements_to_fit, out_fitp);
        switch (outcome)
        {
        case fitting::optimizers::OPTIMIZER_OUTCOME::CONVERGED:
            // if we have a good fit, update our fit parameters so we are closer for the next fit
            params_override->fit_params.append_and_update(&out_fitp);
            ret_val = true;
            break;
        case fitting::optimizers::OPTIMIZER_OUTCOME::EXHAUSTED:
        case fitting::optimizers::OPTIMIZER_OUTCOME::F_TOL_LT_TOL:
        case fitting::optimizers::OPTIMIZER_OUTCOME::X_TOL_LT_TOL:
        case fitting::optimizers::OPTIMIZER_OUTCOME::G_TOL_LT_TOL:
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
        io::save_optimized_fit_params(dataset_directory, dataset_filename, detector_num, &out_fitp, &int_spectra, &(params_override->elements_to_fit));
    }
    
    return ret_val;

}

// ----------------------------------------------------------------------------

void generate_optimal_params(data_struct::Analysis_Job* analysis_job)
{
    std::unordered_map<int, data_struct::Fit_Parameters> fit_params_avgs;
    std::unordered_map<int, data_struct::Params_Override*> params;
    std::unordered_map<int, float> detector_file_cnt;
    data_struct::Params_Override* params_override = nullptr;

    for (size_t detector_num : analysis_job->detector_num_arr)
    {
        detector_file_cnt[detector_num] = 0.0;
    }


    for(auto &itr : analysis_job->optimize_dataset_files)
    {
        for(size_t detector_num : analysis_job->detector_num_arr)
        {
            //reuse previous param override if it exists
            if (params.count(detector_num) > 0)
            {
                params_override = params[detector_num];
            }
            else
            {
                params_override = new data_struct::Params_Override();
                //load override parameters
                if (false == io::load_override_params(analysis_job->dataset_directory, detector_num, params_override))
                {
                    if (false == io::load_override_params(analysis_job->dataset_directory, -1, params_override))
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

            data_struct::Fit_Parameters out_fitp;
            if (optimize_integrated_fit_params(analysis_job->dataset_directory, itr, detector_num, params_override, analysis_job->optimize_fit_params_preset, analysis_job->optimizer(), out_fitp))
            {
                detector_file_cnt[detector_num] += 1.0;
                if (fit_params_avgs.count(detector_num) > 0)
                {
                    fit_params_avgs[detector_num].sum_values(params_override->fit_params);
                }
                else
                {
                    fit_params_avgs[detector_num] = params_override->fit_params;
                }
            }
        }
    }
    for(size_t detector_num : analysis_job->detector_num_arr)
    {
        if (detector_file_cnt[detector_num] > 0.)
        {
            fit_params_avgs[detector_num].divide_fit_values_by(detector_file_cnt[detector_num]);
        }
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

    io::save_averaged_fit_params(analysis_job->dataset_directory, fit_params_avgs, analysis_job->detector_num_arr);

}

// ----------------------------------------------------------------------------

void proc_spectra(data_struct::Spectra_Volume* spectra_volume,
                  data_struct::Detector * detector,
                  ThreadPool* tp,
                  bool save_spec_vol,
                  Callback_Func_Status_Def* status_callback)
{
    if (detector == nullptr)
    {
        logE << "Detector meta information not loaded. Cannot process!\n";
        return;
    }

    if (spectra_volume == nullptr)
    {
        logE << "Spectra Volume not loaded. Cannot process!\n";
        return;
    }

    data_struct::Params_Override * override_params = &(detector->fit_params_override_dict);

    //Range of energy in spectra to fit
    fitting::models::Range energy_range = data_struct::get_energy_range(spectra_volume->samples_size(), &(detector->fit_params_override_dict.fit_params) );

    std::chrono::time_point<std::chrono::system_clock> start, end;

    for(auto &itr : detector->fit_routines)
    {
        fitting::routines::Base_Fit_Routine *fit_routine = itr.second;

        logI << "Processing  "<< fit_routine->get_name()<<"\n";

        start = std::chrono::system_clock::now();
        if (override_params->elements_to_fit.size() < 1)
        {
            logE<<"No elements to fit. Check  maps_fit_parameters_override.txt0 - 3 exist"<<"\n";
            continue;
        }

		//Fit job queue
		std::queue<std::future<bool> >* fit_job_queue = new std::queue<std::future<bool> >();

        //Allocate memeory to save fit counts
        data_struct::Fit_Count_Dict  *element_fit_count_dict = generate_fit_count_dict(&override_params->elements_to_fit, spectra_volume->rows(), spectra_volume->cols(), true);

        for(size_t i=0; i<spectra_volume->rows(); i++)
        {
            for(size_t j=0; j<spectra_volume->cols(); j++)
            {
                //logD<< i<<" "<<j<<"\n";
                fit_job_queue->emplace( tp->enqueue(fit_single_spectra, fit_routine, detector->model, &(*spectra_volume)[i][j], &override_params->elements_to_fit, element_fit_count_dict, i, j) );
            }
        }

        size_t total_blocks = (spectra_volume->rows() * spectra_volume->cols()) - 1;
        size_t cur_block = 0;
        //wait for queue to finish processing
        while(!fit_job_queue->empty())
        {
            auto ret = std::move(fit_job_queue->front());
            fit_job_queue->pop();
            ret.get();
            if (status_callback != nullptr)
            {
                (*status_callback)(cur_block, total_blocks);
            }
            cur_block++;
        }

        std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end-start;
        logI << "Fitting [ "<< fit_routine->get_name() <<" ] elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

        io::file::HDF5_IO::inst()->save_element_fits(fit_routine->get_name(), element_fit_count_dict);

        if(itr.first == data_struct::Fitting_Routines::GAUSS_MATRIX || itr.first == data_struct::Fitting_Routines::NNLS)
        {
            fitting::routines::Matrix_Optimized_Fit_Routine* matrix_fit = (fitting::routines::Matrix_Optimized_Fit_Routine*)fit_routine;
            io::file::HDF5_IO::inst()->save_fitted_int_spectra( fit_routine->get_name(),
																matrix_fit->fitted_integrated_spectra(),
																matrix_fit->energy_range(),
                                                                matrix_fit->fitted_integrated_background(),
																(*spectra_volume)[0][0].size());
        }
		if (itr.first == data_struct::Fitting_Routines::GAUSS_MATRIX)
		{
			fitting::routines::Matrix_Optimized_Fit_Routine* matrix_fit = (fitting::routines::Matrix_Optimized_Fit_Routine*)fit_routine;
			io::file::HDF5_IO::inst()->save_max_10_spectra(fit_routine->get_name(),
																matrix_fit->energy_range(),
																matrix_fit->max_integrated_spectra(),
																matrix_fit->max_10_integrated_spectra(),
                                                                matrix_fit->fitted_integrated_background());
		}

        delete fit_job_queue;
        element_fit_count_dict->clear();
        delete element_fit_count_dict;
    }

    real_t energy_offset = 0.0;
    real_t energy_slope = 0.0;
    real_t energy_quad = 0.0;
    data_struct::Fit_Parameters fit_params = detector->model->fit_parameters();
    if(fit_params.contains(STR_ENERGY_OFFSET))
    {
        energy_offset = fit_params[STR_ENERGY_OFFSET].value;
    }
    if(fit_params.contains(STR_ENERGY_SLOPE))
    {
        energy_slope = fit_params[STR_ENERGY_SLOPE].value;
    }
    if(fit_params.contains(STR_ENERGY_QUADRATIC))
    {
        energy_quad = fit_params[STR_ENERGY_QUADRATIC].value;
    }

    if(save_spec_vol)
    {
        io::file::HDF5_IO::inst()->save_quantification(detector);
        io::save_volume(spectra_volume, energy_offset, energy_slope, energy_quad);
        io::file::HDF5_IO::inst()->end_save_seq();
    }
    else
    {
        io::file::HDF5_IO::inst()->save_quantification(detector);
        io::file::HDF5_IO::inst()->end_save_seq();
    }

}

// ----------------------------------------------------------------------------

void process_dataset_files(data_struct::Analysis_Job* analysis_job, Callback_Func_Status_Def* status_callback)
{
    ThreadPool tp(analysis_job->num_threads);

    for(auto &dataset_file : analysis_job->dataset_files)
    {
        //if quick and dirty then sum all detectors to 1 spectra volume and process it
        if(analysis_job->quick_and_dirty)
        {
            process_dataset_files_quick_and_dirty(dataset_file, analysis_job, tp);
        }
        //otherwise process each detector separately
        else
        {
            for(size_t detector_num : analysis_job->detector_num_arr)
            {

                data_struct::Detector* detector = analysis_job->get_detector(detector_num);

                //Spectra volume data
                data_struct::Spectra_Volume* spectra_volume = new data_struct::Spectra_Volume();

                std::string fullpath;
                size_t dlen = dataset_file.length();
                if (dataset_file[dlen - 4] == '.' && dataset_file[dlen - 3] == 'm' && dataset_file[dlen - 2] == 'd' && dataset_file[dlen - 1] == 'a')
                {
                    std::string str_detector_num = "";
                    if (detector_num != -1)
                    {
                        str_detector_num = std::to_string(detector_num);
                    }
                    std::string full_save_path = analysis_job->dataset_directory + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + dataset_file + ".h5" + str_detector_num;
                    io::file::HDF5_IO::inst()->set_filename(full_save_path);
                }
                else
                {
                    std::string full_save_path = analysis_job->dataset_directory + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + dataset_file;
                    io::file::HDF5_IO::inst()->set_filename(full_save_path);
                }
                
                bool loaded_from_analyzed_hdf5 = false;
                //load spectra volume
                if (false == io::load_spectra_volume(analysis_job->dataset_directory, dataset_file, detector_num, spectra_volume, &detector->fit_params_override_dict, &loaded_from_analyzed_hdf5, true) )
                {
                    logW<<"Skipping detector "<<detector_num<<"\n";
                    delete spectra_volume;
                    if (status_callback != nullptr)
                    {
                        (*status_callback)(0, 1);
                    }
                    continue;
                }

                analysis_job->init_fit_routines(spectra_volume->samples_size(), true);
                proc_spectra(spectra_volume, detector, &tp, !loaded_from_analyzed_hdf5, status_callback);
				delete spectra_volume;
            }
        }
    }
}

// ----------------------------------------------------------------------------

void process_dataset_files_quick_and_dirty(std::string dataset_file, data_struct::Analysis_Job* analysis_job, ThreadPool &tp, Callback_Func_Status_Def* status_callback)
{
    std::string full_save_path = analysis_job->dataset_directory + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + dataset_file + ".h5";

    data_struct::Detector* detector = analysis_job->get_detector(0);
    //Spectra volume data
    data_struct::Spectra_Volume* spectra_volume = new data_struct::Spectra_Volume();
    data_struct::Spectra_Volume* tmp_spectra_volume = new data_struct::Spectra_Volume();

    io::file::HDF5_IO::inst()->start_save_seq(full_save_path, true); // force to create new file for quick and dirty

    //load the first one
    size_t detector_num = analysis_job->detector_num_arr[0];
    bool is_loaded_from_analyzed_h5 = false;
    if (false == io::load_spectra_volume(analysis_job->dataset_directory, dataset_file, detector_num, spectra_volume, &detector->fit_params_override_dict, &is_loaded_from_analyzed_h5, true))
    {
        logE << "Loading all detectors for " << analysis_job->dataset_directory << DIR_END_CHAR << dataset_file << "\n";
        delete spectra_volume;
        delete tmp_spectra_volume;
        if (status_callback != nullptr)
        {
            (*status_callback)(0, 1);
        }
        return;
    }

    //load spectra volume
    for (int i = 1; i < analysis_job->detector_num_arr.size(); i++)
    {
        if (false == io::load_spectra_volume(analysis_job->dataset_directory, dataset_file, analysis_job->detector_num_arr[i], tmp_spectra_volume, &detector->fit_params_override_dict, &is_loaded_from_analyzed_h5, false))
        {
            logE << "Loading all detectors for " << analysis_job->dataset_directory << DIR_END_CHAR << dataset_file << "\n";
            delete spectra_volume;
            delete tmp_spectra_volume;
            if (status_callback != nullptr)
            {
                (*status_callback)(0, 1);
            }
            return;
        }
        //add all detectors up
        for (size_t j = 0; j < spectra_volume->rows(); j++)
        {
            for (size_t k = 0; k < spectra_volume->cols(); k++)
            {
                real_t elapsed_livetime = (*spectra_volume)[j][k].elapsed_livetime();
                real_t elapsed_realtime = (*spectra_volume)[j][k].elapsed_realtime();
                real_t input_counts = (*spectra_volume)[j][k].input_counts();
                real_t output_counts = (*spectra_volume)[j][k].output_counts();


                (*spectra_volume)[j][k] += (*tmp_spectra_volume)[j][k];

                elapsed_livetime += (*tmp_spectra_volume)[j][k].elapsed_livetime();
                elapsed_realtime += (*tmp_spectra_volume)[j][k].elapsed_realtime();
                input_counts += (*tmp_spectra_volume)[j][k].input_counts();
                output_counts += (*tmp_spectra_volume)[j][k].output_counts();

                (*spectra_volume)[j][k].elapsed_livetime(elapsed_livetime);
                (*spectra_volume)[j][k].elapsed_realtime(elapsed_realtime);
                (*spectra_volume)[j][k].input_counts(input_counts);
                (*spectra_volume)[j][k].output_counts(output_counts);
            }
        }
    }
    delete tmp_spectra_volume;

    analysis_job->init_fit_routines(spectra_volume->samples_size(), true);
	
    proc_spectra(spectra_volume, detector, &tp, !is_loaded_from_analyzed_h5, status_callback);
    delete spectra_volume;
}

// ----------------------------------------------------------------------------

void find_quantifier_scalers(data_struct::Params_Override * override_params, unordered_map<string, string> &pv_map, Quantification_Standard* quantification_standard)
{
    std::string quant_scalers_names[] = {STR_US_IC, STR_DS_IC, "SRCURRENT"};
    real_t *pointer_arr[] = {&(quantification_standard->US_IC),&(quantification_standard->DS_IC), &(quantification_standard->sr_current)};
	real_t scaler_clock = 1.0;
	if (override_params->time_scaler_clock.length() > 0)
	{
		scaler_clock = std::stof(override_params->time_scaler_clock);
	}
    int i =0;
    for(auto &itr : quant_scalers_names)
    {

        Summed_Scaler* sscaler = nullptr;
        for(auto & ssItr : override_params->summed_scalers)
        {
            if (ssItr.scaler_name == itr)
            {
                sscaler = &(ssItr);
                break;
            }
        }
        if(sscaler != nullptr)
        {
            *(pointer_arr[i]) = 0.0;
            for (auto &jitr : sscaler->scalers_to_sum)
            {
                if(override_params->time_normalized_scalers.count(jitr)
               && pv_map.count(override_params->time_normalized_scalers[jitr])
               && pv_map.count(override_params->time_scaler))
                {
                    real_t val = std::stof(pv_map[override_params->time_normalized_scalers[jitr]]);
                    real_t det_time = std::stof(pv_map[override_params->time_scaler]);
                    det_time /= scaler_clock;
                    val /= det_time;
                    *(pointer_arr[i]) += val;
                }
                else if(override_params->scaler_pvs.count(jitr) && pv_map.count(override_params->scaler_pvs[jitr]) > 0)
                {
                    *(pointer_arr[i]) += std::stof(pv_map[override_params->scaler_pvs[jitr]]);
                }
            }
        }
        else if(override_params->time_normalized_scalers.count(itr) && pv_map.count(override_params->time_normalized_scalers.at(itr)))
        {
            *(pointer_arr[i]) = std::stof(pv_map[override_params->time_normalized_scalers[itr]]);
        }
        else if(override_params->scaler_pvs.count(itr) && pv_map.count(override_params->scaler_pvs.at(itr)))
        {
            *(pointer_arr[i]) = std::stof(pv_map[override_params->scaler_pvs[itr]]);
        }
        i++;
    }
}

// ----------------------------------------------------------------------------

void load_and_fit_quatification_datasets(data_struct::Analysis_Job* analysis_job, size_t detector_num)
{
    fitting::models::Gaussian_Model model;
    quantification::models::Quantification_Model quantification_model;

    data_struct::Detector* detector = analysis_job->get_detector(detector_num);
    data_struct::Params_Override* override_params = &(detector->fit_params_override_dict);

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    for (Quantification_Standard& standard_itr : analysis_job->standard_element_weights)
    {
        // detecotr_struct descructor will delete this memory
        detector->quantification_standards[standard_itr.standard_filename] = Quantification_Standard(standard_itr.standard_filename, standard_itr.element_standard_weights);
        Quantification_Standard* quantification_standard = &(detector->quantification_standards[standard_itr.standard_filename]);

        //Output of fits for elements specified
        std::unordered_map<std::string, data_struct::Fit_Element_Map*> elements_to_fit;
        for (auto& itr : standard_itr.element_standard_weights)
        {
            data_struct::Element_Info* e_info = data_struct::Element_Info_Map::inst()->get_element(itr.first);
            elements_to_fit[itr.first] = new data_struct::Fit_Element_Map(itr.first, e_info);
            elements_to_fit[itr.first]->init_energy_ratio_for_detector_element(detector->detector_element, standard_itr.disable_Ka_for_quantification, standard_itr.disable_La_for_quantification);
        }

        unordered_map<string, string> pv_map;
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
                        if (false == io::load_and_integrate_spectra_volume(analysis_job->dataset_directory, quantification_standard->standard_filename, detector_num, &quantification_standard->integrated_spectra, override_params))
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
                    find_quantifier_scalers(override_params, pv_map, quantification_standard);
                }
            }
            else
            {
                find_quantifier_scalers(override_params, pv_map, quantification_standard);
            }
        }
        else
        {
            if (false == io::load_and_integrate_spectra_volume(analysis_job->dataset_directory, quantification_standard->standard_filename, detector_num, &quantification_standard->integrated_spectra, override_params))
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

            fitting::routines::Base_Fit_Routine* fit_routine = fit_itr.second;
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
                Fit_Parameters fit_params = model.fit_parameters();
                
                //add elements to fit parameters if they don't exist
                for (auto& itr2 : elements_to_fit)
                {
                    if (false == override_params->fit_params.contains(itr2.first))
                    {
                        fit_params.add_parameter(Fit_Param(itr2.first, quantification_standard->element_counts.at(fit_itr.first).at(itr2.first)));
                    }
                }
                fitting::routines::Matrix_Optimized_Fit_Routine* f_routine = (fitting::routines::Matrix_Optimized_Fit_Routine*)fit_routine;
                real_t energy_offset = fit_params.value(STR_ENERGY_OFFSET);
                real_t energy_slope = fit_params.value(STR_ENERGY_SLOPE);
                real_t energy_quad = fit_params.value(STR_ENERGY_QUADRATIC);

                data_struct::ArrayXr energy = data_struct::ArrayXr::LinSpaced(energy_range.count(), energy_range.min, energy_range.max);
                data_struct::ArrayXr ev = energy_offset + (energy * energy_slope) + (Eigen::pow(energy, (real_t)2.0) * energy_quad);
                data_struct::ArrayXr sub_spectra = quantification_standard->integrated_spectra.segment(energy_range.min, energy_range.count());

                std::string full_path = analysis_job->dataset_directory + DIR_END_CHAR + "output" + DIR_END_CHAR + "calib_" + fit_routine->get_name() + "_" + standard_itr.standard_filename;
                if (detector_num != -1)
                {
                    full_path += std::to_string(detector_num);
                }
                logI << full_path << "\n";
                #ifdef _BUILD_WITH_QT
                visual::SavePlotSpectrasFromConsole(full_path + ".png", &ev, &sub_spectra, (ArrayXr*)(&f_routine->fitted_integrated_spectra()), (ArrayXr*)(&f_routine->fitted_integrated_background()), true);
                #endif

                io::file::csv::save_fit_and_int_spectra(full_path + ".csv", &ev, &sub_spectra, (ArrayXr*)(&f_routine->fitted_integrated_spectra()), (ArrayXr*)(&f_routine->fitted_integrated_background()));
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
}

// ----------------------------------------------------------------------------

bool perform_quantification(data_struct::Analysis_Job* analysis_job)
{
    quantification::models::Quantification_Model quantification_model;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    logI << "Perform_quantification()"<<"\n";

    vector<string> quant_scaler_name_list = { STR_SR_CURRENT, STR_US_IC, STR_DS_IC };

    if( io::load_quantification_standardinfo(analysis_job->dataset_directory, analysis_job->quantification_standard_filename, analysis_job->standard_element_weights) )
    {
        for(size_t detector_num : analysis_job->detector_num_arr)
        {
            data_struct::Detector* detector = analysis_job->get_detector(detector_num);
            data_struct::Params_Override* override_params = &(detector->fit_params_override_dict);

            
            load_and_fit_quatification_datasets(analysis_job, detector_num);
            detector->generage_avg_quantification_scalers();

            for(auto &fit_itr : detector->fit_routines)
            {
               fitting::optimizers::Optimizer* optimizer = analysis_job->optimizer();
               for (auto& quant_itr : detector->avg_quantification_scaler_map)
               {
                    fitting::routines::Base_Fit_Routine *fit_routine = fit_itr.second;

                    // update e_cal_ratio for elements in standards by average value
                    for (auto& s_itr : detector->quantification_standards)
                    {
                        Quantification_Standard* quantification_standard = &(s_itr.second);
                        detector->update_element_quants(fit_itr.first, quant_itr.first, quantification_standard, &quantification_model, quant_itr.second);
                    }

                    logI << Fitting_Routine_To_Str.at(fit_itr.first) << " "<< quant_itr.first  << "\n";
                    Fit_Parameters fit_params;
                    fit_params.add_parameter(Fit_Param("quantifier", 0.0, std::numeric_limits<real_t>::max(), 1.0, 0.1, E_Bound_Type::FIT));
                    //initial guess: parinfo_value[0] = 100000.0 / factor
                    fit_params["quantifier"].value = (real_t)100000.0 / quant_itr.second;
                    optimizer->minimize_quantification(&fit_params, &detector->all_element_quants[fit_itr.first][quant_itr.first], &quantification_model);
                    real_t val = fit_params["quantifier"].value;

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

        io::save_quantification_plots(analysis_job->dataset_directory, detector);
        }
    }
    else
    {
        logE<<"Loading quantification standard "<<analysis_job->quantification_standard_filename<<"\n";
        return false;
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logI << "quantification elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;

}

// ----------------------------------------------------------------------------

void interate_datasets_and_update(data_struct::Analysis_Job& analysis_job)
{
    for (const auto& dataset_file : analysis_job.dataset_files)
    {
        //average all detectors to one files
        if (analysis_job.generate_average_h5)
        {
            io::generate_h5_averages(analysis_job.dataset_directory, dataset_file, analysis_job.detector_num_arr);
        }

        //generate a list of dataset to update
        std::vector<std::string> hdf5_dataset_list;

        for (size_t detector_num : analysis_job.detector_num_arr)
        {
            if (detector_num != -1)
            {
                hdf5_dataset_list.push_back(analysis_job.dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file + ".h5" + std::to_string(detector_num));
            }
        }
		size_t dlen = dataset_file.length();
        if (dlen > 6 && dataset_file[dlen - 5] == '.' && dataset_file[dlen - 4] == 'h' && dataset_file[dlen - 3] == 'd' && dataset_file[dlen - 2] == 'f' && dataset_file[dlen - 1] == '5')
        {
            hdf5_dataset_list.push_back(analysis_job.dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file);
        }
		else if(dlen > 4 && dataset_file[dlen - 3] == '.' && dataset_file[dlen - 2] == 'h' && dataset_file[dlen - 1] == '5')
		{
			hdf5_dataset_list.push_back(analysis_job.dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file);
		}
		else
		{
			hdf5_dataset_list.push_back(analysis_job.dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file + ".h5");
		}

        

        // Can scan for all hdf5 files in img.dat but will have to do it for multiple ext and filter out by detector num and avg
        //for (auto& itr : io::find_all_dataset_files(analysis_job.dataset_directory + "img.dat", ".hdf5"))
        //{
        //    hdf5_dataset_list.push_back(analysis_job.dataset_directory + "img.dat" + DIR_END_CHAR + itr);
        //}

        for (std::string hdf5_dataset_name : hdf5_dataset_list)
        {
            //export csv
            if (analysis_job.export_int_fitted_to_csv)
            {
                io::file::HDF5_IO::inst()->export_int_fitted_to_csv(hdf5_dataset_name);
            }

            //update theta based on new PV
            if (analysis_job.update_theta_str.length() > 0)
            {
                //data_struct::Params_Override* params_override
                io::file::HDF5_IO::inst()->update_theta(hdf5_dataset_name, analysis_job.update_theta_str);
            }

            //update scalers table in hdf5
            if (analysis_job.update_scalers)
            {
                data_struct::Detector* det = nullptr;
                size_t len = hdf5_dataset_name.length();

                if (len > 4 && hdf5_dataset_name[len - 3] == '.' && hdf5_dataset_name[len - 2] == 'h' && hdf5_dataset_name[len - 1] == '5')
                {
                    det = analysis_job.get_detector(-1);
                }
                else if (len > 4 && hdf5_dataset_name[len - 4] == '.' && hdf5_dataset_name[len - 3] == 'h' && hdf5_dataset_name[len - 2] == '5' && hdf5_dataset_name[len - 1] == '0')
                {
                    det = analysis_job.get_detector(0);
                }
                else if (len > 4 && hdf5_dataset_name[len - 4] == '.' && hdf5_dataset_name[len - 3] == 'h' && hdf5_dataset_name[len - 2] == '5' && hdf5_dataset_name[len - 1] == '1')
                {
                    det = analysis_job.get_detector(1);
                }
                else if (len > 4 && hdf5_dataset_name[len - 4] == '.' && hdf5_dataset_name[len - 3] == 'h' && hdf5_dataset_name[len - 2] == '5' && hdf5_dataset_name[len - 1] == '2')
                {
                    det = analysis_job.get_detector(2);
                }
                else if (len > 4 && hdf5_dataset_name[len - 4] == '.' && hdf5_dataset_name[len - 3] == 'h' && hdf5_dataset_name[len - 2] == '5' && hdf5_dataset_name[len - 1] == '3')
                {
                    det = analysis_job.get_detector(3);
                }

                if (det != nullptr)
                {
                    io::file::HDF5_IO::inst()->update_scalers(hdf5_dataset_name, &det->fit_params_override_dict);
                }
            }


            //add v9 layout soft links
            if (analysis_job.add_v9_layout)
            {
                io::file::HDF5_IO::inst()->add_v9_layout(hdf5_dataset_name);
            }

            //add exchange
            if (analysis_job.add_exchange_layout)
            {
                io::file::HDF5_IO::inst()->add_exchange_layout(hdf5_dataset_name);
            }
        }
    }
}

// ----------------------------------------------------------------------------