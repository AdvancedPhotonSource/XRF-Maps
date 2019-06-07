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
        element_fit_counts_dict->emplace(std::pair<std::string, Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(e_itr.first, Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> ()) );
        element_fit_counts_dict->at(e_itr.first).resize(height, width);
    }

    if (alloc_iter_count)
    {
        //Allocate memeory to save number of fit iterations
        element_fit_counts_dict->emplace(std::pair<std::string, Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(STR_NUM_ITR, Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>() ));
        element_fit_counts_dict->at(STR_NUM_ITR).resize(height, width);
    }

	//  TOTAL_FLUORESCENCE_YIELD
	element_fit_counts_dict->emplace(std::pair<std::string, Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(STR_TOTAL_FLUORESCENCE_YIELD, Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>()));
	element_fit_counts_dict->at(STR_TOTAL_FLUORESCENCE_YIELD).resize(height, width);

	//SUM_ELASTIC_INELASTIC
	element_fit_counts_dict->emplace(std::pair<std::string, Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> >(STR_SUM_ELASTIC_INELASTIC_AMP , Eigen::Matrix<real_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>()));
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
    std::unordered_map<std::string, real_t> counts_dict = fit_routine->fit_spectra(model, spectra, elements_to_fit);
    //save count / sec
    for (auto& el_itr : *elements_to_fit)
    {
        (*out_fit_counts)[el_itr.first](i,j) = counts_dict[el_itr.first] / spectra->elapsed_livetime();
    }
    (*out_fit_counts)[STR_NUM_ITR](i,j) = counts_dict[STR_NUM_ITR];
	// add total fluorescense yield
	if (out_fit_counts->count(STR_TOTAL_FLUORESCENCE_YIELD))
	{
		(*out_fit_counts)[STR_TOTAL_FLUORESCENCE_YIELD](i, j) = spectra->sum();
	}
	// add sum coherent and compton
	if (out_fit_counts->count(STR_SUM_ELASTIC_INELASTIC_AMP) > 0  && counts_dict.count(STR_COHERENT_SCT_AMPLITUDE) > 0 && counts_dict.count(STR_COMPTON_AMPLITUDE) > 0)
	{
		(*out_fit_counts)[STR_SUM_ELASTIC_INELASTIC_AMP](i, j) = counts_dict[STR_COHERENT_SCT_AMPLITUDE] + counts_dict[STR_COMPTON_AMPLITUDE];
	}

    return true;
}

// ----------------------------------------------------------------------------

// struct io::file_name_fit_params* optimize_integrated_fit_params_preloaded(std::string dataset_directory,
//                                                                            std::string  dataset_filename,
//                                                                            size_t detector_num,
//                                                                            data_struct::Spectra * integrated_spectra,
//                                                                            data_struct::Params_Override *params_override,
//                                                                            fitting::models::Fit_Params_Preset optimize_fit_params_preset,
//                                                                            fitting::optimizers::Optimizer *optimizer)
//{

//    //return structure
//    struct io::file_name_fit_params* ret_struct = new struct::io::file_name_fit_params();

//    ret_struct->dataset_dir = dataset_directory;
//    ret_struct->dataset_filename = dataset_filename;
//    ret_struct->detector_num = detector_num;

//    fitting::models::Gaussian_Model model;

//    ret_struct->elements_to_fit = params_override->elements_to_fit;

//    ret_struct->spectra = *integrated_spectra;

//    //Range of energy in spectra to fit
//    fitting::models::Range energy_range = data_struct::get_energy_range(ret_struct->spectra.size(), &(params_override->fit_params));

//    //Fitting routines
//    fitting::routines::Param_Optimized_Fit_Routine fit_routine;
//    fit_routine.set_optimizer(optimizer);

//    //reset model fit parameters to defaults
//    model.reset_to_default_fit_params();
//    //set fixed/fit preset
//    model.set_fit_params_preset(optimize_fit_params_preset);
//    //Update fit parameters by override values
//    model.update_fit_params_values(&(params_override->fit_params));


//    //Initialize the fit routine
//    fit_routine.initialize(&model, &ret_struct->elements_to_fit, energy_range);
//    //Fit the spectra saving the element counts in element_fit_count_dict
//    ret_struct->fit_params = fit_routine.fit_spectra_parameters(&model, &ret_struct->spectra, &ret_struct->elements_to_fit);

//    ret_struct->success = true;

//    delete params_override;
//    delete integrated_spectra;
//    return ret_struct;

//}

// ----------------------------------------------------------------------------

 struct io::file_name_fit_params* optimize_integrated_fit_params(std::string dataset_directory,
                                                                std::string  dataset_filename,
                                                                size_t detector_num,
                                                                fitting::models::Fit_Params_Preset optimize_fit_params_preset,
                                                                fitting::optimizers::Optimizer *optimizer)
{

    //return structure
    struct io::file_name_fit_params* ret_struct = new struct::io::file_name_fit_params();

    ret_struct->dataset_dir = dataset_directory;
    ret_struct->dataset_filename = dataset_filename;
    ret_struct->detector_num = detector_num;

    fitting::models::Gaussian_Model model;

    data_struct::Params_Override params_override;

    //load override parameters
    if(false == io::load_override_params(dataset_directory, detector_num, &params_override) )
    {
        if(false == io::load_override_params(dataset_directory, -1, &params_override))
        {
            logE<<"Loading maps_fit_parameters_override.txt"<<dataset_filename<<" for detector"<<detector_num<<"\n";
            ret_struct->success = false;
            return ret_struct;
        }
    }

    ret_struct->elements_to_fit = params_override.elements_to_fit;

    //load the quantification standard dataset
    if(false == io::load_and_integrate_spectra_volume(dataset_directory, dataset_filename, detector_num, &ret_struct->spectra, &params_override) )
    {
        logE<<"In optimize_integrated_dataset loading dataset"<<dataset_filename<<" for detector"<<detector_num<<"\n";
        ret_struct->success = false;
        return ret_struct;
    }

    //Range of energy in spectra to fit
    fitting::models::Range energy_range = data_struct::get_energy_range(ret_struct->spectra.size(), &(params_override.fit_params));

    //Fitting routines
    fitting::routines::Param_Optimized_Fit_Routine fit_routine;
    fit_routine.set_optimizer(optimizer);

    //reset model fit parameters to defaults
    model.reset_to_default_fit_params();
    //set fixed/fit preset
    model.set_fit_params_preset(optimize_fit_params_preset);
    //Update fit parameters by override values
    model.update_fit_params_values(&(params_override.fit_params));

    
    //Initialize the fit routine
    fit_routine.initialize(&model, &ret_struct->elements_to_fit, energy_range);
    //Fit the spectra saving the element counts in element_fit_count_dict
    ret_struct->fit_params = fit_routine.fit_spectra_parameters(&model, &ret_struct->spectra, &ret_struct->elements_to_fit);

    ret_struct->success = true;

    return ret_struct;

}

// ----------------------------------------------------------------------------

void generate_optimal_params(data_struct::Analysis_Job* analysis_job)
{
    bool first = true;

    std::unordered_map<int, data_struct::Fit_Parameters> fit_params_avgs;
    int file_cnt = 0;
    for(auto &itr : analysis_job->optimize_dataset_files)
    {
        file_cnt += 1.0;
        for(size_t detector_num = analysis_job->detector_num_start; detector_num <= analysis_job->detector_num_end; detector_num++)
        {
            struct io::file_name_fit_params* f_struct = optimize_integrated_fit_params(analysis_job->dataset_directory, itr, detector_num, analysis_job->optimize_fit_params_preset, analysis_job->optimizer());
            if(f_struct->success)
            {
                if(first)
                {
                    fit_params_avgs[f_struct->detector_num] = f_struct->fit_params;
                    first = false;
                }
                else
                {
                    fit_params_avgs[f_struct->detector_num].sum_values(f_struct->fit_params);
                }
                io::save_optimized_fit_params(f_struct);
            }
            delete f_struct;
        }
    }
    for(size_t detector_num = analysis_job->detector_num_start; detector_num <= analysis_job->detector_num_end; detector_num++)
    {
        fit_params_avgs[detector_num].divide_fit_values_by(file_cnt);
    }

    io::save_averaged_fit_params(analysis_job->dataset_directory, fit_params_avgs, analysis_job->detector_num_start, analysis_job->detector_num_end);

}

// ----------------------------------------------------------------------------

//void generate_optimal_params_mp(data_struct::Analysis_Job* analysis_job)
//{
//    bool first = true;
//    std::queue<std::future<struct io::file_name_fit_params*> > job_queue;
//    ThreadPool tp(analysis_job->num_threads);

//    std::unordered_map<int, data_struct::Fit_Parameters> fit_params_avgs;
//    real_t file_cnt = 0.;
//    for(auto &itr : analysis_job->optimize_dataset_files)
//    {
//        file_cnt += 1.0;
//        for(size_t detector_num = analysis_job->detector_num_start; detector_num <= analysis_job->detector_num_end; detector_num++)
//        {
//            data_struct::Params_Override *params_override = new data_struct::Params_Override();
//            data_struct::Spectra * int_spectra = new Spectra();
//            //load override parameters
//            if(false == io::load_override_params(analysis_job->dataset_directory, detector_num, params_override) )
//            {
//                if(false == io::load_override_params(analysis_job->dataset_directory, -1, params_override))
//                {
//                    logE<<"loading maps_fit_parameters_override.txt"<<itr<<" for detector"<<detector_num<<"\n";
//                    continue;
//                }
//            }

//            //load the quantification standard dataset
//            if(false == io::load_and_integrate_spectra_volume(analysis_job->dataset_directory, itr, detector_num, int_spectra, params_override) )
//            {
//                logE<<"loading dataset"<<itr<<" for detector"<<detector_num<<"\n";
//                continue;
//            }

//            job_queue.emplace( tp.enqueue(optimize_integrated_fit_params_preloaded, analysis_job->dataset_directory, itr, detector_num, int_spectra, params_override, analysis_job->optimize_fit_params_preset, analysis_job->optimizer()) );
//        }
//    }

//    while(!job_queue.empty())
//    {
//        auto ret = std::move(job_queue.front());
//        job_queue.pop();
//        struct io::file_name_fit_params* f_struct = ret.get();
//        if(f_struct->success)
//        {
//            if(first)
//            {
//                fit_params_avgs[f_struct->detector_num] = f_struct->fit_params;
//                first = false;
//            }
//            else
//            {
//                fit_params_avgs[f_struct->detector_num].sum_values(f_struct->fit_params);
//            }
//            io::save_optimized_fit_params(f_struct);
//        }
//        delete f_struct;
//    }
//    for(size_t detector_num = analysis_job->detector_num_start; detector_num <= analysis_job->detector_num_end; detector_num++)
//    {
//        fit_params_avgs[detector_num].divide_fit_values_by(file_cnt);
//    }

//    io::save_averaged_fit_params(analysis_job->dataset_directory, fit_params_avgs, analysis_job->detector_num_start, analysis_job->detector_num_end);

//}


// ----------------------------------------------------------------------------

void proc_spectra(data_struct::Spectra_Volume* spectra_volume,
                  data_struct::Detector * detector_struct,
                  ThreadPool* tp,
                  bool save_spec_vol)
{
    data_struct::Params_Override * override_params = &(detector_struct->fit_params_override_dict);

    //Range of energy in spectra to fit
    fitting::models::Range energy_range = data_struct::get_energy_range(spectra_volume->samples_size(), &(detector_struct->fit_params_override_dict.fit_params) );

    std::chrono::time_point<std::chrono::system_clock> start, end;

    for(auto &itr : detector_struct->fit_routines)
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
                fit_job_queue->emplace( tp->enqueue(fit_single_spectra, fit_routine, detector_struct->model, &(*spectra_volume)[i][j], &override_params->elements_to_fit, element_fit_count_dict, i, j) );
            }
        }

        io::save_results( fit_routine->get_name(), element_fit_count_dict, fit_job_queue, start );
    }

    real_t energy_offset = 0.0;
    real_t energy_slope = 0.0;
    real_t energy_quad = 0.0;
    data_struct::Fit_Parameters fit_params = detector_struct->model->fit_parameters();
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
        io::file::HDF5_IO::inst()->save_quantifications(detector_struct->quant_standards);
        io::save_volume(spectra_volume, energy_offset, energy_slope, energy_quad);
        io::file::HDF5_IO::inst()->end_save_seq();
    }
    else
    {
        io::file::HDF5_IO::inst()->save_quantifications(detector_struct->quant_standards);
        io::file::HDF5_IO::inst()->end_save_seq();
    }

}

// ----------------------------------------------------------------------------

void process_dataset_files(data_struct::Analysis_Job* analysis_job)
{
    ThreadPool tp(analysis_job->num_threads);

    for(auto &dataset_file : analysis_job->dataset_files)
    {
        //if quick and dirty then sum all detectors to 1 spectra volume and process it
        if(analysis_job->quick_and_dirty)
        {
            std::string full_save_path = analysis_job->dataset_directory+ DIR_END_CHAR+"img.dat"+ DIR_END_CHAR +dataset_file+".h5";

            data_struct::Detector* detector_struct = analysis_job->get_detector(0);
            //Spectra volume data
            data_struct::Spectra_Volume* spectra_volume = new data_struct::Spectra_Volume();
            data_struct::Spectra_Volume* tmp_spectra_volume = new data_struct::Spectra_Volume();

            io::file::HDF5_IO::inst()->set_filename(full_save_path);

            //load the first one
            size_t detector_num = analysis_job->detector_num_start;
            bool is_loaded_from_analyzed_h5;
            if (false == io::load_spectra_volume(analysis_job->dataset_directory, dataset_file, detector_num, spectra_volume, &detector_struct->fit_params_override_dict, &is_loaded_from_analyzed_h5, true) )
            {
                logE<<"Loading all detectors for "<<analysis_job->dataset_directory<< DIR_END_CHAR <<dataset_file<<"\n";
                delete spectra_volume;
                delete tmp_spectra_volume;
                return;
            }

            //load spectra volume
            for(detector_num = analysis_job->detector_num_start+1; detector_num <= analysis_job->detector_num_end; detector_num++)
            {

                if (false == io::load_spectra_volume(analysis_job->dataset_directory, dataset_file, detector_num, tmp_spectra_volume, &detector_struct->fit_params_override_dict, &is_loaded_from_analyzed_h5, false) )
                {
                    logE<<"Loading all detectors for "<<analysis_job->dataset_directory<< DIR_END_CHAR <<dataset_file<<"\n";
                    delete spectra_volume;
                    delete tmp_spectra_volume;
                    return;
                }
                //add all detectors up
                for(size_t j=0; j<spectra_volume->rows(); j++)
                {
                    for(size_t k=0; k<spectra_volume->cols(); k++)
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

            proc_spectra(spectra_volume, detector_struct, &tp, !is_loaded_from_analyzed_h5);
			delete spectra_volume;
        }
        //otherwise process each detector separately
        else
        {
            for(size_t detector_num = analysis_job->detector_num_start; detector_num <= analysis_job->detector_num_end; detector_num++)
            {

                data_struct::Detector* detector_struct = analysis_job->get_detector(detector_num);

                //Spectra volume data
                data_struct::Spectra_Volume* spectra_volume = new data_struct::Spectra_Volume();

                std::string str_detector_num = std::to_string(detector_num);
                std::string full_save_path = analysis_job->dataset_directory+ DIR_END_CHAR+"img.dat"+ DIR_END_CHAR +dataset_file+".h5"+str_detector_num;
                io::file::HDF5_IO::inst()->set_filename(full_save_path);

                bool loaded_from_analyzed_hdf5 = false;
                //load spectra volume
                if (false == io::load_spectra_volume(analysis_job->dataset_directory, dataset_file, detector_num, spectra_volume, &detector_struct->fit_params_override_dict, &loaded_from_analyzed_hdf5, true) )
                {
                    logW<<"Skipping detector "<<detector_num<<"\n";
                    delete spectra_volume;
                    continue;
                }

                analysis_job->init_fit_routines(spectra_volume->samples_size(), true);
                proc_spectra(spectra_volume, detector_struct, &tp, !loaded_from_analyzed_hdf5);
				delete spectra_volume;
            }
        }
    }
}

// ----------------------------------------------------------------------------

void find_quantifier_scalers(data_struct::Params_Override * override_params, unordered_map<string, string> &pv_map, Quantification_Standard* quantification_standard)
{
    std::string quant_scalers_names[] = {"US_IC", "DS_IC", "SRCURRENT"};
    real_t *pointer_arr[] = {&(quantification_standard->US_IC),&(quantification_standard->DS_IC), &(quantification_standard->sr_current)};
    real_t scaler_clock = std::stof(override_params->time_scaler_clock);
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
                if(override_params->time_normalized_scalers.count(jitr.first)
               && pv_map.count(override_params->time_normalized_scalers[jitr.first])
               && pv_map.count(override_params->time_scaler))
                {
                    real_t val = std::stof(pv_map[override_params->time_normalized_scalers[jitr.first]]);
                    real_t det_time = std::stof(pv_map[override_params->time_scaler]);
                    det_time /= scaler_clock;
                    val /= det_time;
                    *(pointer_arr[i]) += val;
                }
                else if(override_params->scaler_pvs.count(jitr.first) && pv_map.count(override_params->scaler_pvs[jitr.first]) > 0)
                {
                    *(pointer_arr[i]) += std::stof(pv_map[override_params->scaler_pvs[jitr.first]]);
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

bool perform_quantification(data_struct::Analysis_Job* analysis_job)
{
    io::file::MCA_IO mca_io;
    fitting::models::Gaussian_Model model;
    quantification::models::Quantification_Model quantification_model;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    logI << "Perform_quantification()"<<"\n";

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    vector<io::element_weights_struct> standard_element_weights;

    if( io::load_quantification_standardinfo(analysis_job->dataset_directory, analysis_job->quantification_standard_filename, standard_element_weights) )
    {
        for(size_t detector_num = analysis_job->detector_num_start; detector_num <= analysis_job->detector_num_end; detector_num++)
        {
            //                    Proc_type   counts
            //std::unordered_map<std::string, real_t> all_element_counts;
            //  proc_type   quantifier        element     quant_prop
            //map<int, map<int, unordered_map<string, Element_Quant>>> all_element_quants;
            real_t avg_sr_current = 0.0;
            real_t avg_US_IC = 0.0;
            real_t avg_DS_IC = 0.0;
            real_t avg_cnt = 0.0;

            data_struct::Detector* detector_struct = analysis_job->get_detector(detector_num);
            data_struct::Params_Override * override_params = &(detector_struct->fit_params_override_dict);

            for(io::element_weights_struct &standard_itr : standard_element_weights)
            {
                Quantification_Standard* quantification_standard = new Quantification_Standard();
				// detecotr_struct descructor will delete this memory
                detector_struct->quant_standards[standard_itr.standard_file_name] = quantification_standard;

                //Parameters for calibration curve
                if (override_params->detector_element.length() > 0)
                {
                    // Get the element info class                                           // detector element as string "Si" or "Ge" usually
                    quantification_standard->detector_element = (data_struct::Element_Info_Map::inst()->get_element(override_params->detector_element));
                }
                if (override_params->be_window_thickness.length() > 0)
                {
                    quantification_standard->beryllium_window_thickness = (std::stof(override_params->be_window_thickness));
                }
                if (override_params->ge_dead_layer.length() > 0)
                {
                    quantification_standard->germanium_dead_layer = (std::stof(override_params->ge_dead_layer));
                }
                if (override_params->det_chip_thickness.length() > 0)
                {
                    quantification_standard->detector_chip_thickness = (std::stof(override_params->det_chip_thickness));
                }
                if(override_params->fit_params.contains(STR_COHERENT_SCT_ENERGY))
                {
                    quantification_standard->incident_energy = (override_params->fit_params.at(STR_COHERENT_SCT_ENERGY).value);
                }

                quantification_standard->standard_filename = standard_itr.standard_file_name;

                //Output of fits for elements specified
                std::unordered_map<std::string, data_struct::Fit_Element_Map*> elements_to_fit;
                for(auto& itr : standard_itr.element_standard_weights)
                {
                    data_struct::Element_Info* e_info = data_struct::Element_Info_Map::inst()->get_element(itr.first);
                    elements_to_fit[itr.first] = new data_struct::Fit_Element_Map(itr.first, e_info);
                    elements_to_fit[itr.first]->init_energy_ratio_for_detector_element( quantification_standard->detector_element );
                }

                unordered_map<string, string> pv_map;
                //load the quantification standard dataset
                size_t fn_str_len = quantification_standard->standard_filename.length();
                if(fn_str_len > 5 &&
                   quantification_standard->standard_filename[fn_str_len - 4] == '.' &&
                   quantification_standard->standard_filename[fn_str_len - 3] == 'm' &&
                   quantification_standard->standard_filename[fn_str_len - 2] == 'c' &&
                   quantification_standard->standard_filename[fn_str_len - 1] == 'a')
                {
                    //try with adding detector_num on the end for 2ide datasets
                    if(false == mca_io.load_integrated_spectra(analysis_job->dataset_directory + quantification_standard->standard_filename + std::to_string(detector_num), &quantification_standard->integrated_spectra, pv_map) )
                    {
                        //try without detector number on end 2idd
                        if(false == mca_io.load_integrated_spectra(analysis_job->dataset_directory + quantification_standard->standard_filename, &quantification_standard->integrated_spectra, pv_map) )
                        {

                            //legacy code would load mca files, check for mca and replace with mda
                            size_t std_str_len = standard_itr.standard_file_name.length();
                            if(standard_itr.standard_file_name[std_str_len - 4] == '.' && standard_itr.standard_file_name[std_str_len - 3] == 'm' && standard_itr.standard_file_name[std_str_len - 2] == 'c' && standard_itr.standard_file_name[std_str_len - 1] == 'a')
                            {
                                standard_itr.standard_file_name[std_str_len - 2] = 'd';
                                quantification_standard->standard_filename = standard_itr.standard_file_name;
                                if(false == io::load_and_integrate_spectra_volume(analysis_job->dataset_directory, quantification_standard->standard_filename, detector_num, &quantification_standard->integrated_spectra, override_params) )
                                {
                                    logE<<"Could not load file "<< standard_itr.standard_file_name <<" for detector"<<detector_num<<"\n";
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
                                logE<<"Could not load file "<< standard_itr.standard_file_name <<" for detector"<<detector_num<<"\n";
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
                    if(false == io::load_and_integrate_spectra_volume(analysis_job->dataset_directory, quantification_standard->standard_filename, detector_num, &quantification_standard->integrated_spectra, override_params) )
                    {
                        logE<<"Could not load file "<< standard_itr.standard_file_name <<" for detector"<<detector_num<<"\n";
                        continue;
                    }
                    else
                    {
                        quantification_standard->sr_current = override_params->sr_current;
                        quantification_standard->US_IC = override_params->US_IC;
                        quantification_standard->DS_IC = override_params->DS_IC;
                    }
                }
                if(quantification_standard->integrated_spectra.size() == 0)
                {
                    logE<<"Spectra size == 0! Can't process it!\n";
                    continue;
                }

                avg_sr_current += quantification_standard->sr_current;
                avg_US_IC += quantification_standard->US_IC;
                avg_DS_IC += quantification_standard->DS_IC;
                avg_cnt += 1.0;

                analysis_job->init_fit_routines(quantification_standard->integrated_spectra.size(), true);

                //First we integrate the spectra and get the elemental counts
                energy_range.max = quantification_standard->integrated_spectra.size() -1;


                for(auto &itr : detector_struct->fit_routines)
                {
					
                    fitting::routines::Base_Fit_Routine *fit_routine = itr.second;
					for (auto& el_itr : standard_itr.element_standard_weights)
					{
                    //    element_quants[el_itr.first] = Element_Quant(el_itr.second);
                        quantification_standard->append_element(el_itr.first, el_itr.second);
					}

                    //reset model fit parameters to defaults
                    model.reset_to_default_fit_params();
                    //Update fit parameters by override values
                    model.update_fit_params_values(&(override_params->fit_params));
                    //Initialize the fit routine
                    fit_routine->initialize(&model, &elements_to_fit, energy_range);
                    //Fit the spectra
                    std::unordered_map<std::string, real_t>counts_dict = fit_routine->fit_spectra(&model,
                                                                                                  &quantification_standard->integrated_spectra,
                                                                                                  &elements_to_fit);

                    for (auto& itr2 : elements_to_fit)
                    {
                        counts_dict[itr2.first] /= quantification_standard->integrated_spectra.elapsed_livetime();
                    }

                    unordered_map<size_t, real_t*> quant_list =
                    {
                        {Quantifiers::CURRENT, &quantification_standard->sr_current},
                        {Quantifiers::US_IC, &quantification_standard->US_IC},
                        {Quantifiers::DS_IC, &quantification_standard->DS_IC}
                    };

                    for (auto& quant_itr : quant_list)
                    {
                        quantification_standard->init_element_quants(fit_routine->get_name(), &counts_dict, &quantification_model, quant_itr.first, *quant_itr.second);
                        for(auto& q_itr: quantification_standard->element_quants)
                        {
                            detector_struct->all_element_quants[itr.first][quant_itr.first][q_itr.first] = q_itr.second;
                        }
                    }

                }                 

                //cleanup
                for(auto &itr3 : elements_to_fit)
                {
                 delete itr3.second;
                }
				elements_to_fit.clear();
            }

            avg_sr_current /= avg_cnt;
            avg_US_IC /= avg_cnt;
            avg_DS_IC /= avg_cnt;

            unordered_map<size_t, real_t*> quant_list =
            {
                {Quantifiers::CURRENT, &avg_sr_current},
                {Quantifiers::US_IC, &avg_US_IC},
                {Quantifiers::DS_IC, &avg_DS_IC}
            };

            for(auto &itr : detector_struct->fit_routines)
            {
               fitting::optimizers::Optimizer* optimizer = analysis_job->optimizer();
               for (auto& quant_itr : quant_list)
               {
                    fitting::routines::Base_Fit_Routine *fit_routine = itr.second;

                    Fit_Parameters fit_params;
                    fit_params.add_parameter(Fit_Param("quantifier", 0.0, 0.0, 1.0, 0.001, E_Bound_Type::FIT));
                    //initial guess: parinfo_value[0] = 100000.0 / factor
                    fit_params["quantifier"].value = (real_t)100000.0 / (*quant_itr.second);
                    optimizer->minimize_quantification(&fit_params, &detector_struct->all_element_quants[itr.first][quant_itr.first], &quantification_model);
                    real_t val = fit_params["quantifier"].value;

                    if(std::isinf(val))
                    {
                        logW<<"Quantifier Value = Inf. setting it to 0.\n";
                        val = 0;
                    }
                    else
                    {
                        logI<<"Quantifier Value = "<<val<<"\n";
                    }
                    for(auto &dq_itr: detector_struct->quant_standards)
                    {
                        dq_itr.second->generate_calibration_curve(fit_routine->get_name(), quant_itr.first, val);
                    }
                }
            }

        io::save_quantification_plots(analysis_job, &detector_struct->quant_standards, detector_num);
        }
    }
    else
    {
        logE<<"Loading quantification standard "<<analysis_job->quantification_standard_filename<<"\n";
        return false;
    }

    if(analysis_job->quick_and_dirty)
    {
        ////average_quantification(quant_stand_list, analysis_job->detector_num_start(), analysis_job->detector_num_end());
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logI << "quantification elapsed time: " << elapsed_seconds.count() << "s"<<"\n";

    return true;

}

// ----------------------------------------------------------------------------

void average_quantification(std::vector<data_struct::Quantification_Standard>* quant_stand_list,
                            size_t detector_num_start,
                            size_t detector_num_end)
{
    /*
    data_struct::Quantification_Standard q_standard_0;
    q_standard_0.


    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        data_struct::Quantification_Standard * quantification_standard = &(*quant_stand_list)[detector_num];
        quantification_standard->
    }
    */
}

// ----------------------------------------------------------------------------
