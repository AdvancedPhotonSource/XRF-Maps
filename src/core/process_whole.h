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


#ifndef PROCESS_WHOLE
#define PROCESS_WHOLE

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

#include <stdlib.h>

#include "defines.h"

#if defined _WIN32
#include "support/direct/dirent.h"
#else
#include <dirent.h>
#endif

#include "core/defines.h"

#include "workflow/threadpool.h"

#include "io/file/hl_file_io.h"
#include "io/file/mca_io.h"

#include "data_struct/spectra_volume.h"

#include "fitting/models/gaussian_model.h"

#include "data_struct/element_info.h"

#include "io/file/aps/aps_fit_params_import.h"

#include "fitting/routines/roi_fit_routine.h"
#include "fitting/routines/svd_fit_routine.h"
#include "fitting/routines/nnls_fit_routine.h"
#include "fitting/routines/hybrid_param_nnls_fit_routine.h"

#include "fitting/optimizers/optimizer.h"

#include "data_struct/fit_element_map.h"
#include "data_struct/params_override.h"

#include "data_struct/quantification_standard.h"

#include "io/file/aps/aps_roi.h"


using namespace std::placeholders; //for _1, _2,


// ----------------------------------------------------------------------------

DLL_EXPORT bool perform_quantification(data_struct::Analysis_Job<double>* analysis_job, bool save_when_done);

DLL_EXPORT bool optimize_integrated_fit_params(data_struct::Analysis_Job<double>* analysis_job,
                                                std::string  dataset_filename,
                                                size_t detector_num,
                                                data_struct::Params_Override<double>* params_override,
                                                data_struct::Fit_Parameters<double>& out_fitp,
                                                Callback_Func_Status_Def* status_callback = nullptr);

DLL_EXPORT void generate_optimal_params(data_struct::Analysis_Job<double>* analysis_job);

void load_and_fit_quatification_datasets(data_struct::Analysis_Job<double>* analysis_job, size_t detector_num, std::vector<Quantification_Standard<double>>& standard_element_weights, std::unordered_map<size_t, double>& quant_map);

DLL_EXPORT void optimize_single_roi(data_struct::Analysis_Job<double>& analysis_job, std::string roi_file_name, std::map<int, std::map<std::string, data_struct::Fit_Parameters<double>>>& out_roi_fit_params, Callback_Func_Status_Def* status_callback = nullptr);

DLL_EXPORT void optimize_rois(data_struct::Analysis_Job<double>& analysis_job);

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT data_struct::Fit_Count_Dict<T_real>* generate_fit_count_dict(const Fit_Element_Map_Dict<T_real> *elements_to_fit, size_t height, size_t width, bool alloc_iter_count)
{
    data_struct::Fit_Count_Dict<T_real>* element_fit_counts_dict = new data_struct::Fit_Count_Dict<T_real>();
    for (auto& e_itr : *elements_to_fit)
    {
        element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr<T_real> >(e_itr.first, data_struct::ArrayXXr<T_real>()));
        element_fit_counts_dict->at(e_itr.first).resize(height, width);
    }

    if (alloc_iter_count)
    {
        //Allocate memeory to save number of fit iterations
        element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr<T_real> >(STR_NUM_ITR, data_struct::ArrayXXr<T_real>()));
        element_fit_counts_dict->at(STR_NUM_ITR).resize(height, width);
        element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr<T_real> >(STR_RESIDUAL, data_struct::ArrayXXr<T_real>()));
        element_fit_counts_dict->at(STR_RESIDUAL).resize(height, width);
    }

    //  TOTAL_FLUORESCENCE_YIELD
    element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr<T_real> >(STR_TOTAL_FLUORESCENCE_YIELD, data_struct::ArrayXXr<T_real>()));
    element_fit_counts_dict->at(STR_TOTAL_FLUORESCENCE_YIELD).resize(height, width);

    //SUM_ELASTIC_INELASTIC
    element_fit_counts_dict->emplace(std::pair<std::string, data_struct::ArrayXXr<T_real> >(STR_SUM_ELASTIC_INELASTIC_AMP, data_struct::ArrayXXr<T_real>()));
    element_fit_counts_dict->at(STR_SUM_ELASTIC_INELASTIC_AMP).resize(height, width);


    return element_fit_counts_dict;
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT bool fit_single_spectra(fitting::routines::Base_Fit_Routine<T_real>* fit_routine,
                        const fitting::models::Base_Model<T_real>* const model,
                        const data_struct::Spectra<T_real>* const spectra,
                        const data_struct::Fit_Element_Map_Dict<T_real>* const elements_to_fit,
                        data_struct::Fit_Count_Dict<T_real>* out_fit_counts,
                        size_t i,
                        size_t j)
{
    std::unordered_map<std::string, T_real> counts_dict;
    fit_routine->fit_spectra(model, spectra, elements_to_fit, counts_dict);
    //save count / sec
    for (auto& el_itr : *elements_to_fit)
    {
        (*out_fit_counts)[el_itr.first](i, j) = counts_dict[el_itr.first] / spectra->elapsed_livetime();
    }
    (*out_fit_counts)[STR_NUM_ITR](i, j) = counts_dict[STR_NUM_ITR];

    (*out_fit_counts)[STR_RESIDUAL](i, j) = counts_dict[STR_RESIDUAL];
    // add total fluorescense yield
    if (out_fit_counts->count(STR_TOTAL_FLUORESCENCE_YIELD))
    {
        (*out_fit_counts)[STR_TOTAL_FLUORESCENCE_YIELD](i, j) = spectra->sum();
    }
    // add sum coherent and compton
    if (out_fit_counts->count(STR_SUM_ELASTIC_INELASTIC_AMP) > 0 && counts_dict.count(STR_COHERENT_SCT_AMPLITUDE) > 0 && counts_dict.count(STR_COMPTON_AMPLITUDE) > 0)
    {
        (*out_fit_counts)[STR_SUM_ELASTIC_INELASTIC_AMP](i, j) = counts_dict[STR_COHERENT_SCT_AMPLITUDE] + counts_dict[STR_COMPTON_AMPLITUDE];
        // add total fluorescense yield
        if (out_fit_counts->count(STR_TOTAL_FLUORESCENCE_YIELD))
        {                   //                                      (sum - (elastic + inelastic)) / live time
            (*out_fit_counts)[STR_TOTAL_FLUORESCENCE_YIELD](i, j) = (spectra->sum() - (*out_fit_counts)[STR_SUM_ELASTIC_INELASTIC_AMP](i, j)) / spectra->elapsed_livetime();
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

template<typename T_real>
DLL_EXPORT void proc_spectra(data_struct::Spectra_Volume<T_real>* spectra_volume,
                             data_struct::Detector<T_real>* detector,
                             ThreadPool* tp,
                             bool save_spec_vol,
                             Callback_Func_Status_Def* status_callback = nullptr)
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

    data_struct::Params_Override<T_real>* override_params = &(detector->fit_params_override_dict);

    //Range of energy in spectra to fit
    fitting::models::Range energy_range = data_struct::get_energy_range(spectra_volume->samples_size(), &(detector->fit_params_override_dict.fit_params));

    std::chrono::time_point<std::chrono::system_clock> start, end;

    for (auto& itr : detector->fit_routines)
    {
        fitting::routines::Base_Fit_Routine<T_real>* fit_routine = itr.second;

        logI << "Processing  " << fit_routine->get_name() << "\n";

        start = std::chrono::system_clock::now();
        if (override_params->elements_to_fit.size() < 1)
        {
            logE << "No elements to fit. Check  maps_fit_parameters_override.txt0 - 3 exist" << "\n";
            continue;
        }

        //Fit job queue
        std::queue<std::future<bool> >* fit_job_queue = new std::queue<std::future<bool> >();

        //Allocate memeory to save fit counts
        data_struct::Fit_Count_Dict<T_real>* element_fit_count_dict = generate_fit_count_dict(&override_params->elements_to_fit, spectra_volume->rows(), spectra_volume->cols(), true);

        for (size_t i = 0; i < spectra_volume->rows(); i++)
        {
            for (size_t j = 0; j < spectra_volume->cols(); j++)
            {
                //logD<< i<<" "<<j<<"\n";
                fit_job_queue->emplace(tp->enqueue(fit_single_spectra<T_real>, fit_routine, detector->model, &(*spectra_volume)[i][j], &override_params->elements_to_fit, element_fit_count_dict, i, j));
            }
        }

        size_t total_blocks = (spectra_volume->rows() * spectra_volume->cols()) - 1;
        size_t cur_block = 0;
        //wait for queue to finish processing
        while (!fit_job_queue->empty())
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
        std::chrono::duration<double> elapsed_seconds = end - start;
        logI << "Fitting [ " << fit_routine->get_name() << " ] elapsed time: " << elapsed_seconds.count() << "s" << "\n";

        io::file::HDF5_IO::inst()->save_element_fits(fit_routine->get_name(), element_fit_count_dict);
        io::file::HDF5_IO::inst()->save_params_override(override_params);

        if (itr.first == data_struct::Fitting_Routines::GAUSS_MATRIX
            || itr.first == data_struct::Fitting_Routines::NNLS
            || itr.first == data_struct::Fitting_Routines::SVD)
        {
            fitting::routines::Matrix_Optimized_Fit_Routine<T_real>* matrix_fit = (fitting::routines::Matrix_Optimized_Fit_Routine<T_real>*)fit_routine;
            if(spectra_volume->samples_size() > 0)
            {
                io::file::HDF5_IO::inst()->save_fitted_int_spectra(fit_routine->get_name(),
                    matrix_fit->fitted_integrated_spectra(),
                    matrix_fit->energy_range(),
                    matrix_fit->fitted_integrated_background(),
                    (*spectra_volume)[0][0].size());
            }
            // save png 
            std::string dataset_fullpath = io::file::HDF5_IO::inst()->get_filename();
            int sidx = dataset_fullpath.find("img.dat");
            if (dataset_fullpath.length() > 0 && sidx > 0) 
            {
                dataset_fullpath.replace(sidx, 7, "output"); // 7 = sizeof("img.dat")
                std::string str_path = dataset_fullpath + "_" + fit_routine->get_name() + ".png";
                data_struct::ArrayTr<T_real> ev = data_struct::generate_energy_array(matrix_fit->energy_range(), &(override_params->fit_params));
                Spectra<T_real> int_spec = spectra_volume->integrate();
                if (matrix_fit->energy_range().count() < int_spec.size())
                {
                    int_spec = int_spec.sub_spectra(matrix_fit->energy_range().min, matrix_fit->energy_range().count());
                    #ifdef _BUILD_WITH_QT
                    visual::SavePlotSpectrasFromConsole(str_path, &ev, &int_spec, (&matrix_fit->fitted_integrated_spectra()), (&matrix_fit->fitted_integrated_background()), true);
                    #endif
                }
            }

        }
        if (itr.first == data_struct::Fitting_Routines::GAUSS_MATRIX)
        {
            fitting::routines::Matrix_Optimized_Fit_Routine<T_real>* matrix_fit = (fitting::routines::Matrix_Optimized_Fit_Routine<T_real>*)fit_routine;
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

    T_real energy_offset = 0.0;
    T_real energy_slope = 0.0;
    T_real energy_quad = 0.0;
    data_struct::Fit_Parameters<T_real> fit_params = detector->model->fit_parameters();
    if (fit_params.contains(STR_ENERGY_OFFSET))
    {
        energy_offset = fit_params[STR_ENERGY_OFFSET].value;
    }
    if (fit_params.contains(STR_ENERGY_SLOPE))
    {
        energy_slope = fit_params[STR_ENERGY_SLOPE].value;
    }
    if (fit_params.contains(STR_ENERGY_QUADRATIC))
    {
        energy_quad = fit_params[STR_ENERGY_QUADRATIC].value;
    }

    io::file::HDF5_IO::inst()->save_energy_calib(spectra_volume->samples_size(), energy_offset, energy_slope, energy_quad);

    if (save_spec_vol)
    {
        io::file::HDF5_IO::inst()->save_spectra_volume("mca_arr", spectra_volume);
    }
    
    io::file::HDF5_IO::inst()->end_save_seq();


}


// ----------------------------------------------------------------------------

template<typename T_real>
void find_quantifier_scalers(std::unordered_map<std::string, double>& pv_map, Quantification_Standard<T_real>* quantification_standard)
{

    // find time scaler
    std::string time_pv = "";
    std::string beamline = "";
    double time_clock = 0.0;
    T_real time_val = 1.0;
    if (data_struct::Scaler_Lookup::inst()->search_for_timing_info(pv_map, time_pv, time_clock, beamline))
    {
        time_val = pv_map.at(time_pv);
        time_val /= time_clock;
    }

    // update pv names to labels
    for (auto& itr : pv_map)
    {
        std::string label = "";
        bool is_time_normalized = false;
        if (data_struct::Scaler_Lookup::inst()->search_pv(itr.first, label, is_time_normalized, beamline))
        {
            if (is_time_normalized)
            {
                itr.second /= time_val;

            }
            pv_map[label] = itr.second;
        }
    }

    // add any summded scalers to pv_map
    const std::vector<struct Summed_Scaler>* summed_scalers = data_struct::Scaler_Lookup::inst()->get_summed_scaler_list(beamline);
    if (summed_scalers != nullptr)
    {
        for (const auto& itr : *summed_scalers)
        {
            T_real summed_val = 0.0;
            for (const auto& sitr : itr.scalers_to_sum)
            {
                if (pv_map.count(sitr) > 0)
                {
                    summed_val += pv_map.at(sitr);
                }
            }
            pv_map[itr.scaler_name] = summed_val;
        }
    }
    // search for sr_current, us_ic, and ds_ic
    for (auto& itr : pv_map)
    {
        if (itr.first == STR_SR_CURRENT)
        {
            quantification_standard->sr_current = itr.second;
        }
        else if (itr.first == STR_US_IC)
        {
            quantification_standard->US_IC = itr.second;
        }
        else if (itr.first == STR_US_FM)
        {
            quantification_standard->US_FM = itr.second;
        }
        else if (itr.first == STR_DS_IC)
        {
            quantification_standard->DS_IC = itr.second;
        }
    }
}

// ----------------------------------------------------------------------------

template<typename T_real>
DLL_EXPORT void process_dataset_files(data_struct::Analysis_Job<T_real>* analysis_job, Callback_Func_Status_Def* status_callback = nullptr)
{
    ThreadPool tp(analysis_job->num_threads);

    for (auto& dataset_file : analysis_job->dataset_files)
    {
        //if quick and dirty then sum all detectors to 1 spectra volume and process it
        if (analysis_job->quick_and_dirty)
        {
            process_dataset_files_quick_and_dirty(dataset_file, analysis_job, tp);
        }
        //otherwise process each detector separately
        else
        {
            for (size_t detector_num : analysis_job->detector_num_arr)
            {

                data_struct::Detector<T_real>* detector = analysis_job->get_detector(detector_num);

                //Spectra volume data
                data_struct::Spectra_Volume<T_real>* spectra_volume = new data_struct::Spectra_Volume<T_real>();

                std::string fullpath;
                size_t dlen = dataset_file.length();
                bool is_mda = (dataset_file[dlen - 4] == '.' && dataset_file[dlen - 3] == 'm' && dataset_file[dlen - 2] == 'd' && dataset_file[dlen - 1] == 'a');
                bool is_mca = (dataset_file[dlen - 4] == '.' && dataset_file[dlen - 3] == 'm' && dataset_file[dlen - 2] == 'c' && dataset_file[dlen - 1] == 'a');
                bool is_mcad = (dataset_file[dlen - 5] == '.' && dataset_file[dlen - 4] == 'm' && dataset_file[dlen - 3] == 'c' && dataset_file[dlen - 2] == 'a');
                //bool is_h5 = (dataset_file[dlen - 3] == '.' && dataset_file[dlen - 2] == 'h' && dataset_file[dlen - 1] == '5');
                //bool is_hdf5 = (dataset_file[dlen - 5] == '.' && dataset_file[dlen - 4] == 'h' && dataset_file[dlen - 3] == 'd' && dataset_file[dlen - 2] == 'f' && dataset_file[dlen - 1] == '5');
                if (is_mda || is_mca || is_mcad)
                {
                    std::string str_detector_num = "";
                    if (detector_num != -1)
                    {
                        str_detector_num = std::to_string(detector_num);
                    }
                    std::string full_save_path = analysis_job->output_dir + "img.dat" + DIR_END_CHAR + dataset_file + ".h5" + str_detector_num;
                    io::file::HDF5_IO::inst()->set_filename(full_save_path);
                }
                else if(detector_num == (size_t)-1)
                {
                    std::string full_save_path = analysis_job->output_dir + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + dataset_file + ".h5";
                    io::file::HDF5_IO::inst()->set_filename(full_save_path);
                }
                else
                {
                    std::string full_save_path = analysis_job->output_dir + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + dataset_file + ".h5" + std::to_string(detector_num);
                    io::file::HDF5_IO::inst()->set_filename(full_save_path);
                }

                bool loaded_from_analyzed_hdf5 = false;
                //load spectra volume
                if (false == io::file::load_spectra_volume(analysis_job->dataset_directory, dataset_file, detector_num, spectra_volume, &detector->fit_params_override_dict, &loaded_from_analyzed_hdf5, true))
                {
                    logW << "Skipping detector " << detector_num << "\n";
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

template<typename T_real>
DLL_EXPORT void process_dataset_files_quick_and_dirty(std::string dataset_file, data_struct::Analysis_Job<T_real>* analysis_job, ThreadPool& tp, Callback_Func_Status_Def* status_callback = nullptr)
{
    std::string full_save_path = analysis_job->output_dir + DIR_END_CHAR + "img.dat" + DIR_END_CHAR + dataset_file + ".h5";

    data_struct::Detector<T_real>* detector = analysis_job->get_detector(0);
    //Spectra volume data
    data_struct::Spectra_Volume<T_real>* spectra_volume = new data_struct::Spectra_Volume<T_real>();
    data_struct::Spectra_Volume<T_real>* tmp_spectra_volume = new data_struct::Spectra_Volume<T_real>();

    io::file::HDF5_IO::inst()->start_save_seq(full_save_path, true); // force to create new file for quick and dirty

    //load the first one
    size_t detector_num = analysis_job->detector_num_arr[0];
    bool is_loaded_from_analyzed_h5 = false;
    if (false == io::file::load_spectra_volume(analysis_job->dataset_directory, dataset_file, detector_num, spectra_volume, &detector->fit_params_override_dict, &is_loaded_from_analyzed_h5, true))
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
        if (false == io::file::load_spectra_volume(analysis_job->dataset_directory, dataset_file, analysis_job->detector_num_arr[i], tmp_spectra_volume, &detector->fit_params_override_dict, &is_loaded_from_analyzed_h5, false))
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
                T_real elapsed_livetime = (*spectra_volume)[j][k].elapsed_livetime();
                T_real elapsed_realtime = (*spectra_volume)[j][k].elapsed_realtime();
                T_real input_counts = (*spectra_volume)[j][k].input_counts();
                T_real output_counts = (*spectra_volume)[j][k].output_counts();


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

template<typename T_real>
DLL_EXPORT void iterate_datasets_and_update(data_struct::Analysis_Job<T_real>& analysis_job)
{
    for (auto dataset_file : analysis_job.dataset_files)
    {
        //average all detectors to one files
        if (analysis_job.generate_average_h5)
        {
            io::file::generate_h5_averages(analysis_job.output_dir, dataset_file, analysis_job.detector_num_arr);
        }

        if (analysis_job.add_background)
        {
            data_struct::Detector<T_real>* detector = analysis_job.get_detector(0);

            if (detector != nullptr)
            {
                io::file::HDF5_IO::inst()->add_background(analysis_job.output_dir, dataset_file, detector->fit_params_override_dict);
            }
            else
            {
                logW << "Detector == nullptr for add_background\n";
            }
        }

        //generate a list of dataset to update
        std::vector<std::string> hdf5_dataset_list;

        for (size_t detector_num : analysis_job.detector_num_arr)
        {
            if (detector_num != -1)
            {
                hdf5_dataset_list.push_back(analysis_job.output_dir + "img.dat" + DIR_END_CHAR + dataset_file + ".h5" + std::to_string(detector_num));
            }
        }
        size_t dlen = dataset_file.length();
        if (dlen > 6 && dataset_file[dlen - 5] == '.' && dataset_file[dlen - 4] == 'h' && dataset_file[dlen - 3] == 'd' && dataset_file[dlen - 2] == 'f' && dataset_file[dlen - 1] == '5')
        {
            hdf5_dataset_list.push_back(analysis_job.output_dir + "img.dat" + DIR_END_CHAR + dataset_file);
        }
        else if (dlen > 4 && dataset_file[dlen - 3] == '.' && dataset_file[dlen - 2] == 'h' && dataset_file[dlen - 1] == '5')
        {
            hdf5_dataset_list.push_back(analysis_job.output_dir + "img.dat" + DIR_END_CHAR + dataset_file);
        }
        else
        {
            hdf5_dataset_list.push_back(analysis_job.output_dir + "img.dat" + DIR_END_CHAR + dataset_file + ".h5");
        }


        // Can scan for all hdf5 files in img.dat but will have to do it for multiple ext and filter out by detector num and avg
        //for (auto& itr : io::file::File_Scan::inst()->find_all_dataset_files(analysis_job.output_dir + "img.dat", ".hdf5"))
        //{
        //    hdf5_dataset_list.push_back(analysis_job.output_dir + "img.dat" + DIR_END_CHAR + itr);
        //}

        for (std::string hdf5_dataset_name : hdf5_dataset_list)
        {
            //export csv
            if (analysis_job.export_int_fitted_to_csv)
            {
                io::file::HDF5_IO::inst()->export_int_fitted_to_csv<T_real>(hdf5_dataset_name);
            }

            //update theta based on new PV
            if (analysis_job.update_theta_str.length() > 0)
            {
                io::file::HDF5_IO::inst()->update_theta(hdf5_dataset_name, analysis_job.update_theta_str);
            }

            //update upstream and downstream amps 
            if (analysis_job.update_us_amps_str.length() > 0 && analysis_job.update_ds_amps_str.length() > 0)
            {
                io::file::HDF5_IO::inst()->update_amps(hdf5_dataset_name, analysis_job.update_us_amps_str, analysis_job.update_ds_amps_str);
            }

            //update quantification upstream and downstream amps
            if (analysis_job.update_quant_ds_amps_str.length() > 0 && analysis_job.update_quant_ds_amps_str.length() > 0)
            {
                io::file::HDF5_IO::inst()->update_quant_amps(hdf5_dataset_name, analysis_job.update_quant_us_amps_str, analysis_job.update_quant_ds_amps_str);
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

#endif
