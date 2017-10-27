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

#include "core/process_streaming.h"

// ----------------------------------------------------------------------------

data_struct::xrf::Stream_Block* proc_spectra_block( data_struct::xrf::Stream_Block* stream_block )
{

    for(auto itr : stream_block->fitting_blocks)
    {
        int i = itr.first;
        std::unordered_map<std::string, real_t> counts_dict = stream_block->fitting_blocks[i].fit_routine->fit_spectra(stream_block->model, stream_block->spectra, stream_block->elements_to_fit);
        //make count / sec
        for (auto& el_itr : *(stream_block->elements_to_fit))
        {
            stream_block->fitting_blocks[i].fit_counts[el_itr.first] = counts_dict[el_itr.first] / stream_block->spectra->elapsed_lifetime();
        }
        stream_block->fitting_blocks[i].fit_counts[data_struct::xrf::STR_NUM_ITR] = counts_dict[data_struct::xrf::STR_NUM_ITR];
    }
    return stream_block;
}

// ----------------------------------------------------------------------------

void run_stream_pipeline(data_struct::xrf::Analysis_Job* job)
{
    workflow::xrf::Spectra_File_Source spectra_stream_producer(job);
    workflow::Distributor<data_struct::xrf::Stream_Block*, data_struct::xrf::Stream_Block*> distributor(job->num_threads());
    workflow::xrf::Spectra_Stream_Saver sink;

    distributor.set_function(proc_spectra_block);
    spectra_stream_producer.connect(&distributor);
    sink.connect(&distributor);

    sink.start();
    spectra_stream_producer.run();
    sink.wait_and_stop();
}

// ----------------------------------------------------------------------------

struct io::file_name_fit_params* optimize_integrated_fit_params( data_struct::xrf::Stream_Block* stream_block )
{
    struct io::file_name_fit_params* ret_struct = new struct io::file_name_fit_params();

    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = ret_struct->spectra.size() -1;

    fitting::models::Gaussian_Model model;
    //Fitting routines
    fitting::routines::Param_Optimized_Fit_Routine fit_routine;
    //TODO:
////    fit_routine.set_optimizer(optimizer);

    //reset model fit parameters to defaults
    model.reset_to_default_fit_params();
    //Update fit parameters by override values
    //TODO:
////    model.update_fit_params_values(params_override.fit_params);
    model.set_fit_params_preset(stream_block->optimize_fit_params_preset);
    //Initialize the fit routine
    fit_routine.initialize(&model, &ret_struct->elements_to_fit, energy_range);
    //Fit the spectra saving the element counts in element_fit_count_dict
    ret_struct->fit_params = fit_routine.fit_spectra_parameters(stream_block->model, stream_block->spectra, stream_block->elements_to_fit);

    ret_struct->success = true;

    delete stream_block;

    return ret_struct;
}

// ----------------------------------------------------------------------------

void save_optimal_params(struct io::file_name_fit_params* f_struct)
{
    static bool first = false;
    if(f_struct->success)
    {
        if(first)
        {
            //TODO:
///            fit_params_avgs[f_struct->detector_num] = f_struct->fit_params;
            first = false;
        }
        else
        {
            //TODO:
///            fit_params_avgs[f_struct->detector_num].moving_average_with(f_struct->fit_params);
        }
        io::save_optimized_fit_params(*f_struct);
    }
    delete f_struct;
}

// ----------------------------------------------------------------------------

void run_optimization_stream_pipeline(data_struct::xrf::Analysis_Job* job)
{
    //TODO: Run only on 8 largest files if no files are specified
    workflow::xrf::Integrated_Spectra_Source spectra_stream_producer(job);
    workflow::Distributor<data_struct::xrf::Stream_Block*, struct io::file_name_fit_params*> distributor(job->num_threads());
    workflow::Sink<struct io::file_name_fit_params*> sink;
    sink.set_function(save_optimal_params);

    distributor.set_function(optimize_integrated_fit_params);
    sink.connect(&distributor);
    spectra_stream_producer.connect(&distributor);


    sink.start();
    spectra_stream_producer.run();
    sink.wait_and_stop();
    /*
    io::save_averaged_fit_params(dataset_directory, fit_params_avgs, detector_num_start, detector_num_end);
    */
}

// ----------------------------------------------------------------------------

void run_quick_n_dirty_pipeline(data_struct::xrf::Analysis_Job* job)
{
    workflow::xrf::Detector_Sum_Spectra_Source sum_detectors_spectra_stream_producer(job);
    workflow::Distributor<data_struct::xrf::Stream_Block*, data_struct::xrf::Stream_Block*> distributor(job->num_threads());
    workflow::xrf::Spectra_Stream_Saver sink;


    distributor.set_function(proc_spectra_block);
    sink.connect(&distributor);
    sum_detectors_spectra_stream_producer.connect(&distributor);

    sink.start();
    sum_detectors_spectra_stream_producer.run();
    sink.wait_and_stop();
}

// ----------------------------------------------------------------------------

bool perform_quantification_streaming(std::string dataset_directory,
                            std::string quantification_info_file,
                            std::vector<data_struct::xrf::Fitting_Routines> proc_types,
                            std::vector<data_struct::xrf::Quantification_Standard>* quant_stand_list,
                            std::unordered_map<int, data_struct::xrf::Params_Override> *fit_params_override_dict,
                            size_t detector_num_start,
                            size_t detector_num_end)
{
/*
    bool air_path = false;
    real_t detector_chip_thickness = 0.0;
    real_t beryllium_window_thickness = 0.0;
    real_t germanium_dead_layer = 0.0;
    real_t incident_energy = 10.0;

    fitting::models::Gaussian_Model model;

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    logit << "Perform_quantification()"<<std::endl;

    data_struct::xrf::Element_Info* detector_element = data_struct::xrf::Element_Info_Map::inst()->get_element("Si");

    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;

    std::string standard_file_name;
    std::unordered_map<std::string, real_t> element_standard_weights;
    //data_struct::xrf::Quantification_Standard* quantification_standard = &(*quant_stand_list)[0];

    if( io::load_quantification_standard(dataset_directory, quantification_info_file, &standard_file_name, &element_standard_weights) )
    {
        for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
        {
            data_struct::xrf::Params_Override * override_params = nullptr;
            if(fit_params_override_dict->count(detector_num) > 0)
            {
               override_params = &fit_params_override_dict->at(detector_num);
            }
            if(override_params == nullptr && fit_params_override_dict->count(-1) > 0)
            {
               override_params = &fit_params_override_dict->at(-1);
            }
            if(override_params == nullptr)
            {
                logit<<"Skipping file "<<dataset_directory << " detector "<<detector_num<<" because could not find maps_fit_params_override.txt"<<std::endl;
                continue;
            }


            data_struct::xrf::Quantification_Standard* quantification_standard = &(*quant_stand_list)[detector_num];
            quantification_standard->standard_filename(standard_file_name);
            for(auto& itr : element_standard_weights)
            {
                quantification_standard->append_element(itr.first, itr.second);
            }

            //Parameters for calibration curve

            if (override_params->detector_element.length() > 0)
            {
                // Get the element info class                                           // detector element as string "Si" or "Ge" usually
                detector_element = data_struct::xrf::Element_Info_Map::inst()->get_element(override_params->detector_element);
            }
            if (override_params->be_window_thickness.length() > 0)
            {
                beryllium_window_thickness = std::stof(override_params->be_window_thickness);
            }
            if (override_params->ge_dead_layer.length() > 0)
            {
                germanium_dead_layer = std::stof(override_params->ge_dead_layer);
            }
            if (override_params->det_chip_thickness.length() > 0)
            {
                detector_chip_thickness = std::stof(override_params->det_chip_thickness);
            }

            if(override_params->fit_params.contains(data_struct::xrf::STR_COHERENT_SCT_ENERGY))
            {
                incident_energy = override_params->fit_params.at(data_struct::xrf::STR_COHERENT_SCT_ENERGY).value;
            }

            //Output of fits for elements specified
            std::unordered_map<std::string, data_struct::xrf::Fit_Element_Map*> elements_to_fit;
            for(auto& itr : element_standard_weights)
            {
                data_struct::xrf::Element_Info* e_info = data_struct::xrf::Element_Info_Map::inst()->get_element(itr.first);
                elements_to_fit[itr.first] = new data_struct::xrf::Fit_Element_Map(itr.first, e_info);
                elements_to_fit[itr.first]->init_energy_ratio_for_detector_element( detector_element );
            }

            data_struct::xrf::Spectra_Volume spectra_volume;
            //load the quantification standard dataset
            if(false == io::load_spectra_volume(dataset_directory, quantification_standard->standard_filename(), &spectra_volume, detector_num, override_params, quantification_standard, false) )
            {
                //legacy code would load mca files, check for mca and replace with mda
                int std_str_len = standard_file_name.length();
                if(standard_file_name[std_str_len - 4] == '.' && standard_file_name[std_str_len - 3] == 'm' && standard_file_name[std_str_len - 2] == 'c' && standard_file_name[std_str_len - 1] == 'a')
                {
                    standard_file_name[std_str_len - 2] = 'd';
                    quantification_standard->standard_filename(standard_file_name);
                    if(false == io::load_spectra_volume(dataset_directory, quantification_standard->standard_filename(), &spectra_volume, detector_num, override_params, quantification_standard, false) )
                    {
                        logit<<"Error perform_quantification() : could not load file "<< standard_file_name <<" for detector"<<detector_num<<std::endl;
                        return false;
                    }
                }
                else
                {
                    logit<<"Error perform_quantification() : could not load file "<< standard_file_name <<" for detector"<<detector_num<<std::endl;
                    return false;
                }
            }

            //First we integrate the spectra and get the elemental counts
            data_struct::xrf::Spectra integrated_spectra = spectra_volume.integrate();
            energy_range.max = integrated_spectra.size() -1;


            for(auto proc_type : proc_types)
            {

                //Fitting routines
                fitting::routines::Base_Fit_Routine *fit_routine = generate_fit_routine(proc_type);

                //reset model fit parameters to defaults
                model.reset_to_default_fit_params();
                //Update fit parameters by override values
                model.update_fit_params_values(override_params->fit_params);
                //Initialize the fit routine
                fit_routine->initialize(&model, &elements_to_fit, energy_range);
                //Fit the spectra
                std::unordered_map<std::string, real_t>counts_dict = fit_routine->fit_spectra(&model,
                                                                                              &integrated_spectra,
                                                                                              &elements_to_fit);

                if (fit_routine != nullptr)
                {
                    delete fit_routine;
                    fit_routine = nullptr;
                }

                for (auto& itr : elements_to_fit)
                {
                    counts_dict[itr.first] /= integrated_spectra.elapsed_lifetime();
                }
                quantification_standard->integrated_spectra(integrated_spectra);

                //save for each proc
                quantification_standard->quantifiy(optimizer,
                                                   fit_routine->get_name(),
                                                  &counts_dict,
                                                  incident_energy,
                                                  detector_element,
                                                  air_path,
                                                  detector_chip_thickness,
                                                  beryllium_window_thickness,
                                                  germanium_dead_layer);

            }

            //cleanup
            for(auto &itr : elements_to_fit)
            {
                delete itr.second;
            }

        }
    }
    else
    {
        logit<<"Error loading quantification standard "<<quantification_info_file<<std::endl;
        return false;
    }

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logit << "quantification elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;
*/
    return true;

}

// ----------------------------------------------------------------------------
