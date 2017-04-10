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

#include "hl_file_io.h"

namespace io
{

std::vector<std::string> netcdf_files;
std::vector<std::string> hdf_files;

// ----------------------------------------------------------------------------

bool compare_file_size (const file_name_size& first, const file_name_size& second)
{
    return ( first.total_rank_size > second.total_rank_size );
}

// ----------------------------------------------------------------------------

void populate_netcdf_hdf5_files(std::string dataset_dir)
{
    //populate netcdf and hdf5 files for fly scans
    netcdf_files = find_all_dataset_files(dataset_dir + "flyXRF/", "_0.nc");
    hdf_files = find_all_dataset_files(dataset_dir + "flyXRF.h5/", "_0.h5");
}

// ----------------------------------------------------------------------------

bool load_element_info(std::string element_henke_filename, std::string element_csv_filename, data_struct::xrf::Element_Info_Map *element_info_map)
{

    if (io::file::load_element_info_from_csv(element_csv_filename, element_info_map) == false)
    {
        logit<<"error loading "<< element_csv_filename<<std::endl;
        return false;
    }

    if (io::file::load_henke_from_xdr(element_henke_filename, element_info_map) == false)
    {
        logit<<"error loading "<< element_henke_filename<<std::endl;
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

bool save_results(std::string save_loc,
                  const data_struct::xrf::Fit_Count_Dict * const element_counts,
                  fitting::routines::Base_Fit_Routine* fit_routine,
                  std::queue<std::future<bool> >* job_queue,
                  std::chrono::time_point<std::chrono::system_clock> start)
{


    //wait for queue to finish processing
    while(!job_queue->empty())
    {
        auto ret = std::move(job_queue->front());
        job_queue->pop();
        ret.get();
    }

    std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    logit << "Fitting [ "<< save_loc <<" ] elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;


    io::file::HDF5_IO::inst()->save_element_fits(save_loc, element_counts);

    delete fit_routine;
    delete job_queue;
    delete element_counts;

    return true;
}

// ----------------------------------------------------------------------------

bool save_volume(data_struct::xrf::Quantification_Standard * quantification_standard,
                 data_struct::xrf::Spectra_Volume *spectra_volume,
                 real_t energy_offset,
                 real_t energy_slope,
                 real_t energy_quad)
                 //std::queue<std::future<bool> >* job_queue)
{
	/*
    //wait for queue to finish processing
    while(!job_queue->empty())
    {
        auto ret = std::move(job_queue->front());
        job_queue->pop();
        ret.get();
    }*/

    io::file::HDF5_IO::inst()->save_quantification(quantification_standard);
    io::file::HDF5_IO::inst()->save_spectra_volume("mca_arr", spectra_volume, energy_offset, energy_slope, energy_quad);
    io::file::HDF5_IO::inst()->end_save_seq();

    //delete job_queue;
    delete spectra_volume;

    return true;
}

// ----------------------------------------------------------------------------

void save_optimized_fit_params(struct file_name_fit_params file_and_fit_params)
{
    io::file::CSV_IO csv_io;
    std::string full_path = file_and_fit_params.dataset_dir+"/output/"+file_and_fit_params.dataset_filename+std::to_string(file_and_fit_params.detector_num)+".csv";
    logit<<full_path<<std::endl;
    csv_io.save_fit_parameters(full_path, file_and_fit_params.fit_params );

#ifdef _BUILD_WITH_VTK
    fitting::models::Gaussian_Model model;
    //Range of energy in spectra to fit
    fitting::models::Range energy_range;
    energy_range.min = 0;
    energy_range.max = file_and_fit_params.spectra.size() -1;

    data_struct::xrf::Spectra model_spectra = model.model_spectrum(&file_and_fit_params.fit_params, &file_and_fit_params.elements_to_fit, energy_range);
    std::string str_path = file_and_fit_params.dataset_dir+"/output/fit_"+file_and_fit_params.dataset_filename+"_det"+std::to_string(file_and_fit_params.detector_num)+".png";
    visual::SavePlotSpectras(str_path, file_and_fit_params.spectra, model_spectra, true);
#endif

    for(auto &itr : file_and_fit_params.elements_to_fit)
    {
        delete itr.second;
    }

}

// ----------------------------------------------------------------------------

void save_averaged_fit_params(std::string dataset_dir, std::vector<data_struct::xrf::Fit_Parameters> fit_params_avgs, int detector_num_start, int detector_num_end)
{
    io::file::aps::APS_Fit_Params_Import aps_io;
    int i =0;
    std::string full_path = dataset_dir+"/maps_fit_parameters_override.txt";
    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        aps_io.save(full_path, fit_params_avgs[i], detector_num );
        i++;
    }
}

// ----------------------------------------------------------------------------

bool load_quantification_standard(std::string dataset_directory,
                                  std::string quantification_info_file,
                                  std::string *standard_file_name,
                                  std::unordered_map<std::string, real_t> *element_standard_weights)
                                  //data_struct::xrf::Quantification_Standard * quantification_standard)
{
    std::string path = dataset_directory + quantification_info_file;
    std::ifstream paramFileStream(path);

    if (paramFileStream.is_open() )
    {
        paramFileStream.exceptions(std::ifstream::failbit);
        bool has_filename = false;
        bool has_elements = false;
        bool has_weights = false;
        //std::string line;
        std::string tag;

        std::vector<std::string> element_names;
        std::vector<real_t> element_weights;

        try
        {

            for (std::string line; std::getline(paramFileStream, line); )
            //while(std::getline(paramFileStream, line))
            {
                std::istringstream strstream(line);
                std::getline(strstream, tag, ':');
                //logit<<"tag : "<<tag<<std::endl;
                if (tag == "FILENAME")
                {
                    std::string standard_filename;
                    logit << line << std::endl;
                    std::getline(strstream, standard_filename, ':');
                    standard_filename.erase(std::remove_if(standard_filename.begin(), standard_filename.end(), ::isspace), standard_filename.end());
                    logit << "Standard file name = "<< standard_filename << std::endl;
                    //quantification_standard->standard_filename(standard_filename);
                    (*standard_file_name) = standard_filename;
                    has_filename = true;
                }
                else if (tag == "ELEMENTS_IN_STANDARD")
                {
                    std::string element_symb;
                    while(std::getline(strstream, element_symb, ','))
                    {
                        element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());
                        logit<<"Element : "<<element_symb<<std::endl;
                        element_names.push_back(element_symb);
                    }
                    has_elements = true;
                }
                else if (tag == "WEIGHT")
                {
                    std::string element_weight_str;
                    while(std::getline(strstream, element_weight_str, ','))
                    {
                        element_weight_str.erase(std::remove_if(element_weight_str.begin(), element_weight_str.end(), ::isspace), element_weight_str.end());
                        logit<<"Element weight: "<<element_weight_str<<std::endl;
                        real_t weight = std::stof(element_weight_str);
                        element_weights.push_back(weight);
                    }
                    has_weights = true;
                }

            }
        }
        catch(std::exception e)
        {
            if (paramFileStream.eof() == 0 && (paramFileStream.bad() || paramFileStream.fail()) )
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << paramFileStream.fail()
                    << "\neofbit: " << paramFileStream.eof()
                    << "\nbadbit: " << paramFileStream.bad() << std::endl;
            }
        }


        paramFileStream.close();
        if(has_filename && has_elements && has_weights)
        {
            if(element_names.size() == element_weights.size())
            {
                for(size_t i=0; i<element_names.size(); i++)
                {
                    (*element_standard_weights)[element_names[i]] = element_weights[i];
                    //quantification_standard->append_element(element_names[i], element_weights[i]);
                }
            }
            else
            {
                logit<<"Error: number of element names ["<<element_names.size()<<"] does not match number of element weights ["<<element_weights.size()<<"]!"<<std::endl;
            }

            return true;
        }

    }
    else
    {
        logit<<"Failed to open file "<<path<<std::endl;
    }
    return false;

}

// ----------------------------------------------------------------------------

bool load_override_params(std::string dataset_directory,
                          int detector_num,
                          data_struct::xrf::Params_Override *params_override)
{
    //Importer for APS datasets
    io::file::aps::APS_Fit_Params_Import fit_param_importer;

    std::string det_num = "";
    if(detector_num > -1)
        det_num = std::to_string(detector_num);

    if(false == fit_param_importer.load(dataset_directory+"maps_fit_parameters_override.txt"+det_num,
                                        data_struct::xrf::Element_Info_Map::inst(),
                                        params_override))
    {
        logit<<"Error loading fit param override file: "<<std::endl;
        return false;
    }
    else
    {
        data_struct::xrf::Element_Info* detector_element;
        if(params_override->detector_element.length() > 0)
        {
            // Get the element info class                                   // detector element as string "Si" or "Ge" usually

            detector_element = data_struct::xrf::Element_Info_Map::inst()->get_element(params_override->detector_element);
        }
        else
        {
         //log error or warning
            logit<<"Error, no detector material defined in maps_fit_parameters_override.txt . Defaulting to Si"<<std::endl;
            detector_element = data_struct::xrf::Element_Info_Map::inst()->get_element("Si");
        }


        //add compton and coherant amp
        if(params_override->elements_to_fit.count(data_struct::xrf::STR_COMPTON_AMPLITUDE) == 0)
        {
            params_override->elements_to_fit.insert(std::pair<std::string, data_struct::xrf::Fit_Element_Map*>(data_struct::xrf::STR_COMPTON_AMPLITUDE, new data_struct::xrf::Fit_Element_Map(data_struct::xrf::STR_COMPTON_AMPLITUDE, nullptr)) );
        }
        if(params_override->elements_to_fit.count(data_struct::xrf::STR_COHERENT_SCT_AMPLITUDE) == 0)
        {
            params_override->elements_to_fit.insert(std::pair<std::string, data_struct::xrf::Fit_Element_Map*>(data_struct::xrf::STR_COHERENT_SCT_AMPLITUDE, new data_struct::xrf::Fit_Element_Map(data_struct::xrf::STR_COHERENT_SCT_AMPLITUDE, nullptr)) );
        }

        logit<<"Elements to fit:  ";
        //Update element ratios by detector element
        for(auto& itr : params_override->elements_to_fit)
        {
            itr.second->init_energy_ratio_for_detector_element(detector_element);
            logit_s<<itr.first<<" ";
        }
        logit_s<<std::endl;

    }

    return true;
}

// ----------------------------------------------------------------------------

bool load_spectra_volume(std::string dataset_directory,
                         std::string dataset_file,
                         data_struct::xrf::Spectra_Volume *spectra_volume,
                         size_t detector_num,
                         data_struct::xrf::Params_Override * params_override,
                         data_struct::xrf::Quantification_Standard * quantification_standard,
                         bool save_scalers)
{

    //Dataset importer
    io::file::MDA_IO mda_io;
    //data_struct::xrf::Detector detector;
    std::string tmp_dataset_file = dataset_file;

    logit<<"Loading dataset "<<dataset_directory+"mda/"+dataset_file<<" detector "<<detector_num<<std::endl;

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    bool hasNetcdf = false;
    bool hasHdf = false;
    std::string file_middle = ""; //_2xfm3_ or dxpM...
    for(auto &itr : netcdf_files)
    {
        if (itr.find(tmp_dataset_file) == 0)
        {
            size_t slen = (itr.length()-4) - tmp_dataset_file.length();
            file_middle = itr.substr(tmp_dataset_file.length(), slen);
            hasNetcdf = true;
            break;
        }
    }
    if (hasNetcdf == false)
    {
        for(auto &itr : hdf_files)
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length()-4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasHdf = true;
                break;
            }
        }
    }

    //load spectra
    if (false == mda_io.load_spectra_volume(dataset_directory+"mda/"+dataset_file, detector_num, spectra_volume, hasNetcdf | hasHdf, params_override, quantification_standard) )
    {
        logit<<"Error load spectra "<<dataset_directory+"mda/"+dataset_file<<std::endl;
        return false;
    }
    else
    {
        if(hasNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for(size_t i=0; i<spectra_volume->rows(); i++)
                {
                    full_filename = dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + std::to_string(i) + ".nc";
                    //todo: add verbose option
                    //logit<<"Loading file "<<full_filename<<std::endl;
                    io::file::NetCDF_IO::inst()->load_spectra_line(full_filename, detector_num, &(*spectra_volume)[i]);
                }
            }
            else
            {
                logit<<"Did not find netcdf files "<<dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc"<<std::endl;
                //return false;
            }
        }
        else if (hasHdf)
        {
            io::file::HDF5_IO::inst()->load_spectra_volume(dataset_directory + "flyXRF.h5/" + tmp_dataset_file + file_middle + "0.h5", detector_num, spectra_volume);
        }

    }

    if(save_scalers)
    {
        io::file::HDF5_IO::inst()->start_save_seq();
        io::file::HDF5_IO::inst()->save_scan_scalers(detector_num, mda_io.get_scan_ptr(), params_override, hasNetcdf | hasHdf);
    }

    mda_io.unload();
    logit<<"Finished Loading dataset "<<dataset_directory+"mda/"+dataset_file<<" detector "<<detector_num<<std::endl;
    return true;
}

// ----------------------------------------------------------------------------

bool load_and_integrate_spectra_volume(std::string dataset_directory,
                                       std::string dataset_file,
                                       data_struct::xrf::Spectra *integrated_spectra,
                                       size_t detector_num,
                                       data_struct::xrf::Params_Override * params_override)
{
    //Dataset importer
    io::file::MDA_IO mda_io;
    //data_struct::xrf::Detector detector;
    std::string tmp_dataset_file = dataset_file;
    bool ret_val = true;

    logit<<"Loading dataset "<<dataset_directory+"mda/"+dataset_file<<std::endl;

    //TODO: check if any analyized mda.h5 files are around to load. They should have integrated spectra saved already.

    //check if we have a netcdf file associated with this dataset.
    tmp_dataset_file = tmp_dataset_file.substr(0, tmp_dataset_file.size()-4);
    bool hasNetcdf = false;
    bool hasHdf = false;
    std::string file_middle = ""; //_2xfm3_ or dxpM...
    for(auto &itr : netcdf_files)
    {
        if (itr.find(tmp_dataset_file) == 0)
        {
            size_t slen = (itr.length()-4) - tmp_dataset_file.length();
            file_middle = itr.substr(tmp_dataset_file.length(), slen);
            hasNetcdf = true;
            break;
        }
    }
    if (hasNetcdf == false)
    {
        for(auto &itr : hdf_files)
        {
            if (itr.find(tmp_dataset_file) == 0)
            {
                size_t slen = (itr.length()-4) - tmp_dataset_file.length();
                file_middle = itr.substr(tmp_dataset_file.length(), slen);
                hasHdf = true;
                break;
            }
        }
    }

    //load spectra
    if (false == hasNetcdf && false == hasHdf)
    {
        data_struct::xrf::Spectra_Volume spectra_volume;
        ret_val = mda_io.load_spectra_volume(dataset_directory+"mda/"+dataset_file, detector_num, &spectra_volume, hasNetcdf | hasHdf, params_override, nullptr);
        if(ret_val)
        {
            *integrated_spectra = spectra_volume.integrate();
        }
        mda_io.unload();
    }
    else
    {
        int rank;
        int dims[10];
        dims[0] = 0;
        rank = mda_io.get_rank_and_dims(dataset_directory + "mda/" + dataset_file, &dims[0]);
        if(rank == 3)
        {
            integrated_spectra->resize(dims[2]);
        }
        else
        {
            integrated_spectra->resize(2048);
        }


        if(hasNetcdf)
        {
            std::ifstream file_io(dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc");
            if(file_io.is_open())
            {
                file_io.close();
                std::string full_filename;
                for(size_t i=0; i<dims[0]; i++)
                {
                    data_struct::xrf::Spectra_Line spectra_line;
                    spectra_line.resize(dims[1], integrated_spectra->size());
                    full_filename = dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + std::to_string(i) + ".nc";
                    //logit<<"Loading file "<<full_filename<<std::endl;
                    if( io::file::NetCDF_IO::inst()->load_spectra_line(full_filename, detector_num, &spectra_line) )
                    {
                        for(int k=0; k<spectra_line.size(); k++)
                        {
                            *integrated_spectra += spectra_line[k];
                        }
                    }
                }
            }
            else
            {
                logit<<"Did not find netcdf files "<<dataset_directory + "flyXRF/" + tmp_dataset_file + file_middle + "0.nc"<<std::endl;
                //return false;
            }
        }
        else if (hasHdf)
        {
            ret_val = io::file::HDF5_IO::inst()->load_and_integrate_spectra_volume(dataset_directory + "flyXRF.h5/" + tmp_dataset_file + file_middle + "0.h5", detector_num, integrated_spectra);
        }

    }

    logit<<"Finished Loading dataset "<<dataset_directory+"mda/"+dataset_file<<" detector "<<detector_num<<std::endl;
    return ret_val;
}

// ----------------------------------------------------------------------------

void generate_h5_averages(std::string dataset_directory,
                          std::string dataset_file,
                          ThreadPool* tp,
                          size_t detector_num_start,
                          size_t detector_num_end)
{
    logit<<std::endl;

    std::vector<std::string> hdf5_filenames;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    if (detector_num_start == detector_num_end)
    {
        logit << "Warning: detector range "<<detector_num_start<<":"<<detector_num_end<<" is only 1 detector. Nothing to avg."<<std::endl;
        return;
    }


    for(size_t detector_num = detector_num_start; detector_num <= detector_num_end; detector_num++)
    {
        hdf5_filenames.push_back(dataset_directory+"img.dat/"+dataset_file+".h5"+std::to_string(detector_num));
    }

    io::file::HDF5_IO::inst()->generate_avg(dataset_directory+"img.dat/"+dataset_file+".h5", hdf5_filenames);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logit << "elapsed time: " << elapsed_seconds.count() << "s"<<std::endl;


}

// ----------------------------------------------------------------------------

std::vector<std::string> find_all_dataset_files(std::string dataset_directory, std::string search_str)
{
    std::vector<std::string> dataset_files;
    logit<<dataset_directory<<" searching for "<<search_str<<std::endl;
    DIR *dir;
    struct dirent *ent;
    size_t search_str_size = search_str.length();
    if ((dir = opendir (dataset_directory.c_str())) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL)
        {
            std::string fname(ent->d_name);
            // check if extension is .mda
            if (fname.size() > 4)
            {
                if (fname.rfind(search_str) == fname.size() -search_str_size)
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
        logit<<"Error: could not open directory "<<dataset_directory<<std::endl;
    }

    logit<<"found "<<dataset_files.size()<<std::endl;
    return dataset_files;
}

// ----------------------------------------------------------------------------

void check_and_create_dirs(std::string dataset_directory)
{

    bool found_img_dat = false;
    logit<<dataset_directory<<std::endl;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dataset_directory.c_str())) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir (dir)) != NULL)
        {
            if( strcmp(ent->d_name , "img.dat") == 0)
            {
                found_img_dat = true;
                break;
            }
        }
        closedir (dir);
    }
    else
    {
        /* could not open directory */
        logit<<"Error: could not open directory "<<dataset_directory<<std::endl;
    }

    if (false == found_img_dat)
    {
        std::string cmd = "mkdir "+dataset_directory+"img.dat";
        system(cmd.c_str());
    }
    logit<<"done"<<std::endl;

}

// ----------------------------------------------------------------------------

void sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string> *dataset_files)
{

    io::file::MDA_IO mda_io;
    logit<<dataset_directory<<" "<<dataset_files->size()<<" files"<<std::endl;
    std::list<file_name_size> f_list;

    for (auto &itr : *dataset_files)
    {
        std::string full_path = dataset_directory + "/mda/"+itr;
        int fsize = mda_io.get_multiplied_dims(full_path);
        f_list.push_back(file_name_size(itr, fsize));
    }

    f_list.sort(compare_file_size);

    dataset_files->clear();

    for(auto &itr : f_list)
    {
        dataset_files->push_back(itr.filename);
    }

    logit<<"done"<<std::endl;
}

}// end namespace io
