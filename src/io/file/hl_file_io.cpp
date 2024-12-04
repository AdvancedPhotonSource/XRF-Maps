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

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include "data_struct/scaler_lookup.h"
#include "yaml-cpp/yaml.h"

namespace io
{

namespace file
{
//-----------------------------------------------------------------------------

void parse_scalers(const std::string& beamline, const YAML::Node& node, bool time_normalized)
{
	for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
	{
		switch (it->second.Type())
		{
		case YAML::NodeType::Scalar:
            data_struct::Scaler_Lookup::inst()->add_beamline_scaler(beamline, it->first.as<std::string>(), it->second.as<std::string>(), time_normalized);
			break;
		case YAML::NodeType::Sequence:
            for (YAML::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
            {
                data_struct::Scaler_Lookup::inst()->add_beamline_scaler(beamline, it->first.as<std::string>(), it2->as<std::string>(), time_normalized);
            }
			break;
		case YAML::NodeType::Map:
		case YAML::NodeType::Null:
		case YAML::NodeType::Undefined:
		default:
			break;
		}
	}
}

// ----------------------------------------------------------------------------

void parse_summed_scalers(const std::string& beamline, const YAML::Node& node)
{
	for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
	{
        std::vector<std::string> scaler_list;
		switch (it->second.Type())
		{
		case YAML::NodeType::Sequence:
			for (YAML::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
                scaler_list.push_back(it2->as<std::string>());
			}
            data_struct::Scaler_Lookup::inst()->add_summed_scaler(beamline, it->first.as<std::string>(), scaler_list);
			break;
		case YAML::NodeType::Map:
		case YAML::NodeType::Scalar:
		case YAML::NodeType::Null:
		case YAML::NodeType::Undefined:
		default:
			break;
		}
	}
}

// ----------------------------------------------------------------------------

void parse_time_normalized_scalers(const std::string& beamline, const YAML::Node& node)
{
	for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
	{
        std::string val = it->first.as<std::string>();
        switch (it->second.Type())
        {
        case YAML::NodeType::Sequence:
            if (val == STR_TIMING)
            {
                YAML::Node timing_node = it->second.as<YAML::Node>();
                if (timing_node.size() == 2)
                {
                    data_struct::Scaler_Lookup::inst()->add_timing_info(beamline, timing_node[0].as<std::string>(), timing_node[1].as<double>());
                }
                else
                {
                    logW << " Couldnt parse yaml timing info \n";
                }
            }
            break;
		case YAML::NodeType::Map:
            if (val == STR_SCALERS)
			{
				parse_scalers(beamline, it->second.as<YAML::Node>(), true);
			}
            break;
        case YAML::NodeType::Scalar:
		case YAML::NodeType::Null:
		case YAML::NodeType::Undefined:
		default:
			break;
		}
	}
}

// ----------------------------------------------------------------------------

void parse_beamline(const std::string& beamline, const YAML::Node& node)
{
	for (YAML::const_iterator it = node.begin(); it != node.end(); ++it)
	{
        std::string val = it->first.as<std::string>();
		switch (it->second.Type())
		{
		case YAML::NodeType::Map:
			if (val == STR_SCALERS)
			{
				parse_scalers(beamline, it->second.as<YAML::Node>(), false);
			}
			else if (val == STR_TIME_NORMALIZED_SCALERS)
			{
				parse_time_normalized_scalers(beamline, it->second.as<YAML::Node>());
			}
			else if (val == STR_SUMMED_SCALERS)
			{
				parse_summed_scalers(beamline, it->second.as<YAML::Node>());
			}
			break;
		case YAML::NodeType::Scalar:
		case YAML::NodeType::Sequence:
		case YAML::NodeType::Null:
		case YAML::NodeType::Undefined:
		default:
			break;
		}
	}
}

// ----------------------------------------------------------------------------

bool load_scalers_lookup(const std::string filename)
{
    YAML::Node node = YAML::LoadFile(filename);

	if (node[STR_BEAMLINES])
	{
		YAML::Node beamlines_node = node[STR_BEAMLINES];

		for (YAML::const_iterator it = beamlines_node.begin(); it != beamlines_node.end(); ++it)
		{
			parse_beamline(it->first.as<std::string>(), it->second.as<YAML::Node>());
		}
		return true;
	}  
    return false;
}

// ----------------------------------------------------------------------------

void save_optimized_fit_params(std::string dataset_dir, std::string dataset_filename, int detector_num, std::string result, data_struct::Fit_Parameters<double> *fit_params, const data_struct::Spectra<double>* const spectra, const data_struct::Fit_Element_Map_Dict<double>* const elements_to_fit)
{
    std::string full_path = dataset_dir + DIR_END_CHAR + "output" + DIR_END_CHAR + STR_FIT_SPEC_DIR + DIR_END_CHAR + dataset_filename;
    std::string mca_full_path = dataset_dir + DIR_END_CHAR + "output" + DIR_END_CHAR + STR_INT_SPEC_DIR + DIR_END_CHAR + "intspec_" + dataset_filename;
    std::string fp_full_path = dataset_dir + DIR_END_CHAR + "output" + DIR_END_CHAR + STR_FIT_PARAM_DIR + DIR_END_CHAR + "fit_param_" + dataset_filename;
    
    if (detector_num != -1)
    {
        full_path += std::to_string(detector_num) + ".csv";
        mca_full_path += std::to_string(detector_num) + ".mca";
        fp_full_path += std::to_string(detector_num) + ".csv";
    }
    else
    {
        full_path += ".csv";
        mca_full_path += ".txt";
        fp_full_path += ".txt";
    }
    logI<<full_path<<"\n";

    if (fit_params == nullptr)
    {
        logE << "Fit Parameters == nullptr. Can not save!\n";
		return;
    }
	else
	{
		if (fit_params->size() == 0)
		{
			logE << "Fit Parameters size = 0. Can not save!\n";
			return;
		}
	}

    if (spectra == nullptr)
    {
        logE << "int Spectra == nullptr. Can not save!\n";
		return;
    }

    if (elements_to_fit == nullptr)
    {
        logE << "Elements to Fit == nullptr. Can not save!\n";
		return;
    }

    fitting::models::Gaussian_Model<double> model;
    //Range of energy in spectra to fit
    fitting::models::Range energy_range = data_struct::get_energy_range(spectra->size(), fit_params);
    data_struct::Spectra<double> snip_spectra = spectra->sub_spectra(energy_range.min, energy_range.count());

    std::unordered_map<std::string, ArrayTr<double>> labeled_spectras;
    data_struct::Spectra<double> model_spectra = model.model_spectrum(fit_params, elements_to_fit, &labeled_spectras, energy_range);
    
    data_struct::ArrayTr<double> background;

    double energy_offset = fit_params->value(STR_ENERGY_OFFSET);
    double energy_slope = fit_params->value(STR_ENERGY_SLOPE);
    double energy_quad = fit_params->value(STR_ENERGY_QUADRATIC);
    const data_struct::ArrayTr<double> ev = data_struct::generate_energy_array(energy_range, energy_offset, energy_slope, energy_quad);


    if (fit_params->contains(STR_SNIP_WIDTH))
	{
        data_struct::ArrayTr<double> s_background = data_struct::snip_background<double>(spectra,
                                                                        fit_params->value(STR_ENERGY_OFFSET),
                                                                        fit_params->value(STR_ENERGY_SLOPE),
                                                                        fit_params->value(STR_ENERGY_QUADRATIC),
                                                                        fit_params->value(STR_SNIP_WIDTH),
                                                                        energy_range.min,
                                                                        energy_range.max);
        if (s_background.size() >= energy_range.count())
        {
            background = s_background.segment(energy_range.min, energy_range.count());
            model_spectra += background;
        }
        if (background.size() == 0)
        {
            background.resize(energy_range.count());
            background.setZero();
        }
	}
    else
    {
        background.resize(energy_range.count());
        background.setZero();
    }

    std::string str_path = dataset_dir + "/output/" +  STR_FIT_SPEC_DIR + "fit_" + dataset_filename + "_det";
    if (detector_num != -1)
    {
        str_path += std::to_string(detector_num) + ".png";
    }
    else
    {
        str_path += ".png";
    }
    #ifdef _BUILD_WITH_QT
    visual::SavePlotSpectrasFromConsole(str_path, &ev, &snip_spectra, &model_spectra, &background, true);
    #endif
    for (const auto& eitr : *elements_to_fit)
    {
        if (fit_params->contains(eitr.first))
        {
            double val = (*fit_params)[eitr.first].value;
            val = pow(10.0, val);
            (*fit_params)[eitr.first].value = val;
        }
    }

    io::file::csv::save_fit_and_int_spectra(full_path, &ev, &snip_spectra, &model_spectra, &background, &labeled_spectras);
    io::file::aps::save_fit_parameters_override(fp_full_path, *fit_params, result);
    std::unordered_map<std::string, double> scaler_map;
    scaler_map[STR_ENERGY_OFFSET] = energy_offset;
    scaler_map[STR_ENERGY_SLOPE] = energy_slope;
    scaler_map[STR_ENERGY_QUADRATIC] = energy_quad;
    io::file::mca::save_integrated_spectra(mca_full_path, spectra, scaler_map);

}

// ----------------------------------------------------------------------------

DLL_EXPORT bool load_quantification_standardinfo(std::string dataset_directory,
                                                 std::string quantification_info_file,
                                                std::vector<Quantification_Standard<double>> &standards)
{
    std::string path = dataset_directory + quantification_info_file;
    std::ifstream paramFileStream(path);

    if (paramFileStream.is_open() )
    {
        paramFileStream.exceptions(std::ifstream::failbit);
        bool has_filename = false;
        bool has_elements = false;
        bool disable_Ka_quant = false;
        bool disable_La_quant = false;
        std::string tag;

        std::vector<std::string> element_names;
        std::vector<double> element_weights;

        try
        {
            std::string standard_filename;
            for (std::string line; std::getline(paramFileStream, line); )
            {
                std::istringstream strstream(line);
                std::getline(strstream, tag, ':');
                //logD<<"tag : "<<tag<<"\n";
                if (tag == "FILENAME")
                {
                    standard_filename="";
                    logI << line << "\n";
                    std::getline(strstream, standard_filename, ':');
                    standard_filename.erase(std::remove_if(standard_filename.begin(), standard_filename.end(), ::isspace), standard_filename.end());
                    logI << "Standard file name = "<< standard_filename << "\n";
                    has_filename = true;
                }
                else if (tag == "DISABLE_KA")
                {
                    std::string tmp = "";
                    logI << line << "\n";
                    std::getline(strstream, tmp, ':');
                    tmp.erase(std::remove_if(tmp.begin(), tmp.end(), ::isspace), tmp.end());
                    if (tmp == "1")
                    {
                        disable_Ka_quant = true;
                    }
                    logI << "Disable Ka = " << disable_Ka_quant << "\n";
                }
                else if (tag == "DISABLE_LA")
                {
                    std::string tmp = "";
                    logI << line << "\n";
                    std::getline(strstream, tmp, ':');
                    tmp.erase(std::remove_if(tmp.begin(), tmp.end(), ::isspace), tmp.end());
                    if (tmp == "1")
                    {
                        disable_La_quant = true;
                    }
                    logI << "Disable La = " << disable_La_quant << "\n";
                }
                else if (tag == "ELEMENTS_IN_STANDARD")
                {
                    std::string element_symb;
                    while(std::getline(strstream, element_symb, ','))
                    {
                        element_symb.erase(std::remove_if(element_symb.begin(), element_symb.end(), ::isspace), element_symb.end());
                        logI<<"Element : "<<element_symb<<"\n";
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
                        logI<<"Element weight: "<<element_weight_str<<"\n";
                        double weight = parse_input_real<double>(element_weight_str);
                        element_weights.push_back(weight);
                    }
                    //has_weights = true;
                    if(has_elements && has_filename)
                    {
                        if(element_names.size() == element_weights.size())
                        {
                            standards.emplace_back(Quantification_Standard<double>(standard_filename, element_names, element_weights, disable_Ka_quant, disable_La_quant));
                        }
                        else
                        {
                            logE<<"Number of element names ["<<element_names.size()<<"] does not match number of element weights ["<<element_weights.size()<<"]!"<<"\n";
                        }
                    }
                    standard_filename = "";
                    element_names.clear();
                    element_weights.clear();
                    has_filename = false;
                    has_elements = false;
                    disable_Ka_quant = false;
                    disable_La_quant = false;
                }
            }
        }
        catch(std::exception& e)
        {
            if (paramFileStream.eof() == 0 && (paramFileStream.bad() || paramFileStream.fail()) )
            {
                std::cerr << "ios Exception happened: " << e.what() << "\n"
                    << "Error bits are: "
                    << "\nfailbit: " << paramFileStream.fail()
                    << "\neofbit: " << paramFileStream.eof()
                    << "\nbadbit: " << paramFileStream.bad() << "\n";
            }
        }


        paramFileStream.close();
    }
    else
    {
        logE<<"Failed to open file "<<path<<"\n";
        return false;
    }

    return true;
}

// ----------------------------------------------------------------------------

void generate_h5_averages(std::string dataset_directory,
                          std::string dataset_file,
					      const std::vector<size_t>& detector_num_arr,
                          bool append_h5_with_num)
{
    std::vector<std::string> hdf5_filenames;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    if (detector_num_arr.size() <= 1)
    {
        logW << "Warning: need more than 1 detector. Nothing to avg."<<"\n";
        return;
    }

    int len = dataset_file.length();
    if (dataset_file[len - 3] == '.' && dataset_file[len - 2] == 'h' && dataset_file[len - 1] == '5' && false == append_h5_with_num)
    {
        for (size_t detector_num : detector_num_arr)
        {
            hdf5_filenames.push_back(dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file + std::to_string(detector_num));
        }
    }
    else
    {
        for (size_t detector_num : detector_num_arr)
        {
            hdf5_filenames.push_back(dataset_directory + "img.dat" + DIR_END_CHAR + dataset_file + ".h5" + std::to_string(detector_num));
        }
    }
    io::file::HDF5_IO::inst()->generate_avg(dataset_directory+"img.dat"+ DIR_END_CHAR +dataset_file+".h5", hdf5_filenames);

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    logI << "elapsed time: " << elapsed_seconds.count() << "s"<<"\n";


}

// ----------------------------------------------------------------------------
/*
std::vector<std::string> find_all_dataset_files(std::string dataset_directory, std::string search_str)
{
    std::vector<std::string> dataset_files;
    logI<<dataset_directory<<" searching for "<<search_str<<"\n";
    DIR *dir;
    struct dirent *ent;
    size_t search_str_size = search_str.length();
    if ((dir = opendir (dataset_directory.c_str())) != NULL)
    {
        // print all the files and directories within directory 
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
        // could not open directory 
        logW<<"Could not open directory "<<dataset_directory<<" using search string "<<search_str<<"\n";
    }

    logI<<"found "<<dataset_files.size()<<"\n";
    return dataset_files;
}
*/
// ----------------------------------------------------------------------------

void check_and_create_dirs(std::string dataset_directory)
{

    bool found_img_dat = false;
    bool found_output = false;
    logI<<dataset_directory<<"\n";
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
            }
            if( strcmp(ent->d_name , "output") == 0)
            {
                found_output = true;
            }
			if(found_img_dat && found_output)
			{
                break;
			}
        }
        closedir (dir);
    }
    else
    {
        /* could not open directory */
        logW<<"Could not open directory "<<dataset_directory<<"\n";
    }

    if (false == found_img_dat)
    {
		int retval = system(nullptr);
        std::string cmd = "mkdir -p "+dataset_directory+"img.dat";
        logI << cmd << "\n";
        retval = system(cmd.c_str());
		if (retval != 0)
		{
			logE << "Could not create directory: " << cmd << " . May not be able to save results!\n";
		}
    }
    if (false == found_output)
    {
		int retval = system(nullptr);
        std::string cmd = "mkdir -p "+dataset_directory+"output";
        logI << cmd << "\n";
        retval = system(cmd.c_str());
		if (retval != 0)
		{
			logE << "Could not create directory: " << cmd << " . May not be able to save results!\n";
		}
    }

    int retval = system(nullptr);
    std::string cmd1 = "mkdir -p "+dataset_directory+"output"+DIR_END_CHAR+STR_FIT_SPEC_DIR;
    std::string cmd2 = "mkdir -p "+dataset_directory+"output"+DIR_END_CHAR+STR_INT_SPEC_DIR;
    std::string cmd3 = "mkdir -p "+dataset_directory+"output"+DIR_END_CHAR+STR_FIT_PARAM_DIR;
    
    retval = system(cmd1.c_str());
    logI << cmd1 << " = "<<retval<< "\n";
    retval = system(cmd2.c_str());
    logI << cmd2 << " = "<<retval<< "\n";
    retval = system(cmd3.c_str());
    logI << cmd3 << " = "<<retval<< "\n";
    /*
	if (retval != 0)
	{
		logE << "Could not create directory: " << cmd << " . May not be able to save results!\n";
	}
    */
    logI<<"done"<<"\n";

}

// ----------------------------------------------------------------------------
/*
void sort_dataset_files_by_size(std::string dataset_directory, std::vector<std::string> *dataset_files)
{
    // only supports soring mda files
    std::string ending = ".mda";
    io::file::MDA_IO<float> mda_io;
    logI<<dataset_directory<<" "<<dataset_files->size()<<" files"<<"\n";
    std::list<file_name_size> f_list;

    for (const auto &itr : *dataset_files)
    {
        //check if file ends with .mda
        if (itr.compare(itr.length() - ending.length(), ending.length(), ending) == 0 )
        {
            std::string full_path = dataset_directory + DIR_END_CHAR+"mda"+ DIR_END_CHAR +itr;
            long fsize = file::mda_get_multiplied_dims(full_path);
            f_list.push_back(file_name_size(itr, fsize));
        }
        else
        {
            std::string full_path = dataset_directory + DIR_END_CHAR + itr;
            std::ifstream in(full_path.c_str(), std::ifstream::ate | std::ifstream::binary);
            long fsize = in.tellg();
            if (fsize > -1)
            {
                f_list.push_back(file_name_size(itr, fsize));
            }
        }
    }
    
    f_list.sort(compare_file_size);
    if (f_list.size() > 0)
    {
        dataset_files->clear();

        for (auto& itr : f_list)
        {
            dataset_files->push_back(itr.filename);
        }
    }
    logI<<"done"<<"\n";
}
*/
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

}// end namespace file
}// end namespace io
