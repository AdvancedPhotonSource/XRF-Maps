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


#ifndef HDF5_IO_H
#define HDF5_IO_H


#include <list>
#include <mutex>
#include <queue>
#include <future>
#include <stack>
#include <type_traits>
#include "hdf5.h"
#include "data_struct/spectra_volume.h"
#include "data_struct/fit_element_map.h"
#include "data_struct/detector.h"
#include "data_struct/params_override.h"
#include "data_struct/scan_info.h"

#include "core/mem_info.h"
#include "core/defines.h"

namespace io
{
namespace file
{
enum H5_OBJECTS{H5O_FILE, H5O_GROUP, H5O_DATASPACE, H5O_DATASET, H5O_ATTRIBUTE, H5O_PROPERTY};

enum H5_SPECTRA_LAYOUTS {MAPS_RAW, MAPS_V9, MAPS_V10, XSPRESS, APS_SEC20};

enum GSE_CARS_SAVE_VER {UNKNOWN, XRFMAP, XRMMAP};

class DLL_EXPORT HDF5_IO
{
public:

    static HDF5_IO* inst();

    ~HDF5_IO();

    template<typename T_real>
    bool load_spectra_volume(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol);

    template<typename T_real>
    bool load_spectra_volume_with_callback(std::string path,
											const std::vector<size_t>& detector_num_arr,
										   data_struct::IO_Callback_Func_Def<T_real> callback_func,
                                           void* user_data);

    template<typename T_real>
	bool load_spectra_volume_emd_with_callback(std::string path,
												const std::vector<size_t>& detector_num_arr,
												data_struct::IO_Callback_Func_Def<T_real> callback_func,
												void* user_data);

    template<typename T_real>
    bool load_spectra_volume_emd(std::string path,
                                 size_t frame_num,
                                 data_struct::Spectra_Volume<T_real> *spec_vol,
                                 bool logerr = true);

    template<typename T_real>
    bool load_spectra_line_xspress3(std::string path, size_t detector_num, data_struct::Spectra_Line<T_real>* spec_row);

    template<typename T_real>
    bool load_spectra_volume_confocal(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol, bool log_error=true);

    template<typename T_real>
	bool load_spectra_volume_gsecars(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol, bool log_error = true);

    template<typename T_real>
    bool load_spectra_volume_bnl(std::string path, size_t detector_num, data_struct::Spectra_Volume<T_real>* spec_vol, bool log_error = true);

    template<typename T_real>
    bool load_integrated_spectra_bnl(std::string path, size_t detector_num, data_struct::Spectra<T_real>* spec, bool log_error);

    template<typename T_real>
    bool load_and_integrate_spectra_volume(std::string path, size_t detector_num, data_struct::Spectra<T_real>* spectra);

    template<typename T_real>
    bool load_spectra_vol_analyzed_h5(std::string path,
                                      data_struct::Spectra_Volume<T_real>* spectra_volume,
                                      int row_idx_start = 0,
                                      int row_idx_end = -1,
                                      int col_idx_start = 0,
                                      int col_idx_end = -1);

    template<typename T_real>
    bool load_integrated_spectra_analyzed_h5(std::string path, data_struct::Spectra<T_real>* spectra, bool log_error=true);

    template<typename T_real>
    bool load_quantification_scalers_analyzed_h5(std::string path, data_struct::Params_Override<T_real> *override_values);

    template<typename T_real>
    bool load_quantification_scalers_gsecars(std::string path, data_struct::Params_Override<T_real> *override_values);

    template<typename T_real>
    bool load_quantification_scalers_BNL(std::string path, data_struct::Params_Override<T_real>* override_values);

    bool generate_avg(std::string avg_filename, std::vector<std::string> files_to_avg);

    bool generate_stream_dataset(std::string dataset_directory,
                                 std::string dataset_name,
                                 int detector_num,
                                 size_t height,
                                 size_t width);

    template<typename T_real>
    bool get_scalers_and_metadata_emd(std::string path, data_struct::Scan_Info<T_real>* scan_info);

    template<typename T_real>
    bool get_scalers_and_metadata_confocal(std::string path, data_struct::Scan_Info<T_real>* scan_info);

    template<typename T_real>
    bool get_scalers_and_metadata_gsecars(std::string path, data_struct::Scan_Info<T_real>* scan_info);

    template<typename T_real>
    bool get_scalers_and_metadata_bnl(std::string path, data_struct::Scan_Info<T_real>* scan_info);

    template<typename T_real>
    bool save_stream_row(size_t d_hash,
                         size_t detector_num,
                         size_t row,
                         std::vector< data_struct::Spectra<T_real>* >  *spectra_row);

    template<typename T_real>
    bool save_itegrade_spectra(data_struct::Spectra<T_real> * spectra);

    bool close_dataset(size_t d_hash);

    bool start_save_seq(const std::string filename, bool force_new_file=false);

    bool start_save_seq(bool force_new_file=false){ return start_save_seq(_cur_filename, force_new_file);}

    void set_filename(std::string fname) {_cur_filename = fname;}

    template<typename T_real>
    bool save_spectra_volume(const std::string path,
                            data_struct::Spectra_Volume<T_real>* spectra_volume,
                             size_t row_idx_start=0,
                             int row_idx_end=-1,
                             size_t col_idx_start=0,
                             int col_idx_end=-1);

    template<typename T_real>
    bool save_energy_calib(int spectra_size, T_real energy_offset, T_real energy_slope, T_real energy_quad);

    template<typename T_real>
    bool save_element_fits(const std::string path,
                           const data_struct::Fit_Count_Dict<T_real>* const element_counts,
                           size_t row_idx_start=0,
                           int row_idx_end=-1,
                           size_t col_idx_start=0,
                           int col_idx_end=-1);

    template<typename T_real>
    bool save_fitted_int_spectra(const std::string path,
                                 const data_struct::Spectra<T_real>& spectra,
                                 const data_struct::Range& range,
                                 const data_struct::Spectra<T_real>& background,
								 const size_t save_spectra_size);
	
    template<typename T_real>
	bool save_max_10_spectra(const std::string path,
							const data_struct::Range& range,
							const data_struct::Spectra<T_real>& max_spectra,
							const data_struct::Spectra<T_real>& max_10_spectra,
                            const data_struct::Spectra<T_real>& fit_int_background);

    template<typename T_real>
    bool save_quantification(data_struct::Detector<T_real>* detector);

    template<typename T_real>
    bool save_scan_scalers(size_t detector_num,
                           data_struct::Scan_Info<T_real>* scan_info,
                           data_struct::Params_Override<T_real> * params_override,
                           size_t row_idx_start=0,
                           int row_idx_end=-1,
                           size_t col_idx_start=0,
                           int col_idx_end=-1);

    template<typename T_real>
    bool save_scan_scalers_confocal(std::string path,
                                    size_t detector_num,
                                    size_t row_idx_start=0,
                                    int row_idx_end=-1,
                                    size_t col_idx_start=0,
                                    int col_idx_end=-1);

    template<typename T_real>
	bool save_scan_scalers_gsecars(std::string path,
								size_t detector_num,
								size_t row_idx_start = 0,
								int row_idx_end = -1,
								size_t col_idx_start = 0,
								int col_idx_end = -1);

    template<typename T_real>
    bool save_scan_scalers_bnl(std::string path,
        size_t detector_num,
        size_t row_idx_start = 0,
        int row_idx_end = -1,
        size_t col_idx_start = 0,
        int col_idx_end = -1);

	// Add links to dataset and set version to 9 so legacy software can load it
    void add_v9_layout(std::string dataset_file);

	// Add exchange layout to be loadable by external software
    void add_exchange_layout(std::string dataset_file);

	// update theta value based on new pv name
    void update_theta(std::string dataset_file, std::string theta_pv_str);

	void update_amps(std::string dataset_file, std::string us_amp_str, std::string ds_amp_str);

	void update_quant_amps(std::string dataset_file, std::string us_amp_str, std::string ds_amp_str);

    //update scalers if maps_fit_parameters_override.txt has changes pv's and you don't want to refit
    //void update_scalers(std::string dataset_file, data_struct::Params_Override<T_real>* params_override);
    
    //export integrated spec, fitted, background into csv
    void export_int_fitted_to_csv(std::string dataset_file);

    bool end_save_seq(bool loginfo=true);

    template<typename T_real>
    bool add_background(std::string directory, std::string filename, data_struct::Params_Override<T_real>& params);

private:

    HDF5_IO();

    static HDF5_IO *_this_inst;

    static std::mutex _mutex;

    template<typename T_real>
	bool _load_integrated_spectra_analyzed_h5(hid_t file_id, data_struct::Spectra<T_real>* spectra);

    template<typename T_real>
    bool _save_scan_meta_data(hid_t scan_grp_id, data_struct::Scan_Meta_Info<T_real>* meta_info);

	bool _save_extras(hid_t scan_grp_id, std::vector<data_struct::Extra_PV>* extra_pvs);

    template<typename T_real>
    bool _save_scalers(hid_t maps_grp_id, std::vector<data_struct::Scaler_Map<T_real>>*scalers_map, T_real us_amps_val, std::string us_amps_unit, T_real ds_amps_val, string ds_amps_unit);
    
    template<typename T_real>
    void _save_amps(hid_t scalers_grp_id, T_real us_amp_sens_num_val, string us_amp_sens_unit_val, T_real ds_amp_sens_num_val, string ds_amp_sens_unit_val);
    
    template<typename T_real>
	bool _save_params_override(hid_t group_id, data_struct::Params_Override<T_real> * params_override);

    void _gen_average(std::string full_hdf5_path, std::string dataset_name, hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids, bool avg=true);
    
    void _generate_avg_analysis(hid_t src_maps_grp_id, hid_t dst_maps_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids);
    
    void _generate_avg_integrated_spectra(hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids);

    void _add_v9_quant(hid_t file_id, hid_t chan_names, hid_t chan_space, int chan_amt, std::string quant_str, std::string new_loc);
    void _add_v9_scalers(hid_t file_id);
    void _add_extra_pvs(hid_t file_id, std::string group_name);

    bool _add_exchange_meta(hid_t file_id, std::string exchange_idx, std::string fits_link, std::string normalize_scaler);
	
    bool _open_h5_object(hid_t &id, H5_OBJECTS obj, std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map, std::string s1, hid_t id2, bool log_error=true, bool close_on_fail=true);
    bool _open_or_create_group(const std::string name, hid_t parent_id, hid_t& out_id, bool log_error = true, bool close_on_fail = true);
    bool _create_memory_space(int rank, const hsize_t* count, hid_t& out_id);
    bool _open_h5_dataset(const std::string& name, hid_t data_type, hid_t parent_id, int dims_size, const hsize_t* dims, const hsize_t* chunk_dims, hid_t& out_id, hid_t& out_dataspece);
    void _close_h5_objects(std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map);



    template<typename T_real>
    herr_t _read_h5d()
    {
        if (std::is_same<T_real, float>::value)
        {
            //return H5Dread(dset_id, H5T_NATIVE_FLOAT, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
        }
        else if (std::is_same<T_real, double>::value)
        {
            //return H5Dread(dset_id, H5T_NATIVE_DOUBLE, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
        }
        
        return -1;
        
        /*
        if constexpr (std::is_same_v<T_real, double>)
        {
            return H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
        }
        else if constexpr (std::is_same_v<T_real, float>)
        {
            return H5Dread(dset_id, H5T_NATIVE_REAL, memoryspace_id, dataspace_id, H5P_DEFAULT, buffer);
        }
        else
        {
            return -1;
        }
        */
    }
    

    hid_t _cur_file_id;
    std::string _cur_filename;
    std::stack<std::pair<hid_t, H5_OBJECTS> > _global_close_map;

};


}// end namespace file
}// end namespace io

#endif // HDF5_IO_H
