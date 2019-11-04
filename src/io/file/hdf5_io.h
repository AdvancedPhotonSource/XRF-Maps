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
#include "hdf5.h"
#include "data_struct/spectra_volume.h"
#include "data_struct/fit_element_map.h"

//Include mda data structures to save scalers
#include "io/file/mda_io.h"

#include "data_struct/quantification_standard.h"
#include "data_struct/params_override.h"

#include "core/mem_info.h"

namespace io
{
namespace file
{


enum H5_OBJECTS{H5O_FILE, H5O_GROUP, H5O_DATASPACE, H5O_DATASET, H5O_ATTRIBUTE};

enum H5_SPECTRA_LAYOUTS {MAPS_RAW, MAPS_V9, MAPS_V10, XSPRESS, APS_SEC20};

enum GSE_CARS_SAVE_VER {UNKNOWN, XRFMAP, XRMMAP};

class DLL_EXPORT HDF5_IO
{
public:

    static HDF5_IO* inst();

    ~HDF5_IO();

    bool load_spectra_volume(std::string path, size_t detector_num, data_struct::Spectra_Volume* spec_vol);

    bool load_spectra_volume_with_callback(std::string path,
											const std::vector<size_t>& detector_num_arr,
										   data_struct::IO_Callback_Func_Def callback_func,
                                           void* user_data);

	bool load_spectra_volume_emd_with_callback(std::string path,
												const std::vector<size_t>& detector_num_arr,
												data_struct::IO_Callback_Func_Def callback_func,
												void* user_data);

    bool load_spectra_volume_emd(std::string path,
                                 size_t frame_num,
                                 data_struct::Spectra_Volume *spec_vol);

    bool load_spectra_line_xspress3(std::string path, size_t detector_num, data_struct::Spectra_Line* spec_row);

    bool load_spectra_volume_confocal(std::string path, size_t detector_num, data_struct::Spectra_Volume* spec_vol, bool log_error=true);

	bool load_spectra_volume_gsecars(std::string path, size_t detector_num, data_struct::Spectra_Volume* spec_vol, bool log_error = true);

    bool load_and_integrate_spectra_volume(std::string path, size_t detector_num, data_struct::Spectra* spectra);

    bool load_spectra_vol_analyzed_h5(std::string path,
                                      data_struct::Spectra_Volume* spectra_volume,
                                      int row_idx_start = 0,
                                      int row_idx_end = -1,
                                      int col_idx_start = 0,
                                      int col_idx_end = -1);

    bool load_integrated_spectra_analyzed_h5(std::string path, data_struct::Spectra* spectra, bool log_error=true);

    bool load_quantification_scalers_analyzed_h5(std::string path, data_struct::Params_Override *override_values);

    bool generate_avg(std::string avg_filename, std::vector<std::string> files_to_avg);

    bool generate_stream_dataset(std::string dataset_directory,
                                 std::string dataset_name,
                                 int detector_num,
                                 size_t height,
                                 size_t width);

    bool save_stream_row(size_t d_hash,
                         size_t detector_num,
                         size_t row,
                         std::vector< data_struct::Spectra* >  *spectra_row);


    bool save_itegrade_spectra(data_struct::Spectra * spectra);

    bool close_dataset(size_t d_hash);

    bool start_save_seq(const std::string filename, bool force_new_file=false);

    bool start_save_seq(bool force_new_file=false){ return start_save_seq(_cur_filename, force_new_file);}

    void set_filename(std::string fname) {_cur_filename = fname;}

    bool save_spectra_volume(const std::string path,
                             data_struct::Spectra_Volume * spectra_volume,
                             real_t energy_offset,
                             real_t energy_slope,
                             real_t energy_quad,
                             size_t row_idx_start=0,
                             int row_idx_end=-1,
                             size_t col_idx_start=0,
                             int col_idx_end=-1);

    bool save_element_fits(const std::string path,
                           const data_struct::Fit_Count_Dict * const element_counts,
                           size_t row_idx_start=0,
                           int row_idx_end=-1,
                           size_t col_idx_start=0,
                           int col_idx_end=-1);

    void save_quantifications(std::map<string, data_struct::Quantification_Standard*> &quants);

    bool save_quantification(data_struct::Quantification_Standard * quantification_standard);

    bool save_scan_scalers(size_t detector_num,
                           struct mda_file *mda_scalers,
                           data_struct::Spectra_Volume * spectra_volume,
                           data_struct::Params_Override * params_override,
                           bool hasNetcdf,
                           size_t row_idx_start=0,
                           int row_idx_end=-1,
                           size_t col_idx_start=0,
                           int col_idx_end=-1);

    bool save_scan_scalers_confocal(std::string path,
                                    size_t detector_num,
                                    size_t row_idx_start=0,
                                    int row_idx_end=-1,
                                    size_t col_idx_start=0,
                                    int col_idx_end=-1);

	bool save_scan_scalers_gsecars(std::string path,
								size_t detector_num,
								size_t row_idx_start = 0,
								int row_idx_end = -1,
								size_t col_idx_start = 0,
								int col_idx_end = -1);

	// Add links to dataset and set version to 9 so legacy software can load it
    void add_v9_layout(std::string dataset_directory,
                       std::string dataset_file,
                       const std::vector<size_t>& detector_num_arr);

	// Add exchange layout to be loadable by external software
	void add_exchange_layout(std::string dataset_directory,
							std::string dataset_file,
							const std::vector<size_t>& detector_num_arr);

    bool end_save_seq(bool loginfo=true);

private:

    HDF5_IO();

    static HDF5_IO *_this_inst;

    static std::mutex _mutex;

	bool _load_integrated_spectra_analyzed_h5(hid_t file_id, data_struct::Spectra* spectra);

    bool _save_scan_meta_data(hid_t scan_grp_id, struct mda_file *mda_scalers, data_struct::Params_Override * params_override);
	bool _save_extras(hid_t scan_grp_id, struct mda_file *mda_scalers);
    bool _save_scalers(hid_t maps_grp_id, struct mda_file *mda_scalers, data_struct::Spectra_Volume * spectra_volume, data_struct::Params_Override * params_override, bool hasNetcdf);
    void _save_amps(hid_t scalers_grp_id, struct mda_file *mda_scalers, data_struct::Params_Override * params_override);
	bool _save_params_override(hid_t group_id, data_struct::Params_Override * params_override);

    void _gen_average(std::string full_hdf5_path, std::string dataset_name, hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids, bool avg=true);
    void _generate_avg_analysis(hid_t src_maps_grp_id, hid_t dst_maps_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids);
    void _generate_avg_integrated_spectra(hid_t src_analyzed_grp_id, hid_t dst_fit_grp_id, std::string group_name, hid_t ocpypl_id, std::vector<hid_t> &hdf5_file_ids);

    void _add_v9_layout(std::string dataset_file);
    void _add_v9_quant(hid_t file_id, hid_t quant_space, hid_t chan_names, hid_t chan_space, int chan_amt, std::string quant_str, std::string new_loc);
    void _add_extra_pvs(hid_t file_id, std::string group_name);

    bool _add_exchange_meta(hid_t file_id, std::string exchange_idx, std::string fits_link, std::string normalize_scaler);
	void _add_exchange_layout(std::string dataset_file);

    bool _open_h5_object(hid_t &id, H5_OBJECTS obj, std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map, std::string s1, hid_t id2, bool log_error=true, bool close_on_fail=true);
    void _close_h5_objects(std::stack<std::pair<hid_t, H5_OBJECTS> > &close_map);

    struct scaler_struct
    {
        scaler_struct(std::string name, std::string units, int mda_idx_, int hdf_idx_, bool normalize_by_time_)
        {
             hdf_name = name;
			 hdf_units = units;
             mda_idx = mda_idx_;
             hdf_idx = hdf_idx_;
             normalize_by_time = normalize_by_time_;
        }
        int mda_idx;
        int hdf_idx;
        std::string hdf_name;
		std::string hdf_units;
        bool normalize_by_time;
    };

    hid_t _cur_file_id;
    std::string _cur_filename;

};

}// end namespace file
}// end namespace io

#endif // HDF5_IO_H
