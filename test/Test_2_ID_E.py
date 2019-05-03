
#Don't forget to append XRF-Maps/bin directory to PYTHONPATH

import pyxrfmaps as px

element_csv_filename = "../reference/xrf_library.csv"
element_henke_filename = "../reference/henke.xdr"

def load_spectra_vol_and_print():
	dataset_dir = '2_ID_E_dataset/'
	dataset_file = 'axo_std.mda'
	detector_num = 0
	sv = px.Spectra_Volume()
	po = px.ParamsOverride()
	if False == px.load_override_params(dataset_dir, detector_num, po):
		# if we can't find maps_fit_parameters_override.txt[detector_num] then try to load the averaged maps_fit_parameters_override.txt
		print(px.load_override_params(dataset_dir, -1, po))
	print (px.load_spectra_volume(dataset_dir, dataset_file, detector_num, sv, po, False, False))
	spectra = sv[0][0]
	for i in spectra:
		print (i)


def run_analysis():
    px.load_element_info(element_henke_filename, element_csv_filename)
    job = px.AnalysisJob()
    job.dataset_directory = '2_ID_E_dataset/'
    job.dataset_files = px.find_all_dataset_files(job.dataset_directory + "mda/", ".mda")
    job.optimize_dataset_files = job.dataset_files
    job.detector_num_start = 0
    job.detector_num_end = 3
    job.quick_and_dirty = False
    job.fitting_routines = [px.FittingRoutines.ROI, px.FittingRoutines.SVD, px.FittingRoutines.MATRIX, px.FittingRoutines.NNLS]
    job.quantification_standard_filename = 'maps_standardinfo.txt'
    px.check_and_create_dirs(job.dataset_directory)
    px.populate_netcdf_hdf5_files(job.dataset_directory)
    px.init_analysis_job_detectors(job)
    px.generate_optimal_params(job)
    px.perform_quantification(job)
    px.process_dataset_files(job)
    print('done')

if __name__ == '__main__':
	run_analysis()
