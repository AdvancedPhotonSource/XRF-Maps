#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/chrono.h>
#include <pybind11/functional.h>
#include <pybind11/eigen.h>

#if defined _WIN32 || defined __CYGWIN__
  #ifdef _BUILD_WITH_ZMQ
    #include "winsock2.h"
  #endif
#endif

#include "fitting/models/gaussian_model.h"
#include "fitting/optimizers/lmfit_optimizer.h"
#include "fitting/optimizers/mpfit_optimizer.h"
#include "fitting/routines/roi_fit_routine.h"
#include "fitting/routines/svd_fit_routine.h"
#include "fitting/routines/param_optimized_fit_routine.h"
#include "fitting/routines/matrix_optimized_fit_routine.h"
#include "fitting/routines/nnls_fit_routine.h"
#include "data_struct/analysis_job.h"
#include "data_struct/stream_block.h"
#include "io/net/basic_serializer.h"

#include "io/file/hl_file_io.h"

//#include "core/process_streaming.h"
#include "core/process_whole.h"
#include "workflow/xrf/spectra_net_streamer.h"
#include "workflow/xrf/spectra_file_source.h"
#include "workflow/xrf/detector_sum_spectra_source.h"
//#include "workflow/xrf/spectra_net_source.h"

namespace py = pybind11;

//PYBIND11_MAKE_OPAQUE(std::vector<int>);

PYBIND11_MODULE(pyxrfmaps, m) {
    m.doc() = R"pbdoc(
        PyXrfMaps
        -----------------------

        .. currentmodule:: pybindings

        .. autosummary::
           :toctree: _generate

    )pbdoc";

    m.attr("have_eigen") = true;
    //sub modules
    //sub module fitting
    py::module fit = m.def_submodule("fitting", "Fitting submodule");
    py::module fm = fit.def_submodule("models", "Fitting models submodule");
    py::module fo = fit.def_submodule("optimizers", "Fitting optimizers submodule");
    py::module fr = fit.def_submodule("routines", "Fitting routines submodule");
    //sub module workflow
    py::module workflow = m.def_submodule("workflow", "Workflow submodule");
    //sub module io
    py::module io = m.def_submodule("io", "IO submodule");
    py::module io_net = io.def_submodule("net", "network submodule");
    py::module io_file = io.def_submodule("file", "file submodule");

    py::enum_<data_struct::Fitting_Routines>(m, "FittingRoutines")
    .value("ROI", data_struct::Fitting_Routines::ROI)
    .value("GAUSS_TAILS", data_struct::Fitting_Routines::GAUSS_TAILS)
    .value("MATRIX", data_struct::Fitting_Routines::GAUSS_MATRIX)
    .value("SVD", data_struct::Fitting_Routines::SVD)
    .value("NNLS", data_struct::Fitting_Routines::NNLS);


    //data structures
    /*
    py::class_<data_struct::Spectra, data_struct::ArrayXr>(m, "Spectra", py::buffer_protocol())
    //py::class_<data_struct::Spectra, data_struct::ArrayXr>(m, "Spectra")
        .def(py::init<size_t>())
        .def("add", &data_struct::Spectra::add)
        .def("recalc_elapsed_livetime", &data_struct::Spectra::recalc_elapsed_livetime)
        .def("set_elapsed_livetime", (void (data_struct::Spectra::*)(real_t)) &data_struct::Spectra::elapsed_livetime )
        .def("get_elapsed_livetime", (const real_t (data_struct::Spectra::*)() const) &data_struct::Spectra::elapsed_livetime)
        .def("set_elapsed_realtime", (void (data_struct::Spectra::*)(real_t)) &data_struct::Spectra::elapsed_realtime )
        .def("get_elapsed_realtime", (const real_t (data_struct::Spectra::*)() const) &data_struct::Spectra::elapsed_realtime)
        .def("set_input_counts", (void (data_struct::Spectra::*)(real_t)) &data_struct::Spectra::input_counts )
        .def("get_input_counts", (const real_t (data_struct::Spectra::*)() const) &data_struct::Spectra::input_counts)
        .def("set_output_counts", (void (data_struct::Spectra::*)(real_t)) &data_struct::Spectra::output_counts )
        .def("get_output_counts", (const real_t (data_struct::Spectra::*)() const) &data_struct::Spectra::output_counts)
        .def("sub_spectra", &data_struct::Spectra::sub_spectra);

        .def_buffer([](data_struct::Spectra &m) -> py::buffer_info {
                return py::buffer_info(
                    m.data(),                               // Pointer to buffer
                    sizeof(real_t),                          // Size of one scalar
                    py::format_descriptor<real_t>::format(), // Python struct-style format descriptor
                    1,                                      // Number of dimensions
                    { m.size() },                 // Buffer dimensions
                    { sizeof(real_t) }             // Strides (in bytes) for each index
                );
            });
*/

    py::class_<data_struct::Spectra_Line>(m, "Spectra_Line", py::buffer_protocol())
        .def(py::init<>())
        .def("__getitem__", [](const data_struct::Spectra_Line &s, size_t i) {
        if (i >= s.size()) throw py::index_error();
        return s[i];
        })
        .def("resize_and_zero", &data_struct::Spectra_Line::resize_and_zero)
        .def("alloc_row_size", &data_struct::Spectra_Line::alloc_row_size)
        .def("recalc_elapsed_livetime", &data_struct::Spectra_Line::recalc_elapsed_livetime)
        .def("size", &data_struct::Spectra_Line::size);

    py::class_<data_struct::Spectra_Volume>(m, "Spectra_Volume", py::buffer_protocol())
        .def(py::init<>())
        .def("__getitem__", [](const data_struct::Spectra_Volume &s, size_t i) {
        if (i >= s.rows()) throw py::index_error();
        return s[i];
        })
        .def("resize_and_zero", &data_struct::Spectra_Volume::resize_and_zero)
        .def("integrate", &data_struct::Spectra_Volume::integrate)
        .def("generate_scaler_maps", &data_struct::Spectra_Volume::generate_scaler_maps)
        .def("cols", &data_struct::Spectra_Volume::cols)
        .def("rows", &data_struct::Spectra_Volume::rows)
        .def("recalc_elapsed_livetime", &data_struct::Spectra_Volume::recalc_elapsed_livetime)
        .def("samples_size", &data_struct::Spectra_Volume::samples_size)
        .def("rank", &data_struct::Spectra_Volume::rank);

    py::class_<data_struct::Element_Info>(m, "ElementInfo")
    .def(py::init<>())
    .def("init_f_energies", &data_struct::Element_Info::init_f_energies)
    .def("init_extra_energies", &data_struct::Element_Info::init_extra_energies)
    .def("get_energies_between", &data_struct::Element_Info::get_energies_between)
    .def("calc_beta", &data_struct::Element_Info::calc_beta)
    .def_readwrite("number", &data_struct::Element_Info::number)
    .def_readwrite("name", &data_struct::Element_Info::name)
    .def_readwrite("density", &data_struct::Element_Info::density)
    .def_readwrite("mass", &data_struct::Element_Info::mass)
    .def_readwrite("xrf", &data_struct::Element_Info::xrf)
    .def_readwrite("xrf_abs_yield", &data_struct::Element_Info::xrf_abs_yield)
    .def_readwrite("yieldD", &data_struct::Element_Info::yieldD)
    .def_readwrite("bindingE", &data_struct::Element_Info::bindingE)
    .def_readwrite("jump", &data_struct::Element_Info::jump)
    .def_readwrite("f1_atomic_scattering_real", &data_struct::Element_Info::f1_atomic_scattering_real)
    .def_readwrite("f2_atomic_scattering_imaginary", &data_struct::Element_Info::f2_atomic_scattering_imaginary)
    .def_readwrite("energies", &data_struct::Element_Info::energies)
    .def_readwrite("extra_energies", &data_struct::Element_Info::extra_energies)
    .def_readwrite("extra_f1", &data_struct::Element_Info::extra_f1)
    .def_readwrite("extra_f2", &data_struct::Element_Info::extra_f2);

    py::class_<data_struct::Element_Info_Map>(m, "ElementInfoMap")
    .def("inst", &data_struct::Element_Info_Map::inst)
    .def("clear", &data_struct::Element_Info_Map::clear)
    .def("generate_default_elements", &data_struct::Element_Info_Map::generate_default_elements)
    .def("add_element", &data_struct::Element_Info_Map::add_element)
    .def("calc_beta", &data_struct::Element_Info_Map::calc_beta)
    .def("get_element", (data_struct::Element_Info* (data_struct::Element_Info_Map::*)(int)) &data_struct::Element_Info_Map::get_element)
    .def("get_element", (data_struct::Element_Info* (data_struct::Element_Info_Map::*)(std::string)) &data_struct::Element_Info_Map::get_element)
    .def("contains", &data_struct::Element_Info_Map::contains)
    .def_readwrite("_energies", &data_struct::Element_Info_Map::_energies);

    py::class_<data_struct::Element_Quant>(m, "ElementQuant")
    .def(py::init<>())
    .def_readwrite("weight", &data_struct::Element_Quant::weight)
    .def_readwrite("absorption", &data_struct::Element_Quant::absorption)
    .def_readwrite("transmission_Be", &data_struct::Element_Quant::transmission_Be)
    .def_readwrite("transmission_Ge", &data_struct::Element_Quant::transmission_Ge)
    .def_readwrite("yield", &data_struct::Element_Quant::yield)
    .def_readwrite("transmission_through_Si_detector", &data_struct::Element_Quant::transmission_through_Si_detector)
    .def_readwrite("transmission_through_air", &data_struct::Element_Quant::transmission_through_air)
    .def_readwrite("e_cal_ratio", &data_struct::Element_Quant::e_cal_ratio);

    py::class_<data_struct::Params_Override>(m, "ParamsOverride")
    .def(py::init<>())
    .def_readwrite("dataset_directory", &data_struct::Params_Override::dataset_directory)
    .def_readwrite("detector_num", &data_struct::Params_Override::detector_num)
    .def_readwrite("fit_params", &data_struct::Params_Override::fit_params)
    .def_readwrite("elements_to_fit", &data_struct::Params_Override::elements_to_fit)
    .def_readwrite("detector_element", &data_struct::Params_Override::detector_element)
    .def_readwrite("si_escape_factor", &data_struct::Params_Override::si_escape_factor)
    .def_readwrite("ge_escape_factor", &data_struct::Params_Override::ge_escape_factor)
    .def_readwrite("si_escape_enabled", &data_struct::Params_Override::si_escape_enabled)
    .def_readwrite("ge_escape_enabled", &data_struct::Params_Override::ge_escape_enabled)
    .def_readwrite("fit_snip_width", &data_struct::Params_Override::fit_snip_width)
    .def_readwrite("be_window_thickness", &data_struct::Params_Override::be_window_thickness)
    .def_readwrite("det_chip_thickness", &data_struct::Params_Override::det_chip_thickness)
    .def_readwrite("ge_dead_layer", &data_struct::Params_Override::ge_dead_layer)
    .def_readwrite("us_amp_sens_num_pv", &data_struct::Params_Override::us_amp_sens_num_pv)
    .def_readwrite("us_amp_sens_unit_pv", &data_struct::Params_Override::us_amp_sens_unit_pv)
    .def_readwrite("ds_amp_sens_num_pv", &data_struct::Params_Override::ds_amp_sens_num_pv)
    .def_readwrite("ds_amp_sens_unit_pv", &data_struct::Params_Override::ds_amp_sens_unit_pv)
    .def_readwrite("us_amp_sens_num", &data_struct::Params_Override::us_amp_sens_num)
    .def_readwrite("us_amp_sens_unit", &data_struct::Params_Override::us_amp_sens_unit)
    .def_readwrite("ds_amp_sens_num", &data_struct::Params_Override::ds_amp_sens_num)
    .def_readwrite("ds_amp_sens_unit", &data_struct::Params_Override::ds_amp_sens_unit);
    /*
    py::class_<data_struct::Calibration_Curve>(m, "CalibrationCurve")
    .def(py::init<>())
    .def_readwrite("quantifier_name", &data_struct::Calibration_Curve::quantifier_name)
    .def_readwrite("shell_curves", &data_struct::Calibration_Curve::shell_curves)
    .def_readwrite("shell_curves_labels", &data_struct::Calibration_Curve::shell_curves_labels);

    py::class_<data_struct::Quantifiers>(m, "Quantifiers")
    .def(py::init<>())
    .def("resize", &data_struct::Quantifiers::resize)
    .def_readwrite("calib_curves", &data_struct::Quantifiers::calib_curves);

    py::class_<data_struct::Quantification_Standard>(m, "QuantificationStandard")
    .def(py::init<>())
    .def("append_element", &data_struct::Quantification_Standard::append_element)
    .def("processed", &data_struct::Quantification_Standard::processed)
    .def("init_element_quants", &data_struct::Quantification_Standard::init_element_quants)
    .def("generate_calibration_curve", &data_struct::Quantification_Standard::generate_calibration_curve)
    .def_readwrite("quantifier_map", &data_struct::Quantification_Standard::quantifier_map)
    .def_readwrite("element_quants", &data_struct::Quantification_Standard::element_quants)
    .def_readwrite("element_counts", &data_struct::Quantification_Standard::element_counts)
    .def_readwrite("integrated_spectra", &data_struct::Quantification_Standard::integrated_spectra)
    .def_readwrite("standard_filename", &data_struct::Quantification_Standard::standard_filename)
    .def_readwrite("sr_current", &data_struct::Quantification_Standard::sr_current)
    .def_readwrite("US_IC", &data_struct::Quantification_Standard::US_IC)
    .def_readwrite("DS_IC", &data_struct::Quantification_Standard::DS_IC)
    .def_readwrite("beryllium_window_thickness", &data_struct::Quantification_Standard::beryllium_window_thickness)
    .def_readwrite("germanium_dead_layer", &data_struct::Quantification_Standard::germanium_dead_layer)
    .def_readwrite("detector_chip_thickness", &data_struct::Quantification_Standard::detector_chip_thickness)
    .def_readwrite("incident_energy", &data_struct::Quantification_Standard::incident_energy)
    .def_readwrite("airpath", &data_struct::Quantification_Standard::airpath)
    .def_readwrite("detector_element", &data_struct::Quantification_Standard::detector_element);
    */
    py::class_<data_struct::Stream_Fitting_Block>(m, "StreamFittingBlock")
    .def(py::init<>())
    .def_readwrite("fit_routine", &data_struct::Stream_Fitting_Block::fit_routine)
    .def_readwrite("fit_counts", &data_struct::Stream_Fitting_Block::fit_counts);

    py::class_<data_struct::Stream_Block>(m, "StreamBlock")
    .def(py::init<>())
    .def("init_fitting_blocks", &data_struct::Stream_Block::init_fitting_blocks)
    .def("row", &data_struct::Stream_Block::row)
    .def("col", &data_struct::Stream_Block::col)
    .def("height", &data_struct::Stream_Block::height)
    .def("width", &data_struct::Stream_Block::width)
    .def("is_end_of_row", &data_struct::Stream_Block::is_end_of_row)
    .def("dataset_hash", &data_struct::Stream_Block::dataset_hash)
    .def("detector_number", &data_struct::Stream_Block::detector_number)
    .def("is_end_block", &data_struct::Stream_Block::is_end_block)
    .def_readwrite("fitting_blocks", &data_struct::Stream_Block::fitting_blocks)
    .def_readwrite("dataset_directory", &data_struct::Stream_Block::dataset_directory)
    .def_readwrite("dataset_name", &data_struct::Stream_Block::dataset_name)
    .def_readwrite("spectra", &data_struct::Stream_Block::spectra)
    .def_readwrite("elements_to_fit", &data_struct::Stream_Block::elements_to_fit)
    .def_readwrite("optimize_fit_params_preset", &data_struct::Stream_Block::optimize_fit_params_preset)
    .def_readwrite("model", &data_struct::Stream_Block::model)
    .def_readwrite("theta", &data_struct::Stream_Block::theta);
    

    py::class_<data_struct::Detector>(m, "Detector")
    .def(py::init<>())
    .def_readwrite("fit_routines", &data_struct::Detector::fit_routines)
    .def_readwrite("model", &data_struct::Detector::model)
    .def_readwrite("quant_standards", &data_struct::Detector::quantification_standards)
    .def_readwrite("fit_params_override_dict", &data_struct::Detector::fit_params_override_dict);

    py::class_<data_struct::Analysis_Job>(m, "AnalysisJob")
    .def(py::init<>())
    .def("get_first_detector", &data_struct::Analysis_Job::get_first_detector)
    .def("get_detector", &data_struct::Analysis_Job::get_detector)
    .def("set_optimizer", &data_struct::Analysis_Job::set_optimizer)
    .def("get_optimizer", &data_struct::Analysis_Job::optimizer)
    .def("init_fit_routines", &data_struct::Analysis_Job::init_fit_routines)
    .def_readwrite("command_line", &data_struct::Analysis_Job::command_line)
    .def_readwrite("dataset_directory", &data_struct::Analysis_Job::dataset_directory)
    .def_readwrite("quantification_standard_filename", &data_struct::Analysis_Job::quantification_standard_filename)
    .def_readwrite("theta_pv", &data_struct::Analysis_Job::theta_pv)
    .def_readwrite("dataset_files", &data_struct::Analysis_Job::dataset_files)
    .def_readwrite("optimize_dataset_files", &data_struct::Analysis_Job::optimize_dataset_files)
    .def_readwrite("fitting_routines", &data_struct::Analysis_Job::fitting_routines)
    .def_readwrite("detectors_meta_data", &data_struct::Analysis_Job::detectors_meta_data)
    .def_readwrite("optimize_fit_params_preset", &data_struct::Analysis_Job::optimize_fit_params_preset)
    .def_readwrite("detector_num_arr", &data_struct::Analysis_Job::detector_num_arr)
    .def_readwrite("num_threads", &data_struct::Analysis_Job::num_threads)
    .def_readwrite("quick_and_dirty", &data_struct::Analysis_Job::quick_and_dirty)
    .def_readwrite("generate_average_h5", &data_struct::Analysis_Job::generate_average_h5)
    .def_readwrite("is_network_source", &data_struct::Analysis_Job::is_network_source)
    .def_readwrite("stream_over_network", &data_struct::Analysis_Job::stream_over_network);

    //fitting models
    py::class_<fitting::models::Gaussian_Model>(fm, "GaussModel")
    .def(py::init<>())
    .def("model_spectrum", &fitting::models::Gaussian_Model::model_spectrum)
    .def("model_spectrum_element", &fitting::models::Gaussian_Model::model_spectrum_element);


    //fitting optimizers
    py::class_<fitting::optimizers::LMFit_Optimizer>(fo, "lmfit")
    .def(py::init<>())
    .def("minimize", &fitting::optimizers::LMFit_Optimizer::minimize)
    .def("minimize_func", &fitting::optimizers::LMFit_Optimizer::minimize_func)
    .def("minimize_quantification", &fitting::optimizers::LMFit_Optimizer::minimize_quantification);

    py::class_<fitting::optimizers::MPFit_Optimizer>(fo, "mpfit")
    .def(py::init<>())
    .def("minimize", &fitting::optimizers::MPFit_Optimizer::minimize)
    .def("minimize_func", &fitting::optimizers::MPFit_Optimizer::minimize_func)
    .def("minimize_quantification", &fitting::optimizers::MPFit_Optimizer::minimize_quantification);


    //routines
    py::class_<fitting::routines::ROI_Fit_Routine>(fr, "roi")
    .def(py::init<>())
    .def("fit_spectra", &fitting::routines::ROI_Fit_Routine::fit_spectra)
    .def("get_name", &fitting::routines::ROI_Fit_Routine::get_name)
    .def("initialize", &fitting::routines::ROI_Fit_Routine::initialize);

    py::class_<fitting::routines::SVD_Fit_Routine>(fr, "svd")
    .def(py::init<>())
    .def("fit_spectra", &fitting::routines::SVD_Fit_Routine::fit_spectra)
    .def("get_name", &fitting::routines::SVD_Fit_Routine::get_name)
    .def("initialize", &fitting::routines::SVD_Fit_Routine::initialize);

    py::class_<fitting::routines::Matrix_Optimized_Fit_Routine>(fr, "matrix")
    .def(py::init<>())
    .def("fit_spectra", &fitting::routines::Matrix_Optimized_Fit_Routine::fit_spectra)
    .def("get_name", &fitting::routines::Matrix_Optimized_Fit_Routine::get_name)
    .def("initialize", &fitting::routines::Matrix_Optimized_Fit_Routine::initialize);

    py::class_<fitting::routines::NNLS_Fit_Routine>(fr, "nnls")
    .def(py::init<>())
    .def("fit_spectra", &fitting::routines::NNLS_Fit_Routine::fit_spectra)
    .def("get_name", &fitting::routines::NNLS_Fit_Routine::get_name)
    .def("initialize", &fitting::routines::NNLS_Fit_Routine::initialize);

    py::class_<fitting::routines::Param_Optimized_Fit_Routine>(fr, "param_optimized")
    .def(py::init<>())
    .def("fit_spectra", &fitting::routines::Param_Optimized_Fit_Routine::fit_spectra)
    .def("fit_spectra_parameters", &fitting::routines::Param_Optimized_Fit_Routine::fit_spectra_parameters)
    .def("get_name", &fitting::routines::Param_Optimized_Fit_Routine::get_name)
    .def("initialize", &fitting::routines::Param_Optimized_Fit_Routine::initialize)
    .def("set_optimizer", &fitting::routines::Param_Optimized_Fit_Routine::set_optimizer)
    .def("set_update_coherent_amplitude_on_fit", &fitting::routines::Param_Optimized_Fit_Routine::set_update_coherent_amplitude_on_fit);


    ////// IO //////
    //hl_file_io
    m.def("check_and_create_dirs", &io::check_and_create_dirs);
    m.def("compare_file_size", &io::compare_file_size);
    m.def("find_all_dataset_files", &io::find_all_dataset_files);
    m.def("generate_h5_averages", &io::generate_h5_averages);
    m.def("generate_fit_routine", &io::generate_fit_routine);
    m.def("init_analysis_job_detectors", &io::init_analysis_job_detectors);
    m.def("load_element_info", &io::load_element_info);
    m.def("load_and_integrate_spectra_volume", &io::load_and_integrate_spectra_volume);
    m.def("load_override_params", &io::load_override_params);
  ///  m.def("load_quantification_standard", &io::load_quantification_standard);
    m.def("load_spectra_volume", &io::load_spectra_volume);
    m.def("populate_netcdf_hdf5_files", &io::populate_netcdf_hdf5_files);
   // m.def("save_averaged_fit_params", &io::save_averaged_fit_params);
    m.def("save_optimized_fit_params", &io::save_optimized_fit_params);
//    m.def("save_volume", &io::save_volume);
    m.def("sort_dataset_files_by_size", &io::sort_dataset_files_by_size);

    // IO NET
    //basic serializer
    py::class_<io::net::Basic_Serializer>(io_net, "BasicSerializer")
    .def(py::init<>())
    .def("encode_counts", &io::net::Basic_Serializer::encode_counts)
    .def("decode_counts", &io::net::Basic_Serializer::decode_counts)
    .def("encode_spectra", &io::net::Basic_Serializer::encode_spectra)
    .def("decode_spectra", &io::net::Basic_Serializer::decode_spectra);

    // IO FILE
    //mda_io
    py::class_<io::file::MDA_IO>(io_file, "MDA_IO")
    .def(py::init<>())
    .def("unload", &io::file::MDA_IO::unload)
    .def("load_spectra_volume", &io::file::MDA_IO::load_spectra_volume)
    .def("load_spectra_volume_with_callback", &io::file::MDA_IO::load_spectra_volume_with_callback)
    //.def("find_scaler_index", &io::file::MDA_IO::find_scaler_index)
    .def("get_multiplied_dims", &io::file::mda_get_multiplied_dims)
    .def("get_rank_and_dims", &io::file::mda_get_rank_and_dims);

    //NetCDF_IO
    io_file.def("netcdf_load_spectra_line", [](std::string path,
                size_t detector,
                data_struct::Spectra_Line* spec_line)
                {
                    return io::file::NetCDF_IO::inst()->load_spectra_line(path, detector, spec_line);
                });
    //NetCDF_IO
    io_file.def("netcdf_load_spectra_line_with_callback", [](std::string path,
                std::vector<size_t> detector_num_arr,
                int row,
                size_t max_rows,
                size_t max_cols,
                data_struct::IO_Callback_Func_Def callback_fun,
                void* user_data) { return io::file::NetCDF_IO::inst()->load_spectra_line_with_callback(path,
                                                                                  detector_num_arr,
                                                                                  row,
                                                                                  max_rows,
                                                                                  max_cols,
                                                                                  callback_fun,
                                                                                  user_data);
                                 });



    //workflow
    py::class_<workflow::Sink<data_struct::Stream_Block*> >(workflow, "StreamBlockSink")
    .def(py::init<>())
    //.def("connect", &workflow::Sink<data_struct::Stream_Block>::connect)
    .def("set_function", &workflow::Sink<data_struct::Stream_Block*>::set_function)
    .def("start", &workflow::Sink<data_struct::Stream_Block*>::start)
    .def("stop", &workflow::Sink<data_struct::Stream_Block*>::stop)
    .def("wait_and_stop", &workflow::Sink<data_struct::Stream_Block*>::wait_and_stop)
    .def("set_delete_block", &workflow::Sink<data_struct::Stream_Block*>::set_delete_block)
    .def("sink_function", &workflow::Sink<data_struct::Stream_Block*>::sink_function);
#ifdef _BUILD_WITH_ZMQ
    py::class_<workflow::xrf::Spectra_Net_Streamer, workflow::Sink<data_struct::Stream_Block*> >(workflow, "SpectraNetStreamer")
    .def(py::init<std::string>())
    .def("set_send_counts", &workflow::xrf::Spectra_Net_Streamer::set_send_counts)
    .def("set_send_spectra", &workflow::xrf::Spectra_Net_Streamer::set_send_spectra)
    .def("stream", &workflow::xrf::Spectra_Net_Streamer::stream);
#endif

    py::class_<workflow::xrf::Spectra_File_Source>(workflow, "SpectraFileSource")
    .def(py::init<>())
    .def("connect_sink", &workflow::Source<data_struct::Stream_Block*>::connect_sink)
    .def("set_init_fitting_routines", &workflow::xrf::Spectra_File_Source::set_init_fitting_routines)
    .def("load_netcdf_line", &workflow::xrf::Spectra_File_Source::load_netcdf_line)
    .def("run", &workflow::xrf::Spectra_File_Source::run);

    py::class_<workflow::xrf::Detector_Sum_Spectra_Source, workflow::xrf::Spectra_File_Source>(workflow, "DetectorSumSpectraFileSource")
    .def(py::init<>())
    .def("connect_sink", &workflow::Source<data_struct::Stream_Block*>::connect_sink)
    .def("set_init_fitting_routines", &workflow::xrf::Spectra_File_Source::set_init_fitting_routines)
    .def("load_netcdf_line", &workflow::xrf::Spectra_File_Source::load_netcdf_line)
    .def("run", &workflow::xrf::Spectra_File_Source::run);

    //process_streaming
//    m.def("proc_spectra_block", &proc_spectra_block);
//    m.def("run_stream_pipeline", &run_stream_pipeline);

    //process_whole
    //m.def("generate_fit_count_dict", &generate_fit_count_dict<real_t>);
    m.def("fit_single_spectra", &fit_single_spectra);
    m.def("optimize_integrated_fit_params", &optimize_integrated_fit_params);
    m.def("generate_optimal_params", &generate_optimal_params);
   // m.def("generate_optimal_params_mp", &generate_optimal_params_mp);
    m.def("proc_spectra", &proc_spectra);
    m.def("process_dataset_files", &process_dataset_files);
    m.def("perform_quantification", &perform_quantification);
    //m.def("average_quantification", &average_quantification);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}


