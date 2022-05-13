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

auto fit_counts(fitting::routines::Base_Fit_Routine<float>* fit_route,
	const fitting::models::Base_Model<float>* const model,
	const Spectra<float>* const spectra,
	const Fit_Element_Map_Dict<float>* const elements_to_fit)
{
	std::unordered_map<std::string, float> out_counts;
	fit_route->fit_spectra(model, spectra, elements_to_fit, out_counts);
	return out_counts;
}

auto fit_spectra(fitting::routines::Base_Fit_Routine<float>* fit_route,
	fitting::models::Base_Model<float>* const model,
	const Spectra<float>* const spectra,
	const Fit_Element_Map_Dict<float>* const elements_to_fit)
{
	std::unordered_map<std::string, float> out_counts;
	data_struct::Fit_Parameters<float> fit_params;
	fit_params.append_and_update(model->fit_parameters());
	data_struct::Range energy_range = data_struct::get_energy_range(spectra->size(), &fit_params);
	fit_route->fit_spectra(model, spectra, elements_to_fit, out_counts);
	for (auto& itr : out_counts)
	{
		float val = out_counts.at(itr.first);
		if (val != 0.0)
		{
			fit_params[itr.first] = data_struct::Fit_Param<float>(itr.first, log10(val));
		}
		else
		{
			fit_params[itr.first] = data_struct::Fit_Param<float>(itr.first, val);
		}
	}
	/*
	for (auto itr = fit_params.begin(); itr != fit_params.end(); itr++)
	{
		logI << itr->first << " " << itr->second.value << "\n";
	}
	*/
	//return model->model_spectrum_mp(&fit_params, elements_to_fit, energy_range);
	return model->model_spectrum(&fit_params, elements_to_fit, nullptr, energy_range);
}

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

    py::class_<data_struct::Spectra_Line<float>>(m, "Spectra_Line", py::buffer_protocol())
        .def(py::init<>())
        .def("__getitem__", [](const data_struct::Spectra_Line<float>&s, size_t i) {
        if (i >= s.size()) throw py::index_error();
        return s[i];
        })
        .def("resize_and_zero", &data_struct::Spectra_Line<float>::resize_and_zero)
        .def("alloc_row_size", &data_struct::Spectra_Line<float>::alloc_row_size)
        .def("recalc_elapsed_livetime", &data_struct::Spectra_Line<float>::recalc_elapsed_livetime)
        .def("size", &data_struct::Spectra_Line<float>::size);

    py::class_<data_struct::Spectra_Volume<float>>(m, "Spectra_Volume", py::buffer_protocol())
        .def(py::init<>())
        .def("__getitem__", [](const data_struct::Spectra_Volume<float>&s, size_t i) {
        if (i >= s.rows()) throw py::index_error();
        return s[i];
        })
        .def("resize_and_zero", &data_struct::Spectra_Volume<float>::resize_and_zero)
        .def("integrate", &data_struct::Spectra_Volume<float>::integrate)
        .def("generate_scaler_maps", &data_struct::Spectra_Volume<float>::generate_scaler_maps)
        .def("cols", &data_struct::Spectra_Volume<float>::cols)
        .def("rows", &data_struct::Spectra_Volume<float>::rows)
        .def("recalc_elapsed_livetime", &data_struct::Spectra_Volume<float>::recalc_elapsed_livetime)
        .def("samples_size", &data_struct::Spectra_Volume<float>::samples_size)
        .def("rank", &data_struct::Spectra_Volume<float>::rank);

    py::class_<data_struct::Element_Info<float>>(m, "ElementInfo")
    .def(py::init<>())
    .def("init_f_energies", &data_struct::Element_Info<float>::init_f_energies)
    .def("init_extra_energies", &data_struct::Element_Info<float>::init_extra_energies)
    .def("get_energies_between", &data_struct::Element_Info<float>::get_energies_between)
    .def("calc_beta", &data_struct::Element_Info<float>::calc_beta)
    .def_readwrite("number", &data_struct::Element_Info<float>::number)
    .def_readwrite("name", &data_struct::Element_Info<float>::name)
    .def_readwrite("density", &data_struct::Element_Info<float>::density)
    .def_readwrite("mass", &data_struct::Element_Info<float>::mass)
    .def_readwrite("xrf", &data_struct::Element_Info<float>::xrf)
    .def_readwrite("xrf_abs_yield", &data_struct::Element_Info<float>::xrf_abs_yield)
    .def_readwrite("yieldD", &data_struct::Element_Info<float>::yieldD)
    .def_readwrite("bindingE", &data_struct::Element_Info<float>::bindingE)
    .def_readwrite("jump", &data_struct::Element_Info<float>::jump)
    .def_readwrite("f1_atomic_scattering_real", &data_struct::Element_Info<float>::f1_atomic_scattering_real)
    .def_readwrite("f2_atomic_scattering_imaginary", &data_struct::Element_Info<float>::f2_atomic_scattering_imaginary)
    .def_readwrite("energies", &data_struct::Element_Info<float>::energies)
    .def_readwrite("extra_energies", &data_struct::Element_Info<float>::extra_energies)
    .def_readwrite("extra_f1", &data_struct::Element_Info<float>::extra_f1)
    .def_readwrite("extra_f2", &data_struct::Element_Info<float>::extra_f2);

    py::class_<data_struct::Element_Info_Map<float>>(m, "ElementInfoMap")
    .def("inst", &data_struct::Element_Info_Map<float>::inst)
    .def("clear", &data_struct::Element_Info_Map<float>::clear)
    .def("generate_default_elements", &data_struct::Element_Info_Map<float>::generate_default_elements)
    .def("add_element", &data_struct::Element_Info_Map<float>::add_element)
    .def("calc_beta", &data_struct::Element_Info_Map<float>::calc_beta)
    .def("get_element", (data_struct::Element_Info<float>* (data_struct::Element_Info_Map<float>::*)(int)) &data_struct::Element_Info_Map<float>::get_element)
    .def("get_element", (data_struct::Element_Info<float>* (data_struct::Element_Info_Map<float>::*)(std::string)) &data_struct::Element_Info_Map<float>::get_element)
    .def("contains", &data_struct::Element_Info_Map<float>::contains)
    .def_readwrite("_energies", &data_struct::Element_Info_Map<float>::_energies);

    py::class_<data_struct::Element_Quant<float>>(m, "ElementQuant")
    .def(py::init<>())
    .def_readwrite("weight", &data_struct::Element_Quant<float>::weight)
    .def_readwrite("absorption", &data_struct::Element_Quant<float>::absorption)
    .def_readwrite("transmission_Be", &data_struct::Element_Quant<float>::transmission_Be)
    .def_readwrite("transmission_Ge", &data_struct::Element_Quant<float>::transmission_Ge)
    .def_readwrite("yield", &data_struct::Element_Quant<float>::yield)
    .def_readwrite("transmission_through_Si_detector", &data_struct::Element_Quant<float>::transmission_through_Si_detector)
    .def_readwrite("transmission_through_air", &data_struct::Element_Quant<float>::transmission_through_air)
    .def_readwrite("e_cal_ratio", &data_struct::Element_Quant<float>::e_cal_ratio);

	py::class_<data_struct::Fit_Element_Map<float>>(m, "Fit_Element_Map")
	.def(py::init<std::string, Element_Info*>())
	.def("center", &data_struct::Fit_Element_Map<float>::center)
	.def("width", &data_struct::Fit_Element_Map<float>::width)
	.def("set_custom_multiply_ratio", &data_struct::Fit_Element_Map<float>::set_custom_multiply_ratio)
	.def("multiply_custom_multiply_ratio", &data_struct::Fit_Element_Map<float>::multiply_custom_multiply_ratio)
	.def("init_energy_ratio_for_detector_element", &data_struct::Fit_Element_Map<float>::init_energy_ratio_for_detector_element)
	.def("full_name", &data_struct::Fit_Element_Map<float>::full_name)
	.def("symbol", &data_struct::Fit_Element_Map<float>::symbol)
	.def("Z", &data_struct::Fit_Element_Map<float>::Z)
	.def("energy_ratios", &data_struct::Fit_Element_Map<float>::energy_ratios)
	.def("energy_ratio_multipliers", &data_struct::Fit_Element_Map<float>::energy_ratio_multipliers)
	.def("width_multi", &data_struct::Fit_Element_Map<float>::width_multi)
	.def("set_as_pileup", &data_struct::Fit_Element_Map<float>::set_as_pileup)
	.def("pileup_element", &data_struct::Fit_Element_Map<float>::pileup_element)
	.def("shell_type_as_string", &data_struct::Fit_Element_Map<float>::shell_type_as_string)
	.def("check_binding_energy", &data_struct::Fit_Element_Map<float>::check_binding_energy);

    py::class_<data_struct::Params_Override<float>>(m, "ParamsOverride")
    .def(py::init<>())
    .def_readwrite("dataset_directory", &data_struct::Params_Override<float>::dataset_directory)
    .def_readwrite("detector_num", &data_struct::Params_Override<float>::detector_num)
    .def_readwrite("fit_params", &data_struct::Params_Override<float>::fit_params)
    .def_readwrite("elements_to_fit", &data_struct::Params_Override<float>::elements_to_fit)
    .def_readwrite("detector_element", &data_struct::Params_Override<float>::detector_element)
    .def_readwrite("si_escape_factor", &data_struct::Params_Override<float>::si_escape_factor)
    .def_readwrite("ge_escape_factor", &data_struct::Params_Override<float>::ge_escape_factor)
    .def_readwrite("si_escape_enabled", &data_struct::Params_Override<float>::si_escape_enabled)
    .def_readwrite("ge_escape_enabled", &data_struct::Params_Override<float>::ge_escape_enabled)
    .def_readwrite("fit_snip_width", &data_struct::Params_Override<float>::fit_snip_width)
    .def_readwrite("be_window_thickness", &data_struct::Params_Override<float>::be_window_thickness)
    .def_readwrite("det_chip_thickness", &data_struct::Params_Override<float>::det_chip_thickness)
    .def_readwrite("ge_dead_layer", &data_struct::Params_Override<float>::ge_dead_layer)
    .def_readwrite("us_amp_sens_num", &data_struct::Params_Override<float>::us_amp_sens_num)
    .def_readwrite("us_amp_sens_unit", &data_struct::Params_Override<float>::us_amp_sens_unit)
    .def_readwrite("ds_amp_sens_num", &data_struct::Params_Override<float>::ds_amp_sens_num)
    .def_readwrite("ds_amp_sens_unit", &data_struct::Params_Override<float>::ds_amp_sens_unit);
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

	py::class_<data_struct::Range>(m, "Energy_Range")
	.def(py::init<>())
	.def("count", &data_struct::Range::count)
	.def_readwrite("min", &data_struct::Range::min)
	.def_readwrite("max", &data_struct::Range::max);
	
	py::class_<data_struct::Fit_Param>(m, "Fit_Param")
	.def(py::init<>())
	.def(py::init<std::string, real_t>())
	.def("bound_type_str", &data_struct::Fit_Param::bound_type_str)
	.def_readwrite("name", &data_struct::Fit_Param::name)
	.def_readwrite("min_val", &data_struct::Fit_Param::min_val)
	.def_readwrite("max_val", &data_struct::Fit_Param::max_val)
	.def_readwrite("value", &data_struct::Fit_Param::value)
	.def_readwrite("step_size", &data_struct::Fit_Param::step_size)
	.def_readwrite("bound_type", &data_struct::Fit_Param::bound_type)
	.def_readwrite("opt_array_index", &data_struct::Fit_Param::opt_array_index);

	py::class_<data_struct::Fit_Parameters>(m, "Fit_Parameters")
	.def(py::init<>())
	.def("add_parameter", &data_struct::Fit_Parameters::add_parameter)
	.def("append_and_update", &data_struct::Fit_Parameters::append_and_update)
	.def("begin", &data_struct::Fit_Parameters::begin)
	.def("end", &data_struct::Fit_Parameters::end)
	.def("sum_values", &data_struct::Fit_Parameters::sum_values)
	.def("divide_fit_values_by", &data_struct::Fit_Parameters::divide_fit_values_by)
	.def("contains", &data_struct::Fit_Parameters::contains)
	.def("to_array", &data_struct::Fit_Parameters::to_array)
	.def("names_to_array", &data_struct::Fit_Parameters::names_to_array)
	.def("from_array", (void (data_struct::Fit_Parameters::*)(std::vector<real_t>&))  &data_struct::Fit_Parameters::from_array)
	.def("from_array", (void (data_struct::Fit_Parameters::*)(const real_t*, size_t))  &data_struct::Fit_Parameters::from_array)
	.def("set_all_value", &data_struct::Fit_Parameters::set_all_value)
	.def("set_all", &data_struct::Fit_Parameters::set_all)
	.def("update_values", &data_struct::Fit_Parameters::update_values)
	.def("update_and_add_values", &data_struct::Fit_Parameters::update_and_add_values)
	.def("update_and_add_values_gt_zero", &data_struct::Fit_Parameters::update_and_add_values_gt_zero)
	.def("remove", (void (data_struct::Fit_Parameters::*)(Fit_Parameters*)) &data_struct::Fit_Parameters::remove)
	.def("remove", (void (data_struct::Fit_Parameters::*)(std::string)) &data_struct::Fit_Parameters::remove)
	.def("value", &data_struct::Fit_Parameters::value)
	.def("print", &data_struct::Fit_Parameters::print)
	.def("at", &data_struct::Fit_Parameters::at)
	.def("size", &data_struct::Fit_Parameters::size);


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
	py::class_<fitting::models::Base_Model>(fm, "BaseModel");

	py::class_<fitting::models::Gaussian_Model, fitting::models::Base_Model>(fm, "GaussModel")
	.def(py::init<>())
	.def("fit_parameters", &fitting::models::Gaussian_Model::fit_parameters)
	.def("model_spectrum", &fitting::models::Gaussian_Model::model_spectrum)	
	.def("model_spectrum_mp", &fitting::models::Gaussian_Model::model_spectrum_mp)
	.def("model_spectrum_element", &fitting::models::Gaussian_Model::model_spectrum_element)
	.def("set_fit_params_preset", &fitting::models::Gaussian_Model::set_fit_params_preset)
	.def("reset_to_default_fit_params", &fitting::models::Gaussian_Model::reset_to_default_fit_params)
	.def("update_fit_params_values", &fitting::models::Gaussian_Model::update_fit_params_values)
	.def("update_and_add_fit_params_values_gt_zero", &fitting::models::Gaussian_Model::update_and_add_fit_params_values_gt_zero);

    //fitting optimizers
	py::class_<fitting::optimizers::Optimizer>(fo, "Optimizer");

    py::class_<fitting::optimizers::LMFit_Optimizer, fitting::optimizers::Optimizer>(fo, "lmfit")
    .def(py::init<>())
    .def("minimize", &fitting::optimizers::LMFit_Optimizer::minimize)
    .def("minimize_func", &fitting::optimizers::LMFit_Optimizer::minimize_func)
    .def("minimize_quantification", &fitting::optimizers::LMFit_Optimizer::minimize_quantification);

    py::class_<fitting::optimizers::MPFit_Optimizer, fitting::optimizers::Optimizer>(fo, "mpfit")
    .def(py::init<>())
    .def("minimize", &fitting::optimizers::MPFit_Optimizer::minimize)
    .def("minimize_func", &fitting::optimizers::MPFit_Optimizer::minimize_func)
    .def("minimize_quantification", &fitting::optimizers::MPFit_Optimizer::minimize_quantification);


    //routines
	py::class_<fitting::routines::Base_Fit_Routine>(fr, "base_fit_route");

    py::class_<fitting::routines::ROI_Fit_Routine, fitting::routines::Base_Fit_Routine>(fr, "roi")
    .def(py::init<>())
	.def("fit_spectra", [](fitting::routines::ROI_Fit_Routine& self, 
		fitting::models::Base_Model* const model,
		const Spectra* const spectra,
		const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_spectra(&self, model, spectra, elements_to_fit);
	})
	.def("fit_counts", [](fitting::routines::ROI_Fit_Routine& self,
		fitting::models::Base_Model* const model,
		const Spectra* const spectra,
		const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_counts(&self, model, spectra, elements_to_fit);
	})
    .def("get_name", &fitting::routines::ROI_Fit_Routine::get_name)
    .def("initialize", &fitting::routines::ROI_Fit_Routine::initialize);

	py::class_<fitting::routines::Param_Optimized_Fit_Routine, fitting::routines::Base_Fit_Routine>(fr, "param_optimized")
		.def(py::init<>())
		.def("fit_spectra", [](fitting::routines::Param_Optimized_Fit_Routine& self,
			fitting::models::Base_Model* const model,
			const Spectra* const spectra,
			const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_spectra(&self, model, spectra, elements_to_fit);
	})
		.def("fit_counts", [](fitting::routines::Param_Optimized_Fit_Routine& self,
			const fitting::models::Base_Model* const model,
			const Spectra* const spectra,
			const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_counts(&self, model, spectra, elements_to_fit);
	})
		.def("fit_spectra_parameters", &fitting::routines::Param_Optimized_Fit_Routine::fit_spectra_parameters)
		.def("get_name", &fitting::routines::Param_Optimized_Fit_Routine::get_name)
		.def("initialize", &fitting::routines::Param_Optimized_Fit_Routine::initialize)
		.def("set_optimizer", &fitting::routines::Param_Optimized_Fit_Routine::set_optimizer)
		.def("set_update_coherent_amplitude_on_fit", &fitting::routines::Param_Optimized_Fit_Routine::set_update_coherent_amplitude_on_fit);


    py::class_<fitting::routines::Matrix_Optimized_Fit_Routine, fitting::routines::Param_Optimized_Fit_Routine, fitting::routines::Base_Fit_Routine>(fr, "matrix")
    .def(py::init<>())
	.def("fit_spectra", [](fitting::routines::Matrix_Optimized_Fit_Routine& self, 
		fitting::models::Base_Model* const model,
		const Spectra* const spectra,
		const Fit_Element_Map_Dict* const elements_to_fit)
	{
		/*
		std::unordered_map<std::string, real_t> out_counts;
		data_struct::Fit_Parameters fit_params = model->fit_parameters();
		data_struct::Range energy_range = data_struct::get_energy_range(spectra->size(), &fit_params);
		self.fit_spectra(model, spectra, elements_to_fit, out_counts);
		for (const auto& itr : out_counts)
		{
			fit_params.add_parameter(Fit_Param(itr.first, itr.second));
		}
		
		Spectra spec_model;
		spec_model.resize(energy_range.count());
		spec_model.setZero(energy_range.count());
		self.model_spectrum(&fit_params, &energy_range, &spec_model);
		return spec_model;
		*/
		return fit_spectra(&self, model, spectra, elements_to_fit);
	})
	.def("fit_counts", [](fitting::routines::Matrix_Optimized_Fit_Routine& self,
			const fitting::models::Base_Model* const model,
			const Spectra* const spectra,
			const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_counts(&self, model, spectra, elements_to_fit);
	})
    .def("get_name", &fitting::routines::Matrix_Optimized_Fit_Routine::get_name)
    .def("initialize", &fitting::routines::Matrix_Optimized_Fit_Routine::initialize);

    py::class_<fitting::routines::NNLS_Fit_Routine, fitting::routines::Matrix_Optimized_Fit_Routine, fitting::routines::Param_Optimized_Fit_Routine, fitting::routines::Base_Fit_Routine>(fr, "nnls")
    .def(py::init<>())
	.def("fit_spectra", [](fitting::routines::NNLS_Fit_Routine& self, 
		fitting::models::Base_Model* const model,
		const Spectra* const spectra,
		const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_spectra(&self, model, spectra, elements_to_fit);
	})
	.def("fit_counts", [](fitting::routines::NNLS_Fit_Routine& self, 
		fitting::models::Base_Model* const model,
		const Spectra* const spectra,
		const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_counts(&self, model, spectra, elements_to_fit);
	})
    .def("get_name", &fitting::routines::NNLS_Fit_Routine::get_name)
    .def("initialize", &fitting::routines::NNLS_Fit_Routine::initialize);

	py::class_<fitting::routines::SVD_Fit_Routine, fitting::routines::Matrix_Optimized_Fit_Routine, fitting::routines::Param_Optimized_Fit_Routine, fitting::routines::Base_Fit_Routine>(fr, "svd")
		.def(py::init<>())
		.def("fit_spectra", [](fitting::routines::SVD_Fit_Routine& self,
			fitting::models::Base_Model* const model,
			const Spectra* const spectra,
			const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_spectra(&self, model, spectra, elements_to_fit);
	})
		.def("fit_counts", [](fitting::routines::SVD_Fit_Routine& self,
			fitting::models::Base_Model* const model,
			const Spectra* const spectra,
			const Fit_Element_Map_Dict* const elements_to_fit)
	{
		return fit_counts(&self, model, spectra, elements_to_fit);
	})
		.def("get_name", &fitting::routines::SVD_Fit_Routine::get_name)
		.def("initialize", &fitting::routines::SVD_Fit_Routine::initialize);


	py::enum_<fitting::optimizers::OPTIMIZER_OUTCOME>(fr, "OPTIMIZER_OUTCOME")
	.value("FOUND_ZERO", fitting::optimizers::OPTIMIZER_OUTCOME::FOUND_ZERO)
	.value("CONVERGED", fitting::optimizers::OPTIMIZER_OUTCOME::CONVERGED)
	.value("TRAPPED", fitting::optimizers::OPTIMIZER_OUTCOME::TRAPPED)
	.value("EXHAUSTED", fitting::optimizers::OPTIMIZER_OUTCOME::EXHAUSTED)
	.value("FAILED", fitting::optimizers::OPTIMIZER_OUTCOME::FAILED)
	.value("CRASHED", fitting::optimizers::OPTIMIZER_OUTCOME::CRASHED)
	.value("EXPLODED", fitting::optimizers::OPTIMIZER_OUTCOME::EXPLODED)
	.value("STOPPED", fitting::optimizers::OPTIMIZER_OUTCOME::STOPPED)
	.value("FOUND_NAN", fitting::optimizers::OPTIMIZER_OUTCOME::FOUND_NAN)
	.value("F_TOL_LT_TOL", fitting::optimizers::OPTIMIZER_OUTCOME::F_TOL_LT_TOL)
	.value("X_TOL_LT_TOL", fitting::optimizers::OPTIMIZER_OUTCOME::X_TOL_LT_TOL)
	.value("G_TOL_LT_TOL", fitting::optimizers::OPTIMIZER_OUTCOME::G_TOL_LT_TOL);

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
   // m.def("load_override_params", &io::load_override_params);
	m.def("load_override_params", [](std::string dataset_directory,
									int detector_num,
									bool append_file_name)
	{
		data_struct::Params_Override po;
		io::load_override_params(dataset_directory, detector_num, &po, append_file_name);
		return po;
		
	});
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

	m.def("get_energy_range", (data_struct::Range (*)(size_t, Fit_Parameters*)) &data_struct::get_energy_range);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}


