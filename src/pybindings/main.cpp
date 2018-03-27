#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include <pybind11/eigen.h>

#include "fitting/models/gaussian_model.h"
#include "fitting/optimizers/lmfit_optimizer.h"
#include "fitting/optimizers/mpfit_optimizer.h"
#include "fitting/routines/roi_fit_routine.h"
#include "fitting/routines/svd_fit_routine.h"
#include "fitting/routines/param_optimized_fit_routine.h"
#include "fitting/routines/matrix_optimized_fit_routine.h"
#include "fitting/routines/nnls_fit_routine.h"


namespace py = pybind11;

PYBIND11_MODULE(pyxrfmaps, m) {
    m.doc() = R"pbdoc(
        PyXrfMaps
        -----------------------

        .. currentmodule:: pybindings

        .. autosummary::
           :toctree: _generate

    )pbdoc";

    m.attr("have_eigen") = true;

    py::module fit = m.def_submodule("fitting", "Fitting submodule");
    py::module fm = fit.def_submodule("models", "Fitting models submodule");
    py::module fo = fit.def_submodule("optimizers", "Fitting optimizers submodule");
    py::module fr = fit.def_submodule("routines", "Fitting routines submodule");

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


/*
    //spectra class
    py::class_<data_struct::Spectra, data_struct::ArrayXr>(m, "Spectra", py::buffer_protocol())
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
        .def("sub_spectra", &data_struct::Spectra::sub_spectra)
        .def_buffer([](data_struct::Spectra &m) -> py::buffer_info {
                return py::buffer_info(
                    m.data(),                               // Pointer to buffer
                    sizeof(real_t),                          // Size of one scalar
                    py::format_descriptor<real_t>::format(), // Python struct-style format descriptor
                    1,                                      // Number of dimensions
                    { m.cols() },                 // Buffer dimensions
                    { sizeof(real_t) }             // Strides (in bytes) for each index
                );
            });

*/
#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}


