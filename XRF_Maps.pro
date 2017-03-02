#-------------------------------------------------
#
# Project created by QtCreator 2015-11-22T16:42:59
#
#-------------------------------------------------

QT       += core

QT       -= gui

TARGET = XRF_Maps
#CONFIG   += console
CONFIG   -= app_bundle
CONFIG += c++14

TEMPLATE = app
CONFIG += debug
#CONFIG += release
CONFIG += console

QMAKE_EXT_CPP += .c

INCLUDEPATH += src/core
INCLUDEPATH += src/data_struct
INCLUDEPATH += src/data_struct/xrf
INCLUDEPATH += src/data_struct/xrf/aps
INCLUDEPATH += src/quantification/models
INCLUDEPATH += src/workflow
INCLUDEPATH += src/workflow/distributors
INCLUDEPATH += src/io/file
INCLUDEPATH += src/io/file/aps
INCLUDEPATH += src/fitting/models
INCLUDEPATH += src/fitting/routines
INCLUDEPATH += src/fitting/optimizers
INCLUDEPATH += src/support/mda_utils
INCLUDEPATH += src/support/minpack
INCLUDEPATH += src/support/cmpfit-1.3a
INCLUDEPATH += src/support/nnls
INCLUDEPATH += src/visual
INCLUDEPATH += /usr/include/python2.7

#define _REAL_FLOAT for float, _REAL_DOUBLE for double
#DEFINES += _REAL_FLOAT
DEFINES += _REAL_DOUBLE

DESTDIR = ./bin
MOC_DIR += ./generated
OBJECTS_DIR += ./generated

macx {
#VTK
INCLUDEPATH += /Users/aglowacki/lib/vtk_7/include/vtk-7.0
LIBS += -L/Users/aglowacki/lib/vtk_7/lib -lvtkalglib-7.0 -lvtkFiltersFlowPaths-7.0 -lvtkglew-7.0 -lvtkInteractionWidgets-7.0 -lvtkIOVideo-7.0 -lvtkRenderingLabel-7.0 -lvtkChartsCore-7.0 -lvtkFiltersGeneral-7.0 -lvtkGUISupportQt-7.0 -lvtkIOAMR-7.0 -lvtkIOXML-7.0 -lvtkRenderingLOD-7.0 -lvtkCommonColor-7.0 -lvtkFiltersGeneric-7.0 -lvtkGUISupportQtSQL-7.0 -lvtkIOCore-7.0 -lvtkIOXMLParser-7.0 -lvtkRenderingOpenGL2-7.0 -lvtkCommonComputationalGeometry-7.0 -lvtkFiltersGeometry-7.0 -lvtkhdf5-7.0 -lvtkIOEnSight-7.0 -lvtkjpeg-7.0 -lvtkRenderingQt-7.0 -lvtkCommonCore-7.0 -lvtkFiltersHybrid-7.0 -lvtkhdf5_hl-7.0 -lvtkIOExodus-7.0 -lvtkjsoncpp-7.0 -lvtkRenderingVolume-7.0 -lvtkCommonDataModel-7.0 -lvtkFiltersHyperTree-7.0 -lvtkImagingColor-7.0 -lvtkIOExport-7.0 -lvtklibxml2-7.0 -lvtkRenderingVolumeOpenGL2-7.0 -lvtkCommonExecutionModel-7.0 -lvtkFiltersImaging-7.0 -lvtkImagingCore-7.0 -lvtkIOGeometry-7.0 -lvtkmetaio-7.0 -lvtksqlite-7.0 -lvtkCommonMath-7.0 -lvtkFiltersModeling-7.0 -lvtkImagingFourier-7.0 -lvtkIOImage-7.0 -lvtkNetCDF-7.0 -lvtksys-7.0 -lvtkCommonMisc-7.0 -lvtkFiltersParallel-7.0 -lvtkImagingGeneral-7.0 -lvtkIOImport-7.0 -lvtkNetCDF_cxx-7.0 -lvtktiff-7.0 -lvtkCommonSystem-7.0 -lvtkFiltersParallelImaging-7.0 -lvtkImagingHybrid-7.0 -lvtkIOInfovis-7.0 -lvtkoggtheora-7.0 -lvtkverdict-7.0 -lvtkCommonTransforms-7.0 -lvtkFiltersProgrammable-7.0 -lvtkImagingMath-7.0 -lvtkIOLegacy-7.0 -lvtkParallelCore-7.0 -lvtkViewsContext2D-7.0 -lvtkDICOMParser-7.0 -lvtkFiltersSelection-7.0 -lvtkImagingMorphological-7.0 -lvtkIOLSDyna-7.0 -lvtkpng-7.0 -lvtkViewsCore-7.0 -lvtkDomainsChemistry-7.0 -lvtkFiltersSMP-7.0 -lvtkImagingSources-7.0 -lvtkIOMINC-7.0 -lvtkproj4-7.0 -lvtkViewsInfovis-7.0 -lvtkDomainsChemistryOpenGL2-7.0 -lvtkFiltersSources-7.0 -lvtkImagingStatistics-7.0 -lvtkIOMovie-7.0 -lvtkRenderingAnnotation-7.0 -lvtkViewsQt-7.0 -lvtkexoIIc-7.0 -lvtkFiltersStatistics-7.0 -lvtkImagingStencil-7.0 -lvtkIONetCDF-7.0 -lvtkRenderingContext2D-7.0 -lvtkzlib-7.0 -lvtkexpat-7.0 -lvtkFiltersTexture-7.0 -lvtkInfovisCore-7.0 -lvtkIOParallel-7.0 -lvtkRenderingContextOpenGL2-7.0 -lvtkFiltersAMR-7.0 -lvtkFiltersVerdict-7.0 -lvtkInfovisLayout-7.0 -lvtkIOParallelXML-7.0 -lvtkRenderingCore-7.0 -lvtkFiltersCore-7.0 -lvtkfreetype-7.0 -lvtkInteractionImage-7.0 -lvtkIOPLY-7.0 -lvtkRenderingFreeType-7.0 -lvtkFiltersExtraction-7.0 -lvtkGeovisCore-7.0 -lvtkInteractionStyle-7.0 -lvtkIOSQL-7.0 -lvtkRenderingImage-7.0
QMAKE_CXXFLAGS += -DQT_CORE_LIB -DQT_GUI_LIB -DQT_WIDGETS_LIB -Wall

#end VTK

#eigen
INCLUDEPATH += /Users/aglowacki/lib/eigen/include/eigen3

#lmfit
INCLUDEPATH += /Users/aglowacki/lib/lmfit/include
LIBS += -L/Users/aglowacki/lib/lmfit/lib -llmfit


HDF_BASE=/Users/aglowacki/lib/libhdf5
#MPI_BASE=/Users/aglowacki/lib/openmpi
CONFIG += c++14
INCLUDEPATH += $${HDF_BASE}/include
#INCLUDEPATH += $${MPI_BASE}/include
QMAKE_LIBS += -L$${HDF_BASE}/lib -lhdf5
#QMAKE_LIBS += -L$${MPI_BASE}/lib -lmpi
#QMAKE_CXXFLAGS += -std=c++14 -stdlib=libc++
DEFINES += DARWIN __Unix__ __cminpack_double__
}
unix:!macx {
#VTK
INCLUDEPATH += /local/aglowacki/libs/vtk_7/include/vtk-7.0
LIBS += -L/local/aglowacki/libs/vtk_7/lib -lvtkalglib-7.0 -lvtkFiltersFlowPaths-7.0 -lvtkglew-7.0 -lvtkInteractionWidgets-7.0 -lvtkIOVideo-7.0 -lvtkRenderingLabel-7.0 -lvtkChartsCore-7.0 -lvtkFiltersGeneral-7.0 -lvtkGUISupportQt-7.0 -lvtkIOAMR-7.0 -lvtkIOXML-7.0 -lvtkRenderingLOD-7.0 -lvtkCommonColor-7.0 -lvtkFiltersGeneric-7.0 -lvtkGUISupportQtSQL-7.0 -lvtkIOCore-7.0 -lvtkIOXMLParser-7.0 -lvtkRenderingOpenGL2-7.0 -lvtkCommonComputationalGeometry-7.0 -lvtkFiltersGeometry-7.0 -lvtkhdf5-7.0 -lvtkIOEnSight-7.0 -lvtkjpeg-7.0 -lvtkRenderingQt-7.0 -lvtkCommonCore-7.0 -lvtkFiltersHybrid-7.0 -lvtkhdf5_hl-7.0 -lvtkIOExodus-7.0 -lvtkjsoncpp-7.0 -lvtkRenderingVolume-7.0 -lvtkCommonDataModel-7.0 -lvtkFiltersHyperTree-7.0 -lvtkImagingColor-7.0 -lvtkIOExport-7.0 -lvtklibxml2-7.0 -lvtkRenderingVolumeOpenGL2-7.0 -lvtkCommonExecutionModel-7.0 -lvtkFiltersImaging-7.0 -lvtkImagingCore-7.0 -lvtkIOGeometry-7.0 -lvtkmetaio-7.0 -lvtksqlite-7.0 -lvtkCommonMath-7.0 -lvtkFiltersModeling-7.0 -lvtkImagingFourier-7.0 -lvtkIOImage-7.0 -lvtkNetCDF-7.0 -lvtksys-7.0 -lvtkCommonMisc-7.0 -lvtkFiltersParallel-7.0 -lvtkImagingGeneral-7.0 -lvtkIOImport-7.0 -lvtkNetCDF_cxx-7.0 -lvtktiff-7.0 -lvtkCommonSystem-7.0 -lvtkFiltersParallelImaging-7.0 -lvtkImagingHybrid-7.0 -lvtkIOInfovis-7.0 -lvtkoggtheora-7.0 -lvtkverdict-7.0 -lvtkCommonTransforms-7.0 -lvtkFiltersProgrammable-7.0 -lvtkImagingMath-7.0 -lvtkIOLegacy-7.0 -lvtkParallelCore-7.0 -lvtkViewsContext2D-7.0 -lvtkDICOMParser-7.0 -lvtkFiltersSelection-7.0 -lvtkImagingMorphological-7.0 -lvtkIOLSDyna-7.0 -lvtkpng-7.0 -lvtkViewsCore-7.0 -lvtkDomainsChemistry-7.0 -lvtkFiltersSMP-7.0 -lvtkImagingSources-7.0 -lvtkIOMINC-7.0 -lvtkproj4-7.0 -lvtkViewsInfovis-7.0 -lvtkDomainsChemistryOpenGL2-7.0 -lvtkFiltersSources-7.0 -lvtkImagingStatistics-7.0 -lvtkIOMovie-7.0 -lvtkRenderingAnnotation-7.0 -lvtkViewsQt-7.0 -lvtkexoIIc-7.0 -lvtkFiltersStatistics-7.0 -lvtkImagingStencil-7.0 -lvtkIONetCDF-7.0 -lvtkRenderingContext2D-7.0 -lvtkzlib-7.0 -lvtkexpat-7.0 -lvtkFiltersTexture-7.0 -lvtkInfovisCore-7.0 -lvtkIOParallel-7.0 -lvtkRenderingContextOpenGL2-7.0 -lvtkFiltersAMR-7.0 -lvtkFiltersVerdict-7.0 -lvtkInfovisLayout-7.0 -lvtkIOParallelXML-7.0 -lvtkRenderingCore-7.0 -lvtkFiltersCore-7.0 -lvtkfreetype-7.0 -lvtkInteractionImage-7.0 -lvtkIOPLY-7.0 -lvtkRenderingFreeType-7.0 -lvtkFiltersExtraction-7.0 -lvtkGeovisCore-7.0 -lvtkInteractionStyle-7.0 -lvtkIOSQL-7.0 -lvtkRenderingImage-7.0
QMAKE_CXXFLAGS += -DQT_CORE_LIB -DQT_GUI_LIB -DQT_WIDGETS_LIB -D_BUILD_WITH_VTK
#end VTK

#eigen
INCLUDEPATH += /local/aglowacki/libs/eigen3/include/eigen3

#ceres
INCLUDEPATH += /local/aglowacki/libs/ceres/include
LIBS += /local/aglowacki/libs/ceres/lib64/libceres.a
INCLUDEPATH += /local/aglowacki/libs/glog/include
LIBS += /local/aglowacki/libs/glog/lib/libglog.a
INCLUDEPATH += /local/aglowacki/libs/gflags/include
LIBS += /local/aglowacki/libs/gflags/lib/libgflags.a
QMAKE_CXXFLAGS += -fopenmp
LIBS += -lgomp

#lmfit
INCLUDEPATH += /local/aglowacki/libs/lmfit/include
LIBS += -L/local/aglowacki/libs/lmfit/lib -llmfit

INCLUDEPATH += /usr/include/hdf5/serial
QMAKE_CXXFLAGS += -D__Unix__ -std=c++14 -D__cminpack_double__ -g
QMAKE_LIBS += -lhdf5


#netcdf
INCLUDEPATH += /local/aglowacki/libs/netcdf-c/include
LIBS += -L/local/aglowacki/libs/netcdf-c/lib -lnetcdf
}

win32{
INCLUDEPATH += "C:/Program Files/HDF_Group/HDF5/1.10.0/include"
DEFINES += XDR_HACK XDR_LE _USE_MATH_DEFINES H5_BUILT_AS_DYNAMIC_LIB __cminpack_double__
HEADERS += support/mda_utils/xdr_hack.h
SOURCES += support/mda_utils/xdr_hack.c
QMAKE_LIBS += -L"C:/Program Files/HDF_Group/HDF5/1.10.0/lib" -lszip -lzlib -lhdf5 -lhdf5_hl
}

SOURCES += \
    src/support/nnls/nnls.cc \
    src/support/nnls/denseMatrix.cc \
    src/support/mda_utils/mda_loader.c \
    src/support/minpack/minpack.cpp \
    src/support/cmpfit-1.3a/mpfit.c \
    src/data_struct/xrf/element_info.cpp \
    src/data_struct/xrf/fit_parameters.cpp \
    src/data_struct/xrf/fit_element_map.cpp \
    src/data_struct/xrf/quantification_standard.cpp \
    #src/data_struct/xrf/element_quant.cpp \
    src/data_struct/xrf/detector.cpp \
    src/data_struct/xrf/spectra.cpp \
    src/data_struct/xrf/spectra_line.cpp \
    src/data_struct/xrf/spectra_volume.cpp \
    #src/data_struct/xrf/aps/aps_fit_parameters.cpp \
    src/fitting/models/base_model.cpp \
    src/fitting/models/gaussian_model.cpp \
    src/fitting/routines/base_fit_routine.cpp \
    src/fitting/routines/roi_fit_routine.cpp \
    src/fitting/routines/param_optimized_fit_routine.cpp \
    src/fitting/routines/matrix_optimized_fit_routine.cpp \
    src/fitting/routines/svd_fit_routine.cpp \
    src/fitting/routines/nnls_fit_routine.cpp \
    src/fitting/optimizers/optimizer.cpp \
    src/fitting/optimizers/lmfit_optimizer.cpp \
    src/fitting/optimizers/mpfit_optimizer.cpp \
    src/fitting/optimizers/minpack_optimizer.cpp \
    src/quantification/models/quantification_model.cpp \
    #src/fitting/optimizers/ceres_optimizer.cpp \
    #src/workflow/task.cpp \
    #src/workflow/task_queue.cpp \
    #src/workflow/general_function_task.cpp \
    #src/workflow/discretizer.cpp \
    #src/workflow/distributor.cpp \
    #src/workflow/distributors/threadpool_distributor.cpp \
    src/io/file/mda_io.cpp \
    src/io/file/hdf5_io.cpp \
    src/io/file/netcdf_io.cpp \
    src/io/file/csv_io.cpp \
    src/io/file/aps/aps_fit_params_import.cpp \
    #src/io/file/aps/aps_calibration_io.cpp \
    src/visual/vtk_graph.cpp \
    src/core/main.cpp

HEADERS += \
    #pybind11/cast.h \
    #pybind11/common.h \
    #pybind11/operators.h \
    #pybind11/pybind11.h \
    #pybind11/pytypes.h \
    #pybind11/typeid.h \
    #pybind11/numpy.h \
    src/core/defines.h \
    src/support/nnls/nnls.h \
    src/support/nnls/denseMatrix.h \
    src/support/nnls/matrix.h \
    src/support/nnls/vector.h \
    src/support/mda_utils/mda-load.h \
    src/support/minpack/minpack.hpp \
    src/support/cmpfit-1.3a/mpfit.h \
    src/data_struct/xrf/element_info.h \
    src/data_struct/xrf/fit_parameters.h \
    src/data_struct/xrf/fit_element_map.h \
    src/data_struct/xrf/params_override.h \
    src/data_struct/xrf/element_quant.h \
    src/data_struct/xrf/quantification_standard.h \
    src/data_struct/xrf/detector.h \
    src/data_struct/base_dataset.h \
    #src/data_struct/ndarray.h \
    src/data_struct/xrf/spectra.h \
    src/data_struct/xrf/spectra_line.h \
    src/data_struct/xrf/spectra_volume.h \
    #src/data_struct/xrf/aps/aps_fit_parameters.h \
    src/fitting/models/base_model.h \
    src/fitting/models/gaussian_model.h \
    src/fitting/routines/base_fit_routine.h \
    src/fitting/routines/roi_fit_routine.h \
    src/fitting/routines/param_optimized_fit_routine.h \
    src/fitting/routines/matrix_optimized_fit_routine.h \
    src/fitting/routines/svd_fit_routine.h \
    src/fitting/routines/nnls_fit_routine.h \
    src/fitting/optimizers/optimizer.h \
    src/fitting/optimizers/lmfit_optimizer.h \
    src/fitting/optimizers/mpfit_optimizer.h \
    src/fitting/optimizers/minpack_optimizer.h \
    src/quantification/models/quantification_model.h \
    #src/fitting/optimizers/ceres_optimizer.h \
    #src/workflow/task.h \
    #src/workflow/task_queue.h \
    #src/workflow/general_function_task.h \
    #src/workflow/discretizer.h \
    #src/workflow/distributor.h \
    #src/workflow/distributors/threadpool_distributor.h \
    src/workflow/threadpool.h \
    src/io/file/base_file_io.h \
    src/io/file/mda_io.h \
    src/io/file/hdf5_io.h \
    src/io/file/netcdf_io.h \
    src/io/file/csv_io.h \
    src/io/file/aps/aps_fit_params_import.h \
    #src/io/file/aps/aps_calibration_io.h \
    src/visual/vtk_graph.h \
    src/core/command_line_parser.h


