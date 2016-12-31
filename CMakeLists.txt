# CMakeLists.txt -- Build system for the XRF-Maps
#
#Copyright (c) 2016, UChicago Argonne, LLC. All rights reserved.
#
#Copyright 2016. UChicago Argonne, LLC. This software was produced
#under U.S. Government contract DE-AC02-06CH11357 for Argonne National
#Laboratory (ANL), which is operated by UChicago Argonne, LLC for the
#U.S. Department of Energy. The U.S. Government has rights to use,
#reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
#UChicago Argonne, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
#ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is
#modified to produce derivative works, such modified software should
#be clearly marked, so as not to confuse it with the version available
#from ANL.
#
#Additionally, redistribution and use in source and binary forms, with
#or without modification, are permitted provided that the following
#conditions are met:
#
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimer.
#
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the
#      distribution.
#
#    * Neither the name of UChicago Argonne, LLC, Argonne National
#      Laboratory, ANL, the U.S. Government, nor the names of its
#      contributors may be used to endorse or promote products derived
#      from this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY UChicago Argonne, LLC AND CONTRIBUTORS
#"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL UChicago
#Argonne, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#POSSIBILITY OF SUCH DAMAGE.


cmake_minimum_required(VERSION 2.8)

project(XRF_Maps)

#option(PYBIND11_INSTALL "Install pybind11 header files?" ON)

# Add a CMake parameter for choosing a desired Python version
#set(PYBIND11_PYTHON_VERSION "" CACHE STRING "Python version to use for compiling the XRF Maps application")

set(GCC_INCLUDE_DIR "/usr/include" CACHE PATH "Aditional GCC include folder")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" ${PROJECT_SOURCE_DIR}/cmake)

include_directories(GCC_INCLUDE_DIR)

set(EIGEN3_INCLUDE_DIR "/usr/include/eigen3" CACHE PATH "Eigen include folder")
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()

#set(NetCDF_INCLUDE "" CACHE PATH "NetCDF-C include folder")
#set(NetCDF_LIB_DIR "" CACHE PATH "NetCDF-C lib folder")
find_package(NetCDF REQUIRED)

set(SZip_LIB_DIR "" CACHE PATH "SZip lib folder")


find_package (Threads)
find_package(HDF5 REQUIRED COMPONENTS C)


#find_package(VTK REQUIRED)
#include(${VTK_USE_FILE})

link_directories( ${HDF5_LIB} ${SZip_LIB_DIR} ${NetCDF_LIB_DIR} )

# Set a default build configuration if none is specified. 'MinSizeRel' produces the smallest binaries
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'MinSizeRel' as none was specified.")
  set(CMAKE_BUILD_TYPE MinSizeRel CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

#set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6)
#find_package(PythonLibs ${PYBIND11_PYTHON_VERSION} REQUIRED)
#find_package(PythonInterp ${PYBIND11_PYTHON_VERSION} REQUIRED)

include(CheckCXXCompilerFlag)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Enable C++11 mode on C++ / Clang
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -D__Unix__ -D_REAL_DOUBLE -fPIC -g")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__Unix__ -D_REAL_DOUBLE -fPIC")

  # Enable link time optimization and set the default symbol
  # visibility to hidden (very important to obtain small binaries)
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    # Default symbol visibility
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")

    # Check for Link Time Optimization support
    CHECK_CXX_COMPILER_FLAG("-flto" HAS_LTO_FLAG)
    if (HAS_LTO_FLAG)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto")
    endif()
  endif()
endif()

# Compile with compiler warnings turned on
if(MSVC)
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /DHDF5CPP_USEDLL /D_USE_MATH_DEFINES /DH5_BUILT_AS_DYNAMIC_LIB /D_HDF5USEDLL_ /DBUILDING_DLL /DXDR_HACK /DXDR_LE /D_REAL_DOUBLE /DNOMINMAX")
elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# Include path for Python header files
#include_directories(${PYTHON_INCLUDE_DIR})
# Include path for pybind11 header files
#include_directories(pybind11)

#hdf5 include
include_directories(${HDF5_INCLUDE_DIRS})

# Include Eigen 
include_directories(${EIGEN3_INCLUDE_DIR})

# Include NetCDF 
include_directories(${NETCDF_INCLUDES})

# Include libxrf header files
include_directories(src/support/minpack)
include_directories(src/support/mdautils-1.4.1)
include_directories(src/support/faddeeva)
include_directories(src/support/cmpfit-1.3a)
include_directories(src/support/lmfit_6.1)
include_directories(src/support/nnls)
include_directories(src/data_struct)
include_directories(src/data_struct/xrf)
include_directories(src/data_struct/xrf/aps)
include_directories(src/quantification)
include_directories(src/quantification/models)
include_directories(src/fitting)
include_directories(src/fitting/models)
include_directories(src/fitting/routines)
include_directories(src/fitting/optimizers)

# Include libxrf_io Maps header files
include_directories(src/io/file)
include_directories(src/io/file/aps)
include_directories(src/workflow)

# Include XRF-Maps Maps header files
include_directories(src/core)

# Visual Studio specific include and source
if(MSVC)
  set(XDR_HACK src/support/mdautils-1.4.1/xdr_hack.c )
  include_directories(src/support/direct)
else()
  set(XDR_HACK  )
endif()

#set(PYBIND11_HEADERS
#  pybind11/cast.h
#  pybind11/common.h
#  pybind11/operators.h
#  pybind11/pybind11.h
#  pybind11/pytypes.h
#  pybind11/typeid.h
#  pybind11/numpy.h
#)

set(libxrf_fit_HEADERS
    src/core/defines.h
    src/support/faddeeva/Faddeeva.hh
    src/support/minpack/minpack.hpp
    src/support/cmpfit-1.3a/mpfit.h
    src/support/lmfit_6.1/lmstruct.h
    src/support/lmfit_6.1/lmmin.h
    src/support/lmfit_6.1/lmcurve_tyd.h
    src/support/lmfit_6.1/lmcurve.h
    src/support/nnls/nnls.h
    src/support/nnls/vector.h
    src/support/nnls/matrix.h
    src/support/nnls/denseMatrix.h
    src/data_struct/base_dataset.h
    #src/data_struct/ndarray.h
    src/data_struct/xrf/quantification_standard.h
    src/data_struct/xrf/element_quant.h
    src/data_struct/xrf/element_info.h
    src/data_struct/xrf/fit_parameters.h
    src/data_struct/xrf/fit_element_map.h
    src/data_struct/xrf/detector.h
    src/data_struct/xrf/spectra.h
    src/data_struct/xrf/spectra_line.h
    src/data_struct/xrf/spectra_volume.h
    src/quantification/models/quantification_model.h
    src/fitting/models/base_model.h
    src/fitting/models/gaussian_model.h
    src/fitting/routines/base_fit_routine.h
    src/fitting/routines/param_optimized_fit_routine.h
    src/fitting/routines/matrix_optimized_fit_routine.h
    src/fitting/routines/roi_fit_routine.h
    src/fitting/routines/svd_fit_routine.h
    src/fitting/routines/nnls_fit_routine.h
    src/fitting/optimizers/optimizer.h
    src/fitting/optimizers/mpfit_optimizer.h
    src/fitting/optimizers/lmfit_optimizer.h
    #src/workflow/task.h
    #src/workflow/task_queue.h
    #src/workflow/general_function_task.h
    src/workflow/threadpool.h
)

# Create the binding library
add_library(libxrf_fit SHARED
    #${PYBIND11_HEADERS}
    ${libxrf_fit_HEADERS}
    src/support/faddeeva/Faddeeva.cc
    src/support/minpack/minpack.cpp
    src/support/cmpfit-1.3a/mpfit.c
    src/support/lmfit_6.1/lmmin.c
    src/support/lmfit_6.1/lmcurve_tyd.c
    src/support/lmfit_6.1/lmcurve.c
    src/support/nnls/nnls.cc
    src/support/nnls/denseMatrix.cc
    src/data_struct/xrf/quantification_standard.cpp
    src/data_struct/xrf/element_info.cpp
    src/data_struct/xrf/fit_parameters.cpp
    src/data_struct/xrf/fit_element_map.cpp
    src/data_struct/xrf/detector.cpp
    src/data_struct/xrf/spectra.cpp
    src/data_struct/xrf/spectra_line.cpp
    src/data_struct/xrf/spectra_volume.cpp
    src/quantification/models/quantification_model.cpp
    src/fitting/models/base_model.cpp
    src/fitting/models/gaussian_model.cpp
    src/fitting/routines/base_fit_routine.cpp
    src/fitting/routines/param_optimized_fit_routine.cpp
    src/fitting/routines/matrix_optimized_fit_routine.cpp
    src/fitting/routines/roi_fit_routine.cpp
    src/fitting/routines/svd_fit_routine.cpp
    src/fitting/routines/nnls_fit_routine.cpp
    src/fitting/optimizers/optimizer.cpp
    src/fitting/optimizers/mpfit_optimizer.cpp
    src/fitting/optimizers/lmfit_optimizer.cpp
    #src/workflow/task.cpp
    #src/workflow/task_queue.cpp
    #src/workflow/general_function_task.cpp
)

#--------------- start xrf io lib -----------------
set(XRF_IO_HEADERS
	src/support/mdautils-1.4.1/mda-load.h
	src/io/file/base_file_io.h
	src/io/file/mda_io.h
	src/io/file/hdf5_io.h
	src/io/file/netcdf_io.h
	src/io/file/csv_io.h
        #src/io/file/aps/aps_calibration_io.h
	src/io/file/aps/aps_fit_params_import.h
	src/core/command_line_parser.h
)

add_library(libxrf_io SHARED
	${XDR_HACK}
        ${XRF_IO_HEADERS}
	src/support/mdautils-1.4.1/mda_loader.c
	src/io/file/mda_io.cpp
	src/io/file/hdf5_io.cpp
	src/io/file/netcdf_io.cpp
	src/io/file/csv_io.cpp
        #src/io/file/aps/aps_calibration_io.cpp
	src/io/file/aps/aps_fit_params_import.cpp				 
	src/core/main.cpp
)

#--------------- start xrf maps exec -----------------
add_executable(xrf_maps
        ${XDR_HACK}
        ${XRF_IO_HEADERS}
        ${libxrf_fit_HEADERS}
        src/core/main.cpp
)

# Don't add a 'lib' prefix to the shared library
set_target_properties(libxrf_fit PROPERTIES PREFIX "")
set_target_properties(libxrf_io PROPERTIES PREFIX "")

# Always write the output file directly into the 'libxrf_fit' directory (even on MSVC)
set(CompilerFlags
  LIBRARY_OUTPUT_DIRECTORY LIBRARY_OUTPUT_DIRECTORY_RELEASE LIBRARY_OUTPUT_DIRECTORY_DEBUG
  LIBRARY_OUTPUT_DIRECTORY_MINSIZEREL LIBRARY_OUTPUT_DIRECTORY_RELWITHDEBINFO
  RUNTIME_OUTPUT_DIRECTORY RUNTIME_OUTPUT_DIRECTORY_RELEASE RUNTIME_OUTPUT_DIRECTORY_DEBUG
  RUNTIME_OUTPUT_DIRECTORY_MINSIZEREL RUNTIME_OUTPUT_DIRECTORY_RELWITHDEBINFO)

foreach(CompilerFlag ${CompilerFlags})
  set_target_properties(libxrf_fit libxrf_io xrf_maps PROPERTIES ${CompilerFlag} ${PROJECT_SOURCE_DIR}/bin)
endforeach()

if (WIN32)
  if (MSVC)
    # Enforce size-based optimization and link time code generation
    # on MSVC (~30% smaller binaries in experiments). /bigobj is needed
    # for bigger binding projects due to the limit to 64k addressable sections
    # /MP enables multithreaded builds (relevant when there are many files).
    set_target_properties(libxrf_fit libxrf_io PROPERTIES COMPILE_FLAGS "/Os /GL /MP /bigobj /DDYNAMIC_LIB")
    set_target_properties(libxrf_fit libxrf_io PROPERTIES LINK_FLAGS "/LTCG")

	set(NETCDF_LIBRARIES netcdf.lib)
	set(SZIP_LIB szip.lib)
	# Link against shared libraries
#    target_link_libraries(libxrf_fit  )
#	target_link_libraries(libxrf_fit  )
  endif()

  # .PYD file extension on Windows
  #set_target_properties(libxrf_fit PROPERTIES SUFFIX ".pyd")
 
elseif (UNIX)
  # It's quite common to have multiple copies of the same Python version
  # installed on one's system. E.g.: one copy from the OS and another copy
  # that's statically linked into an application like Blender or Maya.
  # If we link our plugin library against the OS Python here and import it
  # into Blender or Maya later on, this will cause segfaults when multiple
  # conflicting Python instances are active at the same time.

  # Windows is not affected by this issue since it handles DLL imports 
  # differently. The solution for Linux and Mac OS is simple: we just don't
  # link against the Python library. The resulting shared library will have
  # missing symbols, but that's perfectly fine -- they will be resolved at
  # import time.

  #set(NETCDF_LIB netcdf)
  set(SZIP_LIB sz)

#  target_link_libraries(libxrf_fit optimized hdf5 optimized netcdf_c )
#  target_link_libraries(libxrf_fit debug hdf5 debug netcdf_c )

  # Optimize for a small binary size
  if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    set_target_properties(libxrf_fit  libxrf_io PROPERTIES COMPILE_FLAGS "-Os")
  endif()

  # Strip unnecessary sections of the binary on Linux/Mac OS
  if(APPLE)
    set_target_properties(libxrf_fit libxrf_io PROPERTIES SUFFIX ".dylib")
    set_target_properties(libxrf_fit libxrf_io PROPERTIES COMPILE_FLAGS "-DDARWIN")
    set_target_properties(libxrf_fit libxrf_io PROPERTIES MACOSX_RPATH ".")
    set_target_properties(libxrf_fit libxrf_io PROPERTIES LINK_FLAGS "-undefined dynamic_lookup -dead_strip")
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET libxrf_fit POST_BUILD COMMAND strip -u -r ${PROJECT_SOURCE_DIR}/bin/libxrf_fit.dylib)
      add_custom_command(TARGET libxrf_io POST_BUILD COMMAND strip -u -r ${PROJECT_SOURCE_DIR}/bin/libxrf_io.dylib)
    endif()
  else()
    set_target_properties(libxrf_fit libxrf_io PROPERTIES SUFFIX ".so")
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET libxrf_fit POST_BUILD COMMAND strip ${PROJECT_SOURCE_DIR}/bin/libxrf_fit.so)
      add_custom_command(TARGET libxrf_io POST_BUILD COMMAND strip ${PROJECT_SOURCE_DIR}/bin/libxrf_io.so)
    endif()
  endif()
endif()

#enable_testing()
#set(RUN_TEST ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/bin/run_test.py)
#foreach(i RANGE 1 12)
#  add_test(NAME libxrf_fit${i} COMMAND ${RUN_TEST} libxrf_fit${i})
#endforeach()

#if (PYBIND11_INSTALL)
#    install(FILES ${PYBIND11_HEADERS} DESTINATION include/pybind11)
#endif()

target_link_libraries (xrf_maps LINK_PUBLIC optimized libxrf_fit optimized libxrf_io optimized ${NETCDF_LIBRARIES} optimized ${HDF5_LIBRARIES} optimized ${SZIP_LIB} ${CMAKE_THREAD_LIBS_INIT} )
target_link_libraries (xrf_maps LINK_PUBLIC debug libxrf_fit debug libxrf_io debug ${NETCDF_LIBRARIES} debug ${HDF5_LIBRARIES} debug ${SZIP_LIB} ${CMAKE_THREAD_LIBS_INIT} )
