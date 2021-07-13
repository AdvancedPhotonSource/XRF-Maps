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


cmake_minimum_required(VERSION 3.5)

project(XRF_Maps)

set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/" ${PROJECT_SOURCE_DIR}/../libs/libzmq/build ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${PROJECT_SOURCE_DIR}/../libs/libzmq/build)

option(BUILD_WITH_TIRPC "Build tripc library for libxrpc (RedHat 8)" OFF)

IF(UNIX)
  IF(BUILD_WITH_TIRPC)
    set(GCC_INCLUDE_DIR "/usr/include/tirpc/" CACHE PATH "include folder for xrpc (red hat)") 
    include_directories(${GCC_INCLUDE_DIR})
  ENDIF()
  set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake/linux)
ENDIF()

find_package(NetCDF REQUIRED)
find_package(Threads)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(yaml-cpp REQUIRED)

set(EIGEN3_INCLUDES "${PROJECT_SOURCE_DIR}/src/support/eigen-git-mirror" CACHE PATH "Eigen include folder")

option(BUILD_WITH_PYBIND11 "Build python binding with PyBind11" OFF)
option(BUILD_WITH_ZMQ "Build with ZeroMQ" OFF)
option(BUILD_FOR_PHI "Build for Intel Phi" OFF)
option(BUILD_WITH_QT "Build with QT" OFF)
option(DOUBLE_PREC "Double Precision" OFF)
option(STATIC_BUILD "Static build libxrf_io and libxrf_fit" OFF)
# If compiled on some intel mahcines this causes crashes so let user set it for compile
option(AVX512 "Compule with arch AVX512 on MSVC" OFF)
option(AVX2 "Compule with arch AVX2 on MSVC" OFF)
option(AVX "Compule with arch AVX on MSVC" OFF)
option(SSE2 "Compule with arch SSE2 on MSVC" OFF)

# Opitons for getting YouCompleteMe plugin working with this project
SET( CMAKE_EXPORT_COMPILE_COMMANDS ON )

IF( EXISTS "${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json" )
  EXECUTE_PROCESS( COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${CMAKE_CURRENT_BINARY_DIR}/compile_commands.json
    ${CMAKE_CURRENT_SOURCE_DIR}/compile_commands.json
  )
ENDIF()


IF (BUILD_WITH_ZMQ)
  add_definitions(-D_BUILD_WITH_ZMQ)
  IF( NOT ZeroMQ_INCLUDE_DIR AND NOT ZeroMQ_LIBRARY)
      find_package(ZeroMQ REQUIRED libzmq-static)
  ENDIF()
ENDIF()

IF( NOT EIGEN3_INCLUDES )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDES to the include directory of your Eigen3 installation.")
ENDIF()

IF (DOUBLE_PREC)
  set(PREC _REAL_DOUBLE)
ELSE()
  set(PREC _REAL_FLOAT)
ENDIF()

IF (BUILD_WITH_QT)
  IF(NOT DEFINED CMAKE_PREFIX_PATH)
    set(CMAKE_PREFIX_PATH CACHE PATH "")
    MESSAGE( "IF qt5-config.cmake is not found, set CMAKE_PREFIX_PATH = Qt root directory")
  ENDIF()
  add_definitions(-D_BUILD_WITH_QT)
  set(CMAKE_CXX_STANDARD 14)
  find_package(Qt5 COMPONENTS Core Charts REQUIRED)
  set(VISUAL_INC src/visual/grapher.h src/visual/grapher.cpp )
  include_directories(src/visual/)
  include_directories(${Qt5Charts_INCLUDES})
  set(CMAKE_CXX_FLAGS "${Qt5Charts_EXECUTABLE_COMPILE_FLAGS}")
ELSE()
  set(VISUAL_INC )
ENDIF()

# Set a default build configuration if none is specified. 'Release' produces the smallest binaries
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
ENDIF()

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)

include(CheckCXXCompilerFlag)

IF (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  # Enable C++11 mode on C++ / Clang
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y -D__Unix__ -D${PREC} -fPIC -g -fopenmp")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D__Unix__ -D${PREC} -fPIC")
  IF(BUILD_FOR_PHI AND CMAKE_CXX_COMPILER_ID MATCHES "Intel" )
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xMIC-AVX512")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -xMIC-AVX512")
  ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
  ENDIF()
  #set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Waddress -Warray-bounds=1 -Wbool-compare -Wbool-operation -Wc++11-compat -Wc++14-compat -Wcatch-value -Wchar-subscripts -Wcomment -Wenum-compare -Wformat -Winit-self -Wlogical-not-parentheses -Wmain -Wmaybe-uninitialized -Wmemset-elt-size -Wmemset-transposed-args -Wmisleading-indentation -Wmissing-attributes -Wmissing-braces -Wmultistatement-macros -Wnarrowing -Wnonnull -Wnonnull-compare -Wopenmp-simd -Wparentheses -Wreorder -Wrestrict -Wreturn-type -Wsequence-point -Wsign-compare -Wsizeof-pointer-div -Wsizeof-pointer-memaccess -Wstrict-aliasing -Wstrict-overflow=1 -Wswitch -Wtautological-compare -Wtrigraphs -Wuninitialized -Wunknown-pragmas -Wunused-function -Wunused-label -Wunused-value -Wvolatile-register-var")
ELSEIF(MSVC)
  IF(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  ELSE()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  ENDIF()

  IF (AVX512)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX512")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX512")
  ELSE()
    IF (AVX2)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX2")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2")
    ELSE()
      IF (AVX)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:AVX")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX")
      ELSE()
        IF(SSE2)
          set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /arch:SSE2")
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:SSE2")
        ENDIF()
      ENDIF()
    ENDIF()
  ENDIF()
    
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_USE_MATH_DEFINES /DH5_BUILT_AS_DYNAMIC_LIB /D_HDF5USEDLL_ /DBUILDING_DLL /DXDR_HACK /DXDR_LE /D${PREC} /DNOMINMAX /openmp")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox /Ot /Qvec-report:2 /GL /MP /bigobj")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /LTCG")
ENDIF()

IF(BUILD_WITH_ZMQ)
  #zeromq include
  include_directories(${ZeroMQ_INCLUDE_DIR})
ENDIF()

#hdf5 include
include_directories(${HDF5_INCLUDE_DIRS})

# Include Eigen 
include_directories(${EIGEN3_INCLUDES})

# Include NetCDF 
include_directories(${NETCDF_INCLUDE_DIRS})

# Include libxrf header files
include_directories(src)

# Visual Studio specific include and source
IF(MSVC)
  set(XDR_HACK src/support/mdautils-1.4.1/xdr_hack.c src/support/mdautils-1.4.1/xdr_hack.h )
  include_directories(src/support/direct)
ELSE()
  set(XDR_HACK  )
ENDIF()

#--------------- start xrf lib -----------------
set(libxrf_fit_HEADERS
    src/core/defines.h
    src/support/cmpfit-1.3a/mpfit.hpp
    src/support/lmfit_6.1/lmstruct.hpp
    src/support/lmfit_6.1/lmmin.hpp
    src/support/lmfit_6.1/lmcurve_tyd.hpp
    src/support/lmfit_6.1/lmcurve.hpp
    src/support/nnls/nnls.hpp
    src/data_struct/quantification_standard.h
    src/data_struct/element_quant.h
    src/data_struct/element_info.h
    src/data_struct/scaler_lookup.h
    src/data_struct/fit_parameters.h
    src/data_struct/fit_element_map.h
    src/data_struct/params_override.h
    src/data_struct/scan_info.h
    src/data_struct/spectra.h
    src/data_struct/spectra_line.h
    src/data_struct/spectra_volume.h
    src/data_struct/stream_block.h
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
    src/data_struct/detector.h
    src/data_struct/analysis_job.h
    src/workflow/threadpool.h
)

set(libxrf_fit_SOURCE
    src/data_struct/quantification_standard.cpp
    src/data_struct/element_info.cpp
    src/data_struct/scaler_lookup.cpp
    src/data_struct/fit_parameters.cpp
    src/data_struct/fit_element_map.cpp
    src/data_struct/spectra.cpp
    src/data_struct/spectra_line.cpp
    src/data_struct/spectra_volume.cpp
    src/data_struct/stream_block.cpp
    src/quantification/models/quantification_model.cpp
    src/fitting/models/gaussian_model.cpp
    src/fitting/routines/param_optimized_fit_routine.cpp
    src/fitting/routines/matrix_optimized_fit_routine.cpp
    src/fitting/routines/roi_fit_routine.cpp
    src/fitting/routines/svd_fit_routine.cpp
    src/fitting/routines/nnls_fit_routine.cpp
    src/fitting/optimizers/optimizer.cpp
    src/fitting/optimizers/mpfit_optimizer.cpp
    src/fitting/optimizers/lmfit_optimizer.cpp
    src/data_struct/detector.cpp
    src/data_struct/analysis_job.cpp
)

#--------------- start xrf io lib -----------------
set(XRF_IO_HEADERS
	src/core/mem_info.h
	src/support/mdautils-1.4.1/mda-load.h
	src/io/file/mda_io.h
        src/io/file/mca_io.h
	src/io/file/hdf5_io.h
	src/io/file/netcdf_io.h
	src/io/file/csv_io.h
	src/io/file/aps/aps_fit_params_import.h
	src/io/file/hl_file_io.h
	src/io/net/basic_serializer.h
	src/workflow/source.h
	src/workflow/distributor.h
	src/workflow/sink.h
	src/workflow/xrf/spectra_file_source.h
	src/workflow/xrf/spectra_net_source.h
	src/workflow/xrf/integrated_spectra_source.h
	src/workflow/xrf/detector_sum_spectra_source.h
	src/workflow/xrf/spectra_stream_saver.h
	src/workflow/xrf/spectra_net_streamer.h
    src/core/process_streaming.h
    src/core/process_whole.h
   )

set(XRF_IO_SOURCE
    ${VISUAL_INC}
	src/core/mem_info.cpp
    src/support/zmq/zmq.hpp
    src/support/mdautils-1.4.1/mda_loader.c
    src/io/file/mda_io.cpp
    src/io/file/mca_io.cpp
    src/io/file/hdf5_io.cpp
    src/io/file/netcdf_io.cpp
    src/io/file/csv_io.cpp
    src/io/file/aps/aps_fit_params_import.cpp
    src/io/file/hl_file_io.cpp
    src/io/net/basic_serializer.cpp
    src/workflow/xrf/spectra_file_source.cpp
    src/workflow/xrf/spectra_net_source.cpp
    src/workflow/xrf/integrated_spectra_source.cpp
    src/workflow/xrf/detector_sum_spectra_source.cpp
    src/workflow/xrf/spectra_stream_saver.cpp
    src/workflow/xrf/spectra_net_streamer.cpp
    src/core/process_streaming.cpp
    src/core/process_whole.cpp
    )

IF(BUILD_FOR_PHI)
  add_library(libxrf_fit STATIC
              ${libxrf_fit_HEADERS}
              ${libxrf_fit_SOURCE}
              )
  add_library(libxrf_io STATIC
              ${XRF_IO_HEADERS}
              ${XDR_HACK}
              ${XRF_IO_SOURCE}
              )
ELSE()
  add_library(libxrf_fit SHARED 
              ${libxrf_fit_HEADERS}
              ${libxrf_fit_SOURCE}
              )
  add_library(libxrf_io SHARED 
             ${XRF_IO_HEADERS}
             ${XDR_HACK}
             ${XRF_IO_SOURCE}
             )
ENDIF()

#--------------- start pyxrf-maps lib -----------------

IF (BUILD_WITH_PYBIND11)
    add_subdirectory(src/support/pybind11)
    # .PYD file extension on Windows
    #set_target_properties(pyxrfmaps PROPERTIES SUFFIX ".pyd")

    pybind11_add_module(pyxrfmaps src/pybindings/main.cpp)
    IF (BUILD_WITH_ZMQ)
      target_link_libraries(pyxrfmaps LINK_PUBLIC libxrf_fit libxrf_io ${NETCDF_LIBRARIES} ${HDF5_LIBRARIES} libzmq-static ws2_32.lib rpcrt4.lib iphlpapi.lib)
    ELSE()
      target_link_libraries(pyxrfmaps PUBLIC libxrf_fit PUBLIC libxrf_io PUBLIC ${NETCDF_LIBRARIES} PUBLIC ${HDF5_LIBRARIES})
    ENDIF()
ENDIF()

#--------------- start xrf maps exec -----------------
add_executable(xrf_maps
    src/core/command_line_parser.h
    src/core/main.cpp
)

# Don't add a 'lib' prefix to the shared library
set_target_properties(libxrf_fit PROPERTIES PREFIX "")
set_target_properties(libxrf_io PROPERTIES PREFIX "")

# Always write the output file directly into the 'libxrf_fit' directory (even on MSVC)
set(CompilerFlags
  LIBRARY_OUTPUT_DIRECTORY LIBRARY_OUTPUT_DIRECTORY_RELEASE LIBRARY_OUTPUT_DIRECTORY_DEBUG
  RUNTIME_OUTPUT_DIRECTORY RUNTIME_OUTPUT_DIRECTORY_RELEASE RUNTIME_OUTPUT_DIRECTORY_DEBUG)

IF (BUILD_WITH_PYBIND11)
  foreach(CompilerFlag ${CompilerFlags})
    set_target_properties(pyxrfmaps PROPERTIES ${CompilerFlag} ${PROJECT_SOURCE_DIR}/bin)
  endforeach()
ENDIF()

foreach(CompilerFlag ${CompilerFlags})
  set_target_properties(libxrf_fit libxrf_io xrf_maps PROPERTIES ${CompilerFlag} ${PROJECT_SOURCE_DIR}/bin)
endforeach()


set_property(TARGET libxrf_io libxrf_fit APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src/support/eigen-git-mirror> 
  $<INSTALL_INTERFACE:src>  )

install(TARGETS libxrf_io libxrf_fit EXPORT XRF_MapsConfig 
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
#install(EXPORT XRF_MapsConfig DESTINATION share/XRF_Maps/cmake)
#export(TARGETS libxrf_fit libxrf_io FILE "${CMAKE_CURRENT_BINARY_DIR}/XRF_MapsTargets.cmake")
export(EXPORT XRF_MapsConfig FILE "${CMAKE_CURRENT_BINARY_DIR}/XRF_MapsConfig.cmake")
export(PACKAGE XRF_MapsConfig)

IF (WIN32)
  IF (MSVC)
    # /bigobj is needed for bigger binding projects due to the limit to 64k addressable sections
    # /MP enables multithreaded builds (relevant when there are many files).
    set_target_properties(libxrf_fit libxrf_io PROPERTIES COMPILE_FLAGS "/DDYNAMIC_LIB")
	set_target_properties(xrf_maps libxrf_io PROPERTIES COMPILE_FLAGS "/D_WINSOCKAPI_")
  ENDIF()
  target_link_libraries(libxrf_io LINK_PUBLIC optimized libxrf_fit optimized ${NETCDF_LIBRARIES} optimized ${HDF5_LIBRARIES} yaml-cpp ${CMAKE_THREAD_LIBS_INIT} )
  target_link_libraries(libxrf_io LINK_PUBLIC debug libxrf_fit debug ${NETCDF_LIBRARIES} debug ${HDF5_LIBRARIES} yaml-cpp ${CMAKE_THREAD_LIBS_INIT} )
  
  IF (BUILD_WITH_ZMQ)
    #IF (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
    #  string(REGEX REPLACE "mt-s" "mt" ZeroMQ_LIBRARY "${ZeroMQ_LIBRARY}")
    #ELSE()
    #  string(REGEX REPLACE "sgd" "gd" ZeroMQ_LIBRARY "${ZeroMQ_LIBRARY}")
    #ENDIF()
	target_link_libraries(libxrf_io LINK_PUBLIC libzmq-static ws2_32.lib rpcrt4.lib iphlpapi.lib)
    target_link_libraries(xrf_maps LINK_PUBLIC libzmq-static ws2_32.lib rpcrt4.lib iphlpapi.lib)
  
  #  target_link_libraries(libxrf_io LINK_PUBLIC optimized ${ZeroMQ_LIBRARY_RELEASE} optimized ws2_32.lib )
  #  target_link_libraries(libxrf_io LINK_PUBLIC debug ${ZeroMQ_LIBRARY_DEBUG} debug ws2_32.lib )
  #  target_link_libraries(xrf_maps LINK_PUBLIC optimized ${ZeroMQ_LIBRARY_RELEASE} optimized ws2_32.lib )
  #  target_link_libraries(xrf_maps LINK_PUBLIC debug ${ZeroMQ_LIBRARY_DEBUG} debug ws2_32.lib )
  ENDIF()
ELSEIF (UNIX)
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
  IF (BUILD_WITH_ZMQ)
    target_link_libraries(libxrf_io LINK_PUBLIC optimized ${ZeroMQ_LIBRARY} )
    target_link_libraries(libxrf_io LINK_PUBLIC debug ${ZeroMQ_LIBRARY} )
    target_link_libraries(xrf_maps LINK_PUBLIC optimized ${ZeroMQ_LIBRARY} )
    target_link_libraries(xrf_maps LINK_PUBLIC debug ${ZeroMQ_LIBRARY} )
  ENDIF()

  # Strip unnecessary sections of the binary on Linux/Mac OS
  IF(APPLE)
    set_target_properties(libxrf_fit libxrf_io PROPERTIES SUFFIX ".dylib")
    set_target_properties(libxrf_fit libxrf_io PROPERTIES COMPILE_FLAGS "-DDARWIN")
    set_target_properties(libxrf_fit libxrf_io PROPERTIES MACOSX_RPATH ".")
    set_target_properties(libxrf_fit libxrf_io PROPERTIES LINK_FLAGS "-undefined dynamic_lookup -dead_strip")
    IF (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
      add_custom_command(TARGET libxrf_fit POST_BUILD COMMAND strip -u -r ${PROJECT_SOURCE_DIR}/bin/libxrf_fit.dylib)
      add_custom_command(TARGET libxrf_io POST_BUILD COMMAND strip -u -r ${PROJECT_SOURCE_DIR}/bin/libxrf_io.dylib)
    ENDIF()
  ELSE()
    #if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    IF(NOT BUILD_FOR_PHI)
      IF(NOT STATIC_BUILD)
        set_target_properties(libxrf_fit libxrf_io PROPERTIES SUFFIX ".so")
        IF (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
          add_custom_command(TARGET libxrf_fit POST_BUILD COMMAND strip ${PROJECT_SOURCE_DIR}/bin/libxrf_fit.so)
          add_custom_command(TARGET libxrf_io POST_BUILD COMMAND strip ${PROJECT_SOURCE_DIR}/bin/libxrf_io.so)
        ENDIF()
      ENDIF()
    ENDIF()
  ENDIF()
ENDIF()

set(HDF5_LIB_LEN 0)
list(LENGTH HDF5_LIBRARIES HDF5_LIB_LEN)
# Building on theta doesn't need ${HDF5_LIBRARIES} and it is empty os need a check 
IF(${HDF5_LIB_LEN} LESS 1)
  target_link_libraries (xrf_maps LINK_PUBLIC optimized libxrf_io optimized libxrf_fit optimized ${NETCDF_LIBRARIES} yaml-cpp ${CMAKE_THREAD_LIBS_INIT} )
  target_link_libraries (xrf_maps LINK_PUBLIC debug libxrf_io debug libxrf_fit debug ${NETCDF_LIBRARIES} yaml-cpp ${CMAKE_THREAD_LIBS_INIT} )
ELSE()
  target_link_libraries (xrf_maps LINK_PUBLIC optimized libxrf_io optimized libxrf_fit optimized ${NETCDF_LIBRARIES} optimized ${HDF5_LIBRARIES} yaml-cpp ${CMAKE_THREAD_LIBS_INIT} )
  target_link_libraries (xrf_maps LINK_PUBLIC debug libxrf_io debug libxrf_fit debug ${NETCDF_LIBRARIES} debug ${HDF5_LIBRARIES} yaml-cpp ${CMAKE_THREAD_LIBS_INIT} )
ENDIF()

IF (BUILD_WITH_QT)
  target_link_libraries (libxrf_io LINK_PUBLIC optimized ${Qt5Charts_LIBRARIES} )
  target_link_libraries (libxrf_io LINK_PUBLIC debug ${Qt5Charts_LIBRARIES} )

  target_link_libraries (xrf_maps LINK_PUBLIC optimized ${Qt5Charts_LIBRARIES} )
  target_link_libraries (xrf_maps LINK_PUBLIC debug ${Qt5Charts_LIBRARIES} )

ENDIF()

IF (BUILD_WITH_TIRPC)
  link_directories(AFTER "/lib64" "/usr/lib")
  target_link_libraries (libxrf_io LINK_PUBLIC optimized libtirpc.so  )
  target_link_libraries (libxrf_io LINK_PUBLIC debug libtirpc.so)

  target_link_libraries (xrf_maps LINK_PUBLIC optimized libtirpc.so)
  target_link_libraries (xrf_maps LINK_PUBLIC debug libtirpc.so)
ENDIF()

