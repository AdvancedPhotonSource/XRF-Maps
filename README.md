[![Actions Status](https://github.com/aglowacki/XRF-Maps/workflows/CMake/badge.svg)](https://github.com/aglowacki/XRF-Maps/actions)
[![Build Status](https://dev.azure.com/aglow/XRF-Maps/_apis/build/status/aglowacki.XRF-Maps?branchName=master)](https://dev.azure.com/aglow/XRF-Maps/_build/latest?definitionId=2&branchName=master)

# XRF-Maps

X-ray fluorescence (XRF) imaging typically involves the creation and analysis of 3D data sets, where at each scan position a full energy dispersive x-ray spectrum is recorded. This allows one to later process the data in a variety of different approaches, e.g., by spectral region of interest (ROI) summation with or without background subtraction, principal component analysis, or fitting. XRF-Maps is a C++ open source software package that implements these functions to provide a tool set for the analysis of XRF data sets. It is based on the MAPS software http://www.aps.anl.gov/Xray_Science_Division/Xray_Microscopy_and_Imaging/Software_and_Tools/maps.html

# Compiling

## Requires

Visual Studio 2015 or greater (Windows build)
GCC 6.0 or greater (Linux build)
Cmake 3.5 or greater

## Libraries

HDF5 : https://www.hdfgroup.org/downloads/
NetCDF : http://www.unidata.ucar.edu/downloads/netcdf/index.jsp (Download http://www.unidata.ucar.edu/software/netcdf/docs/winbin.html)
Eigen : http://eigen.tuxfamily.org/index.php?title=Main_Page (submodule in src/support)

## Optional Libraries

QT : https://www.qt.io/download
ZeroMQ : http://zeromq.org/area:download

## Compile Default

1) git clone --recurse-submodules https://github.com/AdvancedPhotonSource/XRF-Maps.git
2) cd XRF-Maps
3) mkdir build
4) cd vcpkg
5) vcpkg set Linux
   1) ./bootstrap-vcpkg.sh
   2) ./vcpkg install hdf5 netcdf-c yaml-cpp zeromq
6) vcpkg setup windows
   1) .\bootstrap-vcpkg.bat
   2) .\vcpkg install hdf5 netcdf-c yaml-cpp zeromq --triplet x64-windows
7) cd ../build
8) cmake `-DCMAKE_TOOLCHAIN_FILE=../vcpkg/scripts/buildsystems/vcpkg.cmake -DBUILD_WITH_ZMQ=ON ..`
9) make

## Building optional features

### QT support

-DBUILD_WITH_QT=ON

### Streaming support with ZeroMQ

-DBUILD_WITH_ZMQ=ON -DZeroMQ_INCLUDE_DIR={path to include} -DZeroMQ_LIBRARY_RELEASE={path to zmq.so/zmq.dll}

### Python bindings (NOTE: this may interfere with QT options if QT lib version is different than python qt lib version as with anaconda python)

-DBUILD_WITH_PYBIND11=ON
-DPYTHON_EXECUTABLE={path to python.exe if not found}

### Build for intel phi

-DBUILD_FOR_PHI=ON

### Set double presision instead of single

-DDOUBLE_PREC=ON

### Can't auto find Libraries

#### HDF5

-DHDF5_ROOT={path to hdf5 root dir}

#### NetCDF

-DNETCDF_INCLUDE_DIRS={path} -DNETCDF_LIBRARIES={path including .so/.dll}

Libraries and executable stored in bin directory

# Data Analysis

## General
Run from bin directory. The software looks for references one directory up ( ../references ). 
```bash
Help:
Usage: xrf_maps [Options] --dir [dataset directory]

Options:
--nthreads : <int> number of threads to use (default is all system threads)
--quantify-with : <standard.txt> File to use as quantification standard
--detectors : <int,..> Detectors to process, Defaults to 0,1,2,3 for 4 detector
--generate-avg-h5 : Generate .h5 file which is the average of all detectors .h50 - h.53 or range specified.
--add-v9layout : Generate .h5 file which has v9 layout able to open in IDL MAPS software.
--add-exchange : Add exchange group into hdf5 file with normalized data.
--export-csv : Export Integrated spec, fitted, background to csv file.
--update-theta : <theta_pv_string> Update the theta dataset value using theta_pv_string as new pv string ref.
--update-amps <us_amp>,<ds_amp>: Updates upstream and downstream amps if they changed inbetween scans.
--update-quant-amps <us_amp>,<ds_amp>: Updates upstream and downstream amps for quantification if they changed inbetween scans.
--quick-and-dirty : Integrate the detector range into 1 spectra.
--optimize-fit-override-params : <int> Integrate the 8 largest mda datasets and fit with multiple params.
  1 = matrix batch fit
  2 = batch fit without tails
  3 = batch fit with tails
  4 = batch fit with free E, everything else fixed
--optimizer <lmfit, mpfit> : Choose which optimizer to use for --optimize-fit-override-params or matrix fit routine
Fitting Routines:
--fit <routines,> comma seperated
  roi : element energy region of interest
  roi_plus : SVD method
  nnls : Non-Negative Least Squares
  tails : Fit with multiple parameters
  matrix : Fit with locked parameters

Dataset:
--dir : Dataset directory
--files : Dataset files: comma (',') separated if multiple
Network:
--streamin [source ip] : Accept a ZMQ stream of spectra to process. Source ip defaults to localhost (must compile with -DBUILD_WITH_ZMQ option)
--streamout : Streams the analysis counts over a ZMQ stream (must compile with -DBUILD_WITH_ZMQ option)

Examples:
   Perform roi and matrix analysis on the directory /data/dataset1
xrf_maps --fit roi,matrix --dir /data/dataset1
   Perform roi and matrix analysis on the directory /data/dataset1 but only process scan1 and scan2
xrf_maps --fit roi,matrix --dir /data/dataset1 --files scan1.mda,scan2.mda
   Perform roi, matrix, and nnls  analysis on the directory /data/dataset1, use maps_standard.txt information for quantification
xrf_maps --fit roi,matrix,nnls --quantify-with maps_standard.txt --dir /data/dataset1
```