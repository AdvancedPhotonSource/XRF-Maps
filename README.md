[![Build Status](https://travis-ci.org/aglowacki/XRF-Maps.svg?branch=master)](https://travis-ci.org/aglowacki/XRF-Maps)
[![Actions Status](https://github.com/aglowacki/XRF-Maps/workflows/cmake/badge.svg)](https://github.com/aglowacki/XRF-Maps/actions)
[![Build Status](https://dev.azure.com/aglow/XRF-Maps/_apis/build/status/aglowacki.XRF-Maps?branchName=master)](https://dev.azure.com/aglow/XRF-Maps/_build/latest?definitionId=2&branchName=master)

# XRF-Maps

X-ray fluorescence (XRF) imaging typically involves the creation and analysis of 3D data sets, where at each scan position a full energy dispersive x-ray spectrum is recorded. This allows one to later process the data in a variety of different approaches, e.g., by spectral region of interest (ROI) summation with or without background subtraction, principal component analysis, or fitting. XRF-Maps is a C++ open source software package that implements these functions to provide a tool set for the analysis of XRF data sets. It is based on the MAPS software http://www.aps.anl.gov/Xray_Science_Division/Xray_Microscopy_and_Imaging/Software_and_Tools/maps.html

# Prebuilt Binaries
## Windows 64 (Visual C++ 2015 redis)
### Single precision with QT, and ZeroMQ
 https://jenkins.aps.anl.gov/view/XRF_Maps/job/XRF_Maps_Win64_QT/lastSuccessfulBuild/artifact/*zip*/archive.zip

## RedHat Linux 7 (gcc 6.2)
### Single precision with ZeroMQ, and Anaconda Python 3.6 bindings
 https://jenkins.aps.anl.gov/view/XRF_Maps/job/XRF_Maps_Linux64/lastSuccessfulBuild/artifact/*zip*/archive.zip

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
 3) mkdir build ; cd build
 4) cmake ../
 5) make

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
