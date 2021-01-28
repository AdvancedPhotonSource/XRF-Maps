# XRF-Maps Compiling
## Windows
### MSVC

#### Prepare environment
* git clone --recurse-submodules https://github.com/advancedphotonsource/XRF-Maps.git
* download and install cmake https://cmake.org/download/ 
* download and install netcdf 3 library https://cmake.org/download/
* download and install hdf5 library https://support.hdfgroup.org/HDF5/release/obtain5.html
* (optional) git clone https://github.com/zeromq/libzmq.git


## Linux
### No Streaming
* mkdir build
* cd build
* CC=icc CXX=icpc cmake ../
* cmake --build . --config Release

### Streaming
* mkdir build
* cd build
* CC=icc CXX=icpc cmake -DBUILD_WITH_ZMQ -DZeroMQ_INCLUDE_DIR=/libs/zmq/include/ -DZeroMQ_STATIC_LIBRARY=/libs/zmq/lib/libzmq.a -DEIGEN3_INCLUDES=/usr/common/software/eigen3/3.3.3/include/eigen3/ ../
* cmake --build . --config Release

## Linux HPC
### Nersc Cori
#### Intel Phi
* module swap craype-haswell craype-mic-knl
* module add gcc/9.3.0
* module add netcdf/4.6.1
* module add hdf5/1.10.1
* mkdir build
* cd build
* CC=icc CXX=icpc cmake -DBUILD_FOR_PHI=ON -DEIGEN3_INCLUDES=/usr/common/software/eigen3/3.3.3/include/eigen3/ ../
* cmake --build . --config Release
### ALCF Theta
#### Intel Phi
* module add cmake/3.18.0
* module add cray-hdf5-parallel/1.10.6.0
* module add cray-netcdf-hdf5parallel/4.7.3.3
* mkdir build
* cd build
* cmake -DBUILD_FOR_PHI=ON ../
* cmake --build . --config Release
