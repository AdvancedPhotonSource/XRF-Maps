# XRF-Maps Compiling
## Windows
### MSVC

#### Prepare environment
* git clone https://github.com/advancedphotonsource/XRF-Maps.git
* download and install cmake https://cmake.org/download/ 
* download and install netcdf 3 library https://cmake.org/download/
* download and install hdf5 library https://support.hdfgroup.org/HDF5/release/obtain5.html
* download eigen 3 library http://eigen.tuxfamily.org/index.php?title=Main_Page
* (optional) git clone https://github.com/zeromq/libzmq.git


## Linux
### No Streaming
* mkdir build
* cd build
* CC=icc CXX=icpc cmake -DEIGEN3_INCLUDES=/usr/common/software/eigen3/3.3.3/include/eigen3/ ../
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
* module add gcc/6.3.0
* module add eigen3/3.3.3
* module add netcdf/4.4.1
* export PATH=/usr/common/software/hdf5-serial/1.8.16/hsw/intel/bin/:$PATH
* mkdir build
* cd build
* CC=icc CXX=icpc cmake -DBUILD_FOR_PHI=ON -DEIGEN3_INCLUDES=/usr/common/software/eigen3/3.3.3/include/eigen3/ ../
* cmake --build . --config Release
