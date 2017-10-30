--CORI Nersc--
---Intel Phi---
--------------
module rm darshan/3.1.4
module swap craype-haswell craype-mic-knl
module add craype-mic-knl
module add gcc/6.3.0
module add eigen3/3.3.3
module add netcdf/4.4.1
export PATH=/usr/common/software/hdf5-serial/1.8.16/hsw/intel/bin/:$PATH
mkdir build
cd build
CC=icc CXX=icpc cmake -DBUILD_FOR_PHI=ON -DEIGEN3_INCLUDES=/usr/common/software/eigen3/3.3.3/include/eigen3/ ../
cmake --build . --config Release
