# CudaLatticeGauge

lattice gauge simulation with CUDA

Testing on GTX-1060 (6G), WIN10

Using GMRES with 10 orthrognal basis, maximum supported (32 x 32 x 32 x 16)
For the release built
One trajectory with 6/48 Multi-rate approximate Force-gradiant steps cost about 230 secs

(The speed is related with kappa and beta, it is tested using kappa=0.1355 and beta=2.5, 
GMRES with 10 orthrognal basis will converge with about 4 restarts)

It is tested on a Notebook with two GPUs, Intel(R) HD Graphic 630 for dispalyer and GTX1060 for calculation.

This is important because when connecting an additional displayer, the speed is significantly slower.

See detailed.pdf for more.

Please report bug to yangjichong@fudan.edu.cn

Thank you ^_^

=================== All ==========================

Note: It is necessary to update the drivers for graphic card to support CUDA 10!

=================== WIN8 =========================

The default Code generation is "compute_61,sm_61";. To test on a GTX970m, please add "compute_52,sm_52;" . (NOTE, to support double float, must use a sm60+.)

Tested on a WIN8 without CUDA installed, and with GTX970m.

NOTE: if CUDA is not installed, the curand64_100.dll file is necessary.

=================== Ubuntu ========================

NOTE: The Linux environment is the run-time environment for us. 
Now the default is DOUBLE-precision float, and architecture compute_61,sm_61.

put

CUDACXX=/usr/local/cuda-10.0/bin/nvcc

into /etc/environment so that it is applied when using sudo. (This is for CMake to recognize CUDA)

put

export PATH=$PATH:/usr/local/cuda-10.0/bin

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.0/lib:/usr/local/lib

export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/cuda-10.0/include

into ~/.bashrc

(if bashrc is changed, reboot is necessary. This is for cuda include dirs, when building .cu, nvcc can add cuda headers automatically, but when building .cpp files, it is necessary)

cd /Code/CMake

(before make, a "make clean" can be excute to do a re-build)

cmake CMakeLists.txt

make

then, in /Bin/UbuntuDebug/

use ./CLGTest to run

Testing using Ubuntu18 + GCC 7.3.0 + CMake 3.10 + Cuda-10.0 installed.

(NOTE: Run directly using built binaries on a machine without Cuda, is not tested. One may need curand.so file)

