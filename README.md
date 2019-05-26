# CudaLatticeGauge

lattice gauge simulation with CUDA</br>

Test on laptops using GMRES with 12 orthrognal basis</br>
One trajectory with 8/48 Multi-rate approximate Force-gradiant steps

+--------------------------+---------------------------+---------------------------+</br>
| WIN8 + i7-4720HQ         | WIN10  + i7-7700HQ        |   Ubuntu 18.04 + i5-8400  |</br>
| + GTX970M (3G)           | + GTX1060 (6G,mobile)     |   + GTX1070 (8G, mobile)  |</br>
| Intel(R) HD Graphic 4600 | Intel(R) HD Graphic 630   |                           |</br>
+--------------------------+---------------------------+---------------------------+</br>
|       (32x32x32x12)      |      (32x32x32x16)        |      (32x32x32x16)        |</br>
+--------------------------+---------------------------+---------------------------+</br>
|      compute_52,sm_52    |     compute_61,sm_61      |     compute_61,sm_61      |</br>
+--------------------------+---------------------------+-------------+-------------+</br>
|      Single float        |        Single float       |    Single   |   Double    |</br>
+--------------------------+---------------------------+-------------+-------------+</br>
|         498 secs         |         286 secs          |    214 secs |   404 secs  |</br>
+--------------------------+---------------------------+-------------+-------------+</br>

For 32x32x32x16, H_diff = 0.32 ( Metropolic accept rate > exp(-0.32) ~ 72.6% )</br>
(The speed is related with kappa and beta, it is tested using kappa=0.1355 and beta=2.5, 
GMRES with 12 orthrognal basis will converge with about 4 restarts)

The 32x32x32x12 is the maximum size of lattice using float and 3GB GPU.</br>
For windows, it is tested on a Notebook with two GPUs, Intel(R) HD Graphic for dispalyer and GTX for calculation.</br>
For 32x32x32x12 or 32x32x32x16, the energy of the Hamitonian is too large, 
and the accuracy of single float can only reach about 1E-7 which is about 2.1 >> 0.3.</br>
So H_diff = 0.32 is tested using double float.

See detailed.pdf for more.

Please report bug to yangjichong@fudan.edu.cn

Thank you ^_^

=================== All ==========================

Note: It is necessary to update the drivers for graphic card to support CUDA 10!

=================== WIN8 =========================

The default Code generation is "compute_61,sm_61";. To test on a GTX970m, please add "compute_52,sm_52;" . (NOTE, to support double float, must use a sm60+.)

Tested on a WIN8 without CUDA installed, and with GTX970m.

NOTE: if CUDA is not installed, the curand64_10.dll file is necessary.

=================== Ubuntu ========================

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

