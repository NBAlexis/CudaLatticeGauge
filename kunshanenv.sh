module rm compiler/rocm/dtk-22.10.1
module load compiler/dtk/25.04
module load compiler/cmake/3.23.3
source /public/software/compiler/rocm/dtk-25.04/cuda/env.sh
export CXX=/public/software/compiler/rocm/dtk-25.04/llvm/bin/clang++
cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 -DCMAKE_CUDA_COMPILER_FORCED=TRUE



source /public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/env.sh
source /public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/cuda/env.sh
cd /public/home/nbalexis/CLGLib/Bin/Ubuntu
./RotationImprovedFermion



source /public/home/nbalexis/dtk-25.04.2-DCC2510-1215-2-centos7/env.sh
source /public/home/nbalexis/dtk-25.04.2-DCC2510-1215-2-centos7/cuda/env.sh
cd /public/home/nbalexis/CLGLib25042/Bin/Ubuntu
./RotationImprovedFermion