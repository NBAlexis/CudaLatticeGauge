# CudaLatticeGauge

***NOTE: This GitHub repository is essentially NO LONGER maintained. The only actively maintained repository is on Gitee: https://gitee.com/NBAlexis/CudaLatticeGauge — please use the Gitee one.***

lattice gauge (lattice QCD) simulation with CUDA and DCU now!</br>

Please report bug to yangjichong@lnnu.edu.cn

Thank you ^_^

---

## Wiki 参考

本项目维护了一份面向 LLM / vibe coding 的代码 wiki，由 Kimi 2.6 自动生成，覆盖项目架构、模块说明、关键 API 和常用模式。

> **快速入口**：[`wiki/home.md`](wiki/home.md)

主要内容包括：
- 项目整体架构和模块关系
- 每个模块的职责、主要类和关键 API
- CUDA Kernel 开发规范（强制规则）
- 设备端数学函数参考
- 如何添加新 Action / Measurement / Gauge Fixing 的指南

### YAML 配置字典

按 CLGManager 初始化顺序组织的**完整 YAML 参数参考**，覆盖所有可配置类及其参数：

> **配置字典入口**：[`wiki/yaml-reference/home.md`](wiki/yaml-reference/home.md)

包含 17 章：全局参数 → 格点 → 规范场 → 费米子场 → 作用量 → 求解器 → 测量器 → 更新器 → 积分器等。每章列出所有 YAML 键名、类型、默认值、来源文件及行号。

## 入门教程

面向新用户的完整入门教程，从零开始逐步讲解四种典型模拟场景：

> **教程入口**：[`Tutorial/Tutorial.md`](Tutorial/Tutorial.md)

涵盖内容：
- 纯 SU(3) 规范场的 HMC 模拟
- 纯 Z2 离散规范场的 Heatbath 模拟
- SU(3) + Wilson Dirac 费米子的 HMC 模拟
- SU(3) + Staggered (KS) 费米子的 HMC 模拟
- 对已保存组态的测量（Polyakov Loop / Wilson Loop / Chiral Condensate）
- 如何将自己的项目加入 CMake 编译

=================== Windows ========================

***CUDA***:

just build with Visual Studio

=================== Linux ========================

cd /Code/CMake

***Note, CMAKE_CUDA_ARCHITECTURES may need version 3.18 or above, for Ubuntu-20, pip install --upgrade cmake=3.18.4 is tested***

(before make, a "make clean" can be excute to do a re-build)

### CUDA:

export PATH=$PATH:/usr/local/cuda-10.0/bin

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-10.0/lib:/usr/local/lib

export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/cuda-10.0/include

(change /usr/local/cuda-10.0 to your path)

cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86

to enable applications, use -DCLG_application=1, where 'application'_is the name, for example,

cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86 -DCLG_ElectricChemical=1 -DCLG_StaggeredRotation=1 -DCLG_StaggeredSpectrum=1 -DCLG_WilsonDiracElectricMagnetic=1

to enable double-presicion float, use -DCLG_DOUBLE=1

### Multi-GPU (one MPI process per GPU)

Build with `-DCLG_MULTI_GPU=1` (requires MPI at configure time; the
single-GPU builds never need MPI). The lattice is split into sub-lattices
(GpuGrid in the yaml; the product of the grid factors MUST equal the number
of MPI processes). All MG tests live in their own group
(`Bin/Debug/TestSuit_MG.yaml`, `TestMG*` entries, n1 baseline with
`./CLGTest MG`); each runs `-n 1` vs `-n N` and the driver compares the printed
values (bitwise for MD5/accept sequences, DOUBLE tolerance for float-sum
reductions). Multi-rank runs go through the test driver (no hand-editing of
yaml grids):

```bash
./runtest.sh TestName GpuCount      # Windows: runtest.bat TestName GpuCount
```

The script reads `rankCount`/`GpuGrid` from the test's yaml metadata
(`MultiGPUTestRankCount`/`MultiGPUTestGrid`) via `CLGTest --mg-config`, then
launches `mpiexec -n rankCount CLGTest TestName --mg-worker --gpu-grid ...
--device-per-node GpuCount`. How to configure the split grid in the yaml
(`GpuGrid`/`HaloWidth`/`DevicePerNode`) is documented in
[wiki/multi-gpu.md](wiki/multi-gpu.md).

- **Windows (CUDA, VS)**: build with the `Release_MG`/`Debug_MG`
  configurations, then `runtest.bat <test-name> 2` (MS-MPI).
- **WSL / Linux (CUDA)**: build with
  `cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86 -DCLG_MULTI_GPU=1`,
  then `./runtest.sh <test-name> 1` (OpenMPI; a single visible GPU means all
  ranks oversubscribe device 0, which is how the 1-vs-N consistency tests run
  on a single card -- the worker report marks it as oversubscription).
- **DCU**: build with the DTK backend plus `-DCLG_MULTI_GPU=1`, then launch
  with `srun --gres=dcu:N -n N ./CLGTest <test-name>` on the cluster
  (GpuGrid product must equal N).

Known limitations: a t-direction Polyakov/Wilson-chain measurement needs a
grid that does NOT split t (use e.g. [1,1,2,1] for -n 2); split-z Z-slice
profiles and split-x/y XY distributions are explicitly rejected
(appCrucial) rather than silently wrong. See
`Docs/MultiGPU-Status-Roadmap.md` for the full status.

### DCU

#### DCU 25.04 

export CXX=path_to_clang_compilier (something like /public/software/compiler/rocm/dtk-25.04/llvm/bin/clang++)

cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 -DCMAKE_CUDA_COMPILER_FORCED=TRUE

make

then, in /Bin/Ubuntu/

use ./CLGTest to run

The '-DCMAKE_CUDA_COMPILER_FORCED=TRUE' is to skip Check for working cuda compilier which might be unnecessary for some version of cmake

some notes:

the '-DCMAKE_CUDA_COMPILER_FORCED=TRUE' work around does not work for dtk-25.04.1

something like this will work for 25.04:

  source /public/software/compiler/rocm/dtk-25.04/cuda/env.sh

  export CXX=/public/software/compiler/rocm/dtk-25.04/llvm/bin/clang++

  cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 -DCMAKE_CUDA_COMPILER_FORCED=TRUE

something like this will work for 25.04.1: (test use cmake 3.23, but should work for 3.18 and newer)

  source /public/home/nbalexis/dtk-25.04.1/env.sh

  source /public/home/nbalexis/dtk-25.04.1/cuda/cuda-12/env.sh

  export CUDAToolkit_ROOT=/public/home/nbalexis/dtk-25.04.1/cuda/cuda-12/bin

  export CXX=/public/home/nbalexis/dtk-25.04.1/llvm/bin/clang++

  cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 -DCMAKE_CUDA_COMPILER_FORCED=TRUE -DCMAKE_CUDA_COMPILER=/public/home/nbalexis/dtk-25.04.1/cuda/cuda-12/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=75

replace '/public/home/nbalexis/dtk-25.04.1' to your path


something like this will work for 25.04.1-2506:

  source /public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/env.sh

  source /public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/cuda/env.sh

  export CXX=/public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/llvm/bin/clang++

  export CUDAToolkit_ROOT=/public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/cuda/bin

  cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 -DCMAKE_CUDA_COMPILER_FORCED=TRUE -DCMAKE_CUDA_COMPILER=/public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=75

replace '/public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7' to your path

#### DCU 25.04.2

for 'dtk-25.04.2-DCC2510-1215-2'

  module load compiler/cmake/3.23.3

  source /public/home/nbalexis/dtk-25.04.2-DCC2510-1215-2-centos7/env.sh

  source /public/home/nbalexis/dtk-25.04.2-DCC2510-1215-2-centos7/cuda/env.sh

  export CXX=/public/home/nbalexis/dtk-25.04.2-DCC2510-1215-2-centos7/llvm/bin/clang++

  export CUDAToolkit_ROOT=/public/home/nbalexis/dtk-25.04.2-DCC2510-1215-2-centos7/cuda/bin

  cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 -DCMAKE_CUDA_COMPILER_FORCED=TRUE -DCMAKE_CUDA_COMPILER=/public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=75 -DCLG_RotationImprovedFermion=1

note: the link cost a lot of time, it can only be done on a node

  salloc --nodes=1 --ntasks=1 --exclusive --time=06:00:00 --partition=some_cpu_node

  echo $SLURM_JOB_NODELIST

  ssh xxx (switch to that node)

  build the code

After swith to node, before build, you still need to set the environment variables (only need to run the two env.sh files if cmake is already done.).

#### DCU 26.04

    module load compiler/cmake/3.23.3 && source /public/home/nbalexis/dtk-26.04/env.sh && source /public/home/nbalexis/dtk-26.04/cuda/env.sh && export CXX=/public/home/nbalexis/dtk-26.04/llvm/bin/clang++ && export CUDAToolkit_ROOT=/public/home/nbalexis/dtk-26.04/cuda/cuda/bin
    cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 -DCMAKE_CUDA_COMPILER_FORCED=TRUE -DCMAKE_CUDA_COMPILER=/public/home/nbalexis/dtk-26.04/cuda/cuda/bin/nvcc -DCMAKE_CUDA_ARCHITECTURES=75 -DCLG_RotationImprovedFermion=1

#### DCU RUN

to run:


  source /public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/env.sh

  source /public/home/nbalexis/dtk-25.04.1-DCC-2506-0730-centos7/cuda/env.sh

  srun --pty -n 1 -p partition_name --gres dcu:1 ./CLGTest

replace './CLGTest' to the path of application you want to run

submit job:

source /public/home/nbalexis/dtk-25.04.2-DCC2510-1215-2-centos7/env.sh
source /public/home/nbalexis/dtk-25.04.2-DCC2510-1215-2-centos7/cuda/env.sh
cd /public/home/nbalexis/CLGLib/Bin/Ubuntu
./RotationImprovedFermion
#SBATCH --gres=dcu:1