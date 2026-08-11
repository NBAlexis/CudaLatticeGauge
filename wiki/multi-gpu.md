# Multi-GPU 配置指南

> 本文档说明如何在 CLGLib 中配置并运行多 GPU（多 MPI 进程）模拟。
> 实现跟踪与任务清单见仓库根目录 [`multi-gpu.md`](../multi-gpu.md)（执行计划）
> 与 [`Docs/MultiGPU-Status-Roadmap.md`](../Docs/MultiGPU-Status-Roadmap.md)（状态路线图）。

## 1. 概述

Multi-GPU 支持把 4D 格点按进程网格切分为多个**子格点**，每个 MPI 进程（rank）持有
一个子格点并驱动一块 GPU。跨进程的邻居数据通过 **halo 交换**（`CLGComm::RefillHalo`）
同步。

- 构建宏：`_CLG_MULTI_GPU`（CMake 传 `-DCLG_MULTI_GPU=1`），要求配置期 MPI 可用；
- **单 GPU 构建（不加该宏）不需要 MPI**，且从不包含任何多 GPU 代码路径；
- 一个进程对应一块 GPU（DCU 集群）或单卡 oversubscribe（一致性测试）；
- 切分方式由 YAML 的 `GpuGrid` 指定，**不是**由 MPI 启动器隐式决定。

## 2. 构建

### Linux / WSL（CUDA）

**前置依赖**：多 GPU 构建在 CMake 配置期要求 MPI（`find_package(MPI REQUIRED)`）
与 `cmake`（≥ 3.17）。单 GPU 构建两者都不需要。

```bash
# Ubuntu / Debian（有 sudo）：
sudo apt install cmake openmpi-bin libopenmpi-dev

# 免 root（集群 / 无 sudo 的 WSL2）：把 cmake 官方二进制与源码编译的
# OpenMPI 装到仓库 Libs/ 目录（已 gitignore），用时加进 PATH：
#   Libs/cmake/bin/cmake, Libs/openmpi/bin/{mpicxx,mpiexec}
```

```bash
cd Code/CMake
# 架构按显卡填：RTX 30 系 86，RTX 40 系 89
cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=89 -DCLG_MULTI_GPU=1
make -j
```

- 不带 `-DCLG_MULTI_GPU=1` 即为普通单 GPU 构建；
- 加 `-DCLG_DOUBLE=1` 为双精度多 GPU 构建（一致性测试常用）。

**WSL2 已知坑（libcuda 失配）**：Windows 侧驱动升级后，运行中的 WSL 会话里
`/usr/lib/wsl/lib/libcuda.so.1` 可能仍是旧版 shim（`nvidia-smi` 头部显示
`CUDA UMD Version: N/A`），CUDA 程序会段错误或报
`cudaErrorInsufficientDriver`。永久修复：Windows PowerShell 执行
`wsl --shutdown` 后重开会话（重新挂载新驱动库）。临时绕过：活动驱动仓库目录
（`/usr/lib/wsl/drivers/nv_dispi.inf_*/`，9p 挂载、实时反映 Windows 侧驱动）
里有完整的 `libcuda.so.1.1`，但缺加载器要找的 `libcuda.so.1` 软链，自建一个
指向它的目录并放到 `LD_LIBRARY_PATH` 最前即可：

```bash
D=$(ls -d /usr/lib/wsl/drivers/nv_dispi.inf_*/ | head -1)   # 活动驱动仓库
mkdir -p Libs/wsl-cuda
ln -sf "${D}libcuda.so.1.1" Libs/wsl-cuda/libcuda.so.1
export LD_LIBRARY_PATH="$PWD/Libs/wsl-cuda:$D:$LD_LIBRARY_PATH"
```

真实多卡运行（2 张 GPU 的机器）：

```bash
./runtest.sh TestName 2     # 单个 MG 测试：每 rank 独占一张物理 GPU
./runall-mg.sh 2            # 批量跑 TestSuit_MG.yaml 全部 MG 测试并汇总 PASS/FAIL
```

`runall-mg.sh`（仓库根目录）从 `Bin/Debug/TestSuit_MG.yaml` 解析每个测试的
`MultiGPUTestRankCount`，逐个调 `runtest.sh`，日志写到
`Bin/Ubuntu/mg_logs_<时间戳>/`，末尾打印 `SUMMARY: PASS=n FAIL=m`；
`GpuCount` 超过某测试 rankCount 时自动按 rankCount 截断（单卡机传 `1` 即全部
oversubscribe）。

### Windows（CUDA + VS）

用 Visual Studio 打开 `Code/CudaLatticeGauge.sln`，配置选择 `Release_MG` / `Debug_MG`
（已带 `_CLG_MULTI_GPU` 宏）。运行需 MS-MPI：

```bat
mpiexec -n 2 CLGTest.exe <test-name>
```

### DCU 集群

```bash
export CXX=path_to_clang_compiler   # DTK 自带 clang++
cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 -DCLG_MULTI_GPU=1 -DCMAKE_CUDA_COMPILER_FORCED=TRUE
make
# 运行时一个 rank 占一块 DCU
srun --gres=dcu:N -n N ./CLGTest <test-name>
```

## 3. YAML 配置

三个键控制多 GPU 行为，均在 YAML **顶层**：

| 键名 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `GpuGrid` | INT[]（4 个） | `[1,1,1,1]` | 进程网格 `[Gx, Gy, Gz, Gt]`，**乘积必须等于 MPI 进程数 N** |
| `HaloWidth` | INT | `2` | 单方向最远 stencil 半径（HISQ Naik 需 2） |
| `DevicePerNode` | INT | 机器 GPU 数 | 每个节点的 GPU 数；`rank % DevicePerNode` 选设备。单卡测试机设为 `1` 让所有 rank 共享设备 0（oversubscribe） |

### 示例 1：8⁴ 格点，2 进程切 t

```yaml
LatticeLength : [8, 8, 8, 8]
GpuGrid : [1, 1, 1, 2]     # 2 个 rank，每 rank 本地 8×8×8×4
HaloWidth : 2
DevicePerNode : 1          # 单卡测试：两个 rank 都跑在设备 0
```

```bash
mpiexec --oversubscribe -n 2 ./CLGTest <test-name>
```

### 示例 2：16³×8 格点，4 进程二维切分

```yaml
LatticeLength : [16, 16, 16, 8]
GpuGrid : [2, 2, 1, 1]     # 乘积 4 == -n 4；本地 8×8×16×8
```

```bash
mpiexec --oversubscribe -n 4 ./CLGTest <test-name>
```

### 示例 3：DCU 集群 2 卡

```yaml
LatticeLength : [16, 16, 16, 16]
GpuGrid : [1, 1, 2, 1]     # 切 z；t 链测量不切 t（见 §6）
```

```bash
srun --gres=dcu:2 -n 2 ./CLGTest <test-name>
```

> **注意**：`GpuGrid` 与 `mpiexec -n N` 必须匹配，否则 `CLGComm::SetGpuGrid`
> 会 `appCrucial` 报错退出（见 §5 错误表）。MG 测试（`TestSuit_MG.yaml`）**不再**
> 手改 yaml 跑多进程：每个开放测试段声明 `MultiGPUTestRankCount` / `MultiGPUTestGrid`
> 元数据，由 `runtest.sh TestName GpuCount` 读取并经 `--mg-worker` 注入（见 §7）；
> 测试 yaml 中也不再固定 `DevicePerNode`（它等于 `runtest` 的第二个参数）。

## 4. 切分语义

- **均匀切分**：每方向本地长度 = 全局长度 / grid 分量，目前只支持整除的均匀切分；
- **rank 编号约定**：`GpuGrid` 按 [x, y, z, t]，rank 编号 t 变化最快（与格点 site 索引
  约定一致，见 `Docs/MultiGPU-Plan.md` §9.2）；
- **全局坐标**：位置相关物理（旋转/加速度/cylinder 作用量、位置相关测量）使用全局
  坐标，经 `_HC_GlobalLx...` + `_HC_OffsetX...` 恢复；`_deviceLocalToGlobalInt4` /
  `_deviceSiteIndexToGlobalInt4` 为 kernel 辅助；
- **halo 布局**：face → edge → corner 三层槽位（`CLGHaloLayout.h`），
  `RefillHalo` 完成跨 rank 交换；单方向切分时 edge/corner 槽位数为 0。

## 5. 校验规则与常见错误

`CLGComm::SetGpuGrid`（`Code/CLGLib/Core/Distributed/CLGComm.cpp`）在初始化时校验，
任一失败都会 `appCrucial` 并终止（宁可失败也不静默算错）：

| 错误信息 | 原因 | 修复 |
|----------|------|------|
| `GpuGrid product %d does not match MPI world size %d` | `GpuGrid` 乘积 ≠ `-n N` | 调整网格使乘积 == N |
| `lattice length %d in direction %d is not divisible by GpuGrid %d` | 该方向全局长度不能被 grid 分量整除 | 换可整除的网格或格点 |
| `local length %d in split direction %d is smaller than halo width %d` | 被切方向本地长度 < HaloWidth | 减小 grid 分量或增大格点 |
| `GpuGrid[%d] = %d must be >= 1` | 分量 < 1 | 使用 ≥ 1 的分量 |

## 6. 已知限制

- **t 方向链测量**（Polyakov / Wilson 链）要求网格**不切 t**（如 `-n 2` 用 `[1,1,2,1]`）；
- **split-z Z-slice 剖面**与 **split-x/y XY 分布**类测量被显式拒绝（`appCrucial`），
  不会静默算错；
- 非规则跨 rank 拼接（如任意路径 Wilson loop）、分布式 FFT 等暂不支持，
  见 `Docs/MultiGPU-Status-Roadmap.md`。
- **Dirichlet 边界 gauge 场暂不可分解**（I10 门槛实测发现）：分解后每 rank 的
  Dirichlet 索引/link 表按本地子格点重新收缩（real link 表 n1 为 0、n4 每 rank
  1024），物理边界平面附近的 force 与 n1 不一致，1-vs-N MD5 无法相同。Torus 与
  ProjectivePlane 均已逐 bit 验证；Dirichlet 的 MG 支持是后续修复项。

## 7. 一致性验证（1-vs-N）

MG 一致性测试（`TestSuit_MG.yaml`，`./CLGTest MG`）的核心判据：
**同一配置在 `mpiexec -n 1`（`GpuGrid [1,1,1,1]` 恒等分解）与 `mpiexec -n N`
下运行，输出标量 / MD5 / 接受序列一致**。

- 整数 / MD5 / 接受序列：逐比特一致；
- 浮点求和（double）：相对容差 1e-10；
- 单卡验证用 `runtest.sh TestName 1`（worker 报告会标记 single-GPU oversubscription）；
- 新建 MG 测试的检查清单见 `Docs/MultiGPU-Status-Roadmap.md` 附录（halo 读集 /
  写集 handle 注册、halo 宽度取 stencil 半径上限）。

### CLGTest 运行模式与 `runtest`

CLGTest 把交互选择器与 MPI worker 严格分开，共四种固定模式：

| 命令 | 语义 |
|------|------|
| `CLGTest`（无参数） | 单进程交互菜单（严格 single-GPU，无 `g` 切换）。`mpiexec -n N`（N>1）启动时所有 rank 在进入 stdin 读取前一致退出（非零），仅 rank 0 报错 |
| `CLGTest TestName` | 单进程 batch，identity `GpuGrid=[1,1,1,1]`，即 n1 baseline |
| `CLGTest --mg-config TestName` | 只读测试元数据：stdout 恰好一行 `rankCount gx gy gz gt`；日志走 stderr；未知 / 非 `_TEST_MULTIGPU` / 元数据非法均非零退出 |
| `CLGTest TestName --mg-worker --gpu-grid gx,gy,gz,gt --device-per-node G` | mpiexec 下的无交互 worker：校验 world size == 网格乘积、`G>=1`、每节点可见 GPU 数 >= G（CUDA API 实测），任一失败 collective 同码退出；测试结束 collective 汇总 errors/fatals/skip，仅 rank 0 输出摘要 |

驱动脚本（仓库根目录，语义一致的两个实现）：

```bash
./runtest.sh TestName GpuCount     # Linux；runtest.bat TestName GpuCount（Windows）
```

- `rankCount` 与 `GpuGrid` 只从 `--mg-config` 读取，脚本不维护 test-name 表；
- 校验 `GpuCount` 为正整数且 `GpuCount <= rankCount`；
- 执行 `mpiexec -n rankCount CLGTest TestName --mg-worker --gpu-grid ... --device-per-node GpuCount`，
  原样传播退出码；`CLGTEST_BIN` / `MPIEXEC` 环境变量可覆盖可执行文件位置；
- 设备映射固定 `deviceId = globalRank % GpuCount`（单节点）；多节点 world 被明确拒绝；
- 只有 1 张可见 GPU 时用 `runtest.sh TestName 1`：全部 rank 映射到设备 0，报告标记
  **oversubscription**（验证分解/halo/collective，不代表真实跨设备验证）；`GpuCount 2`
  在单卡机上会在跑任何 kernel 前非零退出。

批量回归用仓库根目录的 `runall-mg.sh GpuCount`：解析 `TestSuit_MG.yaml` 中每个测试的
`MultiGPUTestRankCount`，逐个经 `runtest.sh` 运行（`GpuCount` 超 rankCount 自动截断），
日志存 `Bin/Ubuntu/mg_logs_<时间戳>/`，末尾汇总 `PASS/FAIL`。真实多卡验证记录：
2026-08-11，2 × RTX 4090（WSL2，CUDA 13.3 + OpenMPI 4.1.6，`CLG_GPU_ARCH=89`），
`./runall-mg.sh 2` 全部 52 个 MG 测试 PASS（47 个 rankCount=2 每 rank 独占一卡、
4 个 rankCount=4、1 个 rankCount=16）。

### 测试注册约定

- 默认所有测试 **single-GPU only**；只有用 `_TEST_MULTIGPU` 显式注册（列表显示 `M`
  前缀）且 yaml 段补齐 `MultiGPUTestRankCount` / `MultiGPUTestGrid` 元数据的测试才会被
  `--mg-config` / `--mg-worker` 接受；
- 元数据约束：`MultiGPUTestGrid` 四个因子均 ≥ 1 且乘积 == `MultiGPUTestRankCount`；
- MG 一致性测试注册在统一测试列表中（无 `#if _CLG_MULTI_GPU` 守卫），单 GPU 构建 /
  普通 batch 下以恒等分解运行。

## 8. Halo guard 自动协议与例外清单

普通 kernel 不需要任何 MG 标注（写法见 [kernel-macros.md](kernel-macros.md)
「Multi-GPU」章）：`_LAUNCH_KERNEL*` 在 MG 构建下由 registry 按浅层 pointer 参数
发现已注册 handle，launch 前统一 `Ensure` 读集 halo、`Commit` invalidate 写集。
以下访问无法被浅层参数检查发现，按下表例外分类显式
处置；普通直接 pointer launch **不在**此清单内。

### 例外清单（最终）

| 类别 | 位置 | 处置 |
|------|------|------|
| 3.4.1 device pointer 数组 | `Data/Field/Staggered/CFieldFermionKST.h:95-99`（rational 场指针数组）；`Update/CStapleCache.cpp:101,118`（staple/plaq 指针数组） | owner 注册 `CHaloBufferSetHandle`，guard 命中外层指针时展开子 handle 集合统一 Ensure/invalidate |
| 3.4.1 审计豁免（不注册） | `Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.cu:1718`；`Data/Field/Staggered/CFieldFermionKST.h:1249` | 指针表仅被 kernel 当 rank 局部读写，无跨 rank 语义 |
| 3.4.2 隐藏/global/constant 指针 | `__boundaryFieldPointers`（`Core/CudaHelper.cu:32,1245`）；gauge fixing 全局 context（`GaugeFixing/CGaugeFixing.cpp:130` + 12 处调用点） | 最接近 launch 的 host helper 显式声明；禁止按 field ID 猜测 |
| 3.4.3 动态 stencil | `Data/Lattice/CIndexSquare.cu:785`（`CIndexSquare::BakeNaikTable`，Naik 3-hop） | cache owner 验证最大绝对坐标位移，超出 `HaloWidth` 在 bake 时 fail-fast |
| 3.4.4 rank 非对称 launch | 无生产调用点 | `_LAUNCH_KERNEL_RANK_LOCAL` 机制就绪（HaloCapable 命中即 fail-fast）；gauge fixing 需求由 3.4.2 吸收 |
| 附加：kernel body MG 分支（审计豁免） | `Tools/Math/DeviceTemplates/DeviceInlineGaugeChair.h:274-286`（1 处 device inline）；边界条件 bake `__global__` kernel 19 处（TorusSquare 7、PeriodicAndDirichlet 6、ProjectivePlane 6） | 普通 kernel body 不得新增 MG 条件分支；既有分支逐项审计豁免 |

## 9. 相关文档

- YAML 参数字典：[yaml-reference/02-lattice.md](yaml-reference/02-lattice.md)（`GpuGrid`/`HaloWidth`）、
  [yaml-reference/01-global.md](yaml-reference/01-global.md)（`DevicePerNode`）
- 执行计划：[multi-gpu.md](../multi-gpu.md)
- 状态路线图：[Docs/MultiGPU-Status-Roadmap.md](../Docs/MultiGPU-Status-Roadmap.md)
- README 构建说明：[README.md](../README.md#multi-gpu-one-mpi-process-per-gpu)
