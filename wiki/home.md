# CLGLib Wiki

> 本 wiki 是 CudaLatticeGauge（CLGLib）的 LLM / 开发者索引，供快速定位代码、理解模块关系、查阅开发规范。
> 项目地址：https://gitee.com/NBAlexis/CudaLatticeGauge
> 维护建议：修改代码时同步更新本 wiki 相关章节。

---

## 快速开始（5 分钟）

1. **读项目定位**：见下方 [项目简介](#项目简介)。
2. **跑最小示例**：
   - Windows：打开 `Code/CudaLatticeGauge.sln`，启动项目选 `CLGExample`。
   - Linux：`cd Code/CMake && cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86 && make`，然后运行 `../../Bin/Ubuntu/CLGTest`。
3. **看 YAML 配置**：`Bin/Debug/CLGExample.yaml` 是最小可运行的参数文件。
4. **学开发规范**：必读 [kernel-macros.md](kernel-macros.md) 的 CUDA Kernel 强制规则。
5. **查 API / 参数**：用本页下方的文件索引和 [yaml-reference/home.md](yaml-reference/home.md)。

---

## 项目简介

CLGLib（CudaLatticeGauge）是一个基于 CUDA/DCU 的格点规范理论（Lattice QCD / Lattice Gauge Theory）模拟库。支持 SU(2)/SU(3)/U(1)/Z_N/D_N 规范场、Wilson-Dirac 费米子、Staggered (KS) 费米子、HMC/RHMC 更新算法、规范固定 / smearing、以及多种物理量测量。

- **后端**：CUDA（NVIDIA）和 DTK（国产 DCU / HIP）
- **维度**：4D 格点（可配置空间/时间维度）
- **精度**：单精度（默认）或双精度（`-DCLG_DOUBLE=1`）
- **构建**：CMake（Linux / DCU）/ Visual Studio（Windows）
- **配置**：运行时 YAML 参数驱动
- **对象创建**：自定义 RTTI + 工厂（`__CLGDECLARE_CLASS` / `__CLGIMPLEMENT_CLASS`）

---

## 学习路径推荐

| 阶段 | 阅读内容 | 目的 |
|------|---------|------|
| 1 | 本页 [快速参考](#快速参考) + [kernel-macros.md](kernel-macros.md) | 了解写 kernel 的强制规则 |
| 2 | [architecture.md](architecture.md) | 理解模块关系、HMC 流程、初始化流程 |
| 3 | [Tutorial/Tutorial.md](../Tutorial/Tutorial.md) | 跟着例子跑模拟 |
| 4 | 按需深入：[data-field-gauge.md](data-field-gauge.md)、[data-field-fermion.md](data-field-fermion.md)、[data-action.md](data-action.md)、[measurement.md](measurement.md) | 理解具体实现 |
| 5 | [yaml-reference/home.md](yaml-reference/home.md) | 查所有 YAML 参数 |
| 6 | 看 `Code/Applications/CLGExample/` | 了解如何写新应用 |

---

## 模块导航

| 模块 | Wiki | 说明 |
|------|------|------|
| Architecture | [architecture.md](architecture.md) | 整体架构、模块关系、初始化流程、跨模块设计模式 |
| Core | [core.md](core.md) | CUDA 基础设施、内存管理、kernel 宏、全局管理器、RTTI |
| Kernel Macros | [kernel-macros.md](kernel-macros.md) | `_LAUNCH_KERNEL`、线程分解、kernel 入口宏大全 |
| Lattice | [data-lattice.md](data-lattice.md) | 格点索引、边界条件、并行分解 |
| Gauge Field | [data-field-gauge.md](data-field-gauge.md) | 规范场（SU2/SU3/U1/Z2/DN）、设备端结构 |
| Fermion Field | [data-field-fermion.md](data-field-fermion.md) | Wilson-Dirac、Staggered/KS 费米子 |
| Boson Field | [data-field-boson.md](data-field-boson.md) | 标量场、辅助场 |
| Action | [data-action.md](data-action.md) | 规范作用量、费米子作用量 |
| Gauge Fixing | [gauge-fixing.md](gauge-fixing.md) | Landau、Coulomb、MAG、MCG |
| Gauge Smearing | [gauge-smearing.md](gauge-smearing.md) | Stout、APE、HISQ、ASQTAD |
| Measurement | [measurement.md](measurement.md) | 物理量测量 |
| Multi-GPU | [multi-gpu.md](multi-gpu.md) | 多进程分解配置（GpuGrid/HaloWidth/DevicePerNode）、构建、运行 |
| Update | [update.md](update.md) | HMC、积分器、热浴 |
| Sparse Linear Algebra | [sparse-linalg.md](sparse-linalg.md) | BiCGStab、GMRES、GCRODR 等求解器 |
| Math Tools | [tools-math.md](tools-math.md) | SU(N)、FFT、随机数、Gamma 矩阵 |
| Platform Tools | [tools-platform.md](tools-platform.md) | 文件、字符串、YAML、Profiler |
| Test Framework | [test-framework.md](test-framework.md) | CLGTest、参数文件、测试注册 |

### YAML 配置参考

按 `CLGLibManager` 初始化顺序组织的完整参数字典（17 章）：

| 章节 | 内容 | 文件 |
|------|------|------|
| 全局参数 → 测量器 | 所有可配置类的 YAML 键名、类型、默认值、代码位置 | [`yaml-reference/home.md`](yaml-reference/home.md) |

---

## 应用与示例

应用层在 `Code/Applications/`，每个目录一个独立可执行文件。

| 应用 | 目录 | 主要用途 |
|------|------|---------|
| CLGExample | `Code/Applications/CLGExample/` | 最小示例，适合学习 |
| HISQ | `Code/Applications/HISQ/` | HISQ 费米子模拟 |
| StaggeredSpectrum | `Code/Applications/StaggeredSpectrum/` | Staggered 能谱 |
| StaggeredRotation | `Code/Applications/StaggeredRotation/` | 旋转系 Staggered |
| RotationImprovedFermion | `Code/Applications/RotationImprovedFermion/` | 旋转改进费米子 |
| WilsonDiracElectricMagnetic | `Code/Applications/WilsonDiracElectricMagnetic/` | Wilson 费米子 + 电磁场 |
| ElectricChemical | `Code/Applications/ElectricChemical/` | 电化学势 |
| ConstAcc / BetaGradient / MatchingRho 等 | 对应目录 | 专题物理研究 |

启用方式：
- Windows：在 VS 解决方案中勾选对应项目。
- Linux：CMake 加 `-DCLG_<AppName>=1`，例如 `-DCLG_HISQ=1 -DCLG_StaggeredSpectrum=1`。

---

## 关键文件索引

### 必读的入门文件

| 文件 | 路径 | 说明 |
|------|------|------|
| `CLGLib.h` | `Code/CLGLib/CLGLib.h` | 主头文件，按模块顺序包含所有公共头 |
| `CLGLib_Private.h` | `Code/CLGLib/CLGLib_Private.h` | 私有超级头文件：kernel 入口宏、线程分解宏、所有内部头 |
| `CLGDefine.h` | `Code/CLGLib/Core/CLGDefine.h` | 核心宏：namespace、`_T()`、`F(x)`、`TMPARG`、安全删除 |
| `CudaHelper.h` | `Code/CLGLib/Core/CudaHelper.h` | CUDA 内存分配、kernel 启动、设备常量、归约 |
| `CBase.h` | `Code/CLGLib/Core/CBase.h` | RTTI 基类系统：`__CLGDECLARE_CLASS`、`__CLGIMPLEMENT_CLASS` |
| `CLGLibManager.h` | `Code/CLGLib/Core/CLGLibManager.h` | 全局管理器单例，持有所有模块实例 |
| `CLatticeData.h` | `Code/CLGLib/Data/Lattice/CLatticeData.h` | 中央数据容器：场、作用量、更新器、测量器、求解器 |
| `CCommonData.h` | `Code/CLGLib/Data/CCommonData.h` | 枚举、常量、`_HC_*` / `_DC_*` 宏定义 |
| `CField.h` | `Code/CLGLib/Data/Field/CField.h` | 场基类：BLAS 接口、`ApplyOperator` |

### 设备端数学

| 文件 | 路径 | 说明 |
|------|------|------|
| `DeviceInlineTemplate.h` | `Code/CLGLib/Tools/Math/DeviceInlineTemplate.h` | 设备端模板数学函数（`_add`、`_mul`、`_dagger` 等） |
| `SUN.h` | `Code/CLGLib/Tools/Math/SUN.h` | 通用 `deviceSUN<N>` 模板 |
| `SU2.h` / `SU3.h` | `Code/CLGLib/Tools/Math/` | SU2 / SU3 设备端结构 |
| `VectorsN.h` | `Code/CLGLib/Tools/Math/VectorsN.h` | 设备端向量模板 |
| `GammaMatrix.h` | `Code/CLGLib/Tools/Math/GammaMatrix.h` | Gamma 矩阵定义 |
| `CLinearAlgebraHelper.h` | `Code/CLGLib/Tools/Math/CLinearAlgebraHelper.h` | 小矩阵运算、QR、Hessenberg、特征值 |

---

## 快速参考

### CUDA Kernel 规则（强制）

1. **禁止裸 `<<< >>>`**，必须使用 `_LAUNCH_KERNEL` 宏：
   ```cpp
   _LAUNCH_KERNEL(func, block, threads, args);
   // 模板 kernel：
   _LAUNCH_KERNEL(_kernelFunc TMPARG(deviceSU3), block, threads, args);
   ```

2. **模板设备代码用全局函数，不用成员函数**：
   ```cpp
   // 错误
   x.Add(y); x.Dagger();
   // 正确
   _add(x, y); _dagger(x);
   ```

3. **Thread 分解用宏**：
   ```cpp
   preparethread;        // 标准格点分解
   preparethreadHalf;    // 半格点（even/odd）
   preparethreadDir;     // per-link 分解
   preparethreadE(4);    // 每 site 4 个元素
   ```

4. **Kernel 索引用宏**：
   ```cpp
   intokernal;           // uiSiteIndex + 边界检查
   intokernalInt4;       // + sSite4 坐标
   intokernalDirInt4;    // + sSite4 + dir
   intokernalEOHalf;     // even/odd + eta table
   ```

详细规则、宏列表、示例见 [kernel-macros.md](kernel-macros.md)。

### 常用宏

| 宏 | 定义位置 | 说明 |
|----|---------|------|
| `__BEGIN_NAMESPACE` / `__END_NAMESPACE` | `CLGDefine.h` | `namespace CLGLib { }` |
| `__CLGDECLARE_CLASS(name)` | `CBase.h` | 类声明宏（RTTI 支持） |
| `__CLGIMPLEMENT_CLASS(name)` | `CBase.h` | 类实现宏（放在 .cpp） |
| `__CLG_REGISTER_HELPER_HEADER(name)` | `CBase.h` | 强制 Linux 静态链接保留构造函数（放在 .h） |
| `_LAUNCH_KERNEL(func, block, threads, args)` | `CudaHelper.h` | 跨平台 kernel 启动 |
| `TMPARG(...)` | `CLGDefine.h` | 模板参数宏：`#define TMPARG(...) <__VA_ARGS__>` |
| `preparethread` | `CLGLib_Private.h` | `block = _HC_DecompBlock, threads = _HC_DecompThread` |
| `intokernal` | `CLGLib_Private.h` | 计算 `uiSiteIndex` 并做边界检查 |
| `_HC_Volume` | `CCommonData.h` | host 端格点体积 |
| `_DC_Volume` | `CCommonData.h` | device 端格点体积（常量内存） |
| `F(x)` | `CLGDefine.h` | 字面量宏，单精度时 `x##f`，双精度时 `x` |
| `_T("...")` | 平台相关 | 跨平台字符串字面量 |

### 全局单例访问

| 函数 | 返回类型 | 说明 |
|------|---------|------|
| `appGetLattice()` | `CLatticeData*` | 格点数据（场、作用量、测量器） |
| `appGetCudaHelper()` | `CCudaHelper*` | CUDA 辅助（内存、stream、常量） |
| `appGetRandom()` | `CRandom*` | 随机数生成器 |
| `appGetFermionSolver(BYTE id)` | `CSLASolver*` | 按 field id 取单偏移求解器 |
| `appGetMultiShiftSolver(BYTE id)` | `CMultiShiftSolver*` | 按 field id 取多偏移求解器 |

### 场类型枚举（EFieldType）

| 枚举值 | 说明 |
|--------|------|
| `EFT_GaugeSU3` | SU(3) 规范场 |
| `EFT_GaugeSU3_12` | SU(3) 12 参数表示 |
| `EFT_GaugeU1` | U(1) 规范场 |
| `EFT_GaugeReal` | 实规范场 |
| `EFT_GaugeZ2` ~ `EFT_GaugeZ6` | Z_N 离散规范场 |
| `EFT_GaugeD3` / `EFT_GaugeD4` / `EFT_GaugeD8` | D_N 二元多面体群 |
| `EFT_GaugeSU2` / `EFT_GaugeSU4` ~ `EFT_GaugeSU8` / `EFT_GaugeSUN` | SU(N) 规范场 |
| `EFT_FermionWilsonSquareSU3` | Wilson-Dirac 费米子 |
| `EFT_FermionStaggeredSU3` / `EFT_FermionStaggeredU1` / `EFT_FermionStaggeredSU2` / ... | Staggered 费米子 |
| `EFT_BosonComplex` / `EFT_BosonComplexVectorN` / `EFT_BosonReal` | 玻色子场 |

完整定义见 `Code/CLGLib/Data/CCommonData.h`。

---

## 构建速查

### Windows

```bash
# 直接用 Visual Studio 打开
Code/CudaLatticeGauge.sln
# 输出目录：Bin/Debug/ 或 Bin/Release/
```

### Linux / CUDA

```bash
cd Code/CMake
cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86 \
    -DCLG_HISQ=1 -DCLG_StaggeredSpectrum=1
make
# 输出目录：Bin/Ubuntu/
# 运行测试：../../Bin/Ubuntu/CLGTest
```

### Linux / DCU（DTK）

```bash
cd Code/CMake
export CXX=/path/to/dtk/llvm/bin/clang++
cmake CMakeLists.txt -DCLG_BACKEND=DTK -DCLG_GPU_ARCH=gfx906 \
    -DCMAKE_CUDA_COMPILER_FORCED=TRUE
make
```

常用选项：
- `-DCLG_DOUBLE=1`：双精度
- `-DCLG_DEBUG=1`：输出到 `Bin/UbuntuDebug/`
- `-DCLG_PROFILER=1`：启用 `_RECORD` 性能分析（默认关闭；开启后每次记录都会同步设备以测量 kernel 耗时，开销大）
- `-DCLG_CHECKSYNCHRONIZE=1`：每次 kernel launch 后同步设备并尽早报错（默认关闭；调试用，会使所有 kernel 串行化）
- `-DCLG_<AppName>=1`：启用对应应用

详细说明见 [README.md](../README.md) 和 [Tutorial/07-cmake.md](../Tutorial/07-cmake.md)。

---

## 如何添加新功能

所有 YAML 可实例化的对象都通过 RTTI 工厂自动注册，**不需要手动修改 `CLGLibManager` 的 switch-case**。注册分两步：

1. 在 `.h` 中：
   ```cpp
   __CLG_REGISTER_HELPER_HEADER(CMyClass)
   class CLGAPI CMyClass : public CAction   // 或 CMeasure / CField / ...
   {
       __CLGDECLARE_CLASS(CMyClass)
       // ...
   };
   ```

2. 在 `.cpp` 中：
   ```cpp
   __CLGIMPLEMENT_CLASS(CMyClass)
   ```

3. 在 YAML 中通过类名字符串使用：
   ```yaml
   Action1:
       ActionName: CMyClass
       Beta: 5.5
   ```

### 添加一个新的 Action

1. 继承 `CAction`（`Data/Action/CAction.h`）。
2. 实现：
   - `EnergySingleField(...)` — 能量计算
   - `CalculateForceOnGaugeSingleField(...)` — 力计算
   - `PrepareForHMCSingleField(...)` — HMC 前准备
3. 头文件加 `__CLG_REGISTER_HELPER_HEADER(CMyClass)`，cpp 加 `__CLGIMPLEMENT_CLASS(CMyClass)`。
4. 设备端 kernel 放在 `.cu` 文件，使用 `_LAUNCH_KERNEL` 启动；模板化以支持多规范群。
5. 参考 [data-action.md](data-action.md)。

### 添加一个新的 Measurement

1. 继承 `CMeasure`（或 `CMeasureStochastic` 用于随机估计）。
2. 实现回调：
   - `OnConfigurationAcceptedSingleField(...)` — 规范/玻色子测量
   - `OnConfigurationAcceptedZ4SingleField(...)` — 费米子随机估计
   - `SourceSanningSingleField(...)` — 源扫描模式
3. 实现 `Report()` 和 `Reset()`。
4. 使用 `__CLG_REGISTER_HELPER_HEADER` / `__CLGIMPLEMENT_CLASS` 注册。
5. 参考 [measurement.md](measurement.md)。

### 添加一个新的 Gauge Fixing 算法

1. 继承 `CGaugeFixing`（`GaugeFixing/CGaugeFixing.h`）。
2. 实现 `GaugeFixing(pResGauge)` 和 `CheckRes(pGauge)`。
3. kernel 使用 `_LAUNCH_KERNEL`，设备代码用 `DeviceInlineTemplate.h` 中的全局函数。
4. 使用 `__CLG_REGISTER_HELPER_HEADER` / `__CLGIMPLEMENT_CLASS` 注册。
5. 参考 [gauge-fixing.md](gauge-fixing.md)。

### 添加一个新的 CUDA Kernel

1. 声明：`__global__ void _CLG_LAUNCH_BOUND _kernelFunc(...)`。
2. 入口用 `intokernal*` 系列宏。
3. 启动用 `_LAUNCH_KERNEL(func, block, threads, args)`。
4. 设备端数学用 `DeviceInlineTemplate.h` 中的全局模板函数。
5. 参考 [core.md](core.md) 和 [kernel-macros.md](kernel-macros.md)。

---

## 常见任务索引

| 任务 | 关键文件 / 文档 |
|------|----------------|
| 修改或新增 YAML 参数 | `Code/CLGLib/Tools/CYAMLParser.h`，[yaml-reference/home.md](yaml-reference/home.md) |
| 调试 kernel | `Code/CLGLib/Core/CudaHelper.h` 中的 `checkCudaErrors`，`appGeneral` 日志 |
| 添加新的物理量输出 | `Measurement/CMeasure*.h` 模式，[measurement.md](measurement.md) |
| 修改 HMC 积分器 | `Code/CLGLib/Update/Continous/CIntegrator*.h`，[update.md](update.md) |
| 修改或新增求解器 | `Code/CLGLib/SparseLinearAlgebra/`，[sparse-linalg.md](sparse-linalg.md) |
| 支持新规范群 | `Code/CLGLib/Tools/Math/SUN.h`、`Data/Field/Gauge/CFieldGaugeLink.h` |
| 在 VS 和 CMake 间同步文件列表 | `Code/Tools/CLGMakeWriter/` |
| 跑回归测试 | `Code/CLGTest/`，[test-framework.md](test-framework.md) |

---

## 贡献与维护

- 修改了 `AGENTS.md` 中提及的文件、风格、结构、工作流时，同步更新对应文档。
- 新增类或 YAML 参数时，同步更新：
  - 对应模块 wiki（如 `data-action.md`、`measurement.md`）
  - [yaml-reference/](yaml-reference/) 相关章节
  - 本页的关键文件索引和场类型枚举
