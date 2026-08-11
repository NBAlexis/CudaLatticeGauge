> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 02. 格点常数

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::InitialLatticeAndConstant()` (lines 58-433)

这些参数定义格点的几何结构、CUDA 线程分解、随机数种子等全局常数，存储在 `CCommonData` 中供所有 kernel 访问。

---

## 格点几何

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Dim` | INT | `4` | 否 | `CLGLibManager.cpp` ~L100 | 时空维度。代码中写入 `_HC_Dim` |
| `Dir` | INT | `4` | 否 | `CLGLibManager.cpp` ~L102 | 方向数（通常 = Dim）。代码中写入 `_HC_Dir` |
| `LatticeLength` | INT[] | `[8,8,8,8]` | 否 | `CLGLibManager.cpp` ~L108 | 各方向格点数 `[Lx, Ly, Lz, Lt]`。代码中写入 `_HC_Lx`, `_HC_Ly`, `_HC_Lz`, `_HC_Lt` |
| `LatticeLengthDebug` | INT[] | （无） | 否 | `CLGLibManager.cpp` ~L110 | Debug 模式下覆盖 `LatticeLength`。仅当 `_CLG_DEBUG` 宏定义时生效 |
| `Center` | INT[] | `[Lx/2, Ly/2, Lz/2, Lt/2]` | 否 | `CLGLibManager.cpp` ~L128 | 中心点坐标（整数除法）。代码中写入 `_HC_Centerx`, `_HC_Centery`, `_HC_Centerz`, `_HC_Centert` |

## 多 GPU 进程网格（Multi-GPU）

> 以下参数仅在 `-DCLG_MULTI_GPU=1` 构建下生效（此时框架要求 MPI 可用）。
> 单 GPU 构建（`_CLG_MULTI_GPU=0`）下 `GpuGrid` 恒为 `[1,1,1,1]`、`HaloWidth` 读取但不影响分解，
> 这些配置对数值结果无影响。详细指南见 [multi-gpu.md](../multi-gpu.md)。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `GpuGrid` | INT[] | `[1,1,1,1]` | 否 | `CLGLibManager.cpp` ~L153-176 | 进程网格 `[Gx, Gy, Gz, Gt]`。**各方向乘积必须等于 MPI 进程数**。格点被均分为每 rank 一个子格点。代码中写入 `_HC_GpuGridX/Y/Z/T` |
| `HaloWidth` | INT | `2` | 否 | `CLGLibManager.cpp` ~L163-169 | 单方向最远 stencil 半径（HISQ Naik 需要 2）。被切分方向的本地长度必须 ≥ 该值。代码中写入 `_HC_HaloWidth` |

**校验规则**（`CLGComm::SetGpuGrid`，`Code/CLGLib/Core/Distributed/CLGComm.cpp` ~L101-167）：

1. 每个分量 ≥ 1；
2. `GpuGrid` 乘积 == MPI world size（`mpiexec -n N` 的 N）；
3. 每个方向的全局格点长度必须能被对应 grid 分量整除（目前只支持均匀切分）；
4. 被切分方向（grid 分量 > 1）的本地长度 ≥ `HaloWidth`。

任一规则不满足都会以 `appCrucial` 报错并终止，不会静默算错。

**MG 测试驱动元数据**（仅 `CLGTest` 的 `--mg-config` / `runtest` 读取，库不消费；
`TestSuit_MG.yaml` 中每个 `_TEST_MULTIGPU` 开放测试段必须补齐）：

| 键名 | 类型 | 说明 |
|------|------|------|
| `MultiGPUTestRankCount` | INT | 该测试要求的 MPI rank 数（≥ 1） |
| `MultiGPUTestGrid` | INT[]（4 个） | 该测试要求的进程网格；四因子均 ≥ 1 且乘积 == `MultiGPUTestRankCount` |

测试 yaml 中**不**固定 `GpuGrid`（worker 经 `--gpu-grid` 注入元数据网格）与
`DevicePerNode`（`runtest` 第二参数经 `--device-per-node` 注入）。

**使用示例**（8⁴ 格点拆 2 个进程，切 t 方向）：

```yaml
Dim : 4
Dir : 4
LatticeLength : [8, 8, 8, 8]
GpuGrid : [1, 1, 1, 2]   # 乘积 2 == mpiexec -n 2；本地子格点 8×8×8×4
HaloWidth : 2            # 默认即可（HISQ Naik）
```

**已写入 CCommonData 的常量**（单 GPU 构建下退化为全局值 / `[1,1,1,1]` / 0，已有代码读取安全）：

| 常量名 | 说明 |
|--------|------|
| `_HC_GpuGridX/Y/Z/T` | 进程网格各方向分量（`ECI_GpuGridX+i`） |
| `_HC_GlobalLx/Ly/Lz/Lt` | **全局**格点长度（`ECI_GlobalLx+i`） |
| `_HC_Lx/Ly/Lz/Lt`、`_HC_Volume*` | **本地**（本 rank）子格点尺寸 |
| `_HC_OffsetX/Y/Z/T` | 本 rank 子格点在全局坐标中的起点（`ECI_GlobalOffsetX+i`） |
| `_HC_HaloWidth` | halo 宽度 |

## CUDA 线程分解

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MaxThreadPerBlock` | INT | `_CLG_LAUNCH_MAX_THREAD` | 否 | `CLGLibManager.cpp` ~L104 | 每个 CUDA block 的最大线程数。代码中写入 `_HC_Thread` |
| `ThreadAutoDecompose` | INT | `1` | 否 | `CLGLibManager.cpp` ~L134 | 是否自动分解 block/thread。`1` 自动，`0` 使用 `ThreadDecompose` |
| `ThreadDecompose` | INT[] | （无） | 否 | `CLGLibManager.cpp` ~L138 | 手动线程分解 `[blocks_x, blocks_y, blocks_z]`。仅当 `ThreadAutoDecompose=0` 时使用 |

## 随机数

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `RandomSeed` | INT | `81192` | 否 | `CLGLibManager.cpp` ~L144 | 随机数种子 |
| `RandomSeedType` | string | （无） | 否 | `CLGLibManager.cpp` ~L148 | 若设为 `ERST_Timestamp`，则用当前时间作为种子，覆盖 `RandomSeed` |
| `RandomType` | string | `ER_Schrage` | 否 | `CLGLibManager.cpp` ~L146 | 随机数生成器类型。见 [03-random.md](03-random.md) |

## 物理常数与精度

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `GaugeMomentumFactor` | Real | `1.0` | 否 | `CLGLibManager.cpp` ~L132 | 规范动量因子 `m_fGaugeMomentumFactor`。代码中写入 `_HC_GaugeMomentum` |
| `ExponentialPrecision` | INT | `8` | 否 | `CLGLibManager.cpp` ~L150 | 指数函数精度。代码中写入 `_HC_ExpPrecision` |
| `StochasticGaussian` | INT | `0` | 否 | `CLGLibManager.cpp` ~L152 | 是否使用高斯随机源。代码中写入 `CCommonData::m_bStochasticGaussian` |

## 其他开关

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ShowParameterContent` | INT | `0` | 否 | `CLGLibManager.cpp` ~L154 | 是否在初始化时打印参数内容。`1` 打印 |

---

## CCommonData 中的常量

上述参数读取后，全部被写入 `CCommonData`（设备端常量内存），供所有 CUDA kernel 访问：

| 常量名 | 来源键 | 类型 |
|--------|--------|------|
| `_HC_Dim` | `Dim` | UINT |
| `_HC_Dir` | `Dir` | UINT |
| `_HC_Lx` / `Ly` / `Lz` / `Lt` | `LatticeLength` | UINT |
| `_HC_Volume` | `Lx * Ly * Lz * Lt` | UINT |
| `_HC_Volume_xyz` | `Lx * Ly * Lz` | UINT |
| `_HC_Centerx` / `y` / `z` / `t` | `Center` | SINT |
| `_HC_Thread` | `MaxThreadPerBlock` | UINT |
| `_HC_GaugeMomentum` | `GaugeMomentumFactor` | Real |
| `_HC_ExpPrecision` | `ExponentialPrecision` | UINT |

---

[< 返回目录](home.md) | [< 上一章：全局参数](01-global.md) | [下一章：随机数 >](03-random.md)
