> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 01. CLGManager 全局参数

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::InitialWithParameter()` (lines 997-1260)

这些参数在初始化流程的最开始读取，控制 CUDA 设备、核心对象分配和整体框架行为。

---

## YAML 格式注意

本框架的 YAML 解析器（`CYAMLParser`）只实现了最基本的键值对语法，**不支持完整 YAML 1.x 规范**。配置时请遵循以下规则：

- **数组必须写成单行方括号形式**：`Key: [a, b, c]`。
- **数组不能跨行**：`[a, b,\n c]` 会因右括号在下一行而解析失败。
- **不支持 `-` 开头的 YAML 列表**：`- a\n- b\n- c` 会被当成普通键解析失败。
- **嵌套块通过缩进表示**：子块比父块多缩进若干空格即可，具体缩进量不影响，只要同一层级保持一致。

---

## 设备与内存

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `DeviceIndex` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1002 | CUDA 设备编号，传给 `cudaSetDevice()`。多 GPU 构建下被 `DevicePerNode` 覆盖 |
| `DevicePerNode` | INT | 机器 GPU 数 | 否 | `CLGLibManager.cpp` ~L1137-1150 | （仅 `_CLG_MULTI_GPU` 构建）每节点 GPU 数；实际设备 = `Rank % DevicePerNode`。单卡测试机设 `1` 让所有 rank 共享设备 0（oversubscribe）。MG 测试 yaml 不固定此键，由 `runtest` 第二参数经 `--device-per-node` 注入。见 [multi-gpu.md](../multi-gpu.md) |
| `AllocateBuffer` | Real | `0.0` | 否 | `CLGLibManager.cpp` ~L1018 | 预分配 GPU 缓冲区大小（GB）。若 `0.1 < value < 32.0` 则分配 |

## 列表长度声明

以下参数必须在 YAML 顶层显式声明，框架根据它们创建对应数量的对象：

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ActionListLength` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1030 | 作用量个数。框架将依次读取 `Action1`, `Action2`, ... |
| `GaugeFieldCount` | INT | `1` | 否 | `CLGLibManager.cpp` ~L1026 | 规范场个数。第一个为 `Gauge`，后续为 `Gauge2`, `Gauge3`, ... |
| `FermionFieldCount` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1027 | 费米子场个数。框架将依次读取 `FermionField1`, `FermionField2`, ... |
| `BosonFieldCount` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1028 | 玻色子场个数。框架将依次读取 `BosonField1`, `BosonField2`, ... |
| `MeasureListLength` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1032 | 测量器个数。框架将依次读取 `Measure1`, `Measure2`, ... |
| `OtherGaugeFieldCount` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1029 | （已废弃）保留兼容性 |

## 调试与输出

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `VerboseLevel` | string | `PARANOIAC` | 否 | `CLGLibManager.cpp` ~L1040 | 日志详细程度。见下表 |
| `VerboseOutput` | string | （无） | 否 | `CLGLibManager.cpp` ~L1042 | 日志输出格式。`datetime` 表示带时间戳 |
| `Profiler` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1056 | 是否启用性能分析器。`1` 启用，`0` 禁用 |
| `ShowParameterContent` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1044 | 是否在初始化时打印所有参数内容。`1` 打印 |

**VerboseLevel 枚举值**（来自 `Tools/Tracer.h`）：

| 枚举值 | 说明 |
|--------|------|
| `CRUCIAL` | 仅关键错误/警告 |
| `WARNING` | 包含警告 |
| `GENERAL` | 一般信息（默认级别） |
| `DETAILED` | 详细信息 |
| `PARANOIAC` | 最详细（调试级） |

## 其他全局开关

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `UseLogADefinition` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1054 | 是否使用 `log(A)` 定义。代码中用于 `CCommonData::m_bUseLogADefinition` |
| `MILC_StaggeredPhase` | INT | `0` | 否 | `CLGLibManager.cpp` ~L1058 | 是否使用 MILC staggered phase 约定。代码中传给 `CIndexData::SetMILCPhase()` |
| `SummationDecompose` | INT | 设备最大线程数 | 否 | `CLGLibManager.cpp` ~L1052 | 求和分解的线程数。默认使用设备属性中的最大线程数 |
| `CacheSolution` | INT | `1` | 否 | `CLGLibManager.cpp` ~L1050 | 是否缓存上一次求解器的解。`1` 缓存，`0` 不缓存 |

---

## 初始化顺序

这些全局参数之后，CLGManager 按以下顺序调用各创建函数：

1. `InitialLatticeAndConstant(params)` — 格点常数 → 见 [02-lattice.md](02-lattice.md)
2. `InitialRandom(params)` — 随机数 → 见 [03-random.md](03-random.md)
3. `CreateIndexAndBoundary(params)` — 索引/边界 → 见 [04-index-boundary.md](04-index-boundary.md)
4. `CreateGaugeFields(params)` — 规范场 → 见 [05-gauge-fields.md](05-gauge-fields.md)
5. `CreateBoundaryFields(params)` — 边界场 → 见 [06-boundary-fields.md](06-boundary-fields.md)
6. `CreateBosonFields(params)` — 玻色子场 → 见 [07-boson-fields.md](07-boson-fields.md)
7. `CreateFermionFields(params)` — 费米子场 → 见 [08-fermion-fields.md](08-fermion-fields.md)
8. `InitialIndexBuffer()` — 索引缓冲烘焙
9. `InitialFieldBuffer()` — 场缓冲初始化
10. `CreateActionList(params)` — 作用量 → 见 [09-actions.md](09-actions.md)
11. `CreateSolver(params)` — 求解器 → 见 [10-solvers.md](10-solvers.md)
12. `CreateMultiShiftSolver(params)` — 多移求解器 → 见 [11-multi-shift-solvers.md](11-multi-shift-solvers.md)
13. `CreateGaugeSmearing(params)` — 规范平滑 → 见 [12-gauge-smearing.md](12-gauge-smearing.md)
14. `CreateGaugeStapleCache(params)` — Staple 缓存 → 见 [13-staple-cache.md](13-staple-cache.md)
15. `CreateGaugeFixing(params)` — 规范固定 → 见 [14-gauge-fixing.md](14-gauge-fixing.md)
16. `CreateUpdator(params)` — 更新器 → 见 [15-updators.md](15-updators.md)
17. `CreateMeasurement(params)` — 测量器 → 见 [17-measurements.md](17-measurements.md)

---

[< 返回目录](home.md) | [下一章：格点常数 >](02-lattice.md)
