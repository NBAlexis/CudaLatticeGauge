> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 11. 多移求解器 (Multi-Shift Solvers)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateMultiShiftSolver()` (lines 854-866)

---

## 配置参数

YAML 键名为 `MSSolver`（第一个）或 `MSSolver2`, `MSSolver3`, ...。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `SolverName` | string | `CMultiShiftGMRES` | 否 | `CLGLibManager.cpp` ~L858 | 多移求解器类名。通过 `appCreate()` 实例化 |
| `SolverForFieldId` | INT | `2` | 否 | `CLGLibManager.cpp` ~L860 | 此求解器服务的费米子场 ID |

## 支持的多移求解器类

| 类名 | 说明 |
|------|------|
| `CMultiShiftGMRES` | 多移 GMRES |
| `CMultiShiftFOM` | 多移 FOM |
| `CMultiShiftBiCGStab` | 多移 BiCGStab |
| `CMultiShiftNested` | 嵌套多移求解器（内层使用普通求解器） |

## 多移 GMRES 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CMultiShiftGMRES.cpp`

**读取函数**：`CMultiShiftGMRES::Configurate()` (lines 42-94)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MaxDim` | INT | `20` | 否 | `CMultiShiftGMRES.cpp` ~L47 | Krylov 子空间维度（限制 5-30） |
| `Restart` | INT | `3` | 否 | `CMultiShiftGMRES.cpp` ~L67 | 重启次数 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CMultiShiftGMRES.cpp` ~L71 | 是否使用绝对精度 |
| `Accuracy` | Real | `0.000001` | 否 | `CMultiShiftGMRES.cpp` ~L75 | 收敛容差 |
| `CheckAddSystem` | INT | `0` | 否 | `CMultiShiftGMRES.cpp` ~L85 | 是否检查附加系统。`1`=启用 |

## 多移 FOM 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CMultiShiftFOM.cpp`

**读取函数**：`CMultiShiftFOM::Configurate()` (lines 37-82)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MaxDim` | INT | `20` | 否 | `CMultiShiftFOM.cpp` ~L42 | Krylov 子空间维度（限制 5-100） |
| `UseCudaForSmallMatrix` | INT | `0` | 否 | `CMultiShiftFOM.cpp` ~L54 | 小矩阵是否用 CUDA 计算。`1`=启用 |
| `Restart` | INT | `3` | 否 | `CMultiShiftFOM.cpp` ~L65 | 重启次数 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CMultiShiftFOM.cpp` ~L69 | 是否使用绝对精度 |
| `Accuracy` | Real | `0.000001` | 否 | `CMultiShiftFOM.cpp` ~L73 | 收敛容差 |

## 多移 BiCGStab 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CMultiShiftBiCGStab.cpp`

**读取函数**：`CMultiShiftBiCGStab::Configurate()` (lines 32-61)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `DiviationStep` | INT | `10` | 否 | `CMultiShiftBiCGStab.cpp` ~L36 | 偏差检查步长间隔 |
| `MaxStep` | INT | `20` | 否 | `CMultiShiftBiCGStab.cpp` ~L40 | 最大迭代步数 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CMultiShiftBiCGStab.cpp` ~L44 | 是否使用绝对精度 |
| `Accuracy` | DOUBLE/Real | `0.000001` | 否 | `CMultiShiftBiCGStab.cpp` ~L50 | 收敛容差 |

## 嵌套多移求解器参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CMultiShiftNested.cpp`

**读取函数**：`CMultiShiftNested::Configurate()` (lines 30-63)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `NestedSolver` | 子参数块 | （无） | **是** | `CMultiShiftNested.cpp` ~L35 | 内层求解器配置子块 |
| `NestedSolver.SolverName` | string | `CSLASolverBiCGStab` | 否 | `CMultiShiftNested.cpp` ~L42 | 内层求解器类名 |
| `NestedSolver.SolverForFieldId` | INT | `2` | 否 | `CMultiShiftNested.cpp` ~L44 | 内层求解器服务的费米子场 ID |

**注意**：`NestedSolver` 是一个子参数块。除 `SolverName` 和 `SolverForFieldId` 外，内层求解器自己的所有参数（见 [10-solvers.md](10-solvers.md) 中对应求解器的参数表）都可以放在 `NestedSolver:` 块内，例如 `NestedSolver.MaxStep`、`NestedSolver.Accuracy` 等。

---

## 与场 / 作用量的依赖关系

| 使用场景 A | 必须配置 B | 说明 |
|---|---|---|
| HMC 中 HISQ 费米子作用量（`CActionFermionHISQCombined`） | `MSSolver`（每个 `FieldIds` 对应场一个） | `CalculateForce` 中调用多移求解器做有理逼近。见 [08b-ks-hisq-fermion.md](08b-ks-hisq-fermion.md)、[09b-fermion-actions.md](09b-fermion-actions.md) |
| HMC 中 improved KS 费米子作用量（`CActionFermionKSImprove` / `KSImproveCombined`） | `MSSolver`（每个相关场一个） | 同上。见 [09b-fermion-actions.md](09b-fermion-actions.md) |
| HMC 中 rational 近似的 Wilson / KS 费米子 | `MSSolver` | `CFieldFermion::CalculateForce` 中 `appGetMultiShiftSolver` 不能为空。见 [08a-wilson-fermion.md](08a-wilson-fermion.md) |
| `CMultiShiftNested` | `NestedSolver` 子块 | 内层求解器配置。参数见 [10-solvers.md](10-solvers.md) |

---

[< 返回目录](home.md) | [< 上一章：求解器](10-solvers.md) | [下一章：规范平滑 >](12-gauge-smearing.md)
