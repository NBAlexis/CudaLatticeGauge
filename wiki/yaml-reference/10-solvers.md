> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 10. 求解器 (Solvers)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateSolver()` (lines 840-852)

---

## 配置参数

YAML 键名为 `Solver`（第一个）或 `Solver2`, `Solver3`, ...。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `SolverName` | string | `CSLASolverBiCGStab` | 否 | `CLGLibManager.cpp` ~L844 | 求解器类名。通过 `appCreate()` 实例化 |
| `SolverForFieldId` | INT | `2` | 否 | `CLGLibManager.cpp` ~L846 | 此求解器服务的费米子场 ID |

## 支持的求解器类

| 类名 | 说明 |
|------|------|
| `CSLASolverBiCGStab` | BiCGStab 稳定双共轭梯度 |
| `CSLASolverGMRES` | GMRES 广义最小残差 |
| `CSLASolverGCR` | GCR 广义共轭残差 |
| `CSLASolverGCRODR` | GCRODR 带收缩 GMRES |
| `CSLASolverGMRESMDR` | GMRES-MDR 多收缩 GMRES |
| `CSolverTFQMR` | TFQMR 转置自由准最小残差 |
| `CSLASolverCG` | CG 共轭梯度（仅 Hermitian 正定 `EFO_F_DDdagger`；`D`/`Ddagger` 通过 `DDdagger` 组合求解） |
| `CSLASolverDeflatedCG` | Deflated CG（Arnoldi 收缩子空间 + 粗网格校正，继承 CSLASolverCG） |

## BiCGStab 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CSolverBiCGstab.cpp`

**读取函数**：`CSolverBiCGstab::Configurate()` (lines 34-94)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `DiviationStep` | INT | `10` | 否 | `CSolverBiCGstab.cpp` ~L38 | 偏差检查步长间隔 |
| `MaxStep` | INT | `20` | 否 | `CSolverBiCGstab.cpp` ~L42 | 最大迭代步数 |
| `Restart` | INT | `3` | 否 | `CSolverBiCGstab.cpp` ~L46 | 重启次数 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CSolverBiCGstab.cpp` ~L50 | 是否使用绝对精度。`1`=绝对精度，`0`=相对精度 |
| `Accuracy` | Real/DOUBLE | `0.000001` | 否 | `CSolverBiCGstab.cpp` ~L56 | 收敛容差 |
| `SmallRho` | Real/DOUBLE | `0.00000001` | 否 | `CSolverBiCGstab.cpp` ~L74 | 小 rho 阈值 |

## GMRES 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CSolverGMRES.cpp`

**读取函数**：`CSolverGMRES::Configurate()` (lines 36-100)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MaxDim` | INT | `20` | 否 | `CSolverGMRES.cpp` ~L41 | Krylov 子空间维度（限制 5-100） |
| `UseCudaForSmallMatrix` | INT | `0` | 否 | `CSolverGMRES.cpp` ~L53 | 小矩阵是否用 CUDA 计算。`1`=启用 |
| `Restart` | INT | `3` | 否 | `CSolverGMRES.cpp` ~L66 | 重启次数 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CSolverGMRES.cpp` ~L71 | 是否使用绝对精度 |
| `Accuracy` | Real/DOUBLE | `0.000001` | 否 | `CSolverGMRES.cpp` ~L78 | 收敛容差 |

## GCR 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CSolverGCR.cpp`

**读取函数**：`CSolverGCR::Configurate()` (lines 34-75)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MaxDim` | INT | `20` | 否 | `CSolverGCR.cpp` ~L39 | Krylov 子空间维度（限制 5-100） |
| `Restart` | INT | `3` | 否 | `CSolverGCR.cpp` ~L50 | 重启次数 |
| `Iterate` | INT | `50` | 否 | `CSolverGCR.cpp` ~L54 | 迭代次数（内层迭代控制） |
| `DiviationStep` | INT | `5` | 否 | `CSolverGCR.cpp` ~L58 | 偏差检查步长间隔 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CSolverGCR.cpp` ~L62 | 是否使用绝对精度 |
| `Accuracy` | Real | `0.000001` | 否 | `CSolverGCR.cpp` ~L66 | 收敛容差 |

## GCRODR 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CSolverGCRODR.cpp`

**读取函数**：`CSolverGCRODR::Configurate()` (lines 69-123)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MDim` | INT | `10` | 否 | `CSolverGCRODR.cpp` ~L74 | Krylov 子空间维度（限制 5-31） |
| `KDim` | INT | `5` | 否 | `CSolverGCRODR.cpp` ~L85 | 收缩子空间维度（限制 2 到 MDim-2） |
| `Restart` | INT | `3` | 否 | `CSolverGCRODR.cpp` ~L98 | 重启次数 |
| `RecalculateR` | INT | `5` | 否 | `CSolverGCRODR.cpp` ~L102 | 重新计算 R 矩阵的频率 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CSolverGCRODR.cpp` ~L106 | 是否使用绝对精度 |
| `Accuracy` | Real | `0.000001` | 否 | `CSolverGCRODR.cpp` ~L110 | 收敛容差 |
| `DeflateSpace` | string | `EEDT_SVD` | 否 | `CSolverGCRODR.cpp` ~L121 | 收缩空间类型。见下表 |

**DeflateSpace 枚举值**：

| 枚举值 | 说明 |
|--------|------|
| `EEDT_SVD` | 使用 SVD 构造收缩空间（默认） |
| `EEDT_REV` | 使用 Ritz 向量（实部） |
| `EEDT_HEV` | 使用 Harmonic Ritz 向量 |

**注意**：`CSLASolverGMRESMDR` 继承自 GCRODR，无额外参数。

## CG 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CSolverCG.cpp`

**读取函数**：`CSLASolverCG::Configurate()` (lines 25-47)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MaxStep` | INT | `5000` | 否 | `CSolverCG.cpp` ~L30 | 最大迭代步数 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CSolverCG.cpp` ~L34 | 是否使用绝对精度。`1`=绝对精度，`0`=相对精度 |
| `Accuracy` | Real/DOUBLE | `0.000001` | 否 | `CSolverCG.cpp` ~L38 | 收敛容差 |

**注意**：CG 只接受 Hermitian 正定算子 `EFO_F_DDdagger`。求解 `D` 时内部组合为 $D^\dagger(DD^\dagger)^{-1}$，求解 `Ddagger` 时组合为 $(DD^\dagger)^{-1}D$，其它算子直接报错。

## DeflatedCG 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CSolverDeflatedCG.cpp`

**读取函数**：`CSLASolverDeflatedCG::Configurate()` (lines 114-141)

继承 CG 全部参数，另有：

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `DeflateDim` | INT | `16` | 否 | `CSolverDeflatedCG.cpp` ~L119 | Arnoldi 步数（Hessenberg 维度，不小于 `Deflate`） |
| `Deflate` | INT | `8` | 否 | `CSolverDeflatedCG.cpp` ~L123 | 收缩向量个数（最小 Ritz 向量数，至少为 1） |
| `ReDeflateInterval` | INT | `0` | 否 | `CSolverDeflatedCG.cpp` ~L135 | 每 N 次求解重建收缩子空间；`0`=只建一次，后续复用 |


## TFQMR 参数

**来源文件**：`Code/CLGLib/SparseLinearAlgebra/CSolverTFQMR.cpp`

**读取函数**：`CSolverTFQMR::Configurate()` (lines 33-62)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `DiviationStep` | INT | `10` | 否 | `CSolverTFQMR.cpp` ~L37 | 偏差检查步长间隔 |
| `MaxStep` | INT | `20` | 否 | `CSolverTFQMR.cpp` ~L41 | 最大迭代步数 |
| `Restart` | INT | `1` | 否 | `CSolverTFQMR.cpp` ~L45 | 重启次数 |
| `AbsoluteAccuracy` | INT | `0` | 否 | `CSolverTFQMR.cpp` ~L49 | 是否使用绝对精度 |
| `Accuracy` | Real | `0.000001` | 否 | `CSolverTFQMR.cpp` ~L53 | 收敛容差 |

---

## 与场 / 作用量的依赖关系

| 使用场景 A | 必须配置 B | 说明 |
|---|---|---|
| 任意费米子测量器（Meson / Chiral / AMomentum / BerryPhase 等） | `Solver`（`SolverForFieldId` 匹配被测场） | 测量中通过 `InverseD` / `InverseDDdagger` 反演。见 [17b-fermion-measurements.md](17b-fermion-measurements.md) |
| HMC 中 `ER_NoRational` 的 Wilson / KS 费米子 | `Solver`（`SolverForFieldId` 匹配） | 单移求解器即可。见 [08a-wilson-fermion.md](08a-wilson-fermion.md)、[08b-ks-hisq-fermion.md](08b-ks-hisq-fermion.md) |
| HMC 中 HISQ / improved KS / rational 费米子 | `MSSolver` | 力计算需要多移求解器。见 [11-multi-shift-solvers.md](11-multi-shift-solvers.md)、[09b-fermion-actions.md](09b-fermion-actions.md) |

---

[< 返回目录](home.md) | [< 上一章：作用量](09-actions.md) | [下一章：多移求解器 >](11-multi-shift-solvers.md)
