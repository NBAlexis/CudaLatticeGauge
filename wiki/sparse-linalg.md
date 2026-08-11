# Sparse Linear Algebra 模块

稀疏线性代数求解器：单移和多移 Krylov 子空间方法，用于费米子行列式计算中的大稀疏矩阵求逆。

## 文件清单

### 基类

| 文件 | 路径 | 说明 |
|------|------|------|
| `CSLASolver.h` | `SparseLinearAlgebra/CSLASolver.h` | 单移求解器抽象基类 |
| `CMultiShiftSolver.h` | `SparseLinearAlgebra/CMultiShiftSolver.h` | 多移求解器抽象基类 |

### 单移求解器

| 文件 | 路径 | 说明 |
|------|------|------|
| `CSLASolverBiCGStab.h` | `SparseLinearAlgebra/` | BiCGStab（双共轭梯度稳定化） |
| `CSLASolverGCR.h` | `SparseLinearAlgebra/` | GCR/ORTHODIR |
| `CSLASolverGMRES.h` | `SparseLinearAlgebra/` | GMRES（广义最小残差） |
| `CSLASolverGCRODR.h` | `SparseLinearAlgebra/` | GCRO-DR（带收缩重启的 GMRES） |
| `CSLASolverGMRESMDR.h` | `SparseLinearAlgebra/` | GMRES-MDR（最小收缩重启变体） |
| `CSolverCG.h` | `SparseLinearAlgebra/` | CG（仅 Hermitian 正定 `EFO_F_DDdagger`；`D`/`Ddagger` 通过 `DDdagger` 组合求解） |
| `CSolverDeflatedCG.h` | `SparseLinearAlgebra/` | Deflated CG（Arnoldi 收缩子空间 + 粗网格校正） |
| `CSolverTFQMR.h` | `SparseLinearAlgebra/` | TFQMR（无转置准最小残差） |

### 多移求解器

| 文件 | 路径 | 说明 |
|------|------|------|
| `CMultiShiftBiCGStab.h` | `SparseLinearAlgebra/` | 多移 BiCGStab |
| `CMultiShiftCG.h` | `SparseLinearAlgebra/` | 多移 CG（仅 Hermitian 正定，`EFO_F_DDdagger`） |
| `CMultiShiftFOM.h` | `SparseLinearAlgebra/` | 多移 FOM（完全正交化方法） |
| `CMultiShiftGMRES.h` | `SparseLinearAlgebra/` | 多移 GMRES |
| `CMultiShiftNested.h` | `SparseLinearAlgebra/` | 嵌套多移求解器 |

## CSLASolver（单移求解器基类）

**文件**: `SparseLinearAlgebra/CSLASolver.h`

所有单移求解器的统一接口：

| 方法 | 说明 |
|------|------|
| `Configurate(param)` | 从参数配置求解器（纯虚） |
| `AllocateBuffers(pField)` | 按场类型分配工作缓冲区（纯虚） |
| `Solve(pFieldX, pFieldB, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, uiM, ePhase, pStart)` | 求解 `x = A^{-1}b`（纯虚） |
| `GetInfos(tab)` | 求解器信息字符串 |
| `IsAbsoluteAccuracy()` | 是否使用绝对精度 |

**参数**：
- `gaugeNum/bosonNum/tensor2Num` 与 `gaugeFields/bosonFields/tensor2Fields` — 环境字段束（规范场、玻色场、tensor2 场；tensor2 暂未被算子使用，为将来依赖 tensor2 的算子预留）
- `uiM` — `EFieldOperator`，指定要应用的算子（`EFO_F_D`, `EFO_F_DD`, `EFO_F_InverseD` 等）
- `ePhase` — `ESolverPhase`（`ESP_Once`, `ESP_InTrajectory`, `ESP_AfterTrajectory`）
- `pStart` — 初始猜测解（可选，利用上次解的热启动）

## CMultiShiftSolver（多移求解器基类）

**文件**: `SparseLinearAlgebra/CMultiShiftSolver.h`

多移求解：对多个位移 `(A + c_n) x_n = b` 同时求解，共享 Krylov 向量。

| 方法 | 说明 |
|------|------|
| `Solve(pFieldX, cn, pFieldB, ...)` | `cn` 为位移数组，`pFieldX` 返回所有解 |

## 单移求解器

### BiCGStab

**文件**: `SparseLinearAlgebra/CSLASolverBiCGStab.h`

双共轭梯度稳定化方法，适合非对称系统。

- `m_uiStepCount` — 最大迭代步数
- `m_uiDevationCheck` — 偏差检查间隔
- `m_uiReTry` — 失败重试次数
- 支持单精度和双精度精度参数

### GCR / ORTHODIR

**文件**: `SparseLinearAlgebra/CSLASolverGCR.h`

GCR（广义共轭残差）实现。注意：实际实现为 ORTHODIR，因为 GCR 和 ORTHOMIN 有振荡问题。

- `m_uiMaxDim` — Krylov 子空间最大维度
- `m_uiReStart` — 重启周期

### GMRES

**文件**: `SparseLinearAlgebra/CSLASolverGMRES.h`

广义最小残差法，构建 Arnoldi 基并在 Hessenberg 矩阵上求解最小二乘问题。

- `m_uiMaxDim` — 子空间维度（默认 `_kMaxStep = 100`）
- `m_uiReStart` — 重启步数
- `RotateH()` — Givens 旋转更新 Hessenberg
- `SolveY()` — 求解上三角系统得到系数

### GCRO-DR

**文件**: `SparseLinearAlgebra/CSLASolverGCRODR.h`

GMRES 带收缩重启（Deflated Restarting），在重启之间回收近似特征向量以加速收敛。

- `m_uiMDim` — GMRES 子空间维度
- `m_uiKDim` — 回收的特征向量数
- `m_eDeflationType` — 收缩类型：`EEDT_REV`（特征值）、`EEDT_HEV`（调和特征值）、`EEDT_SVD`
- `GenerateCUFirstTime()` — 初始特征向量生成（虚函数，MDR 覆盖）
- `GenerateCU()` — 更新回收子空间

### GMRES-MDR

**文件**: `SparseLinearAlgebra/CSLASolverGMRESMDR.h`

GCRO-DR 的最小收缩重启变体。覆盖 `GenerateCUFirstTime()` 使用 QR 分解策略。

### TFQMR

**文件**: `SparseLinearAlgebra/CSolverTFQMR.h`

无转置准最小残差法。类似 BiCGStab 但稍慢。

## 多移求解器

### 多移 BiCGStab

**文件**: `SparseLinearAlgebra/CMultiShiftBiCGStab.h`

最大步数 `_kMaxStep = 100`。

### 多移 FOM

**文件**: `SparseLinearAlgebra/CMultiShiftFOM.h`

完全正交化方法的多移版本。构建 Arnoldi 基后直接求解各位移的 Hessenberg 系统。

- `RotateHSolveY()` — 对每个位移执行旋转和回代
- `RotateH()` / `SolveY()` — 静态辅助函数

### 多移 GMRES

**文件**: `SparseLinearAlgebra/CMultiShiftGMRES.h`

最大步数 `_kMaxStep = 30`。

### 多移 CG

**文件**: `SparseLinearAlgebra/CMultiShiftCG.h`

多移共轭梯度，仅支持 Hermitian 正定算子（`EFO_F_DDdagger`，即 RHMC 的算子）与实位移。所有位移共享同一 Krylov 空间：每次迭代仅一次 `DDdagger` 算子作用，位移解由种子（无位移）CG 递推经 zeta 递推获得（Jegerlehner, hep-lat/9608029）。

- 参数：`MaxStep`（默认 5000）、`Accuracy`、`AbsoluteAccuracy`
- 其他算子或复位移会直接 `appCrucial`

### 嵌套多移

**文件**: `SparseLinearAlgebra/CMultiShiftNested.h`

用于质量预条件系统：当轻夸克被预条件化为 `(x + small/x)^{1/2}` 形式时，通常只需一次求逆即可近似。包装一个内层 `CSLASolver` 而非执行真正的多移求逆。

- `m_pNestedSolver` — 内层单移求解器

## 类层次

```
CBase
├── CSLASolver
│   ├── CSLASolverBiCGStab
│   ├── CSLASolverGCR
│   ├── CSLASolverGMRES
│   ├── CSLASolverGCRODR
│   │   └── CSLASolverGMRESMDR
│   └── CSolverTFQMR
└── CMultiShiftSolver
    ├── CMultiShiftBiCGStab
    ├── CMultiShiftFOM
    ├── CMultiShiftCG
    ├── CMultiShiftGMRES
    └── CMultiShiftNested
```

## 使用场景

| 求解器 | 适用场景 |
|--------|---------|
| BiCGStab | 通用非对称系统，默认首选 |
| GMRES |  BiCGStab 不稳定时 |
| GCRO-DR | 序列相关系统（如 HMC 连续轨迹） |
| TFQMR | BiCGStab 失败时的备选 |
| 多移 BiCGStab | RHMC 有理逼近（多极点同时求逆） |
| 多移 FOM | 多移系统的替代方案 |
| 嵌套多移 | 质量预条件轻夸克 |

## 关键设计模式

- **场抽象**：求解器操作 `CField*` 而非原始向量，自动适应 Wilson、Staggered、Boson 等场类型。
- **热启动**：`pStart` 参数允许传入上次解作为初始猜测，在 HMC 中至关重要（规范场变化缓慢，解也变化缓慢）。
- **精度分离**：单精度构建时 `m_fAccuracy` 为 `DOUBLE`，确保混合精度收敛判断的可靠性。
- **设备端小矩阵**：GMRES/GCRO-DR 可选择用 CUDA 或 host 处理 Hessenberg 矩阵运算。
