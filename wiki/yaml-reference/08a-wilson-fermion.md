> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 08a. Wilson 费米子参数

## Wilson 费米子参数

### CFieldFermionWilsonSquareSU3 基础参数

**来源文件**：`Code/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.cu`

**读取函数**：`CFieldFermionWilsonSquareSU3::InitialOtherParameters()` (lines 638-655)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Hopping` | DOUBLE | `0.125` | **强烈建议** | `CFieldFermionWilsonSquareSU3.cu` ~L642 | 跳迁参数（hopping / kappa）。代码默认 `0.125`，但建议显式设置；若设为 0 会报 `appCrucial` 警告 |

### Wilson + 电磁场参数

**来源文件**：`Code/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3EM.cu`

**读取函数**：`CFieldFermionWilsonSquareSU3EM::InitialOtherParameters()` (lines 579-595)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Charge` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3EM.cu` ~L583 | 费米子电荷 |
| `EMFieldID` | INT | `0` | 否 | `CFieldFermionWilsonSquareSU3EM.cu` ~L589 | U(1) 电磁场 ID |

**适用类**：`CFieldFermionWilsonSquareSU3EM`, `CFieldFermionWilsonSquareCloverEMSU3`

### Wilson + Dirichlet + 旋转参数

**来源文件**：`Code/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DR.cu`

**读取函数**：`CFieldFermionWilsonSquareSU3DR::InitialOtherParameters()` (lines 2081-2102)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Naive` | INT | `1` | 否 | `CFieldFermionWilsonSquareSU3DR.cu` ~L2085 | 使用 naive 旋转实现。`1`=启用 |
| `Exponential` | INT | `1` | 否 | `CFieldFermionWilsonSquareSU3DR.cu` ~L2089 | 使用指数形式旋转。`1`=启用 |
| `ShiftCenter` | INT | `0` | 否 | `CFieldFermionWilsonSquareSU3DR.cu` ~L2093 | 旋转中心偏移。`1`=启用 |
| `Omega` | DOUBLE | `0.0` | 否 | `CFieldFermionWilsonSquareSU3DR.cu` ~L2097 | 旋转角速度 |

**适用类**：`CFieldFermionWilsonSquareSU3DR`, `CFieldFermionWilsonSquareCloverSU3DR`

### Wilson + 刚性加速参数

**来源文件**：`Code/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3DRigidAcc.cu`

**读取函数**：`CFieldFermionWilsonSquareSU3DRigidAcc::InitialOtherParameters()` (lines 445-453)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `G2` | Real | `0.2` | 否 | `CFieldFermionWilsonSquareSU3DRigidAcc.cu` ~L449 | 刚性加速参数 G^2 |

**适用类**：`CFieldFermionWilsonSquareSU3DRigidAcc`

**注意**：类似有限化学势，存在符号问题。

### Wilson + Gamma 矩阵凝聚参数

**来源文件**：`Code/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3Gamma.cu`

**读取函数**：`CFieldFermionWilsonSquareSU3Gamma::InitialOtherParameters()` (lines 679-749)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ExpGamma` | INT | `1` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L683 | 以 `gamma4.Exp(gamma)` 形式施加。`1`=启用 |
| `Gamma1` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L689 | Gamma1 系数 |
| `Gamma2` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L694 | Gamma2 系数 |
| `Gamma3` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L698 | Gamma3 系数 |
| `Gamma4` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L702 | Gamma4 系数 |
| `Gamma5` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L706 | Gamma5 系数 |
| `Gamma51` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L710 | Gamma51 系数 |
| `Gamma52` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L714 | Gamma52 系数 |
| `Gamma53` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L718 | Gamma53 系数 |
| `Gamma54` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L722 | Gamma54 系数 |
| `Sigma12` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L726 | Sigma12 系数 |
| `Sigma13` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L730 | Sigma13 系数 |
| `Sigma14` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L734 | Sigma14 系数 |
| `Sigma23` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L738 | Sigma23 系数 |
| `Sigma24` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L742 | Sigma24 系数 |
| `Sigma34` | Real | `0.0` | 否 | `CFieldFermionWilsonSquareSU3Gamma.cu` ~L746 | Sigma34 系数 |

**适用类**：`CFieldFermionWilsonSquareSU3Gamma`

### Clover 参数

**来源文件**：`Code/CLGLib/Data/Field/WilsonDirac/CFieldFermionWilsonSquareCloverSU3.cpp`

**读取函数**：`CFieldFermionWilsonSquareCloverSU3::InitialOtherParameters()` (lines 115-119)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Csw` | DOUBLE | `1.0` | 否 | `CFieldFermionWilsonSquareCloverSU3.cpp` ~L118 | Clover 系数（Sheikholeslami-Wohlert） |

**适用类**：所有 Clover 类（`CFieldFermionWilsonSquareCloverSU3`, `CFieldFermionWilsonSquareCloverEMSU3`, `CFieldFermionWilsonSquareCloverSU3D`, `CFieldFermionWilsonSquareCloverSU3DR`）

---


---

## 依赖关系

使用 Wilson / Clover 费米子场时，通常需要同时配置以下对象：

| 若使用 A | 则必须同时配置 B | 说明 |
|---|---|---|
| `CFieldFermionWilsonSquareCloverSU3` 及其变体（`...CloverSU3D`、`...CloverSU3DR`、`...CloverEMSU3`） | `StapleCache: CStapleCacheSU3` | Clover 核函数需要 staple cache 提供的 `Fmunu`，否则 `_FAIL_EXIT`。见 [13-staple-cache.md](13-staple-cache.md) |
| `CFieldFermionWilsonSquareSU3EM` / `...CloverEMSU3` | `CFieldGaugeU1Real`（`EMFieldID` 指向的 U(1) 场） | 电磁背景场。见 [05-gauge-fields.md](05-gauge-fields.md) |
| StoutLink Clover Wilson HMC（`...CloverSU3` + `CActionFermionKSImproveCombined`） | `GaugeSmearing: CGaugeSmearingStoutSU3` | 改进作用量通过 smearing 取 effective gauge。见 [09b-fermion-actions.md](09b-fermion-actions.md)、[12-gauge-smearing.md](12-gauge-smearing.md) |
| Wilson 费米子测量 | `Solver`（`SolverForFieldId` 匹配） | 测量反演 `D` / `D†D`。见 [10-solvers.md](10-solvers.md) |
| Wilson 费米子 HMC | `Solver` 或 `MSSolver` | `Rational: ER_NoRational` 用 `Solver`；有理逼近用 `MSSolver`。见 [10-solvers.md](10-solvers.md)、[11-multi-shift-solvers.md](11-multi-shift-solvers.md) |

---

[< 返回 08. 费米子场](08-fermion-fields.md) | [< 返回 yaml-reference 目录](home.md)
