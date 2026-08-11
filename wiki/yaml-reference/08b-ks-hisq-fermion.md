> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 08b. KS / HISQ 费米子参数

## Staggered (KS) 费米子参数

### KS 基础参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKS.cpp`

**读取函数**：`CFieldFermionKS::InitialOtherParameters()` (lines 15-30)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Mass` | Real | `0.01` | **强烈建议** | `CFieldFermionKS.cpp` ~L19 | 裸质量参数 `m_f2am = 2 * a * m`。若为 0 会报 `appCrucial` 警告 |
| `EachSiteEta` | INT | `0` | 否 | `CFieldFermionKS.cpp` ~L25 | 是否使用每点独立的 eta。`1`=启用 |

**适用类**：所有 KS 类（`CFieldFermionKSSU2` ~ `CFieldFermionKSSU4`, `CFieldFermionKSU1`, `CFieldFermionKSU1D`, `CFieldFermionKSSU3D`, 等）

**注意**：`CFieldFermionKS::InitialOtherParameters()` 第 29 行无条件设置 `m_eRational = ER_AllRational`，覆盖基类读取的 `Rational` 值。

### KS + 旋转参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKSTR.cpp`

**读取函数**：`CFieldFermionKSTR::InitialOtherParameters()` (lines 134-156)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `RealRotation` | INT | `0` | 否 | `CFieldFermionKSTR.cpp` ~L139 | 使用实数旋转（vs 虚数）。`1`=启用 |
| `Omega` | DOUBLE | `0.1` | 否 | `CFieldFermionKSTR.h` ~L145 | 旋转角速度。非零旋转必须显式设置 |
| `CachedGauge` | INT | `0` | 否 | `CFieldFermionKSTR.cpp` ~L151 | 使用缓存规范场进行旋转。`1`=启用 |

**适用类**：`CFieldFermionKSU1R`, `CFieldFermionKSSU3R`, `CFieldFermionKSSU3DR`, `CFieldFermionHISQSU3R`

### KS + Gamma 矩阵凝聚参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3Gamma.cpp`

**读取函数**：`CFieldFermionKSSU3Gamma::InitialOtherParameters()` (lines 2037-2112)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Imagine` | INT | `1` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2041 | gamma1-4 以虚数形式施加。`1`=启用 |
| `Gamma1` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2051 | Gamma1 系数 |
| `Gamma2` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2056 | Gamma2 系数 |
| `Gamma3` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2060 | Gamma3 系数 |
| `Gamma4` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2064 | Gamma4 系数 |
| `Gamma5` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2068 | Gamma5 系数 |
| `Gamma51` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2072 | Gamma51 系数 |
| `Gamma52` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2076 | Gamma52 系数 |
| `Gamma53` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2080 | Gamma53 系数 |
| `Gamma54` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2084 | Gamma54 系数 |
| `Sigma12` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2088 | Sigma12 系数 |
| `Sigma13` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2092 | Sigma13 系数 |
| `Sigma14` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2096 | Sigma14 系数 |
| `Sigma23` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2100 | Sigma23 系数 |
| `Sigma24` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2104 | Sigma24 系数 |
| `Sigma34` | Real | `0.0` | 否 | `CFieldFermionKSSU3Gamma.cpp` ~L2108 | Sigma34 系数 |

**适用类**：`CFieldFermionKSSU3Gamma`, `CFieldFermionKSSU3GammaEM`

### KS + Gamma + 电磁场参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3GammaEM.cpp`

**读取函数**：`CFieldFermionKSSU3GammaEM::InitialOtherParameters()` (lines 2598-2612)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Charge` | Real | `0.0` | 否 | `CFieldFermionKSSU3GammaEM.cpp` ~L2602 | 电荷 |
| `EMFieldID` | INT | `0` | 否 | `CFieldFermionKSSU3GammaEM.cpp` ~L2607 | EM 场 ID |

**适用类**：`CFieldFermionKSSU3GammaEM`

### KS + 旋转 + 电磁场参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3REM.cpp`

**读取函数**：`CFieldFermionKSSU3REM::InitialOtherParameters()` (lines 1080-1102)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Charge` | Real | `0.0` | 否 | `CFieldFermionKSSU3REM.cpp` ~L1085 | 电荷 |
| `EMFieldID` | INT | `0` | 否 | `CFieldFermionKSSU3REM.cpp` ~L1090 | EM 场 ID |
| `Omega` | DOUBLE | `0.1` | 否 | `CFieldFermionKSSU3REM.cpp` ~L1096 | 旋转角速度 |

**适用类**：`CFieldFermionKSSU3REM`

**注意**：`m_bEachSiteEta` 被强制设为 `TRUE`（line 1083）。仅支持磁场。

### KS + 刚性加速参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKSSU3RigidAcc.cu`

**读取函数**：`CFieldFermionKSSU3RigidAcc::InitialOtherParameters()` (lines 301-311)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ImaginaryGamma3` | INT | `1` | 否 | `CFieldFermionKSSU3RigidAcc.cu` ~L305 | 使用虚数 gamma3。`1`=启用 |

**适用类**：`CFieldFermionKSSU3RigidAcc`

**注意**：存在符号问题。质量项是对角的，**不支持**嵌套多移求解器或质量预条件子。

---

## HISQ 费米子参数

### HISQ 基础参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKSHISQ.cpp`

**读取函数**：`CFieldFermionHISQSU3::InitialOtherParameters()` (lines 110-134)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Naik` | Real | 自动计算 | 否 | `CFieldFermionKSHISQ.cpp` ~L114 | Naik 项系数 |
| `Epsilon` | Real | 自动计算 | 否 | `CFieldFermionKSHISQ.cpp` ~L115 | HISQ Epsilon 参数 |

**自动计算逻辑**（lines 116-134）：
- 若指定 `Epsilon` 但未指定 `Naik`：`Naik = -(1 + Epsilon) / 24`
- 若两者都未指定：Epsilon 由质量通过级数展开自动计算：
  - `epsilon = -0.16875 * m^2 + 0.018247767857142858 * m^4 - 0.0009072149367559524 * m^6 - 0.00007302123230773134 * m^8`
  - `Naik = -(1 + epsilon) / 24`

**适用类**：`CFieldFermionHISQSU3`, `CFieldFermionHISQSU3R`, `CFieldFermionHISQWithPhaseSU3`

### aHISQ（各向异性 HISQ）参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKSHISQAnisotropic.cpp`

**读取函数**：`CFieldFermionHISQSU3Anisotropic::InitialOtherParameters()`

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `XiF` | Real | `1.0` | 否 | `CFieldFermionKSHISQAnisotropic.cpp` ~L22 | 裸费米子各向异性 $\xi_0^f$，必须为正 |

**适用类**：`CFieldFermionHISQSU3Anisotropic`（继承 `CFieldFermionHISQSU3`，HISQ 基础参数全部适用）

**说明**：
- `XiF` 是**裸的可调输入**，不是测量/重整化各向异性 $\xi=a_\sigma/a_\tau$，二者一般不同；`XiF` 与规范侧 `Xi` 相互独立，不会自动绑定。
- `XiF` 只乘在时间方向（dir=3）的 hopping 贡献上：外层 Fat7 one-link、时间方向 Naik 三链及 epsilon 修正项、以及对应的 RHMC 连接/外积；不乘质量项，不改 Fat3/Fat5/Fat7/Lepage 系数，不动 smearing 缓存（`CGaugeSmearingHISQSU3` 原样复用）。
- 与 `CActionFermionHISQCombined` 配合使用；多个不同 `XiF` 的场可共存于同一 action。

### HISQ + Phase 场参数

**来源文件**：`Code/CLGLib/Data/Field/Staggered/CFieldFermionKSHISQWithPhase.cu`

**读取函数**：`CFieldFermionHISQWithPhaseSU3::InitialOtherParameters()` (lines 335-342)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Charge` | Real | `0.0` | 否 | `CFieldFermionKSHISQWithPhase.cu` ~L338 | 电荷 |
| `U1FieldId` | INT | `0` | 否 | `CFieldFermionKSHISQWithPhase.cu` ~L339 | U(1) Phase 场 ID |

**适用类**：`CFieldFermionHISQWithPhaseSU3`

---


---

## 依赖关系

使用 KS / HISQ 费米子场时，通常需要同时配置以下对象：

| 若使用 A | 则必须同时配置 B | 说明 |
|---|---|---|
| `CFieldFermionHISQSU3` / `...HISQSU3R` / `...HISQWithPhaseSU3` | `GaugeSmearing: CGaugeSmearingHISQSU3` | HISQ 核函数直接通过 `appGetGaugeSmearing` 取 effective gauge / Naik link。见 [12-gauge-smearing.md](12-gauge-smearing.md) |
| `CFieldFermionHISQWithPhaseSU3` | 上述 HISQ smearing + `CFieldGaugeU1Real`（`U1FieldId`） | 需要外部 U(1) phase 场。见 [05-gauge-fields.md](05-gauge-fields.md) |
| `CFieldFermionKSSU3R` 等旋转 KS 场（`CachedGauge: 1`） | `StapleCache: CStapleCacheSU3` | 旋转矩阵缓存取自 staple cache。`CachedGauge` 默认 `0`，但项目示例通常设为 `1`。见 [13-staple-cache.md](13-staple-cache.md) |
| `CFieldFermionKSSU3GammaEM` / `...KSSU3REM` | `CFieldGaugeU1Real`（`EMFieldID`） | 电磁背景场。见 [05-gauge-fields.md](05-gauge-fields.md) |
| HISQ / improved KS HMC | `MSSolver`（每个相关场一个） | `CalculateForce` 调用多移求解器做有理逼近。见 [11-multi-shift-solvers.md](11-multi-shift-solvers.md) |
| HISQ / KS 费米子测量 | `Solver`（`SolverForFieldId` 匹配） | 测量反演 `D` / `D†D`。见 [10-solvers.md](10-solvers.md) |

---

[< 返回 08. 费米子场](08-fermion-fields.md) | [< 返回 yaml-reference 目录](home.md)
