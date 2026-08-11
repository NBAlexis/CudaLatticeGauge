> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 08. 费米子场 (Fermion Fields)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateFermionFields()` (lines 598-663)

---

## 配置参数

YAML 键名为 `FermionField1`, `FermionField2`, ...。通过 `FermionFieldCount` 控制数量。

### 基类参数 (CField)

**来源文件**：`Code/CLGLib/Data/Field/CField.h`

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldName` | string | `CFieldFermionWilsonSquareSU3` | 否 | `CLGLibManager.cpp` ~L618 | 费米子场类名。通过 `appCreate()` 实例化 |
| `FieldInitialType` | string | `EFIT_Random` | 否 | `CLGLibManager.cpp` ~L620 | 初始场类型 |
| `FieldId` | INT | 自动分配 | 否 | `CLGLibManager.cpp` ~L622 | 场唯一标识符。必须 `> 1`。冲突时自动分配 |
| `Period` | INT[] | `[1,1,1,-1]` | 否 | `CField.cpp` ~L111 | 边界条件。默认时间方向反周期 (`-1`) |
| `PoolNumber` | INT | `0` | 否 | `CField.cpp` ~L115 | 辅助场池大小。费米子 HMC 需要较大值 |
| `GaugeFields` | BYTE[] | `[1]` | 否 | `CField.h` ~L114 | 耦合的规范场 ID 列表 |
| `BosonFields` | BYTE[] | `[]` | 否 | `CField.h` ~L115 | 耦合的玻色子场 ID 列表 |
| `Dynamic` | INT | `1` | 否 | `CField.h` ~L122 | 是否为动态场。`1`=动态，`0`=静态 |

### 费米子基类参数 (CFieldFermion)

**来源文件**：`Code/CLGLib/Data/Field/CFieldFermion.cpp`

**读取函数**：`CFieldFermion::InitialOtherParameters()` (lines 379-407)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Even` | INT | `0` | 否 | `CFieldFermion.cpp` ~L383 | 是否使用偶点伪费米子。`1`=启用 |
| `MD` | Real[] | `[1.0]` | 否 | `CFieldFermion.cpp` ~L387 | 分子动力学有理近似系数 |
| `MC` | Real[] | `[1.0]` | 否 | `CFieldFermion.cpp` ~L391 | Monte Carlo 有理近似系数 |
| `Rational` | string | `ER_NoRational` | 否 | `CFieldFermion.cpp` ~L395 | 有理近似类型 |

**注意**：对于 KS 费米子，`m_eRational` 在 `CFieldFermionKS::InitialOtherParameters()` 第 29 行被**强制覆盖**为 `ER_AllRational`，YAML 中的 `Rational` 设置无效。

---

## 支持的费米子场类

以下所有类均通过 `appCreate(FieldName)` 工厂实例化。需确认 `Code/CLGLib/Core/CLGSetup.h` 中对应宏已启用。

### Wilson 费米子

| 类名 | 宏开关 | 说明 |
|------|--------|------|
| `CFieldFermionWilsonSquareSU3` | `_CLG_WILSON_DIRAC` | 标准 Wilson Dirac 费米子 |
| `CFieldFermionWilsonSquareSU3Acc` | `_CLG_WILSON_DIRAC` | 加速 Wilson 费米子（Galilean 变换） |
| `CFieldFermionWilsonSquareSU3Boost` | `_CLG_WILSON_DIRAC` | Boosted Wilson 费米子 |
| `CFieldFermionWilsonSquareSU3D` | `_CLG_WILSON_DIRAC` | Dirichlet 边界 Wilson 费米子 |
| `CFieldFermionWilsonSquareSU3DR` | `_CLG_WILSON_DIRAC` | Dirichlet + 旋转 Wilson 费米子 |
| `CFieldFermionWilsonSquareSU3DRigidAcc` | `_CLG_WILSON_DIRAC` | Dirichlet + 刚性加速 Wilson 费米子 |
| `CFieldFermionWilsonSquareSU3EM` | `_CLG_WILSON_DIRAC` | 电磁场耦合 Wilson 费米子 |
| `CFieldFermionWilsonSquareSU3Gamma` | `_CLG_WILSON_DIRAC` | Gamma 矩阵凝聚 Wilson 费米子 |

### Clover 费米子

| 类名 | 宏开关 | 说明 |
|------|--------|------|
| `CFieldFermionWilsonSquareCloverSU3` | `_CLG_WILSON_DIRAC` | Clover 改进 Wilson 费米子 |
| `CFieldFermionWilsonSquareCloverEMSU3` | `_CLG_WILSON_DIRAC` | Clover + 电磁场耦合 |
| `CFieldFermionWilsonSquareCloverSU3D` | `_CLG_WILSON_DIRAC` | Clover + Dirichlet 边界 |
| `CFieldFermionWilsonSquareCloverSU3DR` | `_CLG_WILSON_DIRAC` | Clover + Dirichlet + 旋转 |

### Staggered (KS) 费米子

| 类名 | 宏开关 | 说明 |
|------|--------|------|
| `CFieldFermionKSU1` | `_CLG_STAGGERED_DIRAC` | U(1) KS 费米子 |
| `CFieldFermionKSSU2` | `_CLG_STAGGERED_DIRAC` | SU(2) KS 费米子 |
| `CFieldFermionKSSU3` | `_CLG_STAGGERED_DIRAC` | SU(3) KS 费米子（标准） |
| `CFieldFermionKSSU4` | `_CLG_SU4_KS` | SU(4) KS 费米子 |
| `CFieldFermionKSU1D` | `_CLG_STAGGERED_DIRAC` | U(1) KS + Dirichlet 边界 |
| `CFieldFermionKSSU3D` | `_CLG_STAGGERED_DIRAC` | SU(3) KS + Dirichlet 边界 |
| `CFieldFermionKSU1R` | `_CLG_STAGGERED_DIRAC` | U(1) KS + 旋转 |
| `CFieldFermionKSSU3R` | `_CLG_STAGGERED_DIRAC` | SU(3) KS + 旋转 |
| `CFieldFermionKSSU3DR` | `_CLG_STAGGERED_DIRAC` | SU(3) KS + Dirichlet + 旋转 |
| `CFieldFermionKSSU3Acc` | `_CLG_STAGGERED_DIRAC` | 加速 KS 费米子 |
| `CFieldFermionKSSU3Gamma` | `_CLG_STAGGERED_DIRAC` | KS + Gamma 矩阵凝聚 |
| `CFieldFermionKSSU3GammaEM` | `_CLG_STAGGERED_DIRAC` | KS Gamma + 电磁场耦合 |
| `CFieldFermionKSSU3REM` | `_CLG_STAGGERED_DIRAC` | KS + 旋转 + 电磁场 |
| `CFieldFermionKSSU3RigidAcc` | `_CLG_STAGGERED_DIRAC` | KS + 刚性加速 |

### HISQ 费米子

| 类名 | 宏开关 | 说明 |
|------|--------|------|
| `CFieldFermionHISQSU3` | `_CLG_STAGGERED_DIRAC` | HISQ 改进 KS 费米子 |
| `CFieldFermionHISQSU3R` | `_CLG_STAGGERED_DIRAC` | HISQ + 旋转 |
| `CFieldFermionHISQWithPhaseSU3` | `_CLG_STAGGERED_DIRAC` | HISQ + U(1) Phase 场耦合 |

---

## 参数继承关系速查

### Wilson 家族参数继承

```
CFieldFermionWilsonSquareSU3
  ├── Hopping
  ├── CFieldFermionWilsonSquareSU3Acc (无额外参数)
  ├── CFieldFermionWilsonSquareSU3Boost (无额外参数)
  ├── CFieldFermionWilsonSquareSU3EM
  │     ├── Charge, EMFieldID
  │     └── CFieldFermionWilsonSquareCloverEMSU3 (+ Csw)
  ├── CFieldFermionWilsonSquareSU3D
  │     └── CFieldFermionWilsonSquareSU3DR
  │           ├── Naive, Exponential, ShiftCenter, Omega
  │           └── CFieldFermionWilsonSquareCloverSU3DR (+ Csw)
  │     └── CFieldFermionWilsonSquareSU3DRigidAcc
  │           └── G2
  └── CFieldFermionWilsonSquareSU3Gamma
        ├── ExpGamma, Gamma1..Gamma5, Gamma51..Gamma54, Sigma12..Sigma34
        └── (基于 SU3D，但 Gamma 参数在 SU3Gamma 中定义)

CFieldFermionWilsonSquareCloverSU3
  └── Csw (模板参数，所有 Clover 子类都继承)
```

### KS 家族参数继承

```
CFieldFermionKS
  ├── Mass, EachSiteEta (强制 Rational = ER_AllRational)
  ├── CFieldFermionKSU1 / KSSU2 / KSSU3 / KSSU4
  ├── CFieldFermionKSU1D / KSSU3D (Dirichlet)
  ├── CFieldFermionKSU1R / KSSU3R / KSSU3DR (Rotation)
  │     ├── RealRotation, Omega, CachedGauge
  ├── CFieldFermionKSSU3Acc (无额外参数)
  ├── CFieldFermionKSSU3Gamma
  │     ├── Imagine, Gamma1..Gamma5, Gamma51..Gamma54, Sigma12..Sigma34
  │     └── CFieldFermionKSSU3GammaEM (+ Charge, EMFieldID)
  ├── CFieldFermionKSSU3REM
  │     ├── Charge, EMFieldID, Omega (EachSiteEta 强制为 TRUE)
  ├── CFieldFermionKSSU3RigidAcc
  │     └── ImaginaryGamma3
  ├── CFieldFermionHISQSU3 / HISQSU3R
  │     ├── Naik, Epsilon (自动计算)
  │     └── CFieldFermionHISQWithPhaseSU3 (+ Charge, U1FieldId)
```

---

## 常见依赖关系

| 若配置此类费米子场 | 必须同时配置 | 说明 |
|---|---|---|
| `CFieldFermionWilsonSquareCloverSU3` 及其变体 | `StapleCache: CStapleCacheSU3` | Clover 需要 `Fmunu`。见 [13-staple-cache.md](13-staple-cache.md)、[08a-wilson-fermion.md](08a-wilson-fermion.md) |
| `CFieldFermionKSSU3R` 等旋转 KS 场（`CachedGauge: 1`） | `StapleCache: CStapleCacheSU3` | 旋转矩阵缓存。见 [13-staple-cache.md](13-staple-cache.md)、[08b-ks-hisq-fermion.md](08b-ks-hisq-fermion.md) |
| `CFieldFermionHISQSU3` / `...HISQSU3R` / `...HISQWithPhaseSU3` | `GaugeSmearing: CGaugeSmearingHISQSU3` | HISQ 直接取 smearing 结果。见 [12-gauge-smearing.md](12-gauge-smearing.md)、[08b-ks-hisq-fermion.md](08b-ks-hisq-fermion.md) |
| `CFieldFermionHISQWithPhaseSU3` | 上述 HISQ smearing + `CFieldGaugeU1Real` | 需要 U(1) phase 场。见 [05-gauge-fields.md](05-gauge-fields.md) |
| `CFieldFermionKSSU3GammaEM` / `...KSSU3REM` / Wilson EM 类 | `CFieldGaugeU1Real`（`EMFieldID`） | 电磁背景场。见 [05-gauge-fields.md](05-gauge-fields.md) |
| 任意费米子场参与 HMC | `MSSolver` 或 `Solver` | HISQ / improved KS / rational 用 `MSSolver`；`ER_NoRational` 用 `Solver`。见 [10-solvers.md](10-solvers.md)、[11-multi-shift-solvers.md](11-multi-shift-solvers.md) |
| 任意费米子场参与测量 | `Solver` | 测量反演。见 [10-solvers.md](10-solvers.md) |

---

## 子页面

- [08a. Wilson 费米子参数](08a-wilson-fermion.md)
- [08b. KS / HISQ 费米子参数](08b-ks-hisq-fermion.md)
- [08c. KS/HISQ 质量约定](08c-fermion-mass-convention.md)

[< 返回目录](home.md) | [< 上一章：玻色子场](07-boson-fields.md) | [下一章：作用量 >](09-actions.md)
