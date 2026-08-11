# Fermion Field 模块

费米子场实现：Wilson-Dirac、Staggered/Kogut-Susskind、Clover 改进、HISQ，以及各种变形（旋转、加速、boost、EM 耦合）。

## 文件清单

### 基类

| 文件 | 路径 | 说明 |
|------|------|------|
| `CFieldFermion.h` | `Data/Field/CFieldFermion.h` | 费米子场抽象基类 |

### Wilson Dirac 费米子

| 文件 | 路径 | 说明 |
|------|------|------|
| `CFieldFermionWilsonSquareSU3.h` | `Data/Field/WilsonDirac/` | 标准 Wilson 费米子 |
| `CFieldFermionWilsonSquareSU3D.h` | `Data/Field/WilsonDirac/` | Dirichlet 边界 |
| `CFieldFermionWilsonSquareSU3DR.h` | `Data/Field/WilsonDirac/` | Dirichlet + 旋转 |
| `CFieldFermionWilsonSquareSU3Acc.h` | `Data/Field/WilsonDirac/` | 加速参考系 |
| `CFieldFermionWilsonSquareSU3Boost.h` | `Data/Field/WilsonDirac/` | Boost 变换 |
| `CFieldFermionWilsonSquareSU3DRigidAcc.h` | `Data/Field/WilsonDirac/` | 刚性加速 |
| `CFieldFermionWilsonSquareSU3EM.h` | `Data/Field/WilsonDirac/` | EM (U1) 耦合 |
| `CFieldFermionWilsonSquareSU3Gamma.h` | `Data/Field/WilsonDirac/` | Gamma 矩阵插入 |
| `CFieldFermionWilsonSquareCloverSU3.h` | `Data/Field/WilsonDirac/` | Clover 改进 |
| `CFieldFermionWilsonKernel.h` | `Data/Field/WilsonDirac/` | Clover kernel 静态接口 |

### Staggered / KS 费米子

| 文件 | 路径 | 说明 |
|------|------|------|
| `CFieldFermionKS.h` | `Data/Field/Staggered/` | KS 费米子抽象基类 |
| `CFieldFermionKST.h` | `Data/Field/Staggered/` | 模板实现 + 多群实例 |
| `CFieldFermionKSTD.h` | `Data/Field/Staggered/` | Dirichlet 边界 |
| `CFieldFermionKSTR.h` | `Data/Field/Staggered/` | 旋转项 |
| `CFieldFermionKSSU3Acc.h` | `Data/Field/Staggered/` | 加速 |
| `CFieldFermionKSSU3Gamma.h` | `Data/Field/Staggered/` | Gamma 插入 |
| `CFieldFermionKSSU3GammaEM.h` | `Data/Field/Staggered/` | Gamma + EM |
| `CFieldFermionKSSU3REM.h` | `Data/Field/Staggered/` | 旋转 + EM |
| `CFieldFermionKSSU3RigidAcc.h` | `Data/Field/Staggered/` | 刚性加速 |
| `CFieldFermionKSHISQ.h` | `Data/Field/Staggered/` | HISQ 改进 |
| `CFieldFermionKSHISQWithPhase.h` | `Data/Field/Staggered/` | HISQ + 相位 |
| `CFieldFermionKSTKernel.h` | `Data/Field/Staggered/` | 静态 kernel 接口 |
| `CFieldFermionKSTKernelGamma.h` | `Data/Field/Staggered/` | Gamma kernel |
| `CFieldFermionKSTKernelR.h` | `Data/Field/Staggered/` | 旋转 kernel |

## CFieldFermion（抽象基类）

**文件**: `Data/Field/CFieldFermion.h`

所有费米子场的统一接口：

| 方法 | 说明 |
|------|------|
| `ApplyOperator(uiM, ...)` | 按 `EFieldOperator` 分发的统一接口 |
| `D()` / `Ddagger()` / `DD()` / `DDdagger()` | Dirac 算子（纯虚） |
| `InverseD()` / `InverseDdagger()` / `InverseDD()` | 求解器-based 逆算子 |
| `DWithMass()` / `DDdaggerWithMass()` | 带质量参数的变体 |
| `RationalApproximation()` | 有理逼近（RHMC） |
| `D0()` / `D0OnEvenOrOdd()` | 仅 hopping 的算子（用于有理力） |
| `CalculateForce()` | 费米子力计算 |
| `PrepareForHMC()` | HMC 前准备（随机伪费米子场） |
| `ApplyGamma()` | Gamma 矩阵作用 |
| `Energy()` | 能量计算 |
| `InitialAsSource()` | 初始化为点源/噪声源（纯虚） |
| `GetSourcesAtSiteFromPool()` | 从池中获取点源传播子（纯虚） |
| `CalculateF0AndNaik()` | HISQ Naik 力计算 |

**EFieldOperator 枚举**：`EFO_F_D`, `EFO_F_Ddagger`, `EFO_F_DD`, `EFO_F_DDdagger`, `EFO_F_InverseD`, `EFO_F_InverseDdagger`, `EFO_F_InverseDD`, `EFO_F_InverseDDdagger`, `EFO_F_RationalD` 等，以及带质量变体（`_WithMass` 后缀）。

## Wilson Dirac 费米子

### CFieldFermionWilsonSquareSU3

**文件**: `Data/Field/WilsonDirac/CFieldFermionWilsonSquareSU3.h`

标准 Wilson 费米子，每格点存储 `deviceWilsonVectorSU3`（4 自旋 × 3 颜色 = 24 复数）。

**关键方法**：
- `DS()`, `DdaggerS()`, `DDS()`, `DDdaggerS()` — Wilson Dirac 算子
- `DOperator()`, `DerivateDOperator()` — 底层 kernel 分发
- `SetKai()` / `GetKai()` — Wilson hopping 参数 κ
- `PrepareForHMCS()`, `CalculateForceS()`, `EnergyS()` — HMC 生命周期
- `ApplyGammaS()` — Gamma 矩阵作用
- `GetSourcesAtSiteFromPool()`, `InitialAsSource()` — 传播子源
- `TestGamma5Hermitian()` — γ5-厄米性测试

### 变形 Wilson 费米子

通过继承 `CFieldFermionWilsonSquareSU3` 并重写 `DOperator()` 和 `DerivateDOperator()` 实现：

| 类 | 变形内容 |
|----|---------|
| `CFieldFermionWilsonSquareSU3D` | Dirichlet 边界（`FixBoundary()`） |
| `CFieldFermionWilsonSquareSU3DR` | Dirichlet + 旋转（Omega 项，支持 naive/exponential/shift-center 模式） |
| `CFieldFermionWilsonSquareSU3Acc` | 匀加速参考系 |
| `CFieldFermionWilsonSquareSU3Boost` | 常数 Lorentz boost |
| `CFieldFermionWilsonSquareSU3DRigidAcc` | 刚性加速（Rindler），有符号问题 |
| `CFieldFermionWilsonSquareSU3EM` | 外电磁场 U(1) 耦合（`m_fCharge`, `m_byU1FieldId`） |
| `CFieldFermionWilsonSquareSU3Gamma` | 任意 Gamma 矩阵插入（用于凝聚测量，非模拟） |

### Clover 改进

**文件**: `Data/Field/WilsonDirac/CFieldFermionWilsonSquareCloverSU3.h`

模板包装器 `CFieldFermionWilsonSquareClover<CFieldWilson>`：为任意 Wilson 费米子添加 Clover（Sheikholeslami-Wohlert）项。

- `m_fCsw` — Clover 系数
- 使用 staple cache 的 `F_munu` 计算 clover 项
- 实例化：`CloverSU3`, `CloverEMSU3`, `CloverSU3D`, `CloverSU3DR`

## Staggered / KS 费米子

### CFieldFermionKS（抽象基类）

**文件**: `Data/Field/Staggered/CFieldFermionKS.h`

KS 费米子的抽象基类：

| 方法 | 说明 |
|------|------|
| `SetMass()` / `GetMass()` — `m_f2am` | 质量参数（2am，a 为格距） |
| `OnlyMass()` | 纯质量项（纯虚） |
| `OneLink()` / `OneLinkForce()` | 通用路径 hopping（单/多场分发） |
| `OneLinkS()` / `OneLinkForceS()` | 单场实现（纯虚） |
| `DOperatorKS()` / `DOperatorKSOnEvenOrOdd()` | KS Dirac 算子（纯虚） |
| `DerivateD0()` | hopping 力（纯虚） |
| `CalculateForceEvenOddS()` | even-odd 力 |
| `PrepareForHMCOnlyRandomize()` / `PrepareForHMCNotRandomize()` | HMC 准备（纯虚） |
| `TestAntiHermitian()` | 反厄米性测试 |

### CFieldFermionKST（模板实现）

**文件**: `Data/Field/Staggered/CFieldFermionKST.h`

模板类 `CFieldFermionKST<deviceVector, deviceGauge, vectorN>`：

- **设备向量类型**：`CLGComplex` (U1), `deviceSU2Vector` (SU2), `deviceSU3Vector` (SU3), `deviceSU4Vector` (SU4)
- **注册**：`CFieldFermionKSU1`, `CFieldFermionKSSU2`, `CFieldFermionKSSU3`, `CFieldFermionKSSU4`

**关键方法**：
- 完整线性代数（`Axpy`, `Dot`, `Mul`, `Dagger`, `ScalarMultply`）；even 伪费米子（`m_bEvenPseudofermion`）下 `Axpy` 系列只更新 even 站点（奇站点恒为零，见 `CCommonKernelField::Axpy*EvenOdd`）
- `DS()`, `DdaggerS()`, `DDS()`, `DDdaggerS()` — KS Dirac 算子
- `DOperatorKS()`, `DerivateD0()` — 委托给 `CFieldFermionKSTKernel`
- `PrepareForHMC()` — 支持 even-odd 伪费米子
- `CalculateForceS()`, `CalculateForceEvenOddS()` — 通过多移求解器计算力
- `ApplyGammaS()` — taste 重构的 Gamma 应用（通过 `CFieldFermionKSTKernelGamma`）
- `Connection()`, `ConnectionSelf()` — 双线性连接

### 变形 KS 费米子

| 类 | 变形内容 |
|----|---------|
| `CFieldFermionKSTD` | Dirichlet 边界（`_D` kernel 后缀） |
| `CFieldFermionKSTR` | 旋转项（实/虚旋转，支持缓存 gauge）；虚旋转在 Torus/Dirichlet 下支持 even 伪费米子（EO D 算子见 `DOperatorKSOnEvenOrOdd_R_ImaginaryRotation`，EO 旋转力经 `CalculateForceEvenOddS_SingleTermOfRationalR` 钩子复用全校 kernel；projective plane 与实旋转不支持） |
| `CFieldFermionKSSU3Acc` | 匀加速 |
| `CFieldFermionKSSU3Gamma` | Gamma 插入（taste 重构） |
| `CFieldFermionKSSU3GammaEM` | Gamma + EM 耦合 |
| `CFieldFermionKSSU3REM` | 旋转 + EM（仅磁场） |
| `CFieldFermionKSSU3RigidAcc` | 刚性加速，有符号问题，不支持嵌套多移 |

### HISQ

**文件**: `Data/Field/Staggered/CFieldFermionKSHISQ.h`

模板 `CFieldFermionHISQT<CFieldKS>`：Highly Improved Staggered Quark 作用量。

- 添加 Naik link（3-hop）和 Lepage 项
- 仅支持 even-odd
- `CalculateF0AndNaik()` — 分别计算 f0（标准力）、epsilon 项、Naik 力
- `m_fNaik`, `m_fEpsilon` — 系数（可从质量自动计算）

## Kernel 静态接口

### Wilson Kernel

**文件**: `Data/Field/WilsonDirac/CFieldFermionWilsonKernel.h`

- `DOperatorClover()` — clover 项应用
- `PrepareSigmaMunu()` — 构建 sigma_munu 双线性
- `CloverForce()` — clover 力

### KS Kernel

**文件**: `Data/Field/Staggered/CFieldFermionKSTKernel.h`

- `DOperatorKS()`, `DOperatorKS_D()`, `DOperatorKSOnEvenOrOdd()` — Dirac 算子
- `DerivateD0()`, `DerivateD0_D()` — hopping 力
- `OnlyMass()` — 质量项
- `OneLinkS()`, `OneLinkForceS()` — 通用路径
- `DOperatorEM()`, `KSForceEM()` — EM 耦合
- `DOperatorNaik()`, `NaikConnection()` — HISQ Naik link

**Gamma Kernel**（`CFieldFermionKSTKernelGamma.h`）：
- `appApplyGammaKS()`, `GammaKSForce()` — Gamma 应用和力
- `appApplyGammaKSEM()` — Gamma + EM

**Rotation Kernel**（`CFieldFermionKSTKernelR.h`）：
- `DOperatorKS_R_RealRotation()`, `DOperatorKS_R_ImaginaryRotation()` — 旋转 Dirac 算子
- `DerivateD0_R()`, `DerivateD0_REM()` — 旋转力

## 关键设计模式

- **模板化设备向量**：`CFieldFermionKST` 通过 `deviceVector`/`deviceGauge` 模板支持多规范群。
- **EFieldOperator 分发**：基类 `ApplyOperator` 按枚举分发到具体虚函数，使求解器和作用量无需知道具体场类型。
- **Even-Odd 分解**：KS 费米子使用 even-odd checkerboard 分解，D 算子只在同奇偶性子格点间跳跃。
- **Taste 重构**：staggered 费米子没有显式自旋指标，Gamma 矩阵通过 taste 结构（不同夸克味道的组合）重构。
- **多移求解器集成**：KS 费米子力计算通过 `CFieldMatrixOperationKST` 与多移求解器集成，用于 RHMC。
