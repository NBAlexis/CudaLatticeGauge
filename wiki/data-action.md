# Action 模块

作用量定义：规范作用量、费米子作用量、标量作用量，以及 HMC 能量/力计算接口。

## 文件清单

### 基类

| 文件 | 路径 | 说明 |
|------|------|------|
| `CAction.h` | `Data/Action/CAction.h` | 抽象基类 |

### 规范作用量

| 文件 | 路径 | 说明 |
|------|------|------|
| `CActionGaugePlaquette.h` | `Data/Action/` | 标准 Wilson plaquette |
| `CActionGaugePlaquetteAcceleration.h` | `Data/Action/` | 匀加速参考系 |
| `CActionGaugePlaquetteAtGradient.h` | `Data/Action/` | 空间变化各向异性 Xi(x) |
| `CActionGaugePlaquetteBetaGradient.h` | `Data/Action/` | 空间变化耦合 Beta(x) |
| `CActionGaugePlaquetteBoost.h` | `Data/Action/` | 常数洛伦兹 boost |
| `CActionGaugePlaquettePolar.h` | `Data/Action/` | 极坐标/柱坐标 |
| `CActionGaugePlaquetteRigidAcc.h` | `Data/Action/` | 刚性加速（Rindler） |
| `CActionGaugePlaquetteRotatingT.h` | `Data/Action/` | 旋转参考系（模板） |
| `CActionGaugePlaquetteRotatingT3D.h` | `Data/Action/` | 旋转参考系 3D（模板） |

### 费米子作用量

| 文件 | 路径 | 说明 |
|------|------|------|
| `CActionFermionKS.h` | `Data/Action/` | 单费米子场（Wilson/KS/RHMC） |
| `CActionFermionKSCombined.h` | `Data/Action/` | 多费米子场组合 |
| `CActionFermionKSImprove.h` | `Data/Action/` | 改进 KS（含 HISQ） |
| `CActionFermionKSImproveCombined.h` | `Data/Action/` | 改进 KS 组合 |
| `CActionFermionHISQCombined.h` | `Data/Action/` | HISQ 组合 |

### 标量/玻色子作用量

| 文件 | 路径 | 说明 |
|------|------|------|
| `CActionPhi4.h` | `Data/Action/` | phi^4 标量场 |
| `CActionTemperatureDistribution.h` | `Data/Action/` | 温度梯度（实玻色子场） |
| `CActionTemperatureDistribution3D.h/cu` | `Data/Action/` | 温度梯度 3D（实现同 `CActionTemperatureDistribution`） |

## CAction（抽象基类）

**文件**: `Data/Action/CAction.h`

所有作用量的统一接口：

| 方法 | 说明 |
|------|------|
| `Energy(bBeforeEvolution, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, stapleFields)` | 计算能量；gauge action 可使用预计算 staple；携带 gauge/boson/tensor2 工作场（const） |
| `CalculateForce(...)` | 计算 HMC 力（gauge + boson） |
| `PrepareForHMC(...)` | 轨迹前准备：随机 momentum、预计算 staple |
| `OnFinishTrajectory(bAccepted)` | 接受后缓存能量（非费米子作用量） |
| `OnFinishTrajectory(bWillBeAccept, gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields)` | 在 `CIntegrator::OnFinishTrajectory` 中、被接受场拷贝回 lattice 之前调用，给 action 一个修改即将被接受的场的机会；场可修改，数组本身 const；默认空实现 |
| `Initial(pOwner, params)` | 绑定到格点并解析参数 |
| `IsFermion()` | 是否为费米子作用量（影响能量缓存策略） |
| `GetBetaOverN() / SetBeta(...)` | 耦合参数访问 |

## 规范作用量

### CActionGaugePlaquette — 标准 Wilson plaquette

**文件**: `Data/Action/CActionGaugePlaquette.h`

最基本的规范作用量（实现中参数为 `betaOverN = beta / N`）：
```
S = betaOverN * sum_p (1 - Re Tr[U_p])
```

支持 `CloverEnergy` 标志启用 clover 改进能量计算。力计算通过 staple 实现。

### 变形规范作用量

所有变形作用量继承 `CAction`，通过重写 `EnergySingleField`、`CalculateForceOnGaugeSingleField` 和 `PrepareForHMCSingleField` 实现。

| 类 | 物理场景 | 关键参数 |
|----|---------|---------|
| `CActionGaugePlaquetteAcceleration` | 匀加速 | `G` — 加速度 |
| `CActionGaugePlaquetteAtGradient` | 温度梯度（各向异性） | `Xi` 空间分布列表 |
| `CActionGaugePlaquetteBetaGradient` | 空间变化耦合 | `Beta` 空间分布列表 |
| `CActionGaugePlaquetteCylinder` | 柱坐标（0,1,2,3=r,φ,z,t） | `Beta(r)` 逐 r 层数组、`RStart`、`REnd`、`DeltaR`、`CloverEnergy`；含 φ 平面权重 1/r，其余 r；plaquette 耦合取 4 角平均，推导见 `Docs/Applications/Cylinder.tex` |
| `CActionGaugePlaquetteBoost` | 常数 boost | `G` — boost 强度 |
| `CActionGaugePlaquetteRigidAcc` | 刚性加速 | `G`, `Beta`；空间权重 `1+gz`，时间权重 `1/(1+gz)` |

### 旋转参考系

**文件**: `CActionGaugePlaquetteRotatingT.h`（模板）

模板类 `CActionGaugePlaquetteRotatingT<deviceGauge, matrixN>` 支持多规范群：
- `CActionGaugePlaquetteRotatingU1` — U(1)
- `CActionGaugePlaquetteRotating` — SU(3)
- `CActionGaugePlaquetteRotatingSU2` — SU(2)
- `CActionGaugePlaquetteRotatingSU4` — SU(4)

**关键参数**：
- `Omega` — 角速度
- `Xi` — 各向异性参数
- `CloverEnergy` — 是否使用 clover
- `ShiftHalfCoord` — 半格点偏移
- `Torus` — 环面边界
- `RotationBetaScale` — 旋转项 beta 缩放模式（见下）
- `ScaleFactor` — 仅当 `RotationBetaScale = ERBSF_Custom` 时读取的缩放因子（默认 `1.0`）

支持三种边界条件：Dirichlet、射影平面、环面。每种边界有独立的 energy/force kernel。

**旋转项 beta 缩放（`RotationBetaScale`）**：

旋转作用量分两部分：S0（普通 plaquette/rectangle 项，委托给规范场对象计算）与 S1/S2（依赖 `Omega` 的旋转项，由本 action 的 `_kernel...` 计算）。S1/S2 的推导假设 Wilson 归一化 `beta = 6/g^2`，而改进规范场（`CFieldGaugeSU3TreeImproved`、`CFieldGaugeSU3OneLoopImproved`）将 plaquette 系数归一为 1、采用 `beta = 10/g^2`。因此当旋转 action 接入改进规范场时，直接把 YAML 的 `Beta` 传给 S1/S2 会使旋转耦合被放大 `10/6 = 5/3` 倍。

`RotationBetaScale` 选择一个仅作用于 S1/S2 的缩放因子（S0 仍用原始 `beta`，由规范场按自身约定解释）：

| 取值 | 缩放因子 | 用途 |
| --- | --- | --- |
| `ERBSF_Naive` (默认) | `1.0` | 普通 Wilson 规范场（`beta = 6/g^2`） |
| `ERBSF_TreeImprove` | `0.6`（= 6/10） | tree level improved 规范场（`beta = 10/g^2`） |
| `ERBSF_OneLoopImprove` | `0.6`（= 6/10） | one-loop improved 规范场（同样 `beta = 10/g^2`，plaquette 系数归一为 1） |
| `ERBSF_Custom` | 读取 `ScaleFactor`（默认 `1.0`） | 手动指定 |

缩放同时作用于 S1/S2 的 energy 与 force，故 HMC 的 h_diff 仍然守恒；仅平衡态可观测量（如 plaquette energy）随旋转耦合改变而移动。缩放通过 `GetRotationBetaOverN()`（`= m_fBetaOverN * m_fRotationScaleFactor`）实现。

**注意**：旧版 `CActionGaugePlaquetteRotating.h` 和 `CActionGaugePlaquetteRotating3D.h` 已废弃（被 `#if 0` 包裹）。

## 费米子作用量

### CActionFermionKS — 单费米子场

**文件**: `Data/Action/CActionFermionKS.h`

通用费米子作用量包装器。支持 Wilson-Dirac 和 Staggered/Kogut-Susskind 费米子，通过 RHMC（有理混合蒙特卡洛）实现。

**关键方法**：
- `Energy(...)` — 伪费米子双线性 `phi^+ (D^+ D)^{-1} phi`
- `CalculateForce(...)` — 费米子力（需要调用求解器）
- `PrepareForHMC(...)` — 生成高斯分布的伪费米子场

### CActionFermionKSCombined — 多费米子场

**文件**: `Data/Action/CActionFermionKSCombined.h`

支持多个独立的 even-odd staggered 费米子场（如 Nf=2+1）。遍历所有费米子场计算总能量和力。

### CActionFermionKSImprove — 改进 KS

**文件**: `Data/Action/CActionFermionKSImprove.h`

改进的 staggered 费米子作用量。力的计算机制不同：先计算 `X^+(n+mu) Y(n)` 和 `X^+(n) Y(n+mu)`。

支持 HISQ（Highly Improved Staggered Quark）含 Naik 项。

### CActionFermionHISQCombined — HISQ 组合

**文件**: `Data/Action/CActionFermionHISQCombined.h`

显式支持 even-odd 费米子和带相位的 EO 费米子。使用 `CFieldFermionKSHISQ` 场类型。

## 标量/玻色子作用量

### CActionPhi4 — phi^4 理论

**文件**: `Data/Action/CActionPhi4.h`

标量场作用量：
```
S = sum [ (nabla phi)^2 + m^2 phi^2 + lambda (phi^2)^2 ]
```

**参数**：`M`（质量）、`Lambda`（四次耦合）。

### CActionTemperatureDistribution — 温度梯度

**文件**: `Data/Action/CActionTemperatureDistribution.h`

通过实玻色子场实现温度梯度效应（而非规范场各向异性）。规范场与玻色子场耦合。

**参数**：`BosonKinetic`（玻色子动能系数）、`XiPower`。

## 添加新 Action 的指南

1. **继承 `CAction`**（`Data/Action/CAction.h`）。
2. **重写三个核心方法**：
   - `EnergySingleField(...)` — 能量计算
   - `CalculateForceOnGaugeSingleField(...)` — 力计算
   - `PrepareForHMCSingleField(...)` — HMC 前准备（预计算 staple 等）
3. **注册到 RTTI 工厂**：
   - 头文件加 `__CLG_REGISTER_HELPER_HEADER(CMyAction)`
   - 实现文件加 `__CLGIMPLEMENT_CLASS(CMyAction)`
   - 无需手动修改 `CLGLibManager`
4. **设备端 kernel**：放在 `.cu` 文件中，使用 `_LAUNCH_KERNEL` 宏启动

对于 gauge action，通常需要：
- 设备端 staple/force kernel（模板化以支持多规范群）
- 设备端 energy kernel
- Host 端 `switch (GetFieldType())` 分发到具体模板实例

示例 YAML：
```yaml
Action1:
    ActionName: CActionGaugePlaquette
    Beta: 5.5
    GaugeFields: [1]
```

## 类层次速查

```
CAction
├── CActionGaugePlaquette
│   └── （各种变形通过独立类实现）
├── CActionGaugePlaquetteAcceleration
├── CActionGaugePlaquetteAtGradient
├── CActionGaugePlaquetteBetaGradient
├── CActionGaugePlaquetteBoost
├── CActionGaugePlaquettePolar
├── CActionGaugePlaquetteRigidAcc
├── CActionGaugePlaquetteRotatingT<G,N>
│   └── CActionGaugePlaquetteRotatingT3D<G,N>
├── CActionFermionKS
│   └── CActionFermionKSImprove
├── CActionFermionKSCombined
│   └── CActionFermionKSImproveCombined
├── CActionFermionHISQCombined
├── CActionPhi4
└── CActionTemperatureDistribution
```
