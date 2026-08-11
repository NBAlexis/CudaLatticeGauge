> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 09. 作用量 (Actions)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateActionList()` (lines 727-761)

---

## 配置参数

YAML 键名为 `Action1`, `Action2`, ...。通过 `ActionListLength` 控制数量。

### 基类参数 (CAction)

**来源文件**：`Code/CLGLib/Data/Action/CAction.cpp`

**读取函数**：`CAction::Initial()` (lines 19-37)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ActionName` | string | （无） | **是** | `CLGLibManager.cpp` ~L742 | 作用量类名。通过 `appCreate()` 实例化 |
| `GaugeFields` | BYTE[] | `[1]` | 否 | `CAction.cpp` ~L24 | 此作用量作用的规范场 ID 列表 |
| `BosonFields` | BYTE[] | `[]` | 否 | `CAction.cpp` ~L26 | 此作用量作用的玻色子场 ID 列表 |
| `Beta` | DOUBLE | `0.1` | 否 | `CAction.cpp` ~L33 | 耦合常数。代码中存储为 `m_fBetaOverN = Beta / Nc` |

**注意**：如果 `GaugeFields` 和 `BosonFields` 都未指定或为空，默认设置为 `GaugeFields: [1]`。

### 离散规范群作用量基类 (CActionDiscreteGauge)

**来源文件**：`Code/CLGLib/Data/Action/CActionDiscreteGauge.cpp`

**读取函数**：`CActionDiscreteGauge::Initial()` (lines 27-34)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Beta` | DOUBLE | `1.0` | 否 | `CActionDiscreteGauge.cpp` ~L31 | 离散群耦合常数。覆盖 CAction 的默认值 |

---

## 支持的作用量类

以下所有类均通过 `appCreate(ActionName)` 工厂实例化。需确认 `Code/CLGLib/Core/CLGSetup.h` 中对应宏已启用。

| 类名 | 类型 | 说明 |
|------|------|------|
| `CActionGaugePlaquette` | 规范 | 标准 Plaquette 作用量 |
| `CActionGaugePlaquetteBoost` | 规范 | Boosted Plaquette 作用量 |
| `CActionGaugePlaquetteAcceleration` | 规范 | 加速 Plaquette 作用量 |
| `CActionGaugePlaquetteRigidAcc` | 规范 | 刚性加速 Plaquette 作用量 |
| `CActionGaugePlaquetteAnisotropic` | 规范 | 各向异性 tree-level Symanzik 作用量（裸 `Xi`） |
| `CActionGaugePlaquetteCylinder` | 规范 | 柱坐标 Plaquette 作用量（`Beta(r)`、`RStart`、`DeltaR`） |
| `CActionGaugePlaquetteAtGradient` | 规范 | Xi 梯度 Plaquette 作用量 |
| `CActionGaugePlaquetteGradient` | 规范 | Beta 梯度 Plaquette 作用量 |
| `CActionGaugePlaquetteRotating` | 规范 | 旋转 Plaquette 作用量 (SU3) |
| `CActionGaugePlaquetteRotatingU1` | 规范 | 旋转 Plaquette 作用量 (U1) |
| `CActionGaugePlaquetteRotatingSU2` | 规范 | 旋转 Plaquette 作用量 (SU2) |
| `CActionGaugePlaquetteRotatingSU4` | 规范 | 旋转 Plaquette 作用量 (SU4) |
| `CActionGaugePlaquetteRotating3D` | 规范 | 3D 旋转 Plaquette 作用量 (SU3) |
| `CActionGaugePlaquetteRotatingU1_3D` | 规范 | 3D 旋转 Plaquette 作用量 (U1) |
| `CActionFermionKS` | 费米子 | KS 费米子作用量 |
| `CActionFermionKSCombined` | 费米子 | 组合 KS 费米子作用量 |
| `CActionFermionKSImprove` | 费米子 | 改进 KS 费米子作用量 |
| `CActionFermionKSImproveCombined` | 费米子 | 改进组合 KS 费米子作用量 |
| `CActionFermionHISQCombined` | 费米子 | HISQ 组合费米子作用量 |
| `CActionPhi4` | 标量 | Phi^4 标量场作用量 |
| `CActionTemperatureDistribution` | 温度 | 温度分布作用量 |
| `CActionDiscreteZNPlaquette` | 离散 | Z_N Plaquette 作用量 |
| `CActionDiscreteDNPlaquette` | 离散 | D_N Plaquette 作用量 |

---

## 常见依赖关系

| 若配置此作用量 | 必须同时配置 | 说明 |
|---|---|---|
| `CActionFermionHISQCombined` | `GaugeSmearing: CGaugeSmearingHISQSU3` | HISQ 作用量依赖 HISQ smearing 的 effective gauge。见 [12-gauge-smearing.md](12-gauge-smearing.md) |
| `CActionFermionKSImprove` / `CActionFermionKSImproveCombined` | 任意 `GaugeSmearing`（Stout / ASQTAD / HISQ） | 改进作用量通过 smearing 取 effective gauge。见 [12-gauge-smearing.md](12-gauge-smearing.md) |
| 任意费米子作用量用于 HMC | `MSSolver` 或 `Solver`（每个相关场） | HISQ / improved KS / rational 用 `MSSolver`；`ER_NoRational` 用 `Solver`。见 [10-solvers.md](10-solvers.md)、[11-multi-shift-solvers.md](11-multi-shift-solvers.md) |
| 任意费米子作用量用于测量 | `Solver`（每个被测场） | 测量反演。见 [10-solvers.md](10-solvers.md) |

---


## 子页面

- [09a. 规范作用量](09a-gauge-actions.md) — Plaquette、Boost、加速、极坐标、各向异性、旋转等
- [09b. 费米子作用量](09b-fermion-actions.md) — KS、HISQ 费米子作用量
- [09c. 标量、温度与离散作用量](09c-other-actions.md) — Phi4、温度分布、Z_N / D_N 离散群

## 参数继承关系速查

### 规范作用量参数继承

```
CAction
  ├── GaugeFields, BosonFields, Beta
  ├── CActionGaugePlaquette
  │     └── CloverEnergy
  ├── CActionGaugePlaquetteBoost
  │     └── Boost
  ├── CActionGaugePlaquetteAcceleration
  │     └── AccG
  ├── CActionGaugePlaquetteRigidAcc
  │     ├── AccG, Dirichlet
  ├── CActionGaugePlaquettePolar
  │     ├── BetaList, RIn, ROut, Dirichlet
  ├── CActionGaugePlaquetteAtGradient
  │     └── Xi[]
  ├── CActionGaugePlaquetteGradient
  │     └── Beta[]
  └── CActionGaugePlaquetteRotatingT (模板)
        ├── Omega, Xi, CloverEnergy, ShiftCoord, Torus
        └── 实例化为 RotatingU1, Rotating(SU3), RotatingSU2, RotatingSU4,
             RotatingU1_3D, Rotating3D(SU3)
```

### 费米子作用量参数继承

```
CAction
  ├── GaugeFields, BosonFields, Beta
  ├── CActionFermionKS
  │     └── FieldId
  ├── CActionFermionKSCombined
  │     └── FieldIds (必需)
  ├── CActionFermionKSImprove
  │     ├── FieldId, HISQ
  ├── CActionFermionKSImproveCombined
  │     └── FieldIds (必需)
  └── CActionFermionHISQCombined
        └── FieldIds (必需)
```

### 其他作用量

```
CAction
  ├── GaugeFields, BosonFields, Beta
  ├── CActionPhi4
  │     └── Mass, Lambda
  ├── CActionTemperatureDistribution
  │     └── StaticBosonFieldId, BosonKinetic, XiPower
  ├── CActionTemperatureDistribution3D
  │     └── StaticBosonFieldId, BosonKinetic
  └── CActionDiscreteGauge (抽象基类)
        └── Beta (默认 1.0，覆盖 CAction 的 0.1)
        └── CActionDiscreteZNPlaquette (无额外参数)
        └── CActionDiscreteDNPlaquette (无额外参数)
```

---

[< 返回目录](home.md) | [< 上一章：费米子场](08-fermion-fields.md) | [下一章：求解器 >](10-solvers.md)
