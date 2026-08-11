> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 17. 测量器 (Measurements)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateMeasurement()` (lines 802-838)

---

## 配置参数

YAML 键名为 `Measure1`, `Measure2`, ...。通过 `MeasureListLength` 控制数量。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MeasureName` | string | （无） | **是** | `CLGLibManager.cpp` ~L817 | 测量器类名。通过 `appCreate()` 实例化 |

---

## 基类参数 (CMeasure)

**来源文件**：`Code/CLGLib/Measurement/CMeasure.cu`

**读取函数**：`CMeasure::Initial()` (lines 177-202)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `GaugeSmearing` | INT | `0` | 否 | `CMeasure.cu` ~L184 | 测量前是否做规范平滑。`1`=启用 |
| `FieldId` | INT | `0` | 否 | `CMeasure.cu` ~L188 | 费米子场 ID（用于需随机源的测量） |
| `ShowResult` | INT | `1` | 否 | `CMeasure.cu` ~L192 | 是否输出结果到日志。`1`=输出，`0`=静默 |
| `GaugeFields` | BYTE[] | `[1]` | 否 | `CMeasure.cu` ~L195 | 测量的规范场 ID 列表 |
| `BosonFields` | BYTE[] | `[]` | 否 | `CMeasure.cu` ~L196 | 测量的玻色子场 ID 列表 |

**注意**：如果 `GaugeFields` 和 `BosonFields` 都未指定或为空，默认设置为 `GaugeFields: [1]`。

---

## 随机源测量基类 (CMeasureStochastic)

**来源文件**：`Code/CLGLib/Measurement/CMeasure.h`

**读取函数**：`CMeasureStochastic::Initial()` (lines 570-618)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldCount` | INT | `25` | 否 | `CMeasure.h` ~L584 | Z4/Gaussian 随机源个数。越多统计误差越小 |
| `DebugDivation` | INT | `0` | 否 | `CMeasure.h` ~L588 | 是否启用偏差调试输出。`1`=启用 |

---

## 具体测量器类

各测量器详情见子页面：

- [17a. 规范场测量器](17a-gauge-measurements.md) — 作用量、Plaquette、Polyakov、Wilson loop
- [17b. 费米子测量器](17b-fermion-measurements.md) — 手征凝聚、介子关联、随机源角动量
- [17c. 其他测量器](17c-other-measurements.md) — 角动量 JG、Berry phase、玻色子、拓扑荷等

## 支持的测量器类汇总

| 类名 | 来源文件 | 基类 | 说明 |
|------|----------|------|------|
| `CMeasureAction` | `CMeasureAction.cpp` | `CMeasure` | 作用量期望值 |
| `CMeasurePlaqutteEnergy` | `CMeasurePlaqutteEnergy.cpp` | `CMeasure` | Plaquette 能量密度 |
| `CMeasurePolyakovXY` | `CMeasurePolyakovXY.cu` | `CMeasure` | Polyakov loop（全方向+切片） |
| `CMeasurePolyakovXY3D` | `CMeasurePolyakovXY3D.cu` | `CMeasure` | 3D Polyakov loop |
| `CMeasureWilsonLoop` | `CMeasureWilsonLoop.cu` | `CMeasure` | Wilson loop |
| `CMeasureWilsonLoopXY` | `CMeasureWilsonLoopXY.cu` | `CMeasure` | Wilson loop（XY 平面） |
| `CMeasureWilsonLoopWithPath` | `CMeasureWilsonLoopWithPath.cu` | `CMeasure` | 自定义路径 Wilson loop |
| `CMeasureChiralCondensate` | `CMeasureChiralCondensate.cu` | `CMeasureStochastic` | 手征凝聚（Wilson） |
| `CMeasureChiralCondensateKS` | `CMeasureChiralCondensateKS.cu` | `CMeasureStochastic` | 手征凝聚（KS） |
| `CMeasureConnectedSusceptibilityKS` | `CMeasureConnectedChiralSusceptibilityKS.cu` | `CMeasureStochastic` | 连通手征 susceptibility |
| `CMeasureMesonCorrelator` | `CMeasureMesonCorrelator.cu` | `CMeasureStochastic` | 介子关联函数（Wilson） |
| `CMeasureMesonCorrelatorStaggered` | `CMeasureMesonCorrelatorStaggered.cu` | `CMeasureStochastic` | 介子关联函数（Staggered） |
| `CMeasureMesonCorrelatorStaggeredSimple2` | `CMeasureMesonCorrelatorStaggeredSimple2.cu` | `CMeasureStochastic` | 简化介子关联函数（带电） |
| `CMeasureAMomentumJG` | `CMeasureAMomentumJG.cu` | `CMeasure` | 角动量 JG 分解 |
| `CMeasureAMomentumStochastic` | `CMeasureAMomentumStochastic.cu` | `CMeasureStochastic` | 角动量（随机源） |
| `CMeasureAngularMomentumKS` | `CMeasureAngularMomentumKS.cu` | `CMeasure` | 角动量（KS） |
| `CMeasureAngularMomentumKSREM` | `CMeasureAngularMomentumKSREM.cu` | `CMeasureAngularMomentumKS` | 角动量 KS REM |
| `CMeasureBerryPhase` | `CMeasureBerryPhase.cu` | `CMeasure` | Berry Phase |
| `CMeasureBosonCond` | `CMeasureBosonCond.cpp` | `CMeasure` | 玻色子凝聚 |
| `CMeasureBosonValue` | `CMeasureBosonValue.cu` | `CMeasure` | 玻色子场值（多模板） |
| `CMeasureChargeAndCurrents` | `CMeasureChargeAndCurrents.cu` | `CMeasure` | 电荷与流 |
| `CMeasurePandChiralTalor` | `CMeasurePandChiralTalor.cu` | `CMeasureStochastic` | P 和手征 Taylor |
| `CMeasurePandChiralTalorKS` | `CMeasurePandChiralTalorKS.cu` | `CMeasureStochastic` | P 和手征 Taylor（KS） |
| `CMeasureTopologicChargeXY` | `CMeasureTopologicChargeXY.cu` | `CMeasure` | 拓扑荷（XY） |
| `CMeasureGaugeRotation` | `CMeasureGaugeRotation.cu` | `CMeasure` | 规范转动 |
| `CMeasureWilsonNality` | `CMeasureWilsonNality.cu` | — | **已禁用**（`#if 0`） |

---

## 参数速查表（按测量器）

| 测量器 | 特有参数 | 继承 CMeasureStochastic |
|--------|----------|------------------------|
| `CMeasureAction` | `FermiomFieldCount` | 否 |
| `CMeasurePlaqutteEnergy` | `v0` | 否 |
| `CMeasurePolyakovXY` | `MeasureX/Y/Z`, `X/Y/Z/TSlice`, `MeasureDist`, `Absolute`, `ShiftCenter` | 否 |
| `CMeasurePolyakovXY3D` | `U1`, `ShiftCenter` | 否 |
| `CMeasureWilsonLoop` | （无） | 否 |
| `CMeasureWilsonLoopXY` | （无） | 否 |
| `CMeasureWilsonLoopWithPath` | `Path`, `OnePoint`, `Point` | 否 |
| `CMeasureChiralCondensate` | `MeasureDist` | **是** |
| `CMeasureChiralCondensateKS` | `ShiftCenter`, `MeasureSigma12`, `MeasureConnect`, `X/Y/Z/TSlice` | **是** |
| `CMeasureConnectedSusceptibilityKS` | （无） | **是** |
| `CMeasureMesonCorrelator` | `GammaMatrix` | **是** |
| `CMeasureMesonCorrelatorStaggered` | `GaugeFixing`, `SimpleVersion` | **是** |
| `CMeasureMesonCorrelatorStaggeredSimple2` | `FieldId2`, `Source` | **是** |
| `CMeasureAMomentumJG` | `ShowResult`, `MeasureDist`, `MeasureSpin`, `MeasureApprox`, `ProjectivePlane`, `NaiveNabla`, `BetaOverN` | 否 |
| `CMeasureAMomentumStochastic` | `MeasurePure`, `Exponential`, `Naive` | **是** |
| `CMeasureAngularMomentumKS` | `ShiftCenter`, `ZSlice` | 否 |
| `CMeasureAngularMomentumKSREM` | （无，继承 KS） | 否 |
| `CMeasureBerryPhase` | `WilsonDirac`, `DoGaugeFixing` | 否 |
| `CMeasureBosonCond` | （无） | 否 |
| `CMeasureBosonValue` | （无） | 否 |
| `CMeasureChargeAndCurrents` | （无） | 否 |
| `CMeasurePandChiralTalor` | （无） | **是** |
| `CMeasurePandChiralTalorKS` | （无） | **是** |
| `CMeasureTopologicChargeXY` | （无） | 否 |
| `CMeasureGaugeRotation` | （无） | 否 |


---

## 补充参数（代码中存在但文档中缺失）

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ActionIndex` | INT | （无） | required | Code/CLGLib/Measurement/CMeasureRotatingAction.cpp ~L22 | 指定要测量的旋转作用量索引，对应 `ActionN` 中的 `N` |

---

[< 返回目录](home.md) | [< 上一章：积分器](16-integrators.md)
