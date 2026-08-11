> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 17c. 其他测量器

### CMeasureAMomentumJG — 角动量（Jaffe-Goldman 分解）

**来源文件**：`Code/CLGLib/Measurement/CMeasureAMomentumJG.cu`

**读取函数**：`CMeasureAMomentumJG::InitialOtherParameters()` (lines 740-769)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ShowResult` | INT | `1` | 否 | `CMeasureAMomentumJG.cu` ~L746 | 输出结果到日志。`1`=输出 |
| `MeasureDist` | INT | `1` | 否 | `CMeasureAMomentumJG.cu` ~L750 | 测量 R-分布。`1`=测量 |
| `MeasureSpin` | INT | `0` | 否 | `CMeasureAMomentumJG.cu` ~L754 | 测量自旋分解（需要 E 和 DpureA 场）。`1`=测量 |
| `MeasureApprox` | INT | `0` | 否 | `CMeasureAMomentumJG.cu` ~L758 | 测量近似版本（JGChenApprox, JGChenApprox2）。`1`=测量 |
| `ProjectivePlane` | INT | `0` | 否 | `CMeasureAMomentumJG.cu` ~L762 | 使用射影平面几何。`1`=启用 |
| `NaiveNabla` | INT | `0` | 否 | `CMeasureAMomentumJG.cu` ~L766 | 使用 naive nabla 算子。`1`=启用 |
| `BetaOverN` | DOUBLE | `0.0` | 否 | `CMeasureAMomentumJG.cu` ~L769 | Beta/N 耦合常数 |

---

### CMeasureAMomentumStochastic — 角动量（随机源方法）

**来源文件**：`Code/CLGLib/Measurement/CMeasureAMomentumStochastic.cu`

**继承**：`CMeasureStochastic`

**读取函数**：`CMeasureAMomentumStochastic::InitialOtherParameters()` (lines 865-889)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MeasurePure` | INT | `0` | 否 | `CMeasureAMomentumStochastic.cu` ~L870 | 测量 JLPure 和 JLJM（Chen 和 Jaffe-Manohar 分解）。`1`=测量 |
| `Exponential` | INT | `1` | 否 | `CMeasureAMomentumStochastic.cu` ~L885 | 使用指数离散化。`1`=启用 |
| `Naive` | INT | `1` | 否 | `CMeasureAMomentumStochastic.cu` ~L889 | 使用 naive 导数。`1`=启用 |

**继承参数**（来自 `CMeasureStochastic`）：同上。

---

### CMeasureAngularMomentumKS — 角动量（KS 费米子）

**来源文件**：`Code/CLGLib/Measurement/CMeasureAngularMomentumKS.cu`

**读取函数**：`CMeasureAngularMomentumKS::InitialOtherParameters()` (lines 328-338)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ShiftCenter` | INT | `0` | 否 | `CMeasureAngularMomentumKS.cu` ~L334 | 分布测量时偏移中心。`1`=偏移 |
| `ZSlice` | INT | `0` | 否 | `CMeasureAngularMomentumKS.cu` ~L338 | 测量 Z 切片分布。`1`=测量 |

---

### CMeasureAngularMomentumKSREM — 角动量 KS（REM 版本）

**来源文件**：`Code/CLGLib/Measurement/CMeasureAngularMomentumKSREM.cu`

**继承**：`CMeasureAngularMomentumKS`（无额外参数，Initial 方法被注释掉）

---

### CMeasureBerryPhase — Berry Phase

**来源文件**：`Code/CLGLib/Measurement/CMeasureBerryPhase.cu`

**读取函数**：`CMeasureBerryPhase::InitialOtherParameters()` (lines 418-427)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `WilsonDirac` | INT | `1` | 否 | `CMeasureBerryPhase.cu` ~L423 | 使用 Wilson-Dirac 费米子（vs KS 费米子）。`1`=Wilson |
| `DoGaugeFixing` | INT | `0` | 否 | `CMeasureBerryPhase.cu` ~L427 | 计算前做规范固定。`1`=启用 |

---

### CMeasureBosonCond — 玻色子凝聚

**来源文件**：`Code/CLGLib/Measurement/CMeasureBosonCond.cpp`

**继承**：`CMeasure`（无额外参数，使用 `BosonFields` 指定测量的玻色子场）

---

### CMeasureBosonValue — 玻色子场值

**来源文件**：`Code/CLGLib/Measurement/CMeasureBosonValue.cu`

**模板实例化**：Real, U1, SU2, SU3, SU4-SU8

**继承**：`CMeasure`（无额外参数，使用 `BosonFields` 指定测量的玻色子场）

---

### CMeasureChargeAndCurrents — 电荷与流

**来源文件**：`Code/CLGLib/Measurement/CMeasureChargeAndCurrents.cu`

**继承**：`CMeasure`（无额外参数）

---

### CMeasurePandChiralTalor — P 和手征 Taylor 展开

**来源文件**：`Code/CLGLib/Measurement/CMeasurePandChiralTalor.cu`

**继承**：`CMeasureStochastic`（无额外参数）

---

### CMeasurePandChiralTalorKS — P 和手征 Taylor 展开（KS）

**来源文件**：`Code/CLGLib/Measurement/CMeasurePandChiralTalorKS.cu`

**继承**：`CMeasureStochastic`（无额外参数）

---

### CMeasureTopologicChargeXY — 拓扑荷（XY 平面）

**来源文件**：`Code/CLGLib/Measurement/CMeasureTopologicChargeXY.cu`

**继承**：`CMeasure`（无额外参数）

---

### CMeasureGaugeRotation — 规范转动

**来源文件**：`Code/CLGLib/Measurement/CMeasureGaugeRotation.cu`

**继承**：`CMeasure`（无额外参数）

---

### CMeasureWilsonNality — Wilson Nality（已禁用）

**来源文件**：`Code/CLGLib/Measurement/CMeasureWilsonNality.cu`

**状态**：整个文件被 `#if 0` 包裹，类已废弃/禁用。

---


---

[< 返回 17. 测量器目录](17-measurements.md) | [< 返回 yaml-reference 目录](home.md)
