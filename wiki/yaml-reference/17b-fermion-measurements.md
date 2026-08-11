> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 17b. 费米子测量器

> **依赖提示**：本章所有费米子测量器都需要为被测费米子场配置 `Solver`（`SolverForFieldId` 匹配 `FieldId`）。测量过程会调用 `InverseD` / `InverseDDdagger` 进行反演。详见 [10-solvers.md](10-solvers.md)。

### CMeasureChiralCondensate — 手征凝聚（Wilson 费米子）

**来源文件**：`Code/CLGLib/Measurement/CMeasureChiralCondensate.cu`

**继承**：`CMeasureStochastic`

**读取函数**：`CMeasureChiralCondensate::InitialOtherParameters()` (lines 215-221)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MeasureDist` | INT | `1` | 否 | `CMeasureChiralCondensate.cu` ~L221 | 测量 R-分布。`1`=测量 |

**继承参数**（来自 `CMeasureStochastic`）：

| 键名 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `FieldCount` | INT | `25` | Z4 随机源个数 |
| `DebugDivation` | INT | `0` | 偏差调试 |

---

### CMeasureChiralCondensateKS — 手征凝聚（KS 费米子）

**来源文件**：`Code/CLGLib/Measurement/CMeasureChiralCondensateKS.cu`

**继承**：`CMeasureStochastic`

**读取函数**：`CMeasureChiralCondensateKS::InitialOtherParameters()` (lines 207-234)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ShiftCenter` | INT | `0` | 否 | `CMeasureChiralCondensateKS.cu` ~L213 | 分布测量时偏移中心。`1`=偏移 |
| `MeasureSigma12` | INT | `0` | 否 | `CMeasureChiralCondensateKS.cu` ~L217 | 测量 sigma_12 凝聚。`1`=测量 |
| `MeasureConnect` | INT | `0` | 否 | `CMeasureChiralCondensateKS.cu` ~L221 | 测量连通 susceptibility。`1`=测量 |
| `XSlice` | INT | `0` | 否 | `CMeasureChiralCondensateKS.cu` ~L225 | 测量 X 切片分布。`1`=测量 |
| `YSlice` | INT | `0` | 否 | `CMeasureChiralCondensateKS.cu` ~L228 | 测量 Y 切片分布。`1`=测量 |
| `ZSlice` | INT | `0` | 否 | `CMeasureChiralCondensateKS.cu` ~L231 | 测量 Z 切片分布。`1`=测量 |
| `TSlice` | INT | `0` | 否 | `CMeasureChiralCondensateKS.cu` ~L234 | 测量 T 切片分布。`1`=测量 |

**继承参数**（来自 `CMeasureStochastic`）：同上。

---

### CMeasureConnectedSusceptibilityKS — 连通手征 Susceptibility（KS）

**来源文件**：`Code/CLGLib/Measurement/CMeasureConnectedChiralSusceptibilityKS.cu`

**继承**：`CMeasureStochastic`（无额外参数）

---

### CMeasureMesonCorrelator — 介子关联函数（Wilson 费米子）

**来源文件**：`Code/CLGLib/Measurement/CMeasureMesonCorrelator.cu`

**继承**：`CMeasureStochastic`

**读取函数**：`CMeasureMesonCorrelator::InitialOtherParameters()` (lines 172-178)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `GammaMatrix` | String[] | `["UNITY"]` | 否 | `CMeasureMesonCorrelator.cu` ~L178 | Gamma 矩阵列表，如 `GAMMA1`, `GAMMA2`, `GAMMA3`, `GAMMA4`, `GAMMA5`, `SIGMA12` 等 |

**继承参数**（来自 `CMeasureStochastic`）：同上。

---

### CMeasureMesonCorrelatorStaggered — 介子关联函数（Staggered）

**来源文件**：`Code/CLGLib/Measurement/CMeasureMesonCorrelatorStaggered.cu`

**继承**：`CMeasureStochastic`

**读取函数**：`CMeasureMesonCorrelatorStaggered::InitialOtherParameters()` (lines 570-578)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `GaugeFixing` | INT | `1` | 否 | `CMeasureMesonCorrelatorStaggered.cu` ~L574 | 测量前做规范固定。`1`=启用 |
| `SimpleVersion` | INT | `1` | 否 | `CMeasureMesonCorrelatorStaggered.cu` ~L578 | 使用简化版本算法。`1`=启用 |

**继承参数**（来自 `CMeasureStochastic`）：同上。

---

### CMeasureMesonCorrelatorStaggeredSimple2 — 简化介子关联函数（带电介子）

**来源文件**：`Code/CLGLib/Measurement/CMeasureMesonCorrelatorStaggeredSimple2.cu`

**继承**：`CMeasureStochastic`

**读取函数**：`CMeasureMesonCorrelatorStaggeredSimple2::InitialOtherParameters()` (lines 128-152)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldId2` | INT | `0` | 否 | `CMeasureMesonCorrelatorStaggeredSimple2.cu` ~L132 | 第二个费米子场 ID（用于带电介子，如 `ud`, `du`） |
| `Source` | String | `"EFS_Point"` | 否 | `CMeasureMesonCorrelatorStaggeredSimple2.cu` ~L152 | 源类型：`EFS_Point`, `EFS_Wall` 等 |

**继承参数**（来自 `CMeasureStochastic`）：同上。

---


---

[< 返回 17. 测量器目录](17-measurements.md) | [< 返回 yaml-reference 目录](home.md)
