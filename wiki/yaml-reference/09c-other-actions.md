> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 09c. 标量、温度与离散作用量

## 标量场作用量

### CActionPhi4

**来源文件**：`Code/CLGLib/Data/Action/CActionPhi4.cpp`

**读取函数**：`CActionPhi4::InitialOtherParameters()` (lines 18-24)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Mass` | Real | `1.0` | **是** | `CActionPhi4.cpp` ~L22 | 标量场质量 |
| `Lambda` | Real | `1.0` | **是** | `CActionPhi4.cpp` ~L23 | 自耦合常数 |

---

## 温度分布作用量

### CActionTemperatureDistribution

**来源文件**：`Code/CLGLib/Data/Action/CActionTemperatureDistribution.cu`

**读取函数**：`CActionTemperatureDistribution::InitialOtherParameters()` (lines 457-480)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `StaticBosonFieldId` | INT | `0` | 否 | `CActionTemperatureDistribution.cu` ~L462 | 静态玻色子场 ID。`> 0` 时查找对应场 |
| `BosonKinetic` | DOUBLE | `1.0` | 否 | `CActionTemperatureDistribution.cu` ~L471 | 玻色子动能系数 |
| `XiPower` | INT | `-1` | 否 | `CActionTemperatureDistribution.cu` ~L476 | Xi 幂次 |

### CActionTemperatureDistribution3D

**来源文件**：`Code/CLGLib/Data/Action/CActionTemperatureDistribution3D.cu`

**读取函数**：`CActionTemperatureDistribution3D::InitialOtherParameters()` (lines 305-323)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `StaticBosonFieldId` | INT | `0` | 否 | `CActionTemperatureDistribution3D.cu` ~L310 | 静态玻色子场 ID。`> 0` 时查找对应场 |
| `BosonKinetic` | DOUBLE | `1.0` | 否 | `CActionTemperatureDistribution3D.cu` ~L319 | 玻色子动能系数 |

**注意**：3D 变体**没有** `XiPower` 参数。

---

## 离散规范群作用量

### CActionDiscreteZNPlaquette

**来源文件**：`Code/CLGLib/Data/Action/CActionDiscreteZNPlaquette.h` / `CActionDiscreteGauge.cpp`

继承自 `CActionDiscreteGauge::Initial()`：
- `Beta` (DOUBLE) — 默认值 `1.0`

**说明**：离散群作用量无导数，因此**不支持 HMC** 连续更新。只能与热浴或 Metropolis 更新器一起使用。

### CActionDiscreteDNPlaquette

**来源文件**：`Code/CLGLib/Data/Action/CActionDiscreteDNPlaquette.h` / `CActionDiscreteGauge.cpp`

继承自 `CActionDiscreteGauge::Initial()`：
- `Beta` (DOUBLE) — 默认值 `1.0`

**说明**：同 ZN plaquette，用于 D_N 群。

---


---

[< 返回 09. 作用量](09-actions.md) | [< 返回 yaml-reference 目录](home.md)
