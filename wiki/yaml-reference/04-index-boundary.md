> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 04. 索引与边界条件

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateIndexAndBoundary()` (lines 695-725)

---

## 配置参数

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `LatticeIndex` | string | `CIndexSquare` | 否 | `CLGLibManager.cpp` ~L700 | 索引类名。通过 `appCreate()` 实例化 |
| `LatticeBoundary` | string | `CBoundaryConditionTorusSquare` | 否 | `CLGLibManager.cpp` ~L702 | 边界条件类名。通过 `appCreate()` 实例化 |

## 支持的索引类型

| 类名 | 说明 |
|------|------|
| `CIndexSquare` | 标准四维方格索引 |
| `CIndexSquareAMC` | 带 AMC 改进的方格索引 |
| `CIndexSlab` | Slab 几何索引 |
| `CIndexMultiple` | 多重索引（用于多格点） |

## 支持的边界条件

| 类名 | 说明 |
|------|------|
| `CBoundaryConditionTorusSquare` | 周期边界条件（torus） |
| `CBoundaryConditionSphere` | 球面边界条件 |
| `CBoundaryConditionProjectivePlaneSquare` | 射影平面边界（XY 平面为 RP²，其余方向周期/反周期由 `Period` 控制） |

## 初始化行为

1. 调用 `appCreate(LatticeIndex)` 创建索引对象
2. 调用 `appCreate(LatticeBoundary)` 创建边界条件对象
3. 调用 `CIndex::InitialIndex()` 初始化索引表
4. 调用 `CBoundaryCondition::Initial()` 初始化边界条件

## 射影平面边界条件示例

```yaml
LatticeBoundary: CBoundaryConditionProjectivePlaneSquare
```

- XY 平面使用射影平面（RP²）识别：跨越 X 或 Y 边界时，另一方向坐标取反，同时根据 `Period` 决定是否需要加负号/取厄米共轭。
- Z、T 方向由 `Period` 参数控制（周期 `1` / 反周期 `-1`）。
- 代码实现里，FieldId `0` 与 `1` 默认 T 方向为周期，其它 FieldId 默认 T 方向为反周期（可通过各场子块里的 `Period` 覆盖）。

---

[< 返回目录](home.md) | [< 上一章：随机数](03-random.md) | [下一章：规范场 >](05-gauge-fields.md)
