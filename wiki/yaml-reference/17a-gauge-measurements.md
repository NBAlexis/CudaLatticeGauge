> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 17a. 规范场测量器

### CMeasureAction — 作用量期望值

**来源文件**：`Code/CLGLib/Measurement/CMeasureAction.cpp`

**读取函数**：`CMeasureAction::InitialOtherParameters()` (lines 17-23)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FermiomFieldCount` | INT | `1` | 否 | `CMeasureAction.cpp` ~L21 | 平均的费米子场个数。**注意**：键名拼写为 `Fermiom`（代码中的实际拼写） |

---

### CMeasurePlaqutteEnergy — Plaquette 能量密度

**来源文件**：`Code/CLGLib/Measurement/CMeasurePlaqutteEnergy.cpp`

**读取函数**：`CMeasurePlaqutteEnergy::InitialOtherParameters()` (lines 55-64)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `v0` | INT | `0` | 否 | `CMeasurePlaqutteEnergy.cpp` ~L60 | 是否同时测量 v0（平滑后的 plaquette）。`1`=测量 |

---

### CMeasurePolyakovXY — Polyakov loop（XY 平面及全方向）

**来源文件**：`Code/CLGLib/Measurement/CMeasurePolyakovXY.cu`

**读取函数**：`CMeasurePolyakovXY::InitialOtherParameters()` (lines 367-426)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `MeasureX` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L373 | 测量 X 方向 Polyakov loop。`1`=测量 |
| `MeasureY` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L376 | 测量 Y 方向 Polyakov loop。`1`=测量 |
| `MeasureZ` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L379 | 测量 Z 方向 Polyakov loop。`1`=测量 |
| `XSlice` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L396 | 测量 X 方向切片分布。`1`=测量 |
| `YSlice` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L399 | 测量 Y 方向切片分布。`1`=测量 |
| `ZSlice` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L402 | 测量 Z 方向切片分布。`1`=测量 |
| `TSlice` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L405 | 测量 T 方向切片分布。`1`=测量 |
| `MeasureDist` | INT | `1` | 否 | `CMeasurePolyakovXY.cu` ~L409 | 测量距离分布。`1`=测量 |
| `Absolute` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L412 | 测量绝对值。`1`=测量 |
| `ShiftCenter` | INT | `0` | 否 | `CMeasurePolyakovXY.cu` ~L426 | 分布测量时是否偏移中心。`1`=偏移（仅当 `MeasureDist=1` 时生效） |

**字典输出键**（`AddOneConfigurationResult`，键前缀 `measure{id}.`）：

| 键名 | 类型 | 条件 | 说明 |
|------|------|------|------|
| `PolyakovT` | CLGComplex | 总是 | t 方向 Polyakov loop（逐站点平均） |
| `PolyakovX` / `PolyakovY` / `PolyakovZ` | CLGComplex | `MeasureX/Y/Z=1` | X/Y/Z 方向 Polyakov loop |
| `PolyakovInner`, `PolyakovOverR`（+`Abs`） | CLGComplex / TArray | `MeasureDist=1` | 径向分布 |
| `PolyakovXSlice` / `PolyakovYSlice` / `PolyakovZSlice`（+`Abs`） | TArray | `XSlice/YSlice/ZSlice=1` | t 方向 loop 的逐切片轮廓 |
| `PolyakovX_YSlice` / `PolyakovX_ZSlice` / `PolyakovX_TSlice`（+`Abs`） | TArray | `MeasureX=1` 且对应切片开 | X 方向 loop 的逐切片轮廓 |
| `PolyakovY_XSlice` / `PolyakovY_ZSlice` / `PolyakovY_TSlice`（+`Abs`） | TArray | `MeasureY=1` 且对应切片开 | Y 方向 loop 的逐切片轮廓 |
| `PolyakovZ_XSlice` / `PolyakovZ_YSlice` / `PolyakovZ_TSlice`（+`Abs`） | TArray | `MeasureZ=1` 且对应切片开 | Z 方向 loop 的逐切片轮廓 |

**说明**：
- 所有切片值均归一化为切片内逐站点平均（复数与 `Abs` 通道一致，均为切片体积的倒数加权）。
- 切片轮廓同时写入 CSV 文件（`polyakov_XSlice.csv` 等），字典通道由 `AddOneConfigurationResult` 统一管理。
---

### CMeasurePolyakovXY3D — 3D Polyakov loop

**来源文件**：`Code/CLGLib/Measurement/CMeasurePolyakovXY3D.cu`

**读取函数**：`CMeasurePolyakovXY3D::InitialOtherParameters()` (lines 161-173)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `U1` | INT | `0` | 否 | `CMeasurePolyakovXY3D.cu` ~L167 | 使用 U(1) 规范场代替 SU(3)。`1`=启用 |
| `ShiftCenter` | INT | `0` | 否 | `CMeasurePolyakovXY3D.cu` ~L173 | 分布测量时偏移中心。`1`=偏移 |

---

### CMeasureWilsonLoop — Wilson loop

**来源文件**：`Code/CLGLib/Measurement/CMeasureWilsonLoop.cu`

| 键名 | 类型 | 默认值 | 必需 | 说明 |
|------|------|--------|------|------|
| （无特有参数） | — | — | — | `m_uiMaxLengthSq` 由格点尺寸自动计算 |

---

### CMeasureWilsonLoopXY — Wilson loop（XY 平面）

**来源文件**：`Code/CLGLib/Measurement/CMeasureWilsonLoopXY.cu`

| 键名 | 类型 | 默认值 | 必需 | 说明 |
|------|------|--------|------|------|
| （无特有参数） | — | — | — | `m_uiMaxLengthSq` 由格点尺寸自动计算 |

---

### CMeasureWilsonLoopWithPath — 自定义路径 Wilson loop

**来源文件**：`Code/CLGLib/Measurement/CMeasureWilsonLoopWithPath.cu`

**读取函数**：`CMeasureWilsonLoopWithPath::InitialOtherParameters()` (lines 141-164)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Path` | SCHAR[] | （无） | **是** | `CMeasureWilsonLoopWithPath.cu` ~L146 | Wilson loop 路径方向列表，如 `[1,2,-1,-2]` |
| `OnePoint` | INT | `0` | 否 | `CMeasureWilsonLoopWithPath.cu` ~L158 | 仅在单点测量。`1`=单点，`0`=全格点 |
| `Point` | SCHAR[] | （无） | 否 | `CMeasureWilsonLoopWithPath.cu` ~L164 | 单点坐标 `[x,y,z,t]`（仅当 `OnePoint=1` 时生效） |

---


---

[< 返回 17. 测量器目录](17-measurements.md) | [< 返回 yaml-reference 目录](home.md)
