> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 05. 规范场 (Gauge Fields)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateGaugeFields()` (lines 450-527)

---

## 配置参数

YAML 键名为 `Gauge`（第一个）或 `Gauge2`, `Gauge3`, ...（后续）。通过 `GaugeFieldCount` 控制数量。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldName` | string | `CFieldGaugeSU3` | 否 | `CLGLibManager.cpp` ~L470 | 规范场类名。通过 `appCreate()` 实例化 |
| `FieldInitialType` | string | `EFIT_Random` | 否 | `CLGLibManager.cpp` ~L472 | 初始场类型。见下表 |
| `FieldId` | INT | `1` | 否 | `CLGLibManager.cpp` ~L474 | 场唯一标识符（BYTE）。ID 1 保留给主规范场 |
| `GaugeFileType` | string | （无） | 条件 | `CLGLibManager.cpp` ~L476 | 文件格式。当 `FieldInitialType == EFIT_ReadFromFile` 时必需 |
| `GaugeFileName` | string | （无） | 条件 | `CLGLibManager.cpp` ~L478 | 文件路径。当 `FieldInitialType == EFIT_ReadFromFile` 时必需 |
| `Period` | INT[] | `[1,1,1,1]` | 否 | `CField.cpp` ~L111 | 各方向边界条件。`1`=周期，`-1`=反周期 |

## 基类参数 (CField)

**来源文件**：`Code/CLGLib/Data/Field/CField.cpp`

**读取函数**：`CField::InitialOtherParameters()` (lines 111-125)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `GaugeFields` | BYTE[] | `[1]` | 否 | `CField.cpp` ~L117 | 此对象作用的规范场 ID 列表 |
| `BosonFields` | BYTE[] | `[]` | 否 | `CField.cpp` ~L119 | 此对象作用的玻色子场 ID 列表 |
| `Dynamic` | INT | `1` | 否 | `CField.cpp` ~L113 | 场是否随演化更新。`1`=动态，`0`=静态 |

**注意**：如果 `GaugeFields` 和 `BosonFields` 都未指定或为空，默认设置为 `GaugeFields: [1]`。

## 初始场类型 (EFieldInitialType)

| 枚举值 | 说明 |
|--------|------|
| `EFIT_Random` | 随机初始化 |
| `EFIT_Identity` | 单位矩阵（冷启动） |
| `EFIT_ReadFromFile` | 从文件读取 |
| `EFIT_RandomGaussian` | 高斯随机（主要用于费米子场） |
| `EFIT_Zero` | 零初始化 |
| `EFIT_RandomZ4` | Z4 随机（用于随机源） |

## 文件格式 (EGaugeFieldFileType)

| 枚举值 | 说明 |
|--------|------|
| `EFFT_CLGBin` | CLGLib 二进制格式（默认，含 MD5 校验） |
| `EFFT_CLGBinCompressed` | 压缩二进制 |
| `EFFT_CLGBinDouble` | 双精度二进制 |
| `EFFT_CLGBinFloat` | 单精度二进制 |
| `EFFT_BridgePPTXT` | Bridge++ 文本格式 |
| `EFFT_BridgePPBin` | Bridge++ 二进制格式 |

---

## 支持的规范场类

通过 `appCreate(FieldName)` 工厂实例化。需确认 `Code/CLGLib/Core/CLGSetup.h` 中对应宏已启用。

### 连续规范群

| 类名 | 宏开关 | 说明 | 额外参数 |
|------|--------|------|----------|
| `CFieldGaugeSU3` | `_CLG_SU3_GAUGE` | 标准 SU(3) 规范场 | 无 |
| `CFieldGaugeSU3TreeImproved` | `_CLG_SU3_GAUGE` | 树级改进 SU(3) | `RectOverPlaq` |
| `CFieldGaugeSU3OneLoopImproved` | `_CLG_SU3_GAUGE` | 一圈改进 SU(3) | `Cr`, `Ct`, `U0`, `Nf` |
| `CFieldGaugeSU3_12` | `_CLG_SU3_GAUGE` | SU(3) 12-参数优化表示 | 无 |
| `CFieldGaugeSU2` | `_CLG_SU2_GAUGE` | SU(2) 规范场 | 无 |
| `CFieldGaugeU1` | `_CLG_U1_GAUGE` | U(1) 规范场 | 无 |
| `CFieldGaugeU1Real` | `_CLG_U1_REALGAUGE` | U(1) 实规范场（含化学势/电磁场） | `ChemicalType`, `ChemicalValue`, `EzType`, `EzValue`, `BzType`, `BzValue`, `XYShiftCenter`, `CacheFmunu` |
| `CFieldGaugeSU4` | `_CLG_SU4_GAUGE` | SU(4) 规范场 | 无 |
| `CFieldGaugeSU5` | `_CLG_SU5_GAUGE` | SU(5) 规范场 | 无 |
| `CFieldGaugeSU6` | `_CLG_SU6_GAUGE` | SU(6) 规范场 | 无 |
| `CFieldGaugeSU7` | `_CLG_SU7_GAUGE` | SU(7) 规范场 | 无 |
| `CFieldGaugeSU8` | `_CLG_SU8_GAUGE` | SU(8) 规范场 | 无 |

### 离散规范群

| 类名 | 宏开关 | 说明 | 额外参数 |
|------|--------|------|----------|
| `CFieldGaugeZ2` | `_CLG_Z2_GAUGE` | Z₂ 离散规范场 | 无 |
| `CFieldGaugeZ3` | `_CLG_ZN_GAUGE` | Z₃ 离散规范场 | 无 |
| `CFieldGaugeZ4` | `_CLG_ZN_GAUGE` | Z₄ 离散规范场 | 无 |
| `CFieldGaugeZ5` | `_CLG_Z5_GAUGE` | Z₅ 离散规范场 | 无 |
| `CFieldGaugeZ6` | `_CLG_Z6_GAUGE` | Z₆ 离散规范场 | 无 |
| `CFieldGaugeD3` | `_CLG_DN_GAUGE` | D₃ 二面体群 | 无 |
| `CFieldGaugeD4` | `_CLG_D4_GAUGE` | D₄ 二面体群 | 无 |
| `CFieldGaugeD8` | `_CLG_D8_GAUGE` | D₈ 二面体群 | 无 |

### Dirichlet 边界条件变体

以下类在对应基础类名后加 `D`，支持 Dirichlet 边界条件：

`CFieldGaugeU1D`, `CFieldGaugeSU2D`, `CFieldGaugeSU3D`, `CFieldGaugeSU4D` — `SU8D`

这些类本身无额外参数，通过 `CFieldGaugeLinkDirichlet` 模板包装实现。

---

## 规范场特有参数

### SU(3) 树级改进 (`CFieldGaugeSU3TreeImproved`)

**来源文件**：`Code/CLGLib/Data/Field/Gauge/CFieldGaugeSU3TreeImproved.h`

**读取函数**：`CFieldGaugeTreeImproved::InitialOtherParameters()` (line 49)

| 键名 | 类型 | 默认值 | 必需 | 说明 |
|------|------|--------|------|------|
| `RectOverPlaq` | DOUBLE | `-0.05` | 否 | Rectangle 作用量相对于 Plaquette 的权重系数。默认 `-0.05` 对应树级 Symanzik 改进作用量 |

**说明**：`CFieldGaugeSU3TreeImproved` 支持各向异性能量/力虚函数（plaquette+rectangle 均按平面加权），由 `CActionGaugePlaquetteAnisotropic` 驱动；`CFieldGaugeSU3OneLoopImproved` 同样支持（含 twisted loop，按三元组是否含时间方向加权）。Dirichlet×anisotropy 路径未做物理验证。

### SU(3) 一圈改进 (`CFieldGaugeSU3OneLoopImproved`)

**来源文件**：`Code/CLGLib/Data/Field/Gauge/CFieldGaugeSU3OneLoopImproved.h`

**读取函数**：`CFieldGaugeOneLoopImproved::InitialOtherParameters()` (line 49)

| 键名 | 类型 | 默认值 | 必需 | 说明 |
|------|------|--------|------|------|
| `Cr` | DOUBLE | `-0.05` | 否* | Rectangle 系数。若 **同时** 提供 `Cr` 和 `Ct`，则直接使用 |
| `Ct` | DOUBLE | `0.0` | 否* | Twisted-loop 系数。若 **同时** 提供 `Cr` 和 `Ct`，则直接使用 |
| `U0` | DOUBLE | `1.0` | 否 | Mean link (tadpole) 参数。仅在未同时提供 `Cr` 和 `Ct` 时使用 |
| `Nf` | INT | `2` | 否 | 费米子味数。仅在未同时提供 `Cr` 和 `Ct` 时使用 |

**行为说明**：若 YAML 中同时存在 `Cr` 和 `Ct`，则直接使用这两个值（`m_bDirectSetCoefficients = TRUE`）。否则，根据 `U0` 和 `Nf` 通过一圈微扰论公式自动计算（参考 arXiv:1004.0342 Eq. A2）：
- `Cr = (-1/(20*U0^2)) * (1 - (0.6264 - 1.1746*Nf) * log(U0))`
- `Ct = (0.0433 - 0.0156*Nf) * log(U0) / (U0^2)`

### U(1) 实规范场 (`CFieldGaugeU1Real`)

**来源文件**：`Code/CLGLib/Data/Field/Gauge/CFieldGaugeU1Real.cu`

**读取函数**：`CFieldGaugeU1Real::InitialOtherParameters()` (line 1143)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `ChemicalType` | string | `EURT_None` | 否 | `CFieldGaugeU1Real.cu` ~L1158 | 化学势类型。见下表 |
| `ChemicalValue` | Real | `0.0` | 否 | `CFieldGaugeU1Real.cu` ~L1163 | 化学势数值（仅当 `ChemicalType != EURT_None` 时读取） |
| `EzType` | string | `EURT_None` | 否 | `CFieldGaugeU1Real.cu` ~L1168 | Z 方向电场类型 |
| `EzValue` | Real | `0.0` | 否 | `CFieldGaugeU1Real.cu` ~L1173 | Z 方向电场数值（仅当 `EzType != EURT_None` 时读取） |
| `BzType` | string | `EURT_None` | 否 | `CFieldGaugeU1Real.cu` ~L1178 | Z 方向磁场类型 |
| `BzValue` | Real | `0.0` | 否 | `CFieldGaugeU1Real.cu` ~L1183 | Z 方向磁场数值（仅当 `BzType != EURT_None` 时读取） |
| `XYShiftCenter` | INT | `1` | 否 | `CFieldGaugeU1Real.cu` ~L1188 | XY 平面是否偏移中心。`1`=偏移 |
| `CacheFmunu` | INT | `0` | 否 | `CFieldGaugeU1Real.cu` ~L1192 | 是否缓存 F_munu 张量。`1`=缓存 |

**EU1RealType 枚举值**：

| 值 | 说明 |
|----|------|
| `EURT_None` | 无 |
| `EURT_ImagineChemical` | 虚化学势 |
| `EURT_E_t` | 时间方向电场 |
| `EURT_E_z` | Z 方向电场 |
| `EURT_Bp_x` | X 方向磁通 |
| `EURT_Bp_y` | Y 方向磁通 |
| `EURT_Bp_xy` | XY 平面磁通 |
| `EURT_Bp_x_notwist` | X 方向磁通（无 twist） |
| `EURT_Bp_y_notwist` | Y 方向磁通（无 twist） |
| `EURT_Bp_xy_notwist` | XY 平面磁通（无 twist） |
| `EURT_ExBz_ty` | Ex-Bz 组合（ty 方向） |

---

[< 返回目录](home.md) | [< 上一章：索引与边界](04-index-boundary.md) | [下一章：边界场 >](06-boundary-fields.md)
