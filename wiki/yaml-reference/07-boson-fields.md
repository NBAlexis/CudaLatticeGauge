> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 07. 玻色子场 (Boson Fields)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateBosonFields()` (lines 529-596)

---

## 配置参数

YAML 键名为 `BosonField1`, `BosonField2`, ...。通过 `BosonFieldCount` 控制数量。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldName` | string | `CFieldBosonU1` | 否 | `CLGLibManager.cpp` ~L549 | 玻色子场类名。通过 `appCreate()` 实例化 |
| `FieldInitialType` | string | `EFIT_Random` | 否 | `CLGLibManager.cpp` ~L551 | 初始场类型。同规范场 |
| `FieldId` | INT | 自动分配 | 否 | `CLGLibManager.cpp` ~L553 | 场唯一标识符。若未指定或冲突，自动分配 |
| `Period` | INT[] | `[1,1,1,1]` | 否 | `CField.cpp` ~L111 | 各方向边界条件 |

## 基类参数 (CFieldBoson)

**来源文件**：`Code/CLGLib/Data/Field/CFieldBoson.cpp`

**读取函数**：`CFieldBoson::InitialOtherParameters()` (lines 53-68)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Constant` | INT | `0` | 否 | `CFieldBoson.cpp` ~L57 | 场是否为常数（非动态）。`1`=常数 |
| `NoGauge` | INT | `0` | 否 | `CFieldBoson.cpp` ~L59 | 若 `1`，清除 `m_byGaugeFieldIds`（不耦合规范场） |

## 支持的玻色子场类

| 类名 | 说明 |
|------|------|
| `CFieldBosonU1` | U(1) 玻色子场 |
| `CFieldBosonSU3` | SU(3) 玻色子场 |


---

## 补充参数（代码中存在但文档中缺失）

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Omega` | DOUBLE | （无） | required | Code/CLGLib/Data/Field/Boson/CFieldBosonVNRotation.h ~L98 | 旋转坐标系角速度 |
| `ShiftCenter` | INT | （无） | required | Code/CLGLib/Data/Field/Boson/CFieldBosonVNRotation.h ~L104 | 是否将中心平移到原点。`1`=平移 |
| `ValueList` | array<Real> | （无） | required | Code/CLGLib/Data/Field/Boson/CFieldBosonReal3D.cu ~L43 | 与 `ZSliceList` 一一对应的固定 z 切片场值 |
| `ZSliceList` | array<BYTE> | （无） | required | Code/CLGLib/Data/Field/Boson/CFieldBosonReal3D.cu ~L42 | 需要固定场值的 z 切片索引列表（代码另有 1 处读取） |

---

[< 返回目录](home.md) | [< 上一章：边界场](06-boundary-fields.md) | [下一章：费米子场 >](08-fermion-fields.md)
