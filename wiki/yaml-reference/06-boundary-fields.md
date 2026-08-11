> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 06. 边界场 (Boundary Fields)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateBoundaryFields()` (lines 665-693)

---

## 配置参数

边界场用于开放边界条件模拟，在特定边界上引入额外场。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldName` | string | `CFieldBoundaryGaugeSU3` | 否 | `CLGLibManager.cpp` ~L677 | 边界场类名。通过 `appCreate()` 实例化 |
| `FieldId` | INT | `1`（规范）或 `-1` | 否 | `CLGLibManager.cpp` ~L679 | 场唯一标识符 |

## 支持的边界场类

| 类名 | 说明 |
|------|------|
| `CFieldBoundaryGaugeSU3` | SU(3) 规范边界场 |
| `CFieldBoundaryFermionWilsonSquareSU3` | Wilson Dirac 费米子边界场 |

## 初始化行为

边界场的读取逻辑较为复杂：
1. 首先检查 `GaugeBoundary` 子参数
2. 然后检查 `BoundaryFermionField1`, `BoundaryFermionField2`, ... 子参数
3. 通过 `appCreate(FieldName)` 创建实例
4. 调用 `InitialOtherParameters()` 读取额外参数

---

[< 返回目录](home.md) | [< 上一章：规范场](05-gauge-fields.md) | [下一章：玻色子场 >](07-boson-fields.md)
