> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 09b. 费米子作用量

## 费米子作用量

### CActionFermionKS

**来源文件**：`Code/CLGLib/Data/Action/CActionFermionKS.cpp`

**读取函数**：`CActionFermionKS::InitialOtherParameters()` (lines 24-36)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldId` | INT | `-1` | **强烈建议** | `CActionFermionKS.cpp` ~L30 | 关联的费米子场 ID。`-1` 表示未指定 |

### CActionFermionKSCombined

**来源文件**：`Code/CLGLib/Data/Action/CActionFermionKSCombined.cpp`

**读取函数**：`CActionFermionKSCombined::InitialOtherParameters()` (lines 24-32)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldIds` | INT[] | （无） | **是** | `CActionFermionKSCombined.cpp` ~L28 | 关联的费米子场 ID 列表。必须至少有一个元素 |

### CActionFermionKSImprove

**来源文件**：`Code/CLGLib/Data/Action/CActionFermionKSImprove.cpp`

**读取函数**：`CActionFermionKSImprove::InitialOtherParameters()` (lines 31-45)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldId` | INT | `-1` | **强烈建议** | 继承自 `CActionFermionKS` | 关联的费米子场 ID |
| `HISQ` | INT | `0` | 否 | `CActionFermionKSImprove.cpp` ~L41 | 是否使用 HISQ 改进。非零则启用 |

**注意**：若 `m_byGaugeFieldIds.Num() < 1`，自动添加规范场 ID `1`（line 37）。

### CActionFermionKSImproveCombined

**来源文件**：`Code/CLGLib/Data/Action/CActionFermionKSImproveCombined.cpp`

无额外参数。继承自 `CActionFermionKSCombined`：
- `FieldIds` (INT[]) — 必需

### CActionFermionHISQCombined

**来源文件**：`Code/CLGLib/Data/Action/CActionFermionHISQCombined.cpp`

**读取函数**：`CActionFermionHISQCombined::InitialOtherParameters()` (lines 24-32)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldIds` | INT[] | （无） | **是** | `CActionFermionHISQCombined.cpp` ~L28 | 关联的费米子场 ID 列表。必须至少有一个元素 |

---

## 依赖关系

配置费米子作用量时，通常需要同时注册以下配套对象：

| 作用量 A | 必须配套 B | 说明 |
|---|---|---|
| `CActionFermionHISQCombined` | `GaugeSmearing: CGaugeSmearingHISQSU3` | HISQ 作用量直接取 smearing 的 effective gauge。见 [12-gauge-smearing.md](12-gauge-smearing.md) |
| `CActionFermionKSImprove` / `CActionFermionKSImproveCombined` | 任意 `GaugeSmearing`（Stout / ASQTAD / HISQ） | 改进作用量通过 `appGetGaugeSmearing` 取 effective gauge。例如 StoutLink Clover Wilson 配 `CGaugeSmearingStoutSU3`。见 [12-gauge-smearing.md](12-gauge-smearing.md) |
| `CActionFermionKS` / `CActionFermionKSCombined` | 无 smearing 要求 | 直接使用原始规范场 |
| 任意费米子作用量（HMC simulation） | 每个相关场配 `MSSolver` 或 `Solver` | HISQ / improved KS / rational Wilson 用 `MSSolver`；`ER_NoRational` 可用 `Solver`。见 [11-multi-shift-solvers.md](11-multi-shift-solvers.md)、[10-solvers.md](10-solvers.md) |
| 任意费米子作用量（measurement） | 对应场的 `Solver` | 测量中反演 `D` / `D†D`。见 [10-solvers.md](10-solvers.md) |

> **注意**：上述求解器依赖实际来自 `CFieldFermion` 及其子类的 `CalculateForce` / `InverseD` / `InverseDDdagger` 等函数。作用量本身只是把这些调用组织起来。

---


---

[< 返回 09. 作用量](09-actions.md) | [< 返回 yaml-reference 目录](home.md)
