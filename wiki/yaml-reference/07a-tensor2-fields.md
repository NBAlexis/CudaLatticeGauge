> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 07a. Tensor2 场 (Tensor2 / Plaquette Fields)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateTensor2Fields()` (lines 598-666)

---

## 配置参数

YAML 键名为 `Tensor2Field1`, `Tensor2Field2`, ...。通过 `Tensor2FieldCount` 控制数量。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldName` | string | `CFieldTensor2SU3` | 否 | `CLGLibManager.cpp` ~L605 | Tensor2 场类名。通过 `appCreate()` 实例化 |
| `FieldInitialType` | string | `EFIT_Random` | 否 | `CLGLibManager.cpp` ~L607 | 初始场类型。同规范场 |
| `FieldId` | INT | 自动分配 | 否 | `CLGLibManager.cpp` ~L618 | 场唯一标识符。若未指定或冲突，自动分配 |
| `Period` | INT[] | `[1,1,1,1]` | 否 | `CLGLibManager.cpp` ~L636 | 各方向边界条件 |

## 支持的 Tensor2 场类

| 类名 | 元素类型 | 每站点浮点数 | 说明 |
|------|----------|--------------|------|
| `CFieldTensor2Real` | Real | 6 | 实数 plaquette 场 |
| `CFieldTensor2Complex` | CLGComplex | 12 | 复数 plaquette 场 |
| `CFieldTensor2SU2` | deviceSU2 | 48 | SU(2) plaquette 场 |
| `CFieldTensor2SU3` | deviceSU3 | 108 | SU(3) plaquette 场 |

**说明**：
- Tensor2 场即 plaquette 场：每个站点 6 个元素，按 `data[plaquetteIndex * siteCount + siteIndex]`（分量在外，与 `CStapleCache` 的 Fmunu 缓冲布局一致）存储，顺序为 xy, xz, xt, yz, yt, zt（与 `_plaq_idx` 一致）。
- 中间层为 `CFieldTensor2`（`Code/CLGLib/Data/Field/CFieldTensor2.h`），模板层为 `CFieldTensor2T<T>`，BLAS（Axpy 系列、Mul、ScalarMultply、Dagger）与文件 IO 全部复用 `CCommonKernelField<T>`；模板层按惯例 typedef 了 `_FieldKernel`（= `CCommonKernelField<T>`）与 `_Tensor2Kernel`（= `CFieldTensor2Kernel<T>`），子类可直接使用父类 kernel。
- 初始化使用 Tensor2 专用 kernel（随机数种子表每站点只有 4 条流，需按站点顺序抽取 6 个分量）；归约（Dot、GetLength、Sum）直接按 6 个分量段调用 `CCommonKernelField`，无专用 kernel。
- Tensor2 场不参与 hopping，不会烘焙 move index（BakeMoveIndex 只作用于 boson/fermion 场）。
- Tensor2 场登记在 `CLatticeData::m_pTensor2Field`，是 HMC 的**工作场**：`CIntegrator` 为每个 `Dynamic`（默认 1）的 tensor2 场创建副本，轨迹中 action 可通过 `OnFinishTrajectory(bWillBeAccept, ..., tensor2Fields)` 修改它，接受时拷贝回 lattice，拒绝时丢弃。
- `CAction::Energy` 与 `CAction::OnFinishTrajectory` 多参数版本均携带 tensor2 场数组（`tensor2Num`/`tensor2Fields`，位于 boson 参数之后）；`Energy` 中场为 const，`OnFinishTrajectory` 中场可修改（数组本身 const）。
- 环境字段束已统一为 6 元组 `(gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields)`：`_FIELDS` 宏（`CLatticeData.h`）、`CField::ApplyOperator`、`CFieldFermion` 算子族（`D`/`Ddagger`/`DD`/`DDdagger`/`InverseD` 系列/有理逼近/`D_MC`/`D_MD`/`Energy`）、`CSLASolver::Solve`、`CMultiShiftSolver::Solve`、`CMeasure::OnConfigurationAccepted`/`OnConfigurationAcceptedZ4`/`ExportDiagnal` 均按此约定传参；目前 tensor2 在算子内未被使用，为将来依赖 tensor2 的算子与测量预留。
- 冒烟测试：`TestSaveLoadTensor2`（`TestSuit_FileIO.yaml`）。

---

[< 返回目录](home.md) | [< 上一章：玻色子场](07-boson-fields.md) | [下一章：费米子场](08-fermion-fields.md)
