> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 16. 积分器 (Integrators)

**来源文件**：`Code/CLGLib/Update/Continous/CIntegrator.cpp`

**读取函数**：`CIntegrator::Initial()` (lines 52-147)

积分器在更新器初始化时通过 `IntegratorType` 指定并创建，其参数从顶层 `params` 读取（不是子树）。

---

## 基类配置参数

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `IntegratorType` | string | `CIntegratorLeapFrog` | 否 | `CIntegrator.cpp` ~L56 | 积分器类名。通过 `appCreate()` 实例化 |
| `IntegratorStep` | INT | `50` | 否 | `CIntegrator.cpp` ~L61 | 积分步数 |
| `IntegratorStepDebug` | INT | `50` | 否 | `CIntegrator.cpp` ~L63 | Debug 模式下覆盖 `IntegratorStep` |
| `IntegratorStepLength` | Real | `1.0` | 否 | `CIntegrator.cpp` ~L67 | 总轨迹长度（`Tau = Epsilon * Step`） |
| `IntegratorStepWarmup` | INT | `0` | 否 | `CIntegrator.cpp` ~L72 | 热化步数 |
| `DebugForce` | INT | `0` | 否 | `CIntegrator.cpp` ~L76 | 是否调试输出力。`1`=启用 |
| `BindDir` | INT | `0` | 否 | `CIntegrator.cpp` ~L80 | 绑定方向 |
| `BackupFieldTypeName` | string | `Same` | 否 | `CIntegrator.cpp` ~L136 | 备份场类型。`CFieldGaugeSU3_12` 或 `Same`（表示与 gauge 同类型） |

## 支持的积分器类

| 类名 | 说明 |
|------|------|
| `CIntegratorLeapFrog` | 标准 LeapFrog 积分器 |
| `CIntegratorOmelyan` | Omelyan 改进积分器 |
| `CIntegratorForceGradient` | 力梯度积分器 |
| `CIntegratorNestedLeapFrog` | 嵌套 LeapFrog |
| `CIntegratorNestedOmelyan` | 嵌套 Omelyan |
| `CIntegratorNestedForceGradient` | 嵌套力梯度积分器 |
| `CIntegratorNested11Stage` | 嵌套 11-stage 积分器 |
| `CIntegratorMultiLevelOmelyan` | 多层嵌套 Omelyan |
| `CIntegratorMultiLevelNestedForceGradient` | 多层嵌套力梯度 |

## Omelyan 参数

**来源文件**：`Code/CLGLib/Update/Continous/CIntegratorOmelyan.cpp`

**读取函数**：`CIntegratorOmelyan::Initial()` (lines 18-26)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Omelyan2Lambda` | Real | `OmelyanLambda2`（常量） | 否 | `CIntegratorOmelyan.cpp` ~L22 | Omelyan 2lambda 参数 |

## 嵌套积分器参数

**来源文件**：`Code/CLGLib/Update/Continous/CIntegrator.cpp`

**读取函数**：`CNestedIntegrator::Initial()` (lines 413-429)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `NestedStep` | INT | `3` | 否 | `CIntegrator.cpp` ~L418 | 嵌套步数 |
| `InnerLeapfrog` | INT | `0` | 否 | `CIntegrator.cpp` ~L427 | 内层是否使用 LeapFrog。`1`=使用 |

**注意**：以下积分器继承自 `CNestedIntegrator`，均支持上述嵌套参数：
- `CIntegratorNestedLeapFrog`
- `CIntegratorNestedOmelyan`（额外支持 `Omelyan2Lambda`）
- `CIntegratorNestedForceGradient`
- `CIntegratorNested11Stage`

## 嵌套 11-Stage 积分器参数

**来源文件**：`Code/CLGLib/Update/Continous/CIntegratorNested11Stage.cpp`

**读取函数**：`CIntegratorNested11Stage::Initial()` (lines 18-44)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Position` | INT | `0` | 否 | `CIntegratorNested11Stage.cpp` ~L22 | 位置参数。`1`=启用特定位置模式 |
| `Rho` | Real | `_11Stage_Rho` | 否 | `CIntegratorNested11Stage.cpp` ~L25 | Rho 系数 |
| `Theta` | Real | `_11Stage_Theta` | 否 | `CIntegratorNested11Stage.cpp` ~L30 | Theta 系数 |
| `VarTheta` | Real | `_11Stage_VarTheta` | 否 | `CIntegratorNested11Stage.cpp` ~L35 | VarTheta 系数 |
| `Lambda` | Real | `_11Stage_Lambda` | 否 | `CIntegratorNested11Stage.cpp` ~L40 | Lambda 系数 |

## 多层嵌套积分器参数

**来源文件**：`Code/CLGLib/Update/Continous/CIntegrator.cpp`

**读取函数**：`CMultiLevelNestedIntegrator::Initial()` (lines 494-551)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `NestedSteps` | UINT[] | （无） | **是** | `CIntegrator.cpp` ~L499 | 各层嵌套步数列表 |
| `NestedActionList0` | UINT[] | （无） | **是** | `CIntegrator.cpp` ~L530 | 第 0 层作用量索引列表 |
| `NestedActionList1` | UINT[] | （无） | 条件 | `CIntegrator.cpp` ~L530 | 第 1 层作用量索引列表 |
| ... | ... | ... | ... | ... | 依此类推 |
| `InnerLeapfrog` | INT | `0` | 否 | `CIntegrator.cpp` ~L549 | 最内层是否使用 LeapFrog。`1`=使用 |

**注意**：以下积分器继承自 `CMultiLevelNestedIntegrator`：
- `CIntegratorMultiLevelOmelyan`（额外支持 `Omelyan2Lambda`）
- `CIntegratorMultiLevelNestedForceGradient`

---

[< 返回目录](home.md) | [< 上一章：更新器](15-updators.md) | [下一章：测量器 >](17-measurements.md)
