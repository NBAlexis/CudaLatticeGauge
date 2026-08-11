> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 15. 更新器 (Updators)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateUpdator()` (lines 763-800)

---

## 配置参数

YAML 键名为 `Updator`。注意：目前框架只支持**一个**更新器。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `UpdatorType` | string | `CHMC` | 否 | `CLGLibManager.cpp` ~L767 | 更新器类名。通过 `appCreate()` 实例化 |
| `IntegratorType` | string | `CIntegratorLeapFrog` | 否 | `CLGLibManager.cpp` ~L769 | 积分器类名（仅 HMC 需要）。见 [16-integrators.md](16-integrators.md) |
| `Metropolis` | INT | `0` | 否 | `CLGLibManager.cpp` ~L771 | 是否启用 Metropolis 检验。`1`=启用，`0`=禁用 |
| `SaveConfiguration` | INT | `0` | 否 | `CLGLibManager.cpp` ~L773 | 是否保存组态到磁盘。`1`=保存 |
| `ConfigurationFilePrefix` | string | `Untitled` | 否 | `CLGLibManager.cpp` ~L775 | 保存文件的前缀名。实际文件名为 `prefix_1.con`, `prefix_2.con` |
| `Adaptive` | INT | `0` | 否 | `CLGLibManager.cpp` ~L777 | 是否启用自适应步长。`1`=启用 |
| `ReportMeasure` | INT | `1` | 否 | `CLGLibManager.cpp` ~L779 | 更新过程中是否报告测量结果。`1`=报告 |
| `Skip` | INT | `1` | 否 | `CLGLibManager.cpp` ~L781 | 测量间隔（每 Skip 步测量一次） |
| `MinMaxStep` | INT[] | `[5, 100]` | 条件 | `CLGLibManager.cpp` ~L783 | `[minStep, maxStep]`。仅当 `Adaptive=1` 时使用 |
| `GrowReduceThreshold` | Real[] | `[-0.3, 0.03]` | 条件 | `CLGLibManager.cpp` ~L785 | `[growThreshold, reduceThreshold]`。仅当 `Adaptive=1` 时使用 |

## HMC 更新器参数

**来源文件**：`Code/CLGLib/Update/Continous/CHMC.cpp`

**读取函数**：`CHMC::Initial()` (lines 23-78)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Metropolis` | INT | `0` | 否 | `CHMC.cpp` ~L29 | 是否启用 Metropolis 检验。`1`=启用 |
| `SaveConfiguration` | INT | `0` | 否 | `CHMC.cpp` ~L33 | 是否保存组态。`1`=保存 |
| `Adaptive` | INT | `0` | 否 | `CHMC.cpp` ~L37 | 是否启用自适应步长。`1`=启用 |
| `ReportMeasure` | INT | `1` | 否 | `CHMC.cpp` ~L41 | 是否报告测量结果。`1`=报告 |
| `Skip` | INT | `1` | 否 | `CHMC.cpp` ~L45 | 测量间隔 |
| `ConfigurationFilePrefix` | string | `Untitled` | 否 | `CHMC.cpp` ~L50 | 保存文件前缀 |
| `MinMaxStep` | INT[] | `[5, 100]` | 否 | `CHMC.cpp` ~L58 | 自适应步长范围 `[min, max]` |
| `GrowReduceThreshold` | Real[] | `[-0.3, 0.03]` | 否 | `CHMC.cpp` ~L68 | 自适应阈值 `[grow, reduce]` |

## 支持的更新器类

| 类名 | 说明 |
|------|------|
| `CHMC` | Hybrid Monte Carlo（混合蒙特卡洛） |
| `CHeatbath` | 热浴更新器（用于离散规范群） |

## 重要注意

### CHeatbath 的特殊行为

**来源文件**：`Code/CLGLib/Update/Discrete/CHeatbath.cu`

`CHeatbath::Initial()` 当前为**空实现**（参数读取代码被注释掉）。`SaveConfiguration`、`Metropolis` 等参数对 Heatbath **不生效**。如需保存组态，必须在 C++ 代码中手动调用 `SaveToFile()`。详见 [Tutorial/03-pure-z2.md](../../Tutorial/03-pure-z2.md)。

---

[< 返回目录](home.md) | [< 上一章：规范固定](14-gauge-fixing.md) | [下一章：积分器 >](16-integrators.md)
