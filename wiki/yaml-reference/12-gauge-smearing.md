> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 12. 规范平滑 (Gauge Smearing)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateGaugeSmearing()` (lines 868-877)

---

## 配置参数

YAML 键名为 `GaugeSmearing`（第一个）或 `GaugeSmearing2`, `GaugeSmearing3`, ...。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `SmearingName` | string | `CGaugeSmearingAPEStout` | 否 | `CLGLibManager.cpp` ~L872 | 规范平滑类名。通过 `appCreate()` 实例化 |

## 基类参数 (CGaugeSmearing)

**来源文件**：`Code/CLGLib/GaugeSmearing/CGaugeSmearing.cu`

**读取函数**：`CGaugeSmearing::Initial()` (lines 79-108)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `HasT` | INT | `0` | 否 | `CGaugeSmearing.cu` ~L83 | 是否包含 T 方向。`1`=包含 |
| `FieldId` | INT | `1` | 否 | `CGaugeSmearing.cu` ~L87 | 平滑的规范场 ID。`0`=主规范场 |
| `Iterate` | INT | `1` | 否 | `CGaugeSmearing.cu` ~L91 | 平滑迭代次数 |
| `Update` | INT | `0` | 否 | `CGaugeSmearing.cu` ~L106 | 更新模式。代码中用于控制更新行为 |

## 支持的规范平滑类

### CGaugeSmearingAPEStout — APE + Stout 组合平滑

**来源文件**：`Code/CLGLib/GaugeSmearing/CGaugeSmearingAPEStout.cu`

**读取函数**：`CGaugeSmearingAPEStout::InitialOtherParameters()` (lines 46-56)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Rho` | Real | `0.1` | 否 | `CGaugeSmearingAPEStout.cu` ~L50 | Stout 平滑参数 rho |

### CGaugeSmearingStoutSU3 — Stout 平滑

**来源文件**：`Code/CLGLib/GaugeSmearing/CGaugeSmearingStout.cu`

**读取函数**：`CGaugeSmearingStoutSU3::InitialOtherParameters()` (lines 1171-1193)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `rho` | Real | `0.25` | 否 | `CGaugeSmearingStout.cu` ~L1176 | Stout 参数。若指定则同时设置 `rhojk` 和 `rho4mu` |
| `rhojk` | Real | `0.25` | 否 | `CGaugeSmearingStout.cu` ~L1183 | 空间-空间平面 Stout 参数 |
| `rho4mu` | Real | `0.25` | 否 | `CGaugeSmearingStout.cu` ~L1185 | 时间-空间平面 Stout 参数 |
| `Cache` | string | `ESLC_Full` | 否 | `CGaugeSmearingStout.cu` ~L1189 | 缓存模式。见下表 |

**Stout 缓存模式枚举值**：

| 枚举值 | 说明 |
|--------|------|
| `ESLC_Full` | 全缓存 |
| `ESLC_Median` | 中等缓存 |
| `ESLC_Small` | 最小缓存 |

### CGaugeSmearingASQTADSU3 — ASQTAD 改进平滑

**来源文件**：`Code/CLGLib/GaugeSmearing/CGaugeSmearingASQTAD.cu`

**读取函数**：`CGaugeSmearingASQTADSU3::InitialOtherParameters()` (lines 3560-3620)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Proj` | INT | `1` | 否 | `CGaugeSmearingASQTAD.cu` ~L3584 | 是否做投影。`1`=投影到规范群 |
| `ProjSUN` | INT | `1` | 否 | `CGaugeSmearingASQTAD.cu` ~L3588 | 是否投影到 SU(N)。`1`=启用 |
| `UseCaylayHamilton` | INT | `1` | 否 | `CGaugeSmearingASQTAD.cu` ~L3592 | 是否使用 Cayley-Hamilton 方法。`1`=启用 |
| `Origin` | Real | `0.125` | 否 | `CGaugeSmearingASQTAD.cu` ~L3612 | 原始链接权重 |
| `Fat3` | Real | `0.0625` | 否 | `CGaugeSmearingASQTAD.cu` ~L3613 | 3-link staple 权重 |
| `Fat5` | Real | `0.015625` | 否 | `CGaugeSmearingASQTAD.cu` ~L3614 | 5-link staple 权重 |
| `Fat7` | Real | `0.00260417` | 否 | `CGaugeSmearingASQTAD.cu` ~L3615 | 7-link staple 权重 |
| `Lepage` | Real | `0.0` | 否 | `CGaugeSmearingASQTAD.cu` ~L3616 | Lepage 项权重 |

### CGaugeSmearingHISQSU3 — HISQ 改进平滑

**来源文件**：`Code/CLGLib/GaugeSmearing/CGaugeSmearingHISQ.cu`

**读取函数**：`CGaugeSmearingHISQSU3::InitialOtherParameters()` (lines 85-164)

> **注意**：配合 `CFieldFermionHISQSU3` / `CActionFermionHISQCombined` 时必须使用 **`CGaugeSmearingHISQSU3`**，而不是 `CGaugeSmearingASQTADSU3`。
>
> HISQ 实现内部一次性完成 Level 1 + Level 2 smearing，**不读取基类的 `HasT` 和 `Iterate` 参数**。`Iterate` 默认即为 `1`，HISQ 不会重复迭代；`HasT` 对 HISQ 无实际作用。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Proj` | INT | `1` | 否 | `CGaugeSmearingHISQ.cu` ~L96 | 是否做投影。`1`=投影 |
| `ProjSUN` | INT | `1` | 否 | `CGaugeSmearingHISQ.cu` ~L100 | 是否投影到 SU(N)。`1`=启用 |
| `UseCaylayHamilton` | INT | `1` | 否 | `CGaugeSmearingHISQ.cu` ~L104 | 是否使用 Cayley-Hamilton 方法。`1`=启用 |
| `U0` | Real | `1.0` | 否 | `CGaugeSmearingHISQ.cu` ~L125 | 若指定，自动重标度 Level 1 系数 |
| `OriginL1` | Real | `0.125` | 否 | `CGaugeSmearingHISQ.cu` ~L127 | Level 1 原始链接权重 |
| `Fat3L1` | Real | `0.0625` | 否 | `CGaugeSmearingHISQ.cu` ~L128 | Level 1 3-link 权重 |
| `Fat5L1` | Real | `0.015625` | 否 | `CGaugeSmearingHISQ.cu` ~L129 | Level 1 5-link 权重 |
| `Fat7L1` | Real | `0.00260417` | 否 | `CGaugeSmearingHISQ.cu` ~L130 | Level 1 7-link 权重 |
| `LepageL1` | Real | `0.0` | 否 | `CGaugeSmearingHISQ.cu` ~L131 | Level 1 Lepage 权重 |
| `OriginL2` | Real | `1.0` | 否 | `CGaugeSmearingHISQ.cu` ~L151 | Level 2 原始链接权重 |
| `Fat3L2` | Real | `0.0625` | 否 | `CGaugeSmearingHISQ.cu` ~L152 | Level 2 3-link 权重 |
| `Fat5L2` | Real | `0.015625` | 否 | `CGaugeSmearingHISQ.cu` ~L153 | Level 2 5-link 权重 |
| `Fat7L2` | Real | `0.00260417` | 否 | `CGaugeSmearingHISQ.cu` ~L154 | Level 2 7-link 权重 |
| `LepageL2` | Real | `-0.125` | 否 | `CGaugeSmearingHISQ.cu` ~L155 | Level 2 Lepage 权重 |
| `PhaseField` | INT | `-1` | 否 | `CGaugeSmearingHISQ.cu` ~L157 | Phase 场 ID |
| `Cache` | string | `EHLC_Full` | 否 | `CGaugeSmearingHISQ.cu` ~L160 | 缓存模式。见下表 |

**HISQ 缓存模式枚举值**：

| 枚举值 | 说明 |
|--------|------|
| `EHLC_Full` | 全缓存 |
| `EHLC_Median` | 中等缓存 |
| `EHLC_None` | 不缓存 |

### CGaugeSmearingHISQWithPhaseSU3 — HISQ + Phase（测试用）

**来源文件**：`Code/CLGLib/GaugeSmearing/CGaugeSmearingHISQWithPhase.cpp`

**注意**：此类标记为测试/实验性，构造函数会打印警告。

**读取函数**：`CGaugeSmearingHISQWithPhaseSU3::InitialOtherParameters()` (lines 23-32)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `U1FieldId` | INT | `0` | 否 | `CGaugeSmearingHISQWithPhase.cpp` ~L28 | U(1) Phase 场 ID |
| `Charge` | Real | `1.0` | 否 | `CGaugeSmearingHISQWithPhase.cpp` ~L31 | 电荷参数 |

### CGaugeSmearingAPEProj — APE 投影平滑

**来源文件**：`Code/CLGLib/GaugeSmearing/CGaugeSmearingAPEProj.cu`

**读取函数**：`CGaugeSmearingAPEProj::InitialOtherParameters()` (lines 84-114)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `AlphaLeft` | Real | `1.0` | 否 | `CGaugeSmearingAPEProj.cu` ~L87 | APE 左侧权重参数 |
| `AlphaRight` | Real | `0.1` | 否 | `CGaugeSmearingAPEProj.cu` ~L91 | APE 右侧权重参数 |
| `ProjIterate` | INT | `6` | 否 | `CGaugeSmearingAPEProj.cu` ~L97 | 投影迭代次数（限制 4-254） |
| `Cabibbo` | INT | `0` | 否 | `CGaugeSmearingAPEProj.cu` ~L111 | 是否使用 Cabibbo-Marinari 投影。`1`=启用 |

---

## 与场 / 作用量的依赖关系

| Smearing A | 被谁依赖 | 说明 |
|---|---|---|
| `CGaugeSmearingHISQSU3` | `CFieldFermionHISQSU3` / `...HISQSU3R` / `...HISQWithPhaseSU3`、`CActionFermionHISQCombined` | HISQ 场与作用量直接调用。见 [08b-ks-hisq-fermion.md](08b-ks-hisq-fermion.md)、[09b-fermion-actions.md](09b-fermion-actions.md) |
| `CGaugeSmearingStoutSU3` / `CGaugeSmearingASQTADSU3` / `CGaugeSmearingAPEStout` | `CActionFermionKSImprove` / `CActionFermionKSImproveCombined` | 改进作用量通过 smearing 取 effective gauge。例如 StoutLink Clover Wilson 配 `CGaugeSmearingStoutSU3`。见 [09b-fermion-actions.md](09b-fermion-actions.md) |
| 任意 smearing（用于 HMC） | 建议设置 `Update: 1` | HMC 中规范场演化后必须重新计算 fat links。`Update` 默认 `0` |

---

[< 返回目录](home.md) | [< 上一章：多移求解器](11-multi-shift-solvers.md) | [下一章：Staple 缓存 >](13-staple-cache.md)
