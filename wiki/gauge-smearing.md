# Gauge Smearing 模块

规范平滑算法：APE、Stout、ASQTAD、HISQ 等 fat-link 构造，以及力导数计算。

## 文件清单

### 基类

| 文件 | 路径 | 说明 |
|------|------|------|
| `CGaugeSmearing.h` | `GaugeSmearing/CGaugeSmearing.h` | 平滑算法抽象基类 |

### 具体实现

| 文件 | 路径 | 说明 |
|------|------|------|
| `CGaugeSmearingAPEProj.h` | `GaugeSmearing/` | APE 平滑 + 投影 |
| `CGaugeSmearingAPEStout.h` | `GaugeSmearing/` | APE + Stout 混合 |
| `CGaugeSmearingASQTAD.h` | `GaugeSmearing/` | ASQTAD fat-link（模板） |
| `CGaugeSmearingHISQ.h` | `GaugeSmearing/` | HISQ 两级平滑（模板） |
| `CGaugeSmearingHISQWithPhase.h` | `GaugeSmearing/` | HISQ + U(1) 相位（已废弃） |
| `CGaugeSmearingStoutSU3.h` | `GaugeSmearing/` | Stout 平滑（SU3 专用） |

## CGaugeSmearing（抽象基类）

**文件**: `GaugeSmearing/CGaugeSmearing.h`

所有平滑算法的统一接口：

| 方法 | 说明 |
|------|------|
| `Initial(pOwner, params)` | 初始化 |
| `GaugeSmearing(pGauge, pOrignal, pStaple, bProject)` | 核心平滑操作（纯虚） |
| `GaugeSmearingC(pGauge)` | HMC 更新路径包装：缓存有效规范场 |
| `DerivateOnU(...)` | 力导数计算（默认报错，子类覆盖） |
| `GetEffectiveGauge()` | 获取平滑后的有效规范场 |
| `GetEffectiveGaugeLevel1()` | 获取一级 fat-link（HISQ 用） |
| `GetNaikLink()` | 获取 Naik link（HISQ 用） |
| `CalculateSpatialFatLink(...)` | 计算不含 T 方向的 fat-link |

**关键成员**：
- `m_pEffecitveGauge` — 缓存的有效规范场（HMC 中使用）
- `m_bCalledWhenUpdate` — 是否在更新阶段调用
- `m_uiIterate` — 迭代次数

## APE 平滑

### APEProj

**文件**: `GaugeSmearing/CGaugeSmearingAPEProj.h`

标准 APE 平滑：加权平均原链接和 staple，然后投影回规范群。

```cpp
U' = Proj( (1-alpha) * U + alpha * Staple )
```

**参数**：
- `m_fAlphaLeft` / `m_fAlphaRight` — 左右权重
- `m_bCMProj` — 是否使用 Cabibbo-Marinari 投影
- `m_byProjIterate` — 投影迭代次数

### APEStout

**文件**: `GaugeSmearing/CGaugeSmearingAPEStout.h`

APE + Stout 混合：先用 APE 构造 staples，再用 Stout 的指数映射形式更新。

**参数**：`m_fRho` — staple 权重。

## ASQTAD

**文件**: `GaugeSmearing/CGaugeSmearingASQTAD.h`

模板类 `CGaugeSmearingASQTAD<gaugetype, matrixN>`：

ASQTAD fat-link 构造，包含最多 7-link staples 和可选 Lepage 项：

```
U_fat = c1 * U + c3 * fat3 + c5 * fat5 + c7 * fat7 + cLepage * Lepage
```

**投影回 SU(3)**：两种方法
1. **Cayley-Hamilton 解析平方根**：`UR = U (U^+ U)^{-1/2}`
2. **有理逼近**：对 `(U^+ U)^{-1/2}` 的有理函数逼近

**可选 det-投影**：`UR' = UR / det[U]^{1/3}` 确保严格 SU(3)。

**参数**：
- `m_fOriginal` / `m_fFat3` / `m_fFat5` / `m_fFat7` / `m_fLepage` — 各项系数
- `m_bProj` / `m_bProjDet` — 投影开关
- `m_bUseCaylayHamilton` — 使用 Cayley-Hamilton 还是有理逼近

**静态方法**：
- `Fat357Lepage()` — 构造所有 fat-link 项
- `SmearingForce()` — 计算平滑对原始规范场的力
- `ProjectCaylayHamilton()` / `ProjectRationalApproximation()` — 两种投影
- 对应的 `Force` 变体用于 HMC 力链式法则

**注册**：`CGaugeSmearingASQTADSU3`（`EFT_GaugeSU3`）。

## HISQ

**文件**: `GaugeSmearing/CGaugeSmearingHISQ.h`

模板类 `CGaugeSmearingHISQ<gaugetype, matrixN>`：

两级 HISQ 平滑：
- **Level 1**：ASQTAD-like fat links（无 Lepage）
- **Level 2**：含 Lepage 和 Naik 项

```
Level 1: one-link + fat3 + fat5 + fat7
Level 2: one-link + fat3 + fat5 + fat7 + Lepage + Naik
```

**系数默认值**（参考 hep-lat/0710.0737）：
- Level 2 one-link: `1 + epsi/8`
- Level 2 Lepage: `-1/8`
- Level 2 Naik: `-(1+epsi)/24`

**缓存模式**（`EHISQLinkCache`）：
- `EHLC_Full` — 缓存所有中间结果（最快，最耗内存）
- `EHLC_Median` — 缓存部分中间结果
- `EHLC_None` — 不缓存（最慢，最省内存）

**Naik link**：三级跳链接，用于改进 staggered 费米子色散关系。

**U(1) 相位**：支持静态 U(1) 相位场（如化学势、边界相位），缓存 `m_pNaikLinkPhase`。

**注册**：`CGaugeSmearingHISQSU3`（`EFT_GaugeSU3`）。

## Stout

**文件**: `GaugeSmearing/CGaugeSmearingStoutSU3.h`

Stout 平滑（hep-lat/0311018），SU(3) 专用：

```
U' = exp(iQ) * U
Q = -i/2 (Ω - Ω^+),  Ω = ρ * Staple
```

**缓存模式**（`EStoutLinkCache`）：
- `ESLC_Full` — 全缓存
- `ESLC_Median` — 部分缓存
- `ESLC_Small` — 最小缓存

**力导数**：利用 `exp(iQ)` 的导数闭式（特征值分解），需要预计算特征系统相关参数。

**参数**：
- `m_fRhojk` — 空间平面 staple 权重
- `m_fRho4mu` — 时间方向 staple 权重

## 已废弃

### HISQWithPhase

**文件**: `GaugeSmearing/CGaugeSmearingHISQWithPhase.h`

扩展 HISQ 以支持外部 U(1) 相位场。**已废弃，仅供测试。**

## 类层次

```
CGaugeSmearing
├── CGaugeSmearingAPEProj
├── CGaugeSmearingAPEStout
├── CGaugeSmearingASQTAD<gaugetype, matrixN>
│   └── CGaugeSmearingASQTADSU3
├── CGaugeSmearingHISQ<gaugetype, matrixN>
│   ├── CGaugeSmearingHISQSU3
│   └── CGaugeSmearingHISQWithPhase<gaugetype, matrixN> (废弃)
│       └── CGaugeSmearingHISQWithPhaseSU3
└── CGaugeSmearingStoutSU3
```

## HMC 中的使用

1. **测量路径**：`GaugeSmearing(pGauge, pOrignal, pStaple)` 直接平滑
2. **更新路径**：`GaugeSmearingC(pGauge)` 缓存有效规范场到 `m_pEffecitveGauge`
3. **力计算**：`DerivateOnU(pEffectiveGauge, pOrignalGauge, pNaikForce, pf0)` 计算链式法则导数

## 关键设计模式

- **模板化**：ASQTAD 和 HISQ 是模板类，支持多规范群（目前仅 SU3 实例化）。
- **缓存策略**：通过枚举选择内存-速度权衡，HISQ 支持 Full/Median/None 三种模式。
- **投影方法**：ASQTAD/HISQ 支持 Cayley-Hamilton 和有理逼近两种回投影到 SU(3) 的方法。
- **力链式法则**：平滑改变规范场，力计算需要 `dS/dU_original = dS/dU_effective * dU_effective/dU_original`。
