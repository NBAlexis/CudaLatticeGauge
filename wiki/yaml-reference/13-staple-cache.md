> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 13. Staple 缓存 (Staple Cache)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateGaugeStapleCache()` (lines 879-888)

---

## 配置参数

YAML 键名为 `StapleCache`（第一个）或 `StapleCache2`, `StapleCache3`, ...。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `StapleCacheName` | string | `CStapleCacheSU3` | 否 | `CLGLibManager.cpp` ~L883 | Staple 缓存类名。通过 `appCreate()` 实例化 |

## Staple 缓存基类参数 (CStapleCache)

**来源文件**：`Code/CLGLib/Update/CStapleCache.cpp`

**读取函数**：`CStapleCache::Initial()` (lines 16-50)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `FieldId` | INT | `1` | 否 | `CStapleCache.cpp` ~L19 | 缓存的规范场 ID |
| `CacheStaples` | INT | `0` | 否 | `CStapleCache.cpp` ~L24 | 是否缓存 staples。`1`=缓存 |
| `CachePlaquttes` | INT | `0` | 否 | `CStapleCache.cpp` ~L28 | 是否缓存 plaquettes。`1`=缓存 |
| `CacheFmunu` | INT | `1` | 否 | `CStapleCache.cpp` ~L32 | 是否缓存 F_munu。`1`=缓存 |
| `CacheRotationLink` | INT | `0` | 否 | `CStapleCache.cpp` ~L36 | 是否缓存旋转链接。`1`=缓存 |
| `FermionUpdate` | INT | `1` | 否 | `CStapleCache.cpp` ~L40 | 费米子场更新时是否刷新缓存。`1`=刷新 |
| `GaugeUpdate` | INT | `0` | 否 | `CStapleCache.cpp` ~L44 | 规范场更新时是否刷新缓存。`1`=刷新 |
| `UseEffectiveGauge` | INT | `1` | 否 | `CStapleCache.cpp` ~L48 | 是否使用有效规范场。`1`=使用 |

## 支持的 Staple 缓存类

| 类名 | 说明 |
|------|------|
| `CStapleCacheSU3` | SU(3) Staple 缓存（无额外参数） |

---

## 与场 / 作用量的依赖关系

| 若使用 A | 则必须同时配置 B | 说明 |
|---|---|---|
| `CFieldFermionWilsonSquareCloverSU3` 及其变体 | `StapleCache: CStapleCacheSU3` | Clover 核函数需要 staple cache 提供的 `Fmunu`，否则 `_FAIL_EXIT`。见 [08a-wilson-fermion.md](08a-wilson-fermion.md) |
| `CFieldFermionKSSU3R` 等旋转 KS 场（`CachedGauge: 1`） | `StapleCache: CStapleCacheSU3` | 旋转矩阵缓存取自 staple cache。见 [08b-ks-hisq-fermion.md](08b-ks-hisq-fermion.md) |

> **提示**：`CIntegrator` 在 HMC 中会自动使用已注册的 staple cache（如果存在），所以只要有 Clover / 旋转 KS 场，就必须显式写出 `StapleCache` 块。

---

[< 返回目录](home.md) | [< 上一章：规范平滑](12-gauge-smearing.md) | [下一章：规范固定 >](14-gauge-fixing.md)
