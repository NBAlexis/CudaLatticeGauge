> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 03. 随机数生成器

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::InitialRandom()` (lines 435-448)

---

## 配置参数

随机数生成器的参数在 [02-lattice.md](02-lattice.md) 的格点常数阶段读取：

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `RandomType` | string | `ER_Schrage` | 否 | `CLGLibManager.cpp` ~L146 | 随机数生成器类名 |
| `RandomSeed` | INT | `81192` | 否 | `CLGLibManager.cpp` ~L144 | 随机数种子 |
| `RandomSeedType` | string | （无） | 否 | `CLGLibManager.cpp` ~L148 | 若设为 `ERST_Timestamp`，用当前时间覆盖种子 |

## 支持的随机数类型

通过 `appCreate(sRandomType)` 工厂实例化。需确认 `Code/CLGLib/Core/CLGSetup.h` 中对应宏已启用。

| 类名 | 宏开关 | 说明 |
|------|--------|------|
| `CRandomSchrage` | （默认启用） | Schrage 线性同余生成器 |
| `CRandomMTE` | `_CLG_MT19937` | Mersenne Twister MT19937 |
| `CRandomXORWOW` | `_CLG_XORWOW` | NVIDIA CUDA XORWOW |
| `CRandomPHILOX` | `_CLG_PHILOX` | Philox 计数器生成器 |
| `CRandomMRG32K3A` | `_CLG_MRG32K3A` | MRG32k3a 组合生成器 |

## 初始化行为

1. 若 `RandomSeedType == ERST_Timestamp`，用 `time(NULL)` 覆盖 `RandomSeed`
2. 调用 `appCreate(RandomType)` 创建生成器实例
3. 调用 `CRandom::Initial(seed)` 初始化状态
4. 将随机数状态复制到设备端

---

[< 返回目录](home.md) | [< 上一章：格点常数](02-lattice.md) | [下一章：索引与边界 >](04-index-boundary.md)
