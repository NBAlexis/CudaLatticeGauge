> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 14. 规范固定 (Gauge Fixing)

**来源文件**：`Code/CLGLib/Core/CLGLibManager.cpp`

**读取函数**：`CCLGLibManager::CreateGaugeFixing()` (lines 890-899)

---

## 配置参数

YAML 键名为 `GaugeFixing`。

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Name` | string | `CGaugeFixingLandauCornell` | 否 | `CLGLibManager.cpp` ~L894 | 规范固定类名。通过 `appCreate()` 实例化 |

## 支持的规范固定类

| 类名 | 说明 |
|------|------|
| `CGaugeFixingLandauCornell` | Landau 规范固定（Cornell 算法，FFT 加速） |
| `CGaugeFixingCoulombCornell` | Coulomb 规范固定（Cornell 算法，FFT 加速） |
| `CGaugeFixingLandauLosAlamos` | Landau 规范固定（Los Alamos 迭代算法） |
| `CGaugeFixingCoulombLosAlamos` | Coulomb 规范固定（Los Alamos 迭代算法） |
| `CGaugeFixingMAG` | 最大阿贝尔规范 (Maximal Abelian Gauge) |
| `CGaugeFixingMCGDirect` | 最大中心规范 (Maximal Center Gauge) — 直接算法 |
| `CGaugeFixingMCGIndirect` | 最大中心规范 — 间接算法（先 MAG 再 MCG） |
| `CGaugeFixingRandom` | 随机规范变换（无参数） |

## Landau/Coulomb Cornell 算法参数

**来源文件**：`Code/CLGLib/GaugeFixing/CGaugeFixingLandauCornell.cu` / `CGaugeFixingCoulombCornell.cu`

**读取函数**：`CGaugeFixingLandauCornell::InitialOtherParameters()` (lines 343-381) / `CGaugeFixingCoulombCornell::InitialOtherParameters()` (lines 437-476)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Alpha` | DOUBLE/Real | `0.08` | 否 | `CGaugeFixing*Cornell.cu` ~L349 | 迭代步长（松弛参数）。先尝试 `DOUBLE`，失败再尝试 `Real` |
| `Accuracy` | Real | `0.00000000001` | 否 | `CGaugeFixing*Cornell.cu` ~L357 | 收敛精度 |
| `MaxIterate` | INT | `1000000` | 否 | `CGaugeFixing*Cornell.cu` ~L368 | 最大迭代次数 |
| `ShowErrorStep` | INT | `1000` | 否 | `CGaugeFixing*Cornell.cu` ~L375 | 每隔多少步输出一次误差 |
| `FFT` | INT | `1` | 否 | `CGaugeFixing*Cornell.cu` ~L381 | 是否使用 FFT 加速。`1`=使用 FFT |

## Landau/Coulomb Los Alamos 算法参数

**来源文件**：`Code/CLGLib/GaugeFixing/CGaugeFixingLandauLosAlamos.cu` / `CGaugeFixingCoulombLosAlamos.cu`

**读取函数**：`CGaugeFixingLandauLosAlamos::InitialOtherParameters()` (lines 405-437) / `CGaugeFixingCoulombLosAlamos::InitialOtherParameters()` (lines 521-560)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Omega` | Real | `1.0` | 否 | `CGaugeFixing*LosAlamos.cu` ~L410 | 超松弛参数 (SOR) |
| `Accuracy` | Real | `0.00000000001` | 否 | `CGaugeFixing*LosAlamos.cu` ~L415 | 收敛精度 |
| `MaxIterate` | INT | `1000000` | 否 | `CGaugeFixing*LosAlamos.cu` ~L426 | 最大迭代次数 |
| `CheckErrorStep` | INT | `1000` | 否 | `CGaugeFixing*LosAlamos.cu` ~L433 | 每隔多少步检查一次误差 |

## MAG (最大阿贝尔规范) 参数

**来源文件**：`Code/CLGLib/GaugeFixing/CGaugeFixingMAG.cu`

**读取函数**：`CGaugeFixingMAG::InitialOtherParameters()` (lines 1098-1130)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Omega` | Real | `1.5` | 否 | `CGaugeFixingMAG.cu` ~L1102 | 超松弛参数 |
| `Accuracy` | Real | `0.00000000001` | 否 | `CGaugeFixingMAG.cu` ~L1107 | 收敛精度 |
| `MaxIterate` | INT | `1000000` | 否 | `CGaugeFixingMAG.cu` ~L1118 | 最大迭代次数 |
| `CheckErrorStep` | INT | `1000` | 否 | `CGaugeFixingMAG.cu` ~L1125 | 每隔多少步检查一次误差 |

## MCG Direct (直接最大中心规范) 参数

**来源文件**：`Code/CLGLib/GaugeFixing/CGaugeFixingMCGDirect.cu`

**读取函数**：`CGaugeFixingMCGDirect::InitialOtherParameters()` (lines 829-861)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Omega` | Real | `1.5` | 否 | `CGaugeFixingMCGDirect.cu` ~L833 | 超松弛参数 |
| `Accuracy` | Real | `0.00000000001` | 否 | `CGaugeFixingMCGDirect.cu` ~L838 | 收敛精度 |
| `MaxIterate` | INT | `1000000` | 否 | `CGaugeFixingMCGDirect.cu` ~L849 | 最大迭代次数 |
| `CheckErrorStep` | INT | `1000` | 否 | `CGaugeFixingMCGDirect.cu` ~L856 | 每隔多少步检查一次误差 |

## MCG Indirect (间接最大中心规范) 参数

**来源文件**：`Code/CLGLib/GaugeFixing/CGaugeFixingMCGIndirect.cu`

**读取函数**：`CGaugeFixingMCGIndirect::InitialOtherParameters()` (lines 38-87)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Omega` | Real | `1.5` | 否 | `CGaugeFixingMCGIndirect.cu` ~L42 | MCG 阶段超松弛参数 |
| `OmegaStage1` | Real | `1.5` | 否 | `CGaugeFixingMCGIndirect.cu` ~L47 | 第一阶段（MAG）超松弛参数 |
| `OmegaMAG` | Real | `1.5` | 否 | `CGaugeFixingMCGIndirect.cu` ~L50 | `OmegaStage1` 的别名 |
| `Accuracy` | Real | `0.00000000001` | 否 | `CGaugeFixingMCGIndirect.cu` ~L57 | 收敛精度 |
| `MaxIterate` | INT | `1000000` | 否 | `CGaugeFixingMCGIndirect.cu` ~L68 | MCG 阶段最大迭代次数 |
| `Stage1MaxIterate` | INT | `100000` | 否 | `CGaugeFixingMCGIndirect.cu` ~L75 | 第一阶段最大迭代次数 |
| `MAGMaxIterate` | INT | `100000` | 否 | `CGaugeFixingMCGIndirect.cu` ~L78 | `Stage1MaxIterate` 的别名 |
| `CheckErrorStep` | INT | `1000` | 否 | `CGaugeFixingMCGIndirect.cu` ~L86 | 每隔多少步检查一次误差 |

**注意**：间接 MCG 先执行 MAG 固定，再执行 MCG 固定。`OmegaStage1`/`OmegaMAG` 和 `Stage1MaxIterate`/`MAGMaxIterate` 是同一参数的别名。


---

## 补充参数（代码中存在但文档中缺失）

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `IMCGGrid` | INT | `24` | 否 | Code/CLGLib/GaugeFixing/CGaugeFixingMCGIndirect.cu ~L535 | 间接 MCG 的网格数，最小值为 `4` |
| `Mixed` | INT | `0` | 否 | Code/CLGLib/GaugeFixing/CGaugeFixingCoulombLosAlamos.cu ~L575 | 是否使用混合精度更新。`1`=启用 |
| `UseStandardIMCG` | INT | `0` | 否 | Code/CLGLib/GaugeFixing/CGaugeFixingMCGIndirect.cu ~L528 | 是否使用标准间接 MCG。`1`=使用标准算法，`0`=使用 MAG 预处理的直接 MCG |

---

[< 返回目录](home.md) | [< 上一章：Staple 缓存](13-staple-cache.md) | [下一章：更新器 >](15-updators.md)
