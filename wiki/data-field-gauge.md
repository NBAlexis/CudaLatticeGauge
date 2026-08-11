# Gauge Field 模块

规范场实现：SU(2)、SU(3)、U(1)、Z2，以及 Dirichlet 边界和改进作用量变体。

## 文件清单

### 基类与核心模板

| 文件 | 路径 | 说明 |
|------|------|------|
| `CFieldGauge.h` | `Data/Field/CFieldGauge.h` | 抽象基类 |
| `CFieldGaugeLink.h` | `Data/Field/Gauge/CFieldGaugeLink.h` | 链接型规范场模板基类 |
| `CFieldGaugeKernel.h` | `Data/Field/Gauge/CFieldGaugeKernel.h` | 静态 kernel 接口类 |

### 具体实现

| 文件 | 路径 | 说明 |
|------|------|------|
| `CFieldGaugeSU3_12.h` | `Data/Field/Gauge/CFieldGaugeSU3_12.h` | SU(3) 紧凑表示（12 reals） |
| `CFieldGaugeU1Real.h` | `Data/Field/Gauge/CFieldGaugeU1Real.h` | 实值 U(1) 背景场 |
| `CFieldGaugeZ2.h` | `Data/Field/Gauge/CFieldGaugeZ2.h` | Z2 规范场（占位） |
| `CFieldGaugeLinkDirichlet.h` | `Data/Field/Gauge/CFieldGaugeLinkDirichlet.h` | Dirichlet 边界包装模板 |
| `CFieldGaugeSU3OneLoopImproved.h` | `Data/Field/Gauge/CFieldGaugeSU3OneLoopImproved.h` | 一圈 Symanzik 改进 |
| `CFieldGaugeSU3TreeImproved.h` | `Data/Field/Gauge/CFieldGaugeSU3TreeImproved.h` | 树级 Luscher-Weisz 改进 |

## 类层次

```
CBase
└── CField (abstract)
    └── CFieldGauge (abstract)
        ├── CFieldGaugeLink<deviceGauge, matrixN> (template)
        │   ├── CFieldGaugeU1   : CFieldGaugeLink<CLGComplex, 1>
        │   ├── CFieldGaugeSU2  : CFieldGaugeLink<deviceSU2, 2>
        │   ├── CFieldGaugeSU3  : CFieldGaugeLink<deviceSU3, 3>
        │   ├── CFieldGaugeSU4+ : 条件编译（_CLG_SU4_GAUGE 等）
        │   └── CFieldGaugeLinkD<CFeildG> (Dirichlet wrapper)
        │       └── CFieldGaugeSU3D, CFieldGaugeSU2D, ...
        ├── CFieldGaugeSU3_12
        ├── CFieldGaugeU1Real
        ├── CFieldGaugeSU3OneLoopImproved
        └── CFieldGaugeSU3TreeImproved
```

## CFieldGauge（抽象基类）

**文件**: `Data/Field/CFieldGauge.h`

所有规范场的公共接口：

| 方法 | 说明 |
|------|------|
| `CalculateForceAndStaple(pForce, pStaple, betaOverN)` | 计算 HMC 力 + staple（纯虚） |
| `CalculateOnlyStaple(pStaple)` | 仅计算 staple（纯虚） |
| `MakeRandomGenerator()` | 生成随机规范场（纯虚） |
| `CalculatePlaqutteEnergy(betaOverN)` | plaquette 能量（纯虚） |
| `CalculatePlaqutteEnergyUseClover(betaOverN)` | clover 能量（纯虚） |
| `CalculateKinematicEnergy()` | 动能 `-Tr[P^2]`（纯虚） |
| `ExpMult(a, U)` | `U = exp(a * this) * U`（纯虚） |
| `ElementNormalize()` | 重正交化/归一化（纯虚） |
| `MatrixN()` | 返回规范群维度 N（纯虚） |
| `PolyakovOnSpatialSite(buffer, byDir)` | Polyakov 圈（纯虚） |
| `ApplyStaggeredPhase(iType)` | 应用 staggered 相位（标准或 MILC 约定） |
| `TransformToIA()` / `TA()` / `TransformToU()` | 角动量测试转换 |

## CFieldGaugeLink（模板基类）

**文件**: `Data/Field/Gauge/CFieldGaugeLink.h`

链接型规范场的模板实现。`deviceGauge` 是设备端矩阵类型（`deviceSU2`、`deviceSU3`、`CLGComplex` 等），`matrixN` 是群维度。

**核心成员**：`deviceGauge* m_pDeviceData` — 设备端数据数组，大小为 `linkCount`。

**关键方法**：
- `InitialFieldWithFile(sFileName, eFileType)` — 从文件加载（支持 `EFFT_CLGBin`、`EFFT_CLGBinFloat`、`EFFT_CLGBinDouble`、压缩格式）
- `CalculateForceAndStaple(...)` — 委托给 `CFieldGaugeKernel`
- `CalculateAllStaples(ppDeviceStaple)` — 计算所有 staple
- `CalculateFmunu(pDeviceFmunu)` — 计算场强张量
- `Dagger()` / `AxpyPlus(x)` / `Axpy(a, x)` / `Mul(...)` / `LeftMul(...)` — BLAS 操作
- `ExpMult(a, U)` — 指数乘法
- `ElementNormalize()` — 归一化

**注册宏**：`__DEFINE_GAUGE_LINK(classname, deviceGauge, matrixN, fieldTypeEnum)` 用于一次性声明模板特化和注册。

## 具体规范场类

### CFieldGaugeSU3_12 — 紧凑 SU(3)

**文件**: `Data/Field/Gauge/CFieldGaugeSU3_12.h`

使用 `deviceSU3_12`（12 个实数，只存前两列）代替完整的 `deviceSU3`（18 个实数）。第三列通过叉积重构：
```
col2 = conj(col0 × col1)
```

**用途**：作为力梯度积分器的备份缓冲区（`m_pUPrime`），节省 33% 内存。

**限制**：仅对精确 SU(3) 矩阵有效。非 SU(3) 保持的操作会返回完整的 `deviceSU3`。

**静态方法**：
- `CopySU3ToSU3_12(pDest, pSrc, count)` — 压缩
- `CopySU3_12ToSU3(pDest, pSrc, count)` — 解压
- `SU3_12MSE(pCompact, pRef, linkCount)` — 均方误差

### CFieldGaugeU1Real — 背景电磁场

**文件**: `Data/Field/Gauge/CFieldGaugeU1Real.h`

实值 U(1) 场，用于在格点上施加背景场：
- 化学势（时间方向相位）
- 电场（空间方向相位梯度）
- 磁场（空间平面通量）

通过 `InitialU1Real(eChemicalType, eEType, eBType, ...)` 配置。支持检查切片一致性、计算 `F_{mu,nu}`。

### Dirichlet 边界包装

**文件**: `Data/Field/Gauge/CFieldGaugeLinkDirichlet.h`

模板 `CFieldGaugeLinkD<CFeildG>` 包装任意 `CFieldGaugeLink` 子类，重写所有 force/staple/energy 计算为 `_D` 后缀版本（正确处理边界格点）。

通过宏 `_DEFINE_Gauge_Dirichlet(baseclass, classname)` 快速定义 Dirichlet 变体。

### 改进作用量

**一圈 Symanzik**（`CFieldGaugeSU3OneLoopImproved`）：
- 含矩形项（coefficient `Cr`）和 twisted-loop 项（coefficient `Ct`）
- 系数可直接设置，或通过 `U0` 和 `Nf` 按微扰公式计算（arXiv:1004.0342 Eq. A2）

**树级 Luscher-Weisz**（`CFieldGaugeSU3TreeImproved`）：
- 仅矩形项，系数 `RectOverPlaq`（默认 -0.05）
- 支持 Dirichlet 边界（`CFieldGaugeSU3TreeImprovedD`）

## CFieldGaugeKernel（静态 Kernel 接口）

**文件**: `Data/Field/Gauge/CFieldGaugeKernel.h`

所有 gauge field 设备端操作的静态模板接口。方法分为：

| 类别 | 方法示例 |
|------|---------|
| 周期性 BC force/staple | `CalculateForceAndStaple(...)`, `CalculateOnlyStaple(...)` |
| 周期性 BC energy | `CalculatePlaqutteEnergy(...)`, `CalculatePlaqutteEnergyUseClover(...)` |
| Dirichlet BC | `CalculateForceAndStaple_D(...)`, `CalculatePlaqutteEnergy_D(...)` |
| 改进作用量 | `CalculateForceRectangular(...)`, `CalculateTwistedLoopEnergy(...)` |
| 其他 | `CalculateFmunu(...)`, `CacheKSRotationGaugeBuffer(...)`, `NaikForce(...)` |

显式实例化：`CLGComplex/1`、`deviceSU2/2`、`deviceSU3/3`，以及可选的 `deviceSU4/4` 到 `deviceSU8/8`。

## 设备端数据结构

### deviceSU3（`Tools/Math/SU3.h`）

**存储**：`CLGComplex m_me[9]`（行主序 3x3 复矩阵）

**关键方法**：
- 工厂：`makeSU3Zero()`, `makeSU3Id()`, `makeSU3Random(fatIndex)`, `makeSU3TA(...)`, `makeSU3ContractV(...)`
- 运算：`Add()`, `Sub()`, `Mul()`, `MulDagger()`, `DaggerMul()`, `MulReal()`, `MulComp()`, `Dagger()`
- 向量乘法：`MulVector(v)`, `DagMulVector(v)`, `MulWilsonVector(v)`
- 投影：`Norm()` / `Proj(ite=4)`（HYP 投影）, `CabbiboMarinariProj()`
- 代数：`Ta()`（反厄米特无迹部分）, `Th()`（厄米特无迹部分）
- 指数：`Exp(a, precision)`（Taylor 级数）, `QuickExp(a)`（闭式）, `StrictExp()`（特征值分解）, `StrictExpTA(a)`
- 幂/对数：`Power(fPower)`, `Log()`
- 特征系统：`CalculateEigenValues(c1,c2,c3)`, `EigenVectors(c1,c2,c3)`, `Hessenberg()`, `FrancisQRIteration()`
- 行列式/迹：`Det()`, `Tr()`, `ReTr()`, `ImTr()`

### deviceSU2（`Tools/Math/SU2.h`）

**存储**：`CLGComplex m_me[4]`（行主序 2x2）

方法与 `deviceSU3` 类似，但所有操作针对 2x2 矩阵优化：
- `QuickExp(a)` — 使用 Pauli 矩阵闭式
- `Norm()` — 简单 2x2 归一化
- `CalculateEigenValues(c1,c2,vectors)` — 解析特征值

### deviceSU3_12（`Tools/Math/SU3_12.h`）

**存储**：`CLGComplex m_me[6]`（前两列）

**重构**：第三列通过 `conj(col0 × col1)` 实时计算。

**SU(3) 保持的操作**（直接在 12 元素上）：`Mul()`, `MulDagger()`, `DaggerMul()`, `Dagger()`, `Inverse()`

**非 SU(3) 保持的操作**（返回完整 `deviceSU3`）：`AddC()`, `SubC()`, `TaC()`, `Exp()`, `Power()`, `Log()` 等

### deviceSUN<N, NoE>（`Tools/Math/SUN.h`）

通用 SU(N) 模板结构（N=4..8 通过宏 `_TYPEDEFSUN` 生成 typedef：`deviceSU4` 到 `deviceSU8`）。

提供完整的线性代数：QR（Householder）、LU（高斯）、Hessenberg、Francis QR 迭代、特征值分解、矩阵 log/exp/power。

## 关键设计模式

- **模板 + 显式实例化**：`CFieldGaugeLink<deviceGauge, matrixN>` 是模板，但所有实例化在编译期确定，通过 `__DEFINE_GAUGE_LINK` 宏注册到工厂。
- **Kernel 静态接口**：`CFieldGaugeKernel` 将 host 端的 kernel 调用逻辑与具体的 gauge field 类解耦。
- **包装器模式**：`CFieldGaugeLinkD` 和 `CFieldGaugeOneLoopImproved` 通过继承包装基类，只重写需要修改的方法（force/energy）。
- **设备端结构体**：所有矩阵操作在设备端通过结构体方法完成（注意：在模板 kernel 中应使用 `DeviceInlineTemplate.h` 的全局函数，而非成员方法）。

## 常用场类型枚举

| 枚举值 | 对应类 |
|--------|--------|
| `EFT_GaugeSU2` | `CFieldGaugeSU2` |
| `EFT_GaugeSU3` | `CFieldGaugeSU3` |
| `EFT_GaugeSU3D` | `CFieldGaugeSU3D`（Dirichlet） |
| `EFT_GaugeU1` | `CFieldGaugeU1` |
| `EFT_GaugeZ2` | （占位） |
