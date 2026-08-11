# Math Tools 模块

设备端数学函数、SU(N) 群结构、向量、Gamma 矩阵、FFT、随机数生成。

## 文件清单

| 文件 | 路径 | 说明 |
|------|------|------|
| `DeviceInlineTemplate.h` | `Tools/Math/DeviceInlineTemplate.h` | 设备端模板数学函数（`_add`、`_mul`、`_dagger` 等） |
| `DeviceInlineTemplate2.h` | `Tools/Math/DeviceInlineTemplate2.h` | 寄存器优化实验（当前未启用） |
| `SU2.h` | `Tools/Math/SU2.h` | SU(2) 设备端结构 `deviceSU2` |
| `SU3.h` | `Tools/Math/SU3.h` | SU(3) 设备端结构 `deviceSU3` |
| `SU3_12.h` | `Tools/Math/SU3_12.h` | 紧凑 SU(3) 表示 `deviceSU3_12` |
| `SUN.h` | `Tools/Math/SUN.h` | 通用 SU(N) 模板 `deviceSUN` |
| `Vectors.h` | `Tools/Math/Vectors.h` | `deviceSU2Vector`、`deviceSU3Vector`、`deviceWilsonVectorSU3` |
| `VectorsN.h` / `VectorsN2.h` | `Tools/Math/` | 通用 N 维向量模板 |
| `GammaMatrix.h` | `Tools/Math/GammaMatrix.h` | Gamma 矩阵定义与运算 |
| `CudaComplexFunction.h` | `Tools/Math/CudaComplexFunction.h` | 复数辅助函数（幂、对数、平方根） |
| `Random.h` | `Tools/Math/Random.h` | 设备端随机数生成器 `CRandom` |
| `CLGFFT.h` | `Tools/Math/CLGFFT.h` | FFT 包装器 |
| `CLinearAlgebraHelper.h` | `Tools/Math/CLinearAlgebraHelper.h` | Host 端线性代数辅助 |
| `CRationalApproximation.h` | `Tools/Math/CRationalApproximation.h` | 有理逼近（用于 RHMC） |

## DeviceInlineTemplate.h — 模板数学函数

**这是模板 kernel 开发的核心文件。** 所有设备端模板运算通过全局函数实现，支持 `Real`、`CLGComplex`、`deviceSU2`、`deviceSU3`、`deviceSUN`、`deviceSU2Vector`、`deviceSU3Vector`、`deviceWilsonVectorSU3` 等类型。

### 工厂函数

| 函数 | 说明 |
|------|------|
| `_makeId<T>()` | 创建单位元（1、单位矩阵、单位向量） |
| `_makeZero<T>()` | 创建零元 |
| `_makeRandom<T>(fatIdx)` | 均匀随机 |
| `_makeGaussian<T>(fatIdx)` | 高斯随机（用于生成规范场 momentum） |
| `_makeZ4<T>(fatIdx)` | Z4 随机（`±1, ±i`） |
| `_makeSumGenerator<T>(factor)` | 生成元和的倍数 |
| `_makeContract<TMatrix, TVector>(left, right)` | 向量外积（返回矩阵） |
| `_makeColorVector<T>(spin, color)` | 单自旋-颜色基向量 |

### 算术运算

| 函数 | 说明 |
|------|------|
| `_add(left, right)` / `_addC(left, right)` | 加法（in-place / const） |
| `_sub(left, right)` / `_subC(left, right)` | 减法 |
| `_mul(a, b)` / `_mulC(a, b)` | 乘法（矩阵×矩阵、矩阵×标量、向量×标量） |
| `_muldag(a, b)` | `a = a * b†` |
| `_dagmul(a, b)` | `a = a† * b` |
| `_dagger(x)` / `_daggerC(x)` | 厄米共轭（in-place / const） |
| `_oppo(x)` / `_oppoC(x)` | 取反 |
| `_rcpC(x)` | 倒数/逆 |

**关键规则**：在模板设备代码中，必须使用这些全局函数，不能调用成员方法（如 `x.Add(y)`）。这是为了支持模板特化，让同一套 kernel 代码对 `deviceSU2`、`deviceSU3`、`deviceSU4` 等都有效。

## SU(N) 设备端结构

### deviceSU2 / deviceSU3 / deviceSUN

详见 [data-field-gauge.md](data-field-gauge.md) 中的"设备端数据结构"章节。

核心要点：
- `deviceSU2`：2x2 复矩阵，4 个 `CLGComplex`
- `deviceSU3`：3x3 复矩阵，9 个 `CLGComplex`
- `deviceSU3_12`：紧凑表示，6 个 `CLGComplex`（前两列，第三列叉积重构）
- `deviceSUN<N, NoE>`：通用 NxN 模板，N=4..8 通过 `_TYPEDEFSUN` 宏生成 typedef

### 向量结构

**deviceSU3Vector**（`Vectors.h`）：
- 存储：`CLGComplex m_ve[3]`（若 `_CLG_PADDING` 则为 `[4]`）
- 运算：`Add`, `Sub`, `MulReal`, `MulComp`, `MulZ4`
- 内积：`ConjugateDotC(other)` — `Σ v_i* · w_i`
- 归一化：`Norm()` — 除以模长

**deviceWilsonVectorSU3**（`Vectors.h`）：
- 存储：4 个 `deviceSU3Vector`（4 个 Dirac 自旋分量 × 3 个颜色分量）
- 总大小：24 个复数（单精度 192 字节，双精度 384 字节）
- 用途：Wilson 费米子场的基本单元

**deviceSU2Vector**：
- 类似 `deviceSU3Vector`，但 2 个分量

## Gamma 矩阵

**文件**: `Tools/Math/GammaMatrix.h`

Gamma 矩阵有两种实现（由 `__GAMMA_IMPLEMENTATION_SWITCH` 控制）：

| 模式 | 存储 | 特点 |
|------|------|------|
| 0（默认） | `BYTE m_byZ4[4]`, `BYTE m_uiIndex[4]` | 单精度快 8%，行置换 + Z4 相位 |
| 1 | `UINT m_uiValue`（16 bit） | 双精度可能更优，紧凑编码 |

**Z4 相位**：`i^b`，`b=0,1,2,3` → `1, i, -1, -i`

**可用矩阵**（`EGammaMatrix` 枚举）：
- `UNITY`, `GAMMA1`~`GAMMA5`
- `GAMMA51`~`GAMMA54`, `GAMMA15`~`GAMMA45`
- `SIGMA12`, `SIGMA23`, `SIGMA31`（Minkowski sigma）
- `SIGMA12E`, `SIGMA23E`, `SIGMA31E`（Euclidean sigma）
- `SIGMA41`, `SIGMA42`, `SIGMA43`
- `CHARGECONJG`

**初始化**：`gammaMatrixSet::CreateGammaMatrix(eSet, gmarray)` 支持 Dirac 基和 Chiral 基。

**乘法**：`gammaMatrix::MulWilsonC(spinor)` — 返回 `deviceWilsonVectorSU3`。

## 复数辅助函数

**文件**: `Tools/Math/CudaComplexFunction.h`

| 函数 | 说明 |
|------|------|
| `__cuCargf(c)` | 复数幅角 `arg(c)` |
| `__cuCabsSqf(c)` | `|c|^2` |
| `__cuCpowerf(c, p)` | 复数幂 `c^p` |
| `__cuCexpf(c)` | 复指数 `exp(c)` |
| `__cuClogf(c)` | 复对数 `log(c)`，幅角归一化到 `(-π, π]` |
| `__cuCsqrtf(c)` | 复平方根 |
| `cuCaddf_cr(x, y)` | 复数 + 实数 |
| `cuCdivf_cr(x, y)` | 复数 / 实数 |
| `cuCmulf_cr(x, y)` | 复数 × 实数 |

## 随机数生成

**文件**: `Tools/Math/Random.h`

`CRandom` 类封装了多种 CUDA 随机数生成器：

| 类型 | 枚举 | 说明 |
|------|------|------|
| Schrage | `ER_Schrage` | 默认，LCG：`x = 1664525*x + 1013904223`，每个 thread 独立种子表 |
| XORWOW | `ER_XORWOW` | cuRAND 默认伪随机 |
| MRG32K3A | `ER_MRG32K3A` | 条件编译（`_CLG_USE_MRG32K3A`） |
| Philox4 | `ER_PHILOX4_32_10` | 条件编译（`_CLG_USE_PHILOX4`） |
| Scrambled Sobol32 | `ER_SCRAMBLED_SOBOL32` | 准随机，DCU 默认关闭 |

**设备端调用**（通过全局函数）：

```cpp
_deviceRandomF(fatIndex);        // [0, 1) 均匀分布
_deviceRandomGaussF(fatIndex);   // 高斯分布（标准差 1/sqrt(2)）
_deviceRandomGaussFSqrt2(fatIndex); // 高斯分布（标准差 1/2）
_deviceRandomC(fatIndex);        // 复数均匀 [-1,1]×[-1,1]
_deviceRandomGaussC(fatIndex);   // 复高斯
_deviceRandomZ4(fatIndex);       // Z4: {1, i, -1, -i}
```

**Box-Muller 方法**：高斯随机数通过 `sqrt(-2 ln u1) * cos(2π u2)` 生成。

## FFT

**文件**: `Tools/Math/CLGFFT.h`

封装了 cuFFT（CUDA）或 rocFFT（DCU）的 FFT 操作，用于 Cornell 规范固定等需要动量空间运算的场景。

## 有理逼近

**文件**: `Tools/Math/CRationalApproximation.h`

用于 RHMC（Rational Hybrid Monte Carlo）算法，对费米子行列式的有理函数逼近。

## 关键设计模式

- **模板全局函数优先**：所有跨类型的设备端运算通过 `DeviceInlineTemplate.h` 的全局模板函数实现，确保 kernel 代码可复用于多种规范群。
- **C 风格结构体**：设备端结构体（`deviceSU3` 等）使用 C 风格方法（在结构体内定义），但模板 kernel 中不直接调用它们，而是通过 `DeviceInlineTemplate.h` 的包装。
- **精度透明**：所有数学函数通过 `Real`/`CLGComplex` 宏自动适配单/双精度。
- **Gamma 矩阵紧凑编码**：利用 Gamma 矩阵的稀疏性（每行只有一个非零元），只用 `BYTE` 数组存储行置换索引和 Z4 相位，大幅减少寄存器占用。
