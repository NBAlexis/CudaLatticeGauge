# Gauge Fixing 模块

规范固定算法：Landau、Coulomb、MAG、MCG（直接/间接），以及随机规范变换。

## 文件清单

| 文件 | 路径 | 说明 |
|------|------|------|
| `CGaugeFixing.h` | `GaugeFixing/CGaugeFixing.h` | 抽象基类 |
| `CGaugeFixingLandauCornell.h/cu` | `GaugeFixing/` | Landau 规范 — FFT 加速 (Cornell) |
| `CGaugeFixingLandauLosAlamos.h/cu` | `GaugeFixing/` | Landau 规范 — 局部 overrelaxation (Los Alamos) |
| `CGaugeFixingCoulombCornell.h/cu` | `GaugeFixing/` | Coulomb 规范 — FFT 加速 (Cornell) |
| `CGaugeFixingCoulombLosAlamos.h/cu` | `GaugeFixing/` | Coulomb 规范 — 局部 overrelaxation (Los Alamos) |
| `CGaugeFixingMAG.h/cu` | `GaugeFixing/` | 最大阿贝尔规范 (MAG) |
| `CGaugeFixingMCGDirect.h/cu` | `GaugeFixing/` | 直接最大中心规范 (MCG) |
| `CGaugeFixingMCGIndirect.h/cu` | `GaugeFixing/` | 间接最大中心规范 (MAG + MCG) |
| `CGaugeFixingRandom.h/cu` | `GaugeFixing/` | 随机规范变换（测试用） |

## 类层次

```
CBase
└── CGaugeFixing (abstract)
    ├── CGaugeFixingLandauCornell
    ├── CGaugeFixingLandauLosAlamos
    ├── CGaugeFixingCoulombCornell
    ├── CGaugeFixingCoulombLosAlamos
    ├── CGaugeFixingMAG
    ├── CGaugeFixingMCGDirect
    ├── CGaugeFixingMCGIndirect
    └── CGaugeFixingRandom
```

## CGaugeFixing（抽象基类）

**文件**: `GaugeFixing/CGaugeFixing.h`

所有规范固定算法的统一接口：

| 方法 | 说明 |
|------|------|
| `Initial(pOwner, params)` | 从 YAML 参数初始化 |
| `GaugeFixing(pResGauge)` | 执行规范固定（修改传入的场副本） |
| `CheckRes(pGauge)` | 返回残差/目标函数值，用于判断是否收敛 |
| `GetInfos(tab)` | 返回人类可读的状态信息 |

**关键成员**：
- `m_fAccuracy` — 收敛精度（单精度默认 `1e-5`，双精度默认 `1e-11`）
- `m_iMaxIterate` — 最大迭代次数（默认 1,000,000）
- `m_iShowErrorStep` — 每隔多少次迭代打印残差（默认 1000）

## 算法分类

### 1. Landau / Coulomb 规范

**目标**：最大化 `sum_{x,mu} Re Tr[U_mu(x)]`（Landau）或 `sum_{x,i} Re Tr[U_i(x)]`（Coulomb，仅空间方向）。

| 实现 | 方法 | 特点 |
|------|------|------|
| **Cornell** | FFT 加速的 Steepest Descent | 全局更新，收敛快但需要 FFT；`m_fAlpha` 为步长 |
| **Los Alamos** | 局部 overrelaxation | 逐点更新，`m_fOmega` 为过松弛参数；备注：对数定义下似乎无法收敛 |

**Cornell 实现** 使用双精度设备缓冲区存储矩阵的独立分量（`m_pA11`, `m_pA12`, ... `m_pGamma23`），通过 FFT 在动量空间求解变换矩阵 `G(x)`。

**Los Alamos 实现** 使用红黑 checkerboard 并行扫描，每个格点独立求解最优局部变换后应用 overrelaxation。

### 2. MAG（Maximal Abelian Gauge）

**文件**: `GaugeFixing/CGaugeFixingMAG.h`

最大化阿贝尔子群的可见度。对 SU(2) 是最大化 `sigma_3` 分量；对 SU(3) 是使非对角元最小化。

| 规范群 | 局部更新方法 | 投影方法 |
|--------|-------------|---------|
| SU(2) | Cea & Cosmai 算法（hep-lat/9504008） | `MaximalAbelianProjection()` — 提取对角相位 |
| SU(3) | Cabibbo-Marinari 子群分解（hep-lat/0110165） | 同上 — 保留 2x2 对角块 |

**关键参数**：
- `Omega` — overrelaxation 参数（默认 1.5）
- `CheckErrorStep` — 检查收敛的步频（默认 1000）

**收敛判据**：`sum of squared off-diagonal elements < accuracy`

**注意**：MAG 是一种规范变换，不改变物理可观测量。`MaximalAbelianProjection` 是投影操作（破坏规范不变性），用于测试。

### 3. Direct MCG（Maximal Center Gauge）

**文件**: `GaugeFixing/CGaugeFixingMCGDirect.h`

直接最大化中心可见度：
```
R = (1/(N_site * N_dim * N^2)) * sum |Tr U_mu(x)|^2
```

| 规范群 | 局部更新 | 投影 |
|--------|---------|------|
| SU(2) | 构造 trace-weighted staple，投影到 SU(2) | `CenterProjection()` → sign(Tr U) * I (Z_2) |
| SU(3) | Cabibbo-Marinari-Okawa 子群分解 | `CenterProjection()` → 最近 Z_3 元素 |

**关键细节**：
- 使用红黑 checkerboard 并行 sweep
- Overrelaxation：`G_omega = (1-omega)*I + omega*G`，然后重新投影到规范群
- **SU(2) 数值稳定性**：投影后需显式强制 SU(2) 结构（`g.m_me[2] = -g01*`, `g.m_me[3] = g00*`），否则浮点累积会导致偏离 SU(2) 流形
- SU(2) 对 Omega 敏感：Omega=0.1 可能导致发散，Omega=1.0 通常稳定

### 4. Indirect MCG

**文件**: `GaugeFixing/CGaugeFixingMCGIndirect.h`

两步法（Brower et al., hep-lat/9708008）：
1. **Stage 1（MAG）**：先固定到最大阿贝尔规范
2. **Stage 2（MCG）**：在 MAG 固定后的场上执行 Direct MCG

**关键参数**：
- `OmegaStage1` — MAG 阶段的 overrelaxation（默认 1.5）
- `Stage1MaxIterate` — MAG 阶段最大迭代次数（默认 100,000）
- `Omega` — MCG 阶段的 overrelaxation（默认 1.5）

### 5. Random Gauge Transform

**文件**: `GaugeFixing/CGaugeFixingRandom.h`

在每个格点生成随机规范群元素 `g(x)`，执行 `U_mu(x) → g(x) U_mu(x) g†(x+mu)`。用于测试规范不变性。

支持 SU(2) 和 SU(3)。还提供 `AlsoFixingFermion()` 和 `AlsoFixingAphys()` 用于同步变换费米子场和辅助规范场。

## 实现规范

### Kernel 模板化

所有局部更新 kernel 使用 `template<typename deviceGauge>`，支持多规范群：

```cpp
template<typename deviceGauge>
__global__ void _CLG_LAUNCH_BOUND
_kernelLocalUpdate(...)
```

Host 端通过 `switch (GetFieldType())` 显式实例化：

```cpp
switch (pResGauge->GetFieldType())
{
case EFT_GaugeSU2:
    _LAUNCH_KERNEL(_kernelLocalUpdate TMPARG(deviceSU2), ...);
    break;
case EFT_GaugeSU3:
    _LAUNCH_KERNEL(_kernelLocalUpdate TMPARG(deviceSU3), ...);
    break;
}
```

### 红黑 Checkerboard

所有局部迭代算法（Los Alamos、MAG、MCG）使用红黑分解：
1. 先更新所有 "偶" 格点（并行）
2. 再更新所有 "奇" 格点（并行）

这样保证同一轮 sweep 中相邻格点的更新互不干扰。

### 设备端数学

局部更新中的矩阵运算必须使用 `DeviceInlineTemplate.h` 中的全局模板函数：

```cpp
// 正确（模板通用）
_add(x, y); _dagger(x); _mul(x, c);

// 错误（仅对具体类型有效）
x.Add(y); x.Dagger(); x.MulComp(c);
```

## 测试策略

每个规范固定测试遵循三步验证：

1. **投影测试**：将场投影到目标子群（如 Z_2/Z_3/对角），验证 `CheckRes() ≈ 0`
2. **规范不变性**：随机规范变换后能量不变（`|E_before - E_after| < 容差`）
3. **恢复性**：规范固定算法应能将随机变换后的场恢复到投影状态（`CheckRes() ≈ 0`，能量恢复）

**测试文件**: `Code/CLGTest/Tests/TestGaugeFixing.cpp`

**注意**：Indirect MCG 测试依赖于 MAG 和 Direct MCG 各自单独通过。如果单独测试失败，Indirect 的恢复性测试也会失败。
