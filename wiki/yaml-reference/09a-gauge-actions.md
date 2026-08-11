> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

# 09a. 规范作用量

## 规范作用量

### CActionGaugePlaquette

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquette.cpp`

**读取函数**：`CActionGaugePlaquette::InitialOtherParameters()` (lines 42-57)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `CloverEnergy` | INT | `0` | 否 | `CActionGaugePlaquette.cpp` ~L50 | 是否使用 Clover 定义计算 plaquette 能量。`1`=启用 |

### CActionGaugePlaquettePSU3

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquettePSU3.cu`

**读取函数**：`CActionGaugePlaquettePSU3::Initial()` (lines 117–124)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Beta` | DOUBLE | `0.1` | 否 | `CActionGaugePlaquettePSU3.cu` ~L119 | 对应伴随表示耦合常数 $\beta_A$；代码内部使用 $\tilde\beta = \beta_A/8$ |

**说明**：
- 实现 PSU(N) 伴随 Plaquette 作用量 $S_{\rm dyn} = -\tilde\beta \sum_P |\operatorname{Tr} U_P|^2$。
- 仅支持 SU2 与 SU3 规范场。
- 力的计算按每个 surrounding plaquette 单独求迹后再求和，不能先求 staple 和。
- Kernel 只把原始量 $Y_\mu = -\tilde\beta \sum_P \operatorname{Tr}(U_P) \Sigma_P$ 写入 force 场，最终 $U \cdot Y^\dagger + \mathrm{TA}$ 投影由 integrator 统一完成，与其他规范作用量一致。
- 详细推导见 `Docs/Applications/PSU3Adjoint.tex`。

### CActionGaugePlaquetteBoost

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquetteBoost.cu`

**读取函数**：`CActionGaugePlaquetteBoost::InitialOtherParameters()` (lines 352-363)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Boost` | Real | `0.1` | 否 | `CActionGaugePlaquetteBoost.cu` ~L359 | Boost 参数，存入 `CCommonData::m_fG` |

### CActionGaugePlaquetteAcceleration

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquetteAcceleration.cu`

**读取函数**：`CActionGaugePlaquetteAcceleration::InitialOtherParameters()` (lines 335-346)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `AccG` | Real | `0.1` | 否 | `CActionGaugePlaquetteAcceleration.cu` ~L342 | 加速参数 G，存入 `CCommonData::m_fG` |

### CActionGaugePlaquetteRigidAcc

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquetteRigidAcc.cu`

**读取函数**：`CActionGaugePlaquetteRigidAcc::InitialOtherParameters()` (lines 268-299)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `AccG` | Real | `0.1` | 否 | `CActionGaugePlaquetteRigidAcc.cu` ~L275 | 加速参数 G，存入 `CCommonData::m_fG` |
| `Dirichlet` | INT | `1` | 否 | `CActionGaugePlaquetteRigidAcc.cu` ~L295 | 是否使用 Dirichlet 边界。`1`=启用，`0`=禁用 |

### CActionGaugePlaquetteAtGradient

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquetteAtGradient.cu`

**读取函数**：`CActionGaugePlaquetteAtGradient::InitialOtherParameters()` (lines 314-325)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Xi` | DOUBLE[] | `[1.0]` | 否 | `CActionGaugePlaquetteAtGradient.cu` ~L318 | Xi 参数数组。自动扩展至 `_HC_Lz` 个元素（重复最后一个值） |

### CActionGaugePlaquetteGradient (注册名 `CActionGaugePlaquetteGradient`)

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquetteBetaGradient.cu`

**读取函数**：`CActionGaugePlaquetteGradient::InitialOtherParameters()` (lines 291-302)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Beta` | DOUBLE[] | `[5.0]` | 否 | `CActionGaugePlaquetteBetaGradient.cu` ~L295 | 各层 Beta 值数组。每个元素除以 `Nc`。自动扩展至 `_HC_Lz` 个元素 |

### CActionGaugePlaquetteCylinder

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquetteCylinder.cu`

**读取函数**：`CActionGaugePlaquetteCylinder::Initial()` (lines 291-330)

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Beta` | DOUBLE[] | `[5.0]` | 否 | `CActionGaugePlaquetteCylinder.cu` ~L314 | 逐 r 层 Beta 值数组。每个元素除以 `Nc`。自动扩展至 `_HC_Lx` 个元素 |
| `RStart` | Real | `1.0` | 否 | `CActionGaugePlaquetteCylinder.cu` ~L294 | 径向起点，$r(n)=\mathrm{RStart}+\mathrm{DeltaR}\cdot n_0$，必须大于 0（r=0 为奇点） |
| `REnd` | Real | `RStart + DeltaR * Lx` | 否 | `CActionGaugePlaquetteCylinder.cu` ~L296 | 径向终点（仅记录与显示） |
| `DeltaR` | Real | `1.0` | 否 | `CActionGaugePlaquetteCylinder.cu` ~L295 | 径向格距 |
| `CloverEnergy` | INT | `0` | 否 | `CActionGaugePlaquetteCylinder.cu` ~L303 | 能量计算方式。`1`=clover（逐站点），`0`=plaquette（显式 4 角耦合，Dirichlet 边界处与 clover 不同） |

**说明**：
- 柱坐标作用量：方向 0,1,2,3 解释为 $(r,\phi,z,t)$，含 $\phi$ 的平面权重 $1/r$，其余平面权重 $r$；一个 plaquette 的耦合取其 4 个角的平均值。
- r 方向建议 Dirichlet 边界（`CFieldGaugeSU3D` + `Period : [0, 1, 1, 1]`），测试见 `TestUpdatorCylinder`（clover）与 `TestUpdatorCylinderPlaq`（plaquette）。
- 推导见 `Docs/Applications/Cylinder.tex`。

### CActionGaugePlaquetteAnisotropic

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquetteAnisotropic.cpp`

**读取函数**：`CActionGaugePlaquetteAnisotropic::Initial()`

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Xi` | DOUBLE | `1.0` | 否 | `CActionGaugePlaquetteAnisotropic.cpp` ~L34 | 裸规范各向异性 $\xi_0$，必须为正；`Xi=1` 时退化为基类各向同性路径 |

**说明**：
- 各向异性 tree-level Symanzik 作用量 $S=\sum w(\mu,\nu)\beta[c_P P+c_R R]$，$c_P=1$，$c_R=-1/20$（由规范场的 `RectOverPlaq=-0.05` 提供），$\beta=10/g^2$；纯空间平面权重 $1/\xi_0$，含时间方向的平面权重 $\xi_0$。
- `Xi` 是**裸的可调输入**（与其它 action 的 `Xi` 参数惯例一致），与测量/重整化各向异性一般不同，与费米子侧 `XiF` 相互独立。
- 适用规范场：`CFieldGaugeSU3`（仅 plaquette）、`CFieldGaugeSU3TreeImproved`（plaquette+rectangle）、`CFieldGaugeSU3OneLoopImproved`（plaquette+rectangle+twisted loop）。
- twisted loop（Ct 项）的三元组按连续极限匹配加权：弱场展开下一个 twisted loop 同时含三个平面的场强，纯空间三元组 $Q_{123}$ 权重 $2/\xi_0-\xi_0$，含时间三元组 $Q_{\sigma\sigma\tau}$ 权重 $\xi_0$（field kernel 内传 $q=1/\xi_0$，即 $2q-1/q$ 与 $1/q$，能量与力使用完全相同的权重），推导见 `Docs/Applications/AnisotropicSymanzik.tex`。**注意**：该 Ct 项是 continuum-matched 的启发式推广，不是完整的 one-loop 各向异性 Symanzik 改善；严格用法为 tree-level（Ct=0、Cr=-1/20），MILC one-loop 系数只在 $\xi_0=1$ 有效。
- 注意 field 层 kernel 约定与 action 相反，代码内传 $1/\xi_0$；各向异性 staple-energy 不存在，anisotropic 分支收到非空 staple 会报错退出。
- 诊断测试：`TestHmcDiagnostics_TreeImproved_HISQ`（anisotropic action，`LinkTrials: 4`，覆盖含时间方向的 link）、`TestHmcDiagnostics_TreeImproved_AnisotropicHISQ`、`TestHmcDiagnostics_OneLoop_AnisotropicHISQ`（anisotropic action + one-loop 场，FD 覆盖 twisted-loop 分支）、`TestHmcDiagnostics_AnisotropicChecks`。

### CActionGaugePlaquetteRotatingT 模板系列

**来源文件**：`Code/CLGLib/Data/Action/CActionGaugePlaquetteRotatingT.cu`

**读取函数**：`CActionGaugePlaquetteRotatingT::Initial()`

| 键名 | 类型 | 默认值 | 必需 | 代码位置 | 说明 |
|------|------|--------|------|----------|------|
| `Omega` | DOUBLE | `0.1` | 否 | `CActionGaugePlaquetteRotatingT.cu` ~L982 | 旋转角速度 |
| `Xi` | DOUBLE | `1.0` | 否 | `CActionGaugePlaquetteRotatingT.cu` ~L986 | Xi 参数 |
| `CloverEnergy` | INT | `0` | 否 | `CActionGaugePlaquetteRotatingT.cu` ~L990 | 是否使用 Clover 能量。`1`=启用 |
| `ShiftCoord` | INT | `0` | 否 | `CActionGaugePlaquetteRotatingT.cu` ~L999 | 坐标半偏移。非零则启用 |
| `Torus` | INT | `0` | 否 | `CActionGaugePlaquetteRotatingT.cu` ~L1003 | 环面边界。非零则启用 |
| `RotationBetaScale` | STRING | `ERBSF_Naive` | 否 | `CActionGaugePlaquetteRotatingT.cu` ~L1012 | 旋转项（S1/S2）beta 缩放模式，见下表。仅作用于旋转项，S0 plaquette/rectangle 保持原始 beta |
| `ScaleFactor` | DOUBLE | `1.0` | 否 | `CActionGaugePlaquetteRotatingT.cu` ~L1024 | 仅当 `RotationBetaScale = ERBSF_Custom` 时读取 |

**`RotationBetaScale` 取值**（枚举 `ERotationBetaScaleFactor`，定义于 `CActionGaugePlaquetteRotatingT.h`）：

| 取值 | 缩放因子 | 用途 |
|------|----------|------|
| `ERBSF_Naive`（默认） | `1.0` | 普通 Wilson 规范场，`beta = 6/g^2` |
| `ERBSF_TreeImprove` | `0.6`（= 6/10） | tree level improved 规范场，`beta = 10/g^2` |
| `ERBSF_OneLoopImprove` | `0.6`（= 6/10） | one-loop improved 规范场，同为 `beta = 10/g^2`（plaquette 系数归一为 1） |
| `ERBSF_Custom` | 读取 `ScaleFactor`（默认 `1.0`） | 手动指定缩放因子 |

**说明**：旋转项 S1/S2（依赖 `Omega`）的推导假设 Wilson 归一化 `beta = 6/g^2`，而改进规范场（`CFieldGaugeSU3TreeImproved` / `CFieldGaugeSU3OneLoopImproved`）采用 `beta = 10/g^2`。若不缩放，改进场上的旋转耦合会被放大 `10/6 = 5/3` 倍。缩放因子经 `GetRotationBetaOverN()` 同时作用于旋转项的 energy 与 force，故 HMC 的 h_diff 仍守恒；仅平衡态可观测量随之移动。已验证测试：`RQuenchTorusShift`（tree-improved 规范场 + 旋转，`TestSuit_Rotation.yaml`），debug 与 release 均通过（h_diff < 阈值）。

**模板实例化类**：
- `CActionGaugePlaquetteRotatingU1` — U(1) 旋转 Plaquette
- `CActionGaugePlaquetteRotating` — SU(3) 旋转 Plaquette
- `CActionGaugePlaquetteRotatingSU2` — SU(2) 旋转 Plaquette
- `CActionGaugePlaquetteRotatingSU4` — SU(4) 旋转 Plaquette
- `CActionGaugePlaquetteRotatingU1_3D` — 3D U(1) 旋转 Plaquette
- `CActionGaugePlaquetteRotating3D` — 3D SU(3) 旋转 Plaquette

**注意**：旧的非模板版本（`CActionGaugePlaquetteRotating.h` / `CActionGaugePlaquetteRotating3D.h`）被 `#if 0` 包裹，已禁用。活跃版本为模板化的 `CActionGaugePlaquetteRotatingT` 家族。

---


---

[< 返回 09. 作用量](09-actions.md) | [< 返回 yaml-reference 目录](home.md)
