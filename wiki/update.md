# Update 模块

场更新算法：HMC（混合蒙特卡洛）、积分器、热浴、Staple 缓存。

## 文件清单

### 基类

| 文件 | 路径 | 说明 |
|------|------|------|
| `CUpdator.h` | `Update/CUpdator.h` | 更新器抽象基类 |
| `CStapleCache.h` | `Update/CStapleCache.h` | Staple 缓存管理 |

### HMC 与积分器

| 文件 | 路径 | 说明 |
|------|------|------|
| `CHMC.h` | `Update/Continous/` | HMC 更新器 |
| `CIntegrator.h` | `Update/Continous/` | 积分器基类 + 嵌套/多级嵌套 |
| `CIntegratorLeapFrog.h` | `Update/Continous/` | Leapfrog（Verlet）积分器 |
| `CIntegratorOmelyan.h` | `Update/Continous/` | Omelyan 二阶最小范数积分器 |
| `CIntegratorForceGradient.h` | `Update/Continous/` | 力梯度四阶积分器 |
| `CIntegratorNestedLeapFrog.h` | `Update/Continous/` | 嵌套 Leapfrog |
| `CIntegratorNestedOmelyan.h` | `Update/Continous/` | 嵌套 Omelyan |
| `CIntegratorNestedForceGradient.h` | `Update/Continous/` | 嵌套力梯度 |
| `CIntegratorNested11Stage.h` | `Update/Continous/` | 11 阶段高阶积分器 |
| `CIntegratorMultiLevelNestedForceGradient.h` | `Update/Continous/` | 多级嵌套力梯度 |
| `CIntegratorMultiLevelOmelyan.h` | `Update/Continous/` | 多级嵌套 Omelyan |

### 离散更新

| 文件 | 路径 | 说明 |
|------|------|------|
| `CHeatbath.h` | `Update/Discrete/` | 热浴更新（离散群） |

## CUpdator（抽象基类）

**文件**: `Update/CUpdator.h`

所有更新算法的统一接口：

| 方法 | 说明 |
|------|------|
| `Update(iSteps, bMeasure)` | 执行更新（纯虚） |
| `UpdateUntileAccept(iSteps, bMeasure)` | 循环直到接受 |
| `CalculateEnergy()` | 计算能量（纯虚） |
| `Initial(pOwner, params)` | 初始化（纯虚） |
| `SetAutoCorrection(bAutoCorrection)` | 开关自校正（纯虚） |
| `SaveConfiguration(uiUpdateStep)` | 保存规范场到磁盘 |
| `GetHDiff()` / `GetHValue()` | Hamiltonian 差值诊断 |

## HMC 更新

### CHMC

**文件**: `Update/Continous/CHMC.h`

混合蒙特卡洛更新器，用于连续规范群（U(1)、SU(2)、SU(3) 等）。

**流程**：
1. 委托 `CIntegrator` 执行一条分子动力学轨迹
2. 执行 Metropolis 接受/拒绝测试
3. 若拒绝，恢复原场配置

**关键成员**：`m_pIntegrator` — 关联的积分器实例。

### CIntegrator（积分器基类）

**文件**: `Update/Continous/CIntegrator.h`

管理 HMC 轨迹的所有基础设施：

**场管理**：
- 规范场/玻色子场指针列表
- 动量场（`P`）和力场（`F`）
- 备份场（用于 Metropolis 拒绝时恢复）

**核心操作**：
| 方法 | 说明 |
|------|------|
| `Evaluate()` | 执行一条轨迹（纯虚） |
| `UpdateP(fStep, ePhase)` | 从力更新动量：`P += fStep * F` |
| `UpdateU(fStep)` | 从动量更新规范场：`U = exp(fStep * P) * U` |
| `GetEnergy(bBeforeEvolution, actions)` | 计算总能量 |
| `InitialMomentumNoise()` | 抽取随机动量（高斯分布） |
| `PreserveFields()` / `RecoverFields()` | 备份/恢复场配置 |
| `CalcForceOfActions(...)` | 汇总所有作用量的力 |
| `RequireGaugeSmearing()` | 力计算前触发规范平滑 |
| `OnCacheAndSmearing(updateMode)` | 缓存 staple 并应用平滑 |

**自适应步长**：
- `ChangeStepCount(bGrow)` — 倍增或减半步数
- `ChangeStepCountTo(uiStep)` — 设到指定步数

### 积分器变体

| 类 | 父类 | 阶数 | 特点 |
|----|------|------|------|
| `CIntegratorLeapFrog` | `CIntegrator` | 2 | 标准 leapfrog，O(ε²) 误差 |
| `CIntegratorOmelyan` | `CIntegrator` | 2 | 最小范数二阶，系数 λ 优化 |
| `CIntegratorForceGradient` | `CIntegrator` | 4 | 力梯度项，更高精度 |
| `CIntegratorNestedLeapFrog` | `CNestedIntegrator` | 2 | 外层 leapfrog + 内层 leapfrog |
| `CIntegratorNestedOmelyan` | `CNestedIntegrator` | 2 | 外层 Omelyan + 内层 leapfrog/Omelyan |
| `CIntegratorNestedForceGradient` | `CNestedIntegrator` | 4 | 外层力梯度 + 内层 leapfrog/力梯度 |
| `CIntegratorNested11Stage` | `CNestedIntegrator` | 高阶 | 11 阶段速度型积分器 |
| `CIntegratorMultiLevelNestedForceGradient` | `CMultiLevelNestedIntegrator` | 4 | 多级嵌套，最外层力梯度 |
| `CIntegratorMultiLevelOmelyan` | `CMultiLevelNestedIntegrator` | 2 | 多级嵌套，最外层 Omelyan |

**嵌套积分器**（`CNestedIntegrator`）：将力分为 gauge 部分和 fermion 部分，内层用更小的步长积分费米子力，外层用较大步长积分规范力。

**多级嵌套**（`CMultiLevelNestedIntegrator`）：支持多于两级的力拆分，递归执行嵌套积分。

### 积分器参数（YAML）

```yaml
Updator:
    Name: CHMC
    Integrator: CIntegratorOmelyan
    Epsilon: 0.01
    StepCount: 10
    AutoCorrection: 1
    # 嵌套积分器额外参数：
    NestedStepCount: 5
    Lambda: 0.193183
```

## 热浴更新

### CHeatbath

**文件**: `Update/Discrete/CHeatbath.h`

用于离散规范群（Z2、Zn、四面体、八面体、二十面体群）的直接局部更新。无需分子动力学，每个 link 独立按 Boltzmann 权重采样。

**特点**：
- 无 Metropolis 拒绝（每次更新都接受）
- 无自校正（`SetAutoCorrection` 为空操作）
- 适合离散群或强耦合区域

## Staple 缓存

### CStapleCache

**文件**: `Update/CStapleCache.h`

预计算并缓存 gauge field 的 staple、plaquette、Fmunu 和旋转链接缓冲区，避免积分器每一步重复计算。

**模板类** `CStapleCacheT<FieldType>` 按规范场类型参数化：
- `CStapleCacheSU3` — SU(3) staple 缓存

**缓存内容**：
- Staples（用于力计算）
- Plaquettes（用于能量计算）
- Fmunu（场强张量）
- Rotation links（KS 费米子旋转改进）

**关键方法**：
| 方法 | 说明 |
|------|------|
| `InitialBuffers(byFieldId)` | 分配设备缓冲区 |
| `Cache(pGauge, eCC)` | 根据调用阶段计算并缓存 |
| `GetStaples()` / `GetPlaquttes()` / `GetFmunu()` / `GetRotationBuffer()` | 获取缓存指针 |

**调用阶段**（`ECacheCall`）：
- 能量计算前缓存 plaquette
- 力计算前缓存 staple
- HMC 准备阶段缓存所有内容

## 类层次

```
CUpdator
├── CHMC
│   └── uses CIntegrator
└── CHeatbath

CIntegrator
├── CIntegratorLeapFrog
├── CIntegratorOmelyan
├── CIntegratorForceGradient
├── CNestedIntegrator
│   ├── CIntegratorNestedLeapFrog
│   ├── CIntegratorNestedOmelyan
│   ├── CIntegratorNestedForceGradient
│   └── CIntegratorNested11Stage
└── CMultiLevelNestedIntegrator
    ├── CIntegratorMultiLevelNestedForceGradient
    └── CIntegratorMultiLevelOmelyan

CStapleCache
└── CStapleCacheT<FieldType>
    └── CStapleCacheSU3
```

## HMC 轨迹执行流程

1. **Prepare**：
   - 备份当前场配置（`PreserveFields`）
   - 生成随机动量（`InitialMomentumNoise`）
   - 缓存 staple（`OnCacheAndSmearing`）
   - 计算初始能量

2. **Evaluate**（积分器执行轨迹）：
   - 反复调用 `UpdateP` 和 `UpdateU`
   - 每次力计算前刷新 staple 缓存
   - 规范平滑在需要时触发

3. **Metropolis 测试**：
   - 计算最终能量
   - ΔH = H_final - H_initial
   - 以概率 `min(1, exp(-ΔH))` 接受

4. **OnFinishTrajectory**：
   - 若拒绝：`RecoverFields` 恢复备份
   - 若接受：保留新配置，缓存能量
