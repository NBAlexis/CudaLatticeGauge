# AgentDriver

`AgentDriver` 是一个完全由 YAML 驱动的入口程序，用于让外部 agent（或脚本）无需修改代码即可执行 CLGLib 的三大核心任务：

- **Simulate**：热化并生成规范场配置。
- **GaugeFixing**：对已有配置做规范固定。
- **Measure**：加载配置、运行测量器、导出选定的结果。

> 详细架构与模块说明见 [wiki/architecture.md](../../wiki/architecture.md)。

---

## 快速开始

```bash
# 编译（默认 CUDA backend，GPU arch 86）
./build.sh

# 运行任意 YAML 任务
./../../Bin/Ubuntu/AgentDriver YAMLs/t02_simulate.yaml
./../../Bin/Ubuntu/AgentDriver YAMLs/t02_gaugefixing.yaml
./../../Bin/Ubuntu/AgentDriver YAMLs/t02_measure.yaml
```

`YAMLs/` 目录下已经放置了覆盖教程 02–06 的示例 YAML：

| YAML | 任务 | 说明 |
|------|------|------|
| `t02_simulate.yaml` | Simulate | 纯 SU(3) HMC，生成 5 个配置 |
| `t03_simulate.yaml` | Simulate | 纯 Z₂ 热浴，生成 5 个配置 |
| `t04_simulate.yaml` | Simulate | SU(3) + Wilson 费米子 HMC |
| `t05_simulate.yaml` | Simulate | SU(3) + KS 费米子 HMC（参数已调轻，供快速验证） |
| `t07_simulate_only_last.yaml` | Simulate | 只保存末组态并导出 Polyakov loop |
| `t02_gaugefixing.yaml` | GaugeFixing | Landau 规范固定示例 |
| `t02_measure.yaml` | Measure | Plaquette 测量 |
| `t03_measure.yaml` | Measure | Z₂ Plaquette 测量 |
| `t04_measure.yaml` | Measure | SU(3)+Wilson Plaquette 测量 |
| `t05_measure.yaml` | Measure | SU(3)+KS Plaquette 测量 |
| `t06_measure_polyakov.yaml` | Measure | Polyakov loop 测量 |
| `t06_measure_chiral.yaml` | Measure | 手征凝聚测量 |
| `t06_measure_wilson.yaml` | Measure | Wilson loop 测量 |

每个 YAML 运行后都会在同一目录下生成 `.log`，完整记录运行轨迹。

---

## 任务公共参数

所有任务都需要先描述格点、线程、随机数等基础信息：

```yaml
Dim: 4
Dir: 4
LatticeLength: [8, 8, 8, 8]
LatticeIndex: CIndexSquare
LatticeBoundary: CBoundaryConditionTorusSquare
ThreadAutoDecompose: 1
MaxThreadPerBlock: 256
RandomType: ER_Schrage
RandomSeed: 1234567
```

完整参数说明见 wiki：

- [全局参数](../../wiki/yaml-reference/01-global.md)
- [格点参数](../../wiki/yaml-reference/02-lattice.md)
- [随机数](../../wiki/yaml-reference/03-random.md)
- [索引与边界条件](../../wiki/yaml-reference/04-index-boundary.md)

---

## 场的 YAML 配置

CLGLib 的场通过 `FieldName` 指定类名，由 RTTI 工厂自动实例化。

### 规范场（Gauge Fields）

YAML 键名：`Gauge`（第一个）、`Gauge2`、`Gauge3`…，数量由 `GaugeFieldCount` 控制。

常用类名：

| 类名 | 说明 | 详见 |
|------|------|------|
| `CFieldGaugeSU3` | 标准 SU(3) 规范场 | [wiki/data-field-gauge.md](../../wiki/data-field-gauge.md) |
| `CFieldGaugeSU2` | SU(2) 规范场 | [yaml-reference/05-gauge-fields.md](../../wiki/yaml-reference/05-gauge-fields.md) |
| `CFieldGaugeU1` | U(1) 规范场 | [yaml-reference/05-gauge-fields.md](../../wiki/yaml-reference/05-gauge-fields.md) |
| `CFieldGaugeZ2` | Z₂ 离散规范场 | [yaml-reference/05-gauge-fields.md](../../wiki/yaml-reference/05-gauge-fields.md) |
| `CFieldGaugeSU3D` | Dirichlet 边界 SU(3) | [yaml-reference/05-gauge-fields.md](../../wiki/yaml-reference/05-gauge-fields.md) |
| `CFieldGaugeSU3TreeImproved` | 树级 Luscher-Weisz 改进 | [yaml-reference/05-gauge-fields.md](../../wiki/yaml-reference/05-gauge-fields.md) |
| `CFieldGaugeSU3OneLoopImproved` | 一圈 Symanzik 改进 | [yaml-reference/05-gauge-fields.md](../../wiki/yaml-reference/05-gauge-fields.md) |

示例：

```yaml
GaugeFieldCount: 1
Gauge:
    FieldName: CFieldGaugeSU3
    FieldInitialType: EFIT_Random
```

### 费米子场（Fermion Fields）

YAML 键名：`FermionField1`、`FermionField2`…，数量由 `FermionFieldCount` 控制。

常用类名：

| 类名 | 说明 | 详见 |
|------|------|------|
| `CFieldFermionWilsonSquareSU3` | 标准 Wilson Dirac | [wiki/data-field-fermion.md](../../wiki/data-field-fermion.md) |
| `CFieldFermionWilsonSquareCloverSU3` | Clover 改进 Wilson | [yaml-reference/08-fermion-fields.md](../../wiki/yaml-reference/08-fermion-fields.md) |
| `CFieldFermionKSSU3` | Staggered / KS 费米子 | [yaml-reference/08-fermion-fields.md](../../wiki/yaml-reference/08-fermion-fields.md) |
| `CFieldFermionHISQSU3` | HISQ 改进 KS | [yaml-reference/08-fermion-fields.md](../../wiki/yaml-reference/08-fermion-fields.md) |

示例：

```yaml
FermionFieldCount: 1
FermionField1:
    FieldName: CFieldFermionKSSU3
    FieldInitialType: EFIT_RandomGaussian
    Mass: 0.1
    FieldId: 2
```

### 玻色子场（Boson Fields）

YAML 键名：`BosonField1`、`BosonField2`…，数量由 `BosonFieldCount` 控制。

| 类名 | 说明 | 详见 |
|------|------|------|
| `CFieldBosonU1` | U(1) 玻色子场 | [wiki/data-field-boson.md](../../wiki/data-field-boson.md) |
| `CFieldBosonSU3` | SU(3) 玻色子场 | [yaml-reference/07-boson-fields.md](../../wiki/yaml-reference/07-boson-fields.md) |

---

## 作用量（Actions）

YAML 键名：`Action1`、`Action2`…，数量由 `ActionListLength` 控制。

| 类名 | 说明 | 详见 |
|------|------|------|
| `CActionGaugePlaquette` | 标准 Plaquette 作用量 | [wiki/data-action.md](../../wiki/data-action.md) |
| `CActionFermionKS` | KS 费米子作用量 | [yaml-reference/09-actions.md](../../wiki/yaml-reference/09-actions.md) |
| `CActionFermionKSCombined` | 组合 KS 作用量 | [yaml-reference/09-actions.md](../../wiki/yaml-reference/09-actions.md) |
| `CActionFermionHISQCombined` | HISQ 组合作用量 | [yaml-reference/09-actions.md](../../wiki/yaml-reference/09-actions.md) |
| `CActionPhi4` | Φ⁴ 标量场作用量 | [yaml-reference/09-actions.md](../../wiki/yaml-reference/09-actions.md) |
| `CActionDiscreteZNPlaquette` | Z_N 离散作用量 | [yaml-reference/09-actions.md](../../wiki/yaml-reference/09-actions.md) |

示例：

```yaml
ActionListLength: 2
Action1:
    ActionName: CActionGaugePlaquette
    Beta: 5.0
    GaugeFields: {1}

Action2:
    ActionName: CActionFermionKS
    FieldId: 2
    GaugeFields: {1}
```

---

## 测量器（Measurements）

YAML 键名：`Measure1`、`Measure2`…，数量由 `MeasureListLength` 控制。每个测量器会把每配置的结果写入 `CMeasurementManager`，导出 key 的命名规则为：

```
measure<MeasureId>.<Name>
```

例如 `Measure1: CMeasurePlaqutteEnergy` 默认写入的 key 为 `measure1.Plaquette`。

常用测量器：

| 类名 | 说明 | 默认导出 key | 详见 |
|------|------|--------------|------|
| `CMeasurePlaqutteEnergy` | Plaquette 能量 | `measure1.Plaquette` | [yaml-reference/17-measurements.md](../../wiki/yaml-reference/17-measurements.md) |
| `CMeasurePolyakovXY` | Polyakov loop | `measure1.PolyakovT` 等 | [yaml-reference/17-measurements.md](../../wiki/yaml-reference/17-measurements.md) |
| `CMeasureWilsonLoop` | Wilson loop | `measure1.WilsonLoop` | [yaml-reference/17-measurements.md](../../wiki/yaml-reference/17-measurements.md) |
| `CMeasureChiralCondensate` | 手征凝聚（Wilson） | `measure1.ChiralCondensate` | [yaml-reference/17-measurements.md](../../wiki/yaml-reference/17-measurements.md) |
| `CMeasureChiralCondensateKS` | 手征凝聚（KS） | `measure1.ChiralCondensate` | [yaml-reference/17-measurements.md](../../wiki/yaml-reference/17-measurements.md) |
| `CMeasureMesonCorrelator` | 介子关联函数 | 见代码 | [yaml-reference/17-measurements.md](../../wiki/yaml-reference/17-measurements.md) |
| `CMeasureAMomentumJG` | 角动量 JG 分解 | 见代码 | [yaml-reference/17-measurements.md](../../wiki/yaml-reference/17-measurements.md) |

完整测量器列表与参数见：

- [yaml-reference/17-measurements.md](../../wiki/yaml-reference/17-measurements.md)
- [wiki/measurement.md](../../wiki/measurement.md)

## Simulate 任务

`Simulate` 任务热化并生成规范场配置。

### 常用参数

| 键名 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `WarmUp` | INT | `0` | 热化步数。热化期间不保存组态，也不计入 `ConfigurationNumber` |
| `ConfigurationNumber` | INT | `1` | 生产阶段需要生成的被接受组态数 |
| `SaveStartIndex` | INT | `0` | 保存文件的起始编号 |
| `SavePrefix` | string | （必需） | 保存组态文件的前缀，最终文件名为 `<SavePrefix>_<N>.con` |
| `SaveConfiguration` | INT | `1` | `1`=保存每一个被接受组态；`0`=只保存最后一个组态，并建议配合 `Outputs` 导出测量结果 |
| `Outputs` | 块 | （可选） | 生产结束后导出指定测量 key，见 [Outputs 导出](#outputs-导出) |

### 默认保存模式：`SaveConfiguration: 1`

保存每一个被接受的组态，文件名为 `<SavePrefix>_<SaveStartIndex+N>.con`：

```yaml
Task: Simulate
WarmUp: 10
ConfigurationNumber: 5
SaveStartIndex: 0
SavePrefix: t02_su3_conf
SaveConfiguration: 1
```

生成 `t02_su3_conf_1.con` ... `t02_su3_conf_5.con`。

### 轻量保存模式：`SaveConfiguration: 0`

不保存中间组态，只保存最后一个用于衔接 Markov 链，并导出测量结果：

```yaml
Task: Simulate
WarmUp: 10
ConfigurationNumber: 100
SaveStartIndex: 0
SavePrefix: myrun_conf
SaveConfiguration: 0

Outputs:
    Names: [measure1.PolyakovT]
    Files: [myrun_polyakov]
```

运行后只生成：
- `myrun_conf_100.con`（最后一个组态）
- `myrun_polyakov.csv` / `myrun_polyakov.npy`（100 个 Polyakov loop 值）

> **注意**：`SaveConfiguration: 0` 时，如果只设置 `SavePrefix` 而不设置 `Outputs`，则只保存最后一个组态，不会留下中间的测量数据，通常不建议这样使用。

### Outputs 导出

`Outputs` 块与 `Measure` 任务中的用法相同。`Names` 是测量器导出的 key，`Files` 是对应的输出文件名（不带扩展名），会同时生成 `.csv` 和 `.npy`：

```yaml
Outputs:
    Names: [measure1.Plaquette, measure1.PolyakovT]
    Files: [myrun_plaquette, myrun_polyakov_t]
```

可用的 key 由配置的测量器决定，常见 key 见 [测量器](#测量器measurements) 一节。

---

## Measure 任务

`Measure` 任务加载已有配置，运行所有已声明的测量器，并导出选定的 key。

```yaml
Task: Measure
LoadPrefix: t02_su3_conf
StartN: 1
EndN: 5

Outputs:
    Names: [measure1.Plaquette]
    Files: [t02_plaquette]

MeasureListLength: 1
Measure1:
    MeasureName: CMeasurePlaqutteEnergy
    GaugeFields: {1}
```

运行后会生成 `t02_plaquette.csv` 和 `t02_plaquette.npy`。

---

## 更新器与积分器

`Simulate` 任务通过 `Updator` 块指定更新器。

| 键名 | 说明 | 详见 |
|------|------|------|
| `CHMC` | Hybrid Monte Carlo | [yaml-reference/15-updators.md](../../wiki/yaml-reference/15-updators.md) |
| `CHeatbath` | 热浴（主要用于离散群） | [yaml-reference/15-updators.md](../../wiki/yaml-reference/15-updators.md) |
| `IntegratorType` | `CIntegratorLeapFrog`、`CIntegratorOmelyan` 等 | [yaml-reference/16-integrators.md](../../wiki/yaml-reference/16-integrators.md) |

示例：

```yaml
Updator:
    UpdatorType: CHMC
    Metropolis: 1
    IntegratorType: CIntegratorLeapFrog
    IntegratorStepLength: 0.1
    IntegratorStep: 10
```

---

## 规范固定

`GaugeFixing` 任务需要先在场配置中声明 `GaugeFixing` 块：

```yaml
GaugeFixing:
    Name: CGaugeFixingLandauLosAlamos
    Omega: 1.5
    CheckErrorStep: 100
    MaxIterate: 1000
    Accuracy: 0.0001
```

详见：

- [wiki/gauge-fixing.md](../../wiki/gauge-fixing.md)
- [yaml-reference/14-gauge-fixing.md](../../wiki/yaml-reference/14-gauge-fixing.md)

---

## 求解器

涉及费米子时，需要在 YAML 中配置求解器：

```yaml
Solver:
    SolverName: CSLASolverGMRES
    SolverForFieldId: 2
    MaxDim: 20
    Accuracy: 0.0001
    Restart: 15
    AbsoluteAccuracy: 1
```

多移求解器用于 KS/HISQ HMC：

```yaml
MSSolver:
    SolverName: CMultiShiftBiCGStab
    SolverForFieldId: 2
    DiviationStep: 100
    MaxStep: 50
    Accuracy: 0.001
    AbsoluteAccuracy: 1
```

详见：

- [yaml-reference/10-solvers.md](../../wiki/yaml-reference/10-solvers.md)
- [yaml-reference/11-multi-shift-solvers.md](../../wiki/yaml-reference/11-multi-shift-solvers.md)
- [wiki/sparse-linalg.md](../../wiki/sparse-linalg.md)

---

## 输出与日志

- `AgentDriver` 默认使用 `GENERAL` 日志级别输出到 `stdout`。
- 如需写入文件或调整级别，可在 YAML 中设置：

```yaml
VerboseLevel: GENERAL
VerboseOutput: stdout
```

- 测量结果会同时输出 `.csv`（人类可读）和 `.npy`（Python 可读）。
- NumPy 读取示例：

```python
import numpy as np
plaq = np.load('t02_plaquette.npy')
print(plaq.shape, plaq.dtype)
```
