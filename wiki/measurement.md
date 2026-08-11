# Measurement 模块

物理量测量：作用量、拓扑荷、Wilson 圈、Polyakov 圈、介子关联函数、角动量、手征凝聚等。

## 文件清单

### 基类与管理器

| 文件 | 路径 | 说明 |
|------|------|------|
| `CMeasure.h` | `Measurement/CMeasure.h` | 测量基类 + 随机估计基类 |
| `CMeasurementManager.h` | `Measurement/CMeasurementManager.h` | 测量管理器 |

### 基本测量

| 文件 | 路径 | 说明 |
|------|------|------|
| `CMeasureAction.h` | `Measurement/` | 总作用量 |
| `CMeasurePlaqutteEnergy.h` | `Measurement/` | Plaquette 能量 |
| `CMeasurePolyakovXY.h` | `Measurement/` | Polyakov 圈密度（XY 面） |
| `CMeasurePolyakovXY3D.h` | `Measurement/` | Polyakov 圈密度 3D |
| `CMeasureTopologicChargeXY.h` | `Measurement/` | 拓扑荷密度（XY 面） |

### Wilson 圈与势

| 文件 | 路径 | 说明 |
|------|------|------|
| `CMeasureWilsonLoop.h` | `Measurement/` | Wilson 圈 RT 关联 |
| `CMeasureWilsonLoopXY.h` | `Measurement/` | Wilson 圈（XY 面限制） |
| `CMeasureWilsonLoopWithPath.h` | `Measurement/` | 自定义路径 Wilson 圈 |

### 介子关联函数

| 文件 | 路径 | 说明 |
|------|------|------|
| `CMeasureMesonCorrelator.h` | `Measurement/` | Wilson 费米子介子关联 |
| `CMeasureMesonCorrelatorStaggered.h` | `Measurement/` | Staggered 介子关联（20 种） |
| `CMeasureMesonCorrelatorStaggeredSimple2.h` | `Measurement/` | 简化 Staggered 介子关联 |

### 手征与随机估计

| 文件 | 路径 | 说明 |
|------|------|------|
| `CMeasureChiralCondensate.h` | `Measurement/` | 手征凝聚（随机估计） |
| `CMeasureChiralCondensateKS.h` | `Measurement/` | KS 手征凝聚 |
| `CMeasureConnectedSusceptibilityKS.h` | `Measurement/` | KS 连通 susceptibility |

### 角动量

| 文件 | 路径 | 说明 |
|------|------|------|
| `CMeasureAMomentumJG.h` | `Measurement/` | 角动量密度（JG 等） |
| `CMeasureAMomentumStochastic.h` | `Measurement/` | 角动量随机估计 |
| `CMeasureAngularMomentumKS.h` | `Measurement/` | KS 费米子角动量 |
| `CMeasureAngularMomentumKSREM.h` | `Measurement/` | KS 角动量 REM 变体 |

### 其他

| 文件 | 路径 | 说明 |
|------|------|------|
| `CMeasureBerryPhase.h` | `Measurement/` | Berry 相 / 动量空间 phi |
| `CMeasureBosonCond.h` | `Measurement/` | 玻色子凝聚 |
| `CMeasureBosonValue.h` | `Measurement/` | 玻色子场值 |
| `CMeasureChargeAndCurrents.h` | `Measurement/` | 电荷与电流密度（已废弃） |
| `CMeasurePandChiralTalor.h` | `Measurement/` | Polyakov/手征 Taylor 展开 |
| `CMeasurePandChiralTalorKS.h` | `Measurement/` | KS Taylor 展开 |

## CMeasure（抽象基类）

**文件**: `Measurement/CMeasure.h`

所有测量类的统一接口：

| 方法 | 说明 |
|------|------|
| `Initial(pOwner, params)` | 从参数初始化 |
| `OnConfigurationAccepted(gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, staple)` | 每个接受的配置调用 |
| `OnConfigurationAcceptedSingleField(...)` | 单场的配置接受回调 |
| `SourceSanning(...)` | 源扫描模式回调 |
| `OnConfigurationAcceptedZ4(gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, staple, pZ4, pInverseZ4, bStart, bEnd)` | Z4 随机源回调 |
| `Average()` | 对累积结果取平均 |
| `Report()` | 输出结果（纯虚） |
| `Reset()` | 清空累积数据 |
| `IsGaugeOrBosonMeasurement()` | 是否为规范/玻色子测量（纯虚） |
| `IsSourceScanning()` | 是否使用源扫描（纯虚） |
| `NeedGaugeSmearing()` | 是否需要规范平滑 |

**CMeasureStochastic**（随机估计基类）：
- 添加 `m_uiFieldCount` — Z4 随机源的数量
- `OnConfigurationAcceptedZ4SingleField(...)` — 每个 Z4 源的回调

**静态辅助方法**：
- `FillDataWithR` — 径向分布填充
- `ReportDistributionXY` — XY 面分布报告
- `XYDataToRdistri` — XY 数据转径向分布
- `_ZeroXYPlane` / `_AverageXYPlane` — XY 面操作

## CMeasurementManager

**文件**: `Measurement/CMeasurementManager.h`

持有并调度所有测量对象：

| 方法 | 说明 |
|------|------|
| `OnConfigurationAccepted(gaugeNum, bosonNum, tensor2Num, gaugeFields, bosonFields, tensor2Fields, staple)` | 分发到所有注册测量 |
| `OnUpdateFinished(...)` | 调用 `Average()` 和可选 `Report()` |
| `Reset()` | 重置所有测量 |
| `Report()` | 报告所有测量 |
| `GetMeasureById(id)` | 按 ID 查找测量 |

## 测量分类

### 规范场测量

| 类 | 测量内容 | 输出 |
|----|---------|------|
| `CMeasureAction` | 总作用量（gauge + fermion） | 每配置能量值 |
| `CMeasurePlaqutteEnergy` | Plaquette 能量 | 平均 trace，可选 V0 分量 |
| `CMeasurePolyakovXY` | Polyakov 圈密度 | XY/X/Y/Z/T 切片分布，径向分布 |
| `CMeasurePolyakovXY3D` | Polyakov 圈密度（3D） | XY 面分布 |
| `CMeasureTopologicChargeXY` | 拓扑荷密度 | XY 面分布 |
| `CMeasureWilsonLoop` | Wilson 圈 | R×T 关联函数 |
| `CMeasureWilsonLoopXY` | Wilson 圈（XY 面） | R×T 关联函数 |
| `CMeasureWilsonLoopWithPath` | 自定义路径 | 路径 trace |

### 费米子测量

| 类 | 测量内容 | 方法 |
|----|---------|------|
| `CMeasureMesonCorrelator` | Wilson 介子关联 | 可配置 gamma 矩阵插入 |
| `CMeasureMesonCorrelatorStaggered` | Staggered 介子（20 种） | 含 taste 结构 |
| `CMeasureMesonCorrelatorStaggeredSimple2` | 简化 Staggered 介子 | 3 传播子，PRD 38, 2245 |
| `CMeasureChiralCondensate` | 手征凝聚（9 通道） | 随机估计，含径向分布 |
| `CMeasureChiralCondensateKS` | KS 手征凝聚 | 标量、连通 susceptibility、15 gamma/sigma |
| `CMeasureConnectedSusceptibilityKS` | 连通 susceptibility | 零动量源 |

### 角动量测量

| 类 | 内容 | 特点 |
|----|------|------|
| `CMeasureAMomentumJG` | JG, JS, JGChen, JGSurf, JGPot 等 | 旋转参考系，XY 面密度 |
| `CMeasureAMomentumStochastic` | JL, JS, JLPure, JLJM, JPot | Z4 随机估计 |
| `CMeasureAngularMomentumKS` | KS 费米子角动量 | 轨道/自旋/势项，XY/Z 切片 |
| `CMeasureAngularMomentumKSREM` | REM 变体 | 重写矩阵应用方法 |

### 其他测量

| 类 | 内容 |
|----|------|
| `CMeasureBerryPhase` | Berry 相 / 动量空间 phi（Wilson/KS） |
| `CMeasureBosonCond` | 玻色子凝聚期望值 |
| `CMeasureBosonValue*` | 玻色子场值（宏生成：Real/U1/SU2...SU8） |
| `CMeasurePandChiralTalor` | Polyakov/手征 Taylor 展开（旋转系） |
| `CMeasurePandChiralTalorKS` | KS Taylor 展开 |

## 添加新 Measurement 的指南

1. **继承 `CMeasure`**（或 `CMeasureStochastic` 用于随机估计）。
2. **重写回调方法**：
   - `OnConfigurationAcceptedSingleField(...)` — 规范/玻色子测量
   - `OnConfigurationAcceptedZ4SingleField(...)` — 费米子随机估计
   - `SourceSanningSingleField(...)` — 源扫描模式
3. **实现 `Report()`** — 输出累积结果。
4. **实现 `Reset()`** — 清空数据。
5. **注册到 RTTI 工厂**：
   - 头文件加 `__CLG_REGISTER_HELPER_HEADER(CMyMeasure)`
   - 实现文件加 `__CLGIMPLEMENT_CLASS(CMyMeasure)`
   - 无需手动修改 `CMeasurementManager`，`CLGLibManager` 会通过 YAML 中的 `MeasureName` 自动创建

**常见模式**：
- 在 `Initial()` 中分配设备缓冲区（XY 面密度、Z 切片等）
- 在回调中启动 device kernel 计算局部量
- 使用 `UpdateRealResult()` / `UpdateComplexResult()` 累积结果
- 在 `Average()` 中除以配置数
- 在 `Report()` 中打印或导出 CSV

示例 YAML：
```yaml
Measure1:
    MeasureName: CMeasurePlaqutteEnergy
    GaugeFields: [1]
```
