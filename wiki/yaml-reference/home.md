# CLGLib YAML 配置参考

> **警告**：此文档基于代码静态分析生成，**不随代码自动更新**。配置项的默认值、键名、行为请以实际代码为准。

## 概述

CLGLib 的所有物理对象（格点、规范场、费米子场、作用量、测量器、更新器）都通过 YAML 配置文件描述。本参考按 `CLGLibManager::InitialWithParameter()` 的**初始化顺序**组织，便于追踪每个参数在何时、被哪个类读取。

## 初始化流程

```
1. 全局参数 (CLGManager)     → 01-global.md
2. 格点常数 (CCommonData)     → 02-lattice.md
3. 随机数生成器 (CRandom)     → 03-random.md
4. 索引与边界条件            → 04-index-boundary.md
5. 规范场 (Gauge Fields)      → 05-gauge-fields.md
6. 边界场 (Boundary Fields)   → 06-boundary-fields.md
7. 玻色子场 (Boson Fields)    → 07-boson-fields.md
7a. Tensor2 场 (Plaquette)   → 07a-tensor2-fields.md
8. 费米子场 (Fermion Fields)  → 08-fermion-fields.md
9. 作用量 (Actions)           → 09-actions.md
10. 求解器 (Solvers)          → 10-solvers.md
11. 多移求解器 (MSSolvers)    → 11-multi-shift-solvers.md
12. 规范平滑 (Smearing)       → 12-gauge-smearing.md
13. Staple 缓存               → 13-staple-cache.md
14. 规范固定 (Gauge Fixing)   → 14-gauge-fixing.md
15. 更新器 (Updators)         → 15-updators.md
16. 积分器 (Integrators)      → 16-integrators.md
17. 测量器 (Measurements)     → 17-measurements.md
```

## 通用约定

### 参数读取模式

所有参数通过 `CParameters::FetchValueXXX()` 读取：

```cpp
// 基本类型
params.FetchValueINT(_T("KeyName"), iValue);      // 整数
params.FetchValueReal(_T("KeyName"), fValue);     // 实数
params.FetchValueDOUBLE(_T("KeyName"), dValue);   // 双精度
params.FetchValueString(_T("KeyName"), sValue);   // 字符串

// 数组
params.FetchValueArrayINT(_T("KeyName"), intArray);
params.FetchValueArrayReal(_T("KeyName"), realArray);
params.FetchValueArrayDOUBLE(_T("KeyName"), doubleArray);
```

如果键不存在，变量保持调用前的值。框架内部使用 `__FetchIntWithDefault` 等宏来提供默认值。

### 列表索引

当存在多个同类对象时，使用数字后缀：

| 对象类型 | YAML 键名 | 数量键 |
|---------|----------|--------|
| 作用量 | `Action1`, `Action2`, ... | `ActionListLength` |
| 测量器 | `Measure1`, `Measure2`, ... | `MeasureListLength` |
| 规范场 | `Gauge`, `Gauge2`, `Gauge3`, ... | `GaugeFieldCount` |
| 费米子场 | `FermionField1`, `FermionField2`, ... | `FermionFieldCount` |
| 玻色子场 | `BosonField1`, `BosonField2`, ... | `BosonFieldCount` |
| Tensor2 场 | `Tensor2Field1`, `Tensor2Field2`, ... | `Tensor2FieldCount` |
| 求解器 | `Solver`, `Solver2`, ... | (隐式) |
| 多移求解器 | `MSSolver`, `MSSolver2`, ... | (隐式) |
| 规范平滑 | `GaugeSmearing`, `GaugeSmearing2`, ... | (隐式) |
| Staple 缓存 | `StapleCache`, `StapleCache2`, ... | (隐式) |

### 类名注册

所有可通过 YAML 实例化的类都必须注册到工厂：

```cpp
// 头文件中：注册辅助（强制 Linux 静态链接保留构造函数）
__CLG_REGISTER_HELPER_HEADER(CFieldGaugeSU3)

// 类声明中
class CLGAPI CFieldGaugeSU3 : public CFieldGauge
{
    __CLGDECLARE_CLASS(CFieldGaugeSU3)
    // ...
};

// .cpp 文件中
__CLGIMPLEMENT_CLASS(CFieldGaugeSU3)
```

框架通过 `appCreate(_T("ClassName"))` 从字符串类名创建实例。`CLGLibManager` 内部对所有 `ActionName`、`FieldName`、`SolverName`、`MeasureName` 等字段都调用 `appCreate`。如果类名拼写错误、宏未启用，或 Linux 静态链接时注册构造函数被优化掉，返回 `NULL`。

### 危险行为

1. **`GetParameter` 返回 `*this`**：如果请求的键不存在，`CParameters::GetParameter(key)` 不返回错误，而是返回**当前对象自身**。后续读取会从错误的参数树取值。见 [FAQ](../../Tutorial/08-faq.md) Q9。
2. **`FetchValue*` 不区分"未找到"和"值为0"**：如果默认值设为 0 且键不存在，无法判断是否被显式设置。
3. **`FieldId` 冲突**：ID 1 保留给主规范场。费米子/玻色子场 ID 必须 `> 1` 且 `< kMaxFieldCount`（通常 255）。冲突时自动分配下一个可用 ID。

更多问题见 [FAQ](../../Tutorial/08-faq.md)。

---

## 导航

| 章节 | 内容 | 文件 |
|------|------|------|
| 1 | CLGManager 全局参数（含 `DevicePerNode`） | [01-global.md](01-global.md) |
| 2 | 格点常数（含 `GpuGrid` / `HaloWidth`） | [02-lattice.md](02-lattice.md) |
| - | **多 GPU 配置指南**（构建、切分、校验、1-vs-N 验证） | [../multi-gpu.md](../multi-gpu.md) |
| 3 | 随机数生成器 | [03-random.md](03-random.md) |
| 4 | 索引与边界条件 | [04-index-boundary.md](04-index-boundary.md) |
| 5 | 规范场 | [05-gauge-fields.md](05-gauge-fields.md) |
| 6 | 边界场 | [06-boundary-fields.md](06-boundary-fields.md) |
| 7 | 玻色子场 | [07-boson-fields.md](07-boson-fields.md) |
| 7a | Tensor2（plaquette）场 | [07a-tensor2-fields.md](07a-tensor2-fields.md) |
| 8 | 费米子场（索引） | [08-fermion-fields.md](08-fermion-fields.md) |
| 8a | Wilson 费米子参数 | [08a-wilson-fermion.md](08a-wilson-fermion.md) |
| 8b | KS / HISQ 费米子参数 | [08b-ks-hisq-fermion.md](08b-ks-hisq-fermion.md) |
| 8c | KS/HISQ 质量约定 | [08c-fermion-mass-convention.md](08c-fermion-mass-convention.md) |
| 9 | 作用量（索引） | [09-actions.md](09-actions.md) |
| 9a | 规范作用量 | [09a-gauge-actions.md](09a-gauge-actions.md) |
| 9b | 费米子作用量 | [09b-fermion-actions.md](09b-fermion-actions.md) |
| 9c | 标量、温度与离散作用量 | [09c-other-actions.md](09c-other-actions.md) |
| 10 | 求解器 | [10-solvers.md](10-solvers.md) |
| 11 | 多移求解器 | [11-multi-shift-solvers.md](11-multi-shift-solvers.md) |
| 12 | 规范平滑 | [12-gauge-smearing.md](12-gauge-smearing.md) |
| 13 | Staple 缓存 | [13-staple-cache.md](13-staple-cache.md) |
| 14 | 规范固定 | [14-gauge-fixing.md](14-gauge-fixing.md) |
| 15 | 更新器 | [15-updators.md](15-updators.md) |
| 16 | 积分器 | [16-integrators.md](16-integrators.md) |
| 17 | 测量器（索引） | [17-measurements.md](17-measurements.md) |
| 17a | 规范场测量器 | [17a-gauge-measurements.md](17a-gauge-measurements.md) |
| 17b | 费米子测量器 | [17b-fermion-measurements.md](17b-fermion-measurements.md) |
| 17c | 其他测量器 | [17c-other-measurements.md](17c-other-measurements.md) |
