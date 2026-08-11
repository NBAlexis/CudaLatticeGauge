# Boson Field 模块

标量/玻色子场实现：实标量场、复标量场、Dirichlet 边界条件。

## 文件清单

| 文件 | 路径 | 说明 |
|------|------|------|
| `CFieldBoson.h` | `Data/Field/CFieldBoson.h` | 玻色子场抽象基类 |
| `CFieldBosonVN.h` | `Data/Field/Boson/CFieldBosonVN.h` | 模板基类 + 具体实例化 |
| `CFieldBosonVNKernel.h` | `Data/Field/Boson/CFieldBosonVNKernel.h` | 静态 kernel 接口 |
| `CFieldBosonReal.h` | `Data/Field/Boson/CFieldBosonReal.h` | 实值标量场 |

## 类层次

```
CBase
└── CField (abstract)
    └── CFieldBoson (abstract)
        └── CFieldBosonVN<deviceDataBoson, deviceGauge> (template)
            ├── CFieldBosonReal
            └── CFieldBosonRealD_ZSlice (Dirichlet)
```

## CFieldBoson（抽象基类）

**文件**: `Data/Field/CFieldBoson.h`

所有玻色子场的统一接口：

| 方法 | 说明 |
|------|------|
| `D()` / `D_WithMass()` | D 算子（纯虚） |
| `ApplyOperator()` | 按 `EFieldOperator` 分发的统一接口 |
| `CalculateForceOnGauge()` | 玻色子对规范场的力（纯虚） |
| `Energy()` | 能量计算（纯虚） |
| `InitialAsSource()` | 初始化为点源（纯虚） |
| `SetMass()` / `GetMass()` | 质量参数 |
| `CheckHermitian()` | 检查厄米性（纯虚） |
| `AxpyPlus/Minus()`, `ScalarMultply()` | BLAS 操作 |
| `Dot()`, `GetLength()`, `Dagger()` | 内积、长度、共轭 |

**EFieldOperator 分发**：与费米子场类似，`ApplyOperator` 按枚举值调用对应的虚函数（`D`, `Ddagger`, `DD`, `DDdagger`, `InverseD` 等）。

## CFieldBosonVN（模板基类）

**文件**: `Data/Field/Boson/CFieldBosonVN.h`

模板类 `CFieldBosonVN<deviceDataBoson, deviceGauge>`，其中：
- `deviceDataBoson` — 设备端玻色子数据类型（如 `Real`）
- `deviceGauge` — 对应的规范群矩阵类型（如 `deviceSU3`）

**核心成员**：`deviceDataBoson* m_pDeviceData` — 设备端数据数组。

**注册宏**：`__DEFINE_BOSON_FIELD(classname, deviceDataBoson, deviceGauge, fieldTypeEnum)` 一次性声明模板特化和注册。

**已实现的方法**：
- 完整线性代数（`Axpy`, `Dot`, `Mul`, `Dagger`, `ScalarMultply`, `GetLength`）
- `D()` 和 `D_WithMass()` — 委托给 `CFieldBosonVNKernel`
- `Energy()` — 动能 + 质量项 + 相互作用项
- `CalculateForceOnGauge()` — 玻色子对 gauge field 的力
- `InitialAsSource()` — 点源初始化
- `CheckHermitian()` — 厄米性检查

## CFieldBosonReal — 实标量场

**文件**: `Data/Field/Boson/CFieldBosonReal.h`

实值标量场，每个格点存储一个 `Real`：

```cpp
class CFieldBosonReal : public CFieldBosonVN<Real, deviceSU3>
```

**用途**：温度梯度模拟（`CActionTemperatureDistribution`）、实标量场耦合等。

### Dirichlet 变体

`CFieldBosonRealD_ZSlice` — Z 方向切片 Dirichlet 边界条件。通过 `__DEFINE_BOSON_FIELD_D` 宏生成。

## CFieldBosonVNKernel（静态 Kernel 接口）

**文件**: `Data/Field/Boson/CFieldBosonVNKernel.h`

所有玻色子场设备端操作的静态模板接口：

| 类别 | 方法 |
|------|------|
| D 算子 | `DOperator(...)`, `DOperatorWithMass(...)` |
| Dirichlet BC | `DOperator_D(...)`, `DOperatorWithMass_D(...)` |
| 力 | `CalculateForce(...)`, `CalculateForceWithMass(...)` |
| 能量 | `CalculateEnergy(...)`, `CalculateEnergyWithMass(...)` |
| 初始化 | `InitialAsSource(...)` |
| BLAS | `AxpyPlus/Minus(...)`, `ScalarMultply(...)`, `DotReal(...)` |

显式实例化：`Real/deviceSU3`。

## 关键设计模式

- **模板 + 显式实例化**：`CFieldBosonVN` 是模板，但实例化在编译期确定，通过 `__DEFINE_BOSON_FIELD` 宏注册到工厂。
- **静态 kernel 接口**：`CFieldBosonVNKernel` 将 host 端 kernel 调用与具体场类解耦。
- **与规范场耦合**：玻色子场通过 `deviceGauge` 模板参数与规范场耦合，D 算子和力计算都涉及 gauge link。
