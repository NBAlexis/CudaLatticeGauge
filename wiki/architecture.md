# 整体架构

CLGLib 的模块关系、数据流、初始化流程和跨模块设计模式。

## 模块关系图

```
                    ┌─────────────┐
                    │  CLGTest    │
                    │ (测试框架)   │
                    └──────┬──────┘
                           │
                    ┌──────▼──────┐
                    │ CLGLibManager│
                    │  (核心管理器) │
                    └──────┬──────┘
                           │
        ┌──────────────────┼──────────────────┐
        │                  │                  │
   ┌────▼────┐      ┌──────▼──────┐    ┌─────▼─────┐
   │ Lattice │      │    Core     │    │  Platform │
   │ (格点)   │      │ (CUDA基础设施)│    │  (工具库)  │
   └────┬────┘      └─────────────┘    └───────────┘
        │
   ┌────┴────┬──────────┬──────────┬──────────┐
   │         │          │          │          │
┌──▼──┐  ┌──▼──┐   ┌───▼───┐  ┌──▼──┐   ┌───▼───┐
│Gauge│  │Fermion│  │ Boson │  │Action│   │Update │
│(规范场)│  │(费米子)│  │(玻色子) │  │(作用量)│   │(更新器)│
└──┬──┘  └──┬──┘   └───┬───┘  └──┬──┘   └───┬───┘
   │        │          │         │          │
   │   ┌────┘          │         │          │
   │   │               │         │          │
┌──▼───▼───┐      ┌───▼───┐  ┌──▼──┐   ┌───▼────┐
│GaugeFixing│      │Sparse │  │Measure│   │Smearing│
│(规范固定) │      │ LA    │  │(测量) │   │(平滑)  │
│           │      │(求解器)│  │       │   │        │
└───────────┘      └───────┘  └───────┘   └────────┘
        │
   ┌────┘
   │
┌──▼────────┐
│  Math     │
│(SU(N)/Gamma│
│/FFT/Rand) │
└───────────┘
```

## 核心数据流

### HMC 更新流程

```
1. CHMC::Update()
   └── 2. CIntegrator::Evaluate() [轨迹]
       ├── 3. InitialMomentumNoise()  → 随机高斯动量
       ├── 4. OnCacheAndSmearing()    → 缓存 staple / 规范平滑
       ├── 5. 循环 step 次:
       │   ├── UpdateP(fStep)         → P += fStep * F
       │   ├── CalcForceOfActions()   → 汇总所有作用量力
       │   │   ├── CActionGauge::CalculateForce()
       │   │   ├── CActionFermion::CalculateForce() → 调用求解器
       │   │   └── CActionPhi4::CalculateForce()
       │   └── UpdateU(fStep)         → U = exp((fStep / gaugeMomentumFactor) * P) * U
       ├── 6. 计算最终能量
       └── 7. Metropolis 测试
           ├── 接受: 保留新配置
           └── 拒绝: RecoverFields() 恢复备份
```

### 测量流程

```
CMeasurementManager::OnConfigurationAccepted()
├── 对每个注册的 CMeasure:
│   ├── 规范测量: OnConfigurationAcceptedSingleField(pGauge, pBoson)
│   ├── 费米子测量: OnConfigurationAcceptedZ4(pGauge, pFermion)
│   │   └── 生成 Z4 随机源 → 调用求解器 → 累积结果
│   └── 源扫描: SourceScanning(pGauge, pFermion)
└── OnUpdateFinished() → Average() → Report()
```

### 费米子力求解流程

```
CActionFermionKS::CalculateForce()
├── 生成/更新伪费米子场 phi
├── 调用多移求解器: (D^+ D + c_n)^{-1} phi
│   └── CMultiShiftBiCGStab::Solve()
│       └── 每次迭代调用 CFieldFermionKST::ApplyOperator(DD)
└── 从解构造力场
```

## 初始化流程

```
CLGLibManager::Initial()
├── 1. 读取参数文件 (CYAMLParser)
├── 2. 创建 CLatticeData (格点维度、索引、边界条件)
├── 3. 创建 CFieldGauge (规范场)
│   └── 随机初始化或从文件加载
├── 4. 创建 CFieldFermion / CFieldBoson (可选)
├── 5. 创建 CAction 列表
│   └── 每个 Action::Initial() 绑定到场
├── 6. 创建 CUpdator (HMC / Heatbath)
│   └── 创建 CIntegrator (LeapFrog / Omelyan / ...)
│       └── 积分器绑定 Action 列表
├── 7. 创建 CMeasurementManager
│   └── 注册所有 CMeasure 对象
├── 8. 创建 CGaugeSmearing (可选)
├── 9. 创建 CGaugeFixing (可选)
└── 10. 创建 CSLASolver (求解器，用于费米子力)
```

## 跨模块设计模式

### 1. RTTI + 工厂模式

所有可配置对象（Field、Action、Updator、Measure、Solver、GaugeFixing、Smearing）都通过全局工厂创建：

```cpp
// 头文件：声明注册辅助（强制 Linux 静态链接保留构造函数）
__CLG_REGISTER_HELPER_HEADER(CMyClass)

// 头文件：类声明
class CLGAPI CMyClass : public CBase
{
    __CLGDECLARE_CLASS(CMyClass)
    // ...
};

// 实现文件：定义工厂函数并注册到全局类表
__CLGIMPLEMENT_CLASS(CMyClass)

// 创建
CBase* pObj = appCreate(_T("CMyClass"));
```

**机制**：
- `__CLGDECLARE_CLASS` 在类内嵌套定义一个 `classCMyClass` 静态对象，其构造函数调用 `GClassGather.AddClass(this)`，把类名和工厂函数注册到全局链表 `GClassGather.m_pClasses`。
- `__CLG_REGISTER_HELPER_HEADER` 在头文件中定义一个辅助结构体及其全局静态实例，用于解决 Linux 静态链接时未引用对象文件导致注册码被丢弃的问题。
- `CLGLibManager` 通过 `appCreate(sClassName)` 从 `GClassGather` 查找并创建对象。

### 2. 字段池 (Field Pool)

`CLatticeData` 管理同类型场的对象池，避免重复分配 GPU 内存：

```cpp
CField* pCopy = appGetLattice()->GetPooledCopy(pOriginal);
// 使用完后 Return()
pCopy->Return();
```

**用途**：求解器工作向量、HMC 备份场、测量中间场等。

### 3. 模板 Kernel + 显式实例化

所有跨规范群的设备代码通过模板实现：

```cpp
// Kernel 中
template<typename deviceGauge>
__global__ void _kernelFunc(deviceGauge* data) { ... }

// Host 中显式实例化
_LAUNCH_KERNEL(_kernelFunc TMPARG(deviceSU3), block, threads, args);
_LAUNCH_KERNEL(_kernelFunc TMPARG(deviceSU2), block, threads, args);
```

**配合**：`switch (GetFieldType())` 在运行时选择模板实例。

### 4. 静态 Kernel 接口类

将 kernel 调用逻辑从场类中剥离：

```cpp
// CFieldGaugeKernel::CalculateForceAndStaple(...) 
// CFieldFermionKSTKernel::DOperatorKS(...)
// CFieldBosonVNKernel::DOperator(...)
// CFieldGaugeKernel::CalculatePlaqutteEnergy(...)
```

**好处**：场类只负责数据管理和 host 端逻辑，kernel 调度集中在专门的静态类中。

### 5. EFieldOperator 统一分发

所有场（Gauge、Fermion、Boson）的算子应用通过统一枚举分发：

```cpp
enum EFieldOperator {
    EFO_F_D, EFO_F_Ddagger, EFO_F_DD, EFO_F_DDdagger,
    EFO_F_InverseD, EFO_F_InverseDdagger, ...
};

// 场基类
CField::ApplyOperator(uiM, gaugeFields, bosonFields);
// 内部按 uiM 调用 D() / InverseD() / DDdagger() 等
```

**好处**：求解器、作用量、测量器无需知道具体场类型，只需指定算子类型。

### 6. 分层 Action 组合

作用量通过组合而非继承构建：

```cpp
CIntegrator 持有 TArray<CAction*>
├── CActionGaugePlaquette (规范力)
├── CActionFermionKSCombined (多个费米子)
│   └── 内部持有多个 CActionFermionKS
└── CActionPhi4 (玻色子)
```

### 7. 测量回调体系

测量对象按需注册回调类型：

```cpp
// 规范/玻色子测量
OnConfigurationAcceptedSingleField(pGauge, pBoson)

// 费米子随机估计
OnConfigurationAcceptedZ4(pGauge, pFermion)
OnConfigurationAcceptedZ4SingleField(pGauge, pFermion, pZ4)

// 源扫描
SourceScanning(pGauge, pFermion)
```

## 关键依赖关系

| 模块 | 依赖 |
|------|------|
| Core | 无（最底层） |
| Math | Core |
| Platform | Core |
| Lattice | Core + Platform |
| Gauge Field | Core + Math + Lattice |
| Fermion Field | Core + Math + Lattice + Gauge Field |
| Boson Field | Core + Math + Lattice + Gauge Field |
| Action | Core + Lattice + Gauge + Fermion + Boson |
| Update | Core + Action + Gauge + Sparse LA |
| Measurement | Core + Action + Gauge + Fermion + Boson |
| Sparse LA | Core + Fermion + Boson |
| Gauge Fixing | Core + Gauge Field + Math |
| Gauge Smearing | Core + Gauge Field + Math |
| Test | 所有模块 |

## 代码组织原则

1. **头文件即文档**：`.h` 文件包含完整接口，`.cu`/`.cpp` 包含实现
2. **模板声明在 `.h`，显式实例化在 `.cu`**：避免链接问题
3. **设备端代码隔离**：所有 `__global__`/`__device__` 函数在 `.cu` 文件中
4. **宏封装平台差异**：`_LAUNCH_KERNEL`、`preparethread`、`intokernal` 等抽象 CUDA/DCU 差异
5. **精度透明**：`Real` / `CLGComplex` 宏自动映射为 `float`/`cuFloatComplex` 或 `double`/`cuDoubleComplex`
