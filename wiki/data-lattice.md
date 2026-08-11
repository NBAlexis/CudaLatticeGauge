# Lattice 模块

格点索引、边界条件和并行分解基础设施。

## 文件清单

| 文件 | 路径 | 说明 |
|------|------|------|
| `CIndex.h` | `Data/Lattice/CIndex.h` | 抽象格点索引接口 |
| `CIndexSquare.h` | `Data/Lattice/CIndexSquare.h` | 方形格点索引实现 |
| `CIndexData.h` | `Data/Lattice/CIndexData.h` | 索引数据缓存（设备端表） |
| `CLatticeData.h` | `Data/Lattice/CLatticeData.h` | 格点数据管理器（单例） |

## 核心概念

### 索引类型

| 索引 | 类型 | 计算公式 | 说明 |
|------|------|---------|------|
| Site Index | `UINT` | `x*Ly*Lz*Lt + y*Lz*Lt + z*Lt + t` | 格点唯一标识 |
| Link Index | `UINT` | `siteIndex * dir + dir` | 链接唯一标识 |
| Fat Index | `UINT` | `siteIndex*(dir+1) + bSite?0:(dir+1)` | 兼容 site 和 link |
| Data Index | `UINT` | `linkIndex * elementCount + n` | 实际数据偏移（如 SU3: linkIndex*9+n） |
| Big Index | `UINT` | 扩展体积（含边界 ghost 格点） | 用于边界处理 |

### 格点类型

- `EIndexType_Square` — 方形格点（当前唯一实现）
- 未来可能支持三角格点等

## 类层次

```
CBase
└── CIndex (abstract)
    └── CIndexSquare
```

## CIndex（抽象接口）

**文件**: `Data/Lattice/CIndex.h`

所有格点索引实现的基类。核心虚函数：

| 方法 | 说明 |
|------|------|
| `BakeAllIndexBuffer(pData)` | 烘焙所有索引缓冲区 |
| `BakePlaquttes(pData, byFieldId)` | 缓存 plaquette 索引（避免每次 walk） |
| `BakeMoveIndex(pData, byFieldId)` | 缓存 boson/fermion 场邻居索引（只对 boson 和 fermion 场烘焙，gauge/tensor2 等无 hopping 的场不烘焙） |
| `BakeEtaMuTable(pData)` | 预计算 eta_mu 表（Staggered 费米子相位） |
| `BakeNaikTable(pData, byFieldId)` | 缓存 Naik 项索引（改进的费米子作用量） |
| `GetPlaqutteCount(byFieldId)` | 返回 plaquette 数量 |
| `CalculateSiteCount(pData)` | 计算实际格点数（考虑 Dirichlet 边界） |

## CIndexSquare（方形格点实现）

**文件**: `Data/Lattice/CIndexSquare.h`

方形格点的具体实现。支持：
- 4D 方形格点（可配置为 3D=1xLxLxL, 2D=1x1xLxL）
-  even/odd 分解
- 多种边界条件（通过 `CBoundaryCondition`）

**设备端函数**（在 `.cu` 中定义）：

| 函数 | 说明 |
|------|------|
| `_deviceGetBigIndex(sSite, pSmallData)` | 从小坐标计算 big index |
| `_deviceEta2(uiEta, i, j)` | 2D eta 函数（用于 plaquette 相位） |
| `_deviceEta3(sSite, missingDir)` | 3D eta 函数（Staggered 相位） |
| `_deviceEta124(sSite)` | gamma53 等价相位 |

**Even/Odd 枚举**:

```cpp
EIE_All          // 所有格点
EIE_EvenToOdd    // even -> odd
EIE_OddToEven    // odd -> even
EIE_EvenToOddOnEven  // even -> odd，但只在 even 格点上执行
EIE_OddToEvenOnOdd   // odd -> even，但只在 odd 格点上执行
```

## CIndexData（索引数据缓存）

**文件**: `Data/Lattice/CIndexData.h`

存储预计算的索引表（设备端常量内存）：
- `m_pDeviceIndexPositionToSIndex` — 位置到 SIndex 的映射
- `m_pDeviceIndexLinkToSIndex` — 链接索引映射
- `m_pEtaMu` — eta_mu 表
- `m_pPlaqutte` — plaquette 索引表
- `m_pMoveTable` — 费米子移动表
- `m_pNaikTable` — Naik 项表

这些表在初始化时由 `CIndex` 子类烘焙（`BakeAllIndexBuffer` 等），之后设备端 kernel 直接查表。

## CLatticeData（格点数据管理器）

**文件**: `Data/Lattice/CLatticeData.h`

整个模拟的**中心数据容器**，持有所有场、作用量、测量器、求解器等。

**关键成员**:

| 成员 | 类型 | 说明 |
|------|------|------|
| `m_pGaugeField` | `TArray<CFieldGauge*>` | 规范场列表 |
| `m_pBosonField` | `TArray<CFieldBoson*>` | 玻色子场列表 |
| `m_pFermionField` | `TArray<CFieldFermion*>` | 费米子场列表 |
| `m_pActionList` | `TArray<CAction*>` | 作用量列表 |
| `m_pActionMap` | `THashMap<BYTE, CAction*>` | 按 ID 查找作用量 |
| `m_pUpdator` | `CUpdator*` | 更新器（HMC/热浴） |
| `m_pMeasurements` | `CMeasurementManager*` | 测量管理器 |
| `m_pGaugeFixing` | `CGaugeFixing*` | 规范固定器 |
| `m_pGaugeSmearing[]` | `CGaugeSmearing*` | 规范平滑器（按场 ID） |
| `m_pFermionSolver[]` | `CSLASolver*` | 费米子求解器（按场 ID） |
| `m_pStapleCaches[]` | `CStapleCache*` | Staple 缓存（按场 ID） |
| `m_pIndex` | `CIndex*` | 格点索引实现 |
| `m_pIndexCache` | `CIndexData*` | 索引数据缓存 |
| `m_pRandom` | `CRandom*` | 随机数生成器 |

**全局访问**:

```cpp
CLatticeData* appGetLattice();           // 获取格点管理器
CField* appGetLattice()->GetFieldById(id);
CAction* appGetLattice()->GetActionById(id);
```

## 初始化流程

1. `CLGLibManager::Initial()` 读取 YAML 参数
2. 创建 `CIndex` 子类（如 `CIndexSquare`）
3. `CIndex->BakeAllIndexBuffer()` 烘焙索引表
4. 按 YAML 配置创建场（`CFieldGauge`、`CFieldFermion` 等）
5. 按 YAML 配置创建作用量（`CAction` 子类）
6. 创建更新器、测量器、求解器等

## 关键设计模式

- **预计算索引表**: 所有邻居查找、plaquette 遍历在初始化时预计算为设备端查表，避免 kernel 中实时计算坐标
- **Fat Index**: 扩展体积包含 ghost 边界格点，统一处理内部和边界链接
- **Field ID**: 每个场有唯一的 `BYTE byFieldId`，用于在 `THashMap` 中查找

## 边界条件

边界条件由 `CBoundaryCondition` 子类实现，通过 `CIndex::SetBoundaryCondition()` 设置：
- `CBoundaryConditionTorusSquare` — 周期性边界
- `CBoundaryConditionPeriodicAndDirichletSquare` — 混合周期/Dirichlet
- `CBoundaryConditionProjectivePlaneSquare` — 射影平面边界
