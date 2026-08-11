# 4. SU(3) + Wilson Dirac 费米子的 HMC 模拟

加入费米子后，框架需要额外配置：
1. **费米子场** —— Wilson Dirac 费米子
2. **费米子作用量** —— `CActionFermionKS`
3. **线性求解器** —— 用于计算分子动力学中的力
4. **费米子场池（PoolNumber）** —— 费米子算法需要多个辅助场

## 文件组织

```
Code/Applications/MyProject/
└── MySU3Wilson.cpp

Bin/Debug/
└── MySU3Wilson.yaml
```

> **前置检查**：确认 `Code/CLGLib/Core/CLGSetup.h` 中 `_CLG_SU3_GAUGE` 和 `_CLG_WILSON_DIRAC` 宏已设为 `1`。

## 4.1 C++ 代码（`MySU3Wilson.cpp`）

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../Debug/MySU3Wilson.yaml"), params);

    appSetupLog(params);
    appInitialCLG(params);

    // 热化
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(5, FALSE);

    // 正式演化
    appClearProfiler();
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, _T("su3_wilson_conf"));
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->UpdateUntileAccept(20, FALSE);

    appDumpProfiler();
    return 0;
}
```

## 4.2 YAML 配置文件（`Bin/Debug/MySU3Wilson.yaml`）

```yaml
Dim : 4
Dir : 4
LatticeLength : [8, 8, 8, 8]
LatticeIndex : CIndexSquare
LatticeBoundary : CBoundaryConditionTorusSquare
ThreadAutoDecompose : 1
RandomType : ER_Schrage
RandomSeed : 1234567

ActionListLength : 2
GaugeFieldCount : 1
FermionFieldCount : 1
MeasureListLength : 1

Updator:
    UpdatorType : CHMC
    Metropolis : 1
    IntegratorType : CIntegratorLeapFrog
    IntegratorStepLength : 1
    IntegratorStep : 60

Gauge:
    FieldName : CFieldGaugeSU3
    FieldInitialType : EFIT_Random

FermionField1:
    FieldName : CFieldFermionWilsonSquareSU3
    FieldInitialType : EFIT_RandomGaussian
    Hopping : 0.1575       # 跳迁参数 kappa
    FieldId : 2            # 费米子场 ID = 2
    PoolNumber : 26        # 辅助场数量，GMRES 需要较多
    Period : [1, 1, 1, -1] # 时间方向反周期边界（费米子）

Solver:
    SolverName : CSLASolverGMRES
    SolverForFieldId : 2   # 为 FieldId=2 的场求解
    MaxDim : 20
    Accuracy : 0.0001
    Restart : 15
    AbsoluteAccuracy : 1

Action1:
    ActionName : CActionGaugePlaquette
    Beta : 5.5

Action2:
    ActionName : CActionFermionKS
    FieldId : 2

Measure1:
    MeasureName : CMeasurePlaqutteEnergy
```

## 4.3 关键参数详解

- **`Hopping : 0.1575`** —— Wilson 费米子的跳迁参数（`kappa` 或 `hopping`），典型值 0.15 ~ 0.16。
- **`FieldId : 2`** —— 第一个规范场自动占用 ID 1，所以费米子场必须从 2 开始。
- **`PoolNumber : 26`** —— 费米子 HMC 需要大量辅助场用于伪费米子技巧（pseudo-fermion）。对于 GMRES 求解器，经验公式：
  - PoolNumber = 2（phi 左右手） + 3（GMRES 的 x, r, w） + 1（DD^dagger 求解额外需要） + MaxDim
  - 即 `PoolNumber = 6 + MaxDim = 6 + 20 = 26`
- **`Period : [1, 1, 1, -1]`** —— 费米子要求时间方向反周期边界条件（`T` 方向最后一个分量设为 `-1`）。
- **`Solver`** —— 费米子作用量 `CActionFermionKS` 需要求解 `D^dagger D psi = phi`。支持的求解器：
  - `CSLASolverGMRES` —— 通用，Restart 参数影响效率
  - `CSLASolverBiCGstab` —— 双共轭梯度稳定化
  - `CSLASolverGCRODR` —— 带收缩的 GMRES，适合多右端项

---

[< 返回目录](Tutorial.md) | [上一章：纯 Z2 规范场](03-pure-z2.md) | [下一章：SU(3) + Staggered (KS)](05-su3-ks.md)
