# 5. SU(3) + Staggered (KS) 费米子的 HMC 模拟

Staggered 费米子（Kogut-Susskind / KS）是格点 QCD 中常用的费米子离散化方案，相比 Wilson 费米子保留了更多手征对称性。

## 文件组织

```
Code/Applications/MyProject/
└── MySU3KS.cpp

Bin/Debug/
└── MySU3KS.yaml
```

> **前置检查**：确认 `Code/CLGLib/Core/CLGSetup.h` 中 `_CLG_SU3_GAUGE` 和 `_CLG_STAGGERED_DIRAC` 宏已设为 `1`。

## 5.1 C++ 代码（`MySU3KS.cpp`）

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../Debug/MySU3KS.yaml"), params);

    appSetupLog(params);
    appInitialCLG(params);

    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(5, FALSE);

    appClearProfiler();
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, _T("su3_ks_conf"));
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->UpdateUntileAccept(20, FALSE);

    appDumpProfiler();
    return 0;
}
```

## 5.2 YAML 配置文件（`Bin/Debug/MySU3KS.yaml`）

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
    IntegratorStep : 80

Gauge:
    FieldName : CFieldGaugeSU3
    FieldInitialType : EFIT_Random

FermionField1:
    FieldName : CFieldFermionKSSU3
    FieldInitialType : EFIT_RandomGaussian
    Mass : 0.1             # Staggered 费米子用质量参数
    FieldId : 2
    PoolNumber : 15
    Period : [1, 1, 1, -1]
    # 有理近似系数，分别用于不同阶段的分数幂计算
    MC : [1.53128, -0.000947, -0.022930, -1.192485, 0.005144, 0.075515, 1.338794]
    MD : [0.390460, 0.051109, 0.140828, 0.596484, 0.001277, 0.028616, 0.410599]
    EN : [0.653047, 0.008528, 0.051543, 0.458672, 0.002240, 0.039726, 0.583143]

Solver:
    SolverName : CSLASolverGMRES
    SolverForFieldId : 2
    MaxDim : 20
    Accuracy : 0.0001
    Restart : 15
    AbsoluteAccuracy : 1

MSSolver:
    SolverName : CMultiShiftBiCGStab
    SolverForFieldId : 2
    DiviationStep : 100
    MaxStep : 50
    Accuracy : 0.001
    AbsoluteAccuracy : 1

Action1:
    ActionName : CActionGaugePlaquette
    Beta : 5.0

Action2:
    ActionName : CActionFermionKS
    FieldId : 2

Measure1:
    MeasureName : CMeasurePlaqutteEnergy
```

## 5.3 Staggered 特有的参数

- **`FieldName : CFieldFermionKSSU3`** —— 标准 Staggered 费米子。带改进的有一系列变体如 `CFieldFermionKSSU3Acc`、`CFieldFermionKSHISQ` 等。
- **`Mass : 0.1`** —— 裸质量（不是 hopping）。典型值 0.01 ~ 0.5。
- **`MC / MD / EN`** —— **有理近似系数**，分别对应：
  - `MC` —— 用于 MC 步（preparation），计算 `M^(+1/8)` 或 `M^(+Nf/4)`
  - `MD` —— 用于 MD 步（力计算），计算 `M^(-1/4)` 或 `M^(-Nf/4)`
  - `EN` —— 用于能量计算，计算 `M^(-1/4)`

  这些系数需要通过外部有理近似程序生成，匹配你的质量范围和精度要求。
- **`MSSolver`** —— 多移求解器（Multi-Shift Solver），用于 Staggered 费米子有理近似中同时求解多个偏移的线性系统。`CMultiShiftBiCGStab` 是最常用的选择。

---

[< 返回目录](Tutorial.md) | [上一章：SU(3) + Wilson Dirac](04-su3-wilson.md) | [下一章：测量](06-measurements.md)
