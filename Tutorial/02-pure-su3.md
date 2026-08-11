# 2. 纯 SU(3) 规范场的 HMC 模拟

这是最简单的完整模拟：一个 SU(3) 规范场，一个 plaquette 作用量，用 HMC 更新。

## 文件组织

```
Code/Applications/MyProject/
└── MyPureGaugeSU3.cpp

Bin/Debug/
└── MyPureGaugeSU3.yaml
```

> **前置检查**：确认 `Code/CLGLib/Core/CLGSetup.h` 中 `_CLG_SU3_GAUGE` 宏已设为 `1`。否则运行时会报 "Unable to create the gauge field!"

## 2.1 C++ 代码（`MyPureGaugeSU3.cpp`）

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../Debug/MyPureGaugeSU3.yaml"), params);

    appSetupLog(params);
    appInitialCLG(params);

    // 热化 10 步，不保存，不测能量守恒误差
    appGetLattice()->m_pUpdator->SetTestHdiff(FALSE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(10, FALSE);

    // 正式演化 50 步，保存组态，自动修正步长
    // UpdateUntileAccept 的第二个参数 bMeasure=FALSE 表示更新过程中不调用测量器
    // （测量通常在演化完成后单独做，见第 6 章）
    appClearProfiler();
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, _T("su3_conf"));
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->UpdateUntileAccept(50, FALSE);

    appDumpProfiler();
    return 0;
}
```

## 2.2 YAML 配置文件（`Bin/Debug/MyPureGaugeSU3.yaml`）

```yaml
Dim : 4
Dir : 4
LatticeLength : [8, 8, 8, 8]
LatticeIndex : CIndexSquare
LatticeBoundary : CBoundaryConditionTorusSquare
ThreadAutoDecompose : 1
RandomType : ER_Schrage
RandomSeed : 1234567

ActionListLength : 1
GaugeFieldCount : 1
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

Action1:
    ActionName : CActionGaugePlaquette
    Beta : 5.5

Measure1:
    MeasureName : CMeasurePlaqutteEnergy
```

## 2.3 参数解释

- **`Beta : 5.5`** —— 耦合常数 `beta = 6/g^2`。纯规范场常用值范围：5.5 ~ 7.0。
- **`IntegratorStep : 60`** —— 分子动力学轨迹被分为 60 步。步数越多能量守恒越好但越慢。
- **`FieldInitialType : EFIT_Random`** —— 初始化为随机 SU(3) 矩阵（热启动）。也可改为 `EFIT_Identity`（冷启动，全为单位矩阵）。
- **`UpdateUntileAccept(50, FALSE)`** —— 演化直到产生 50 个**被 Metropolis 接受**的组态。实际执行的更新步数可能大于 50（因为有些会被拒绝）。

运行后会在可执行文件同目录下生成 `su3_conf_1.con`、`su3_conf_2.con`... 等组态文件。

---

[< 返回目录](Tutorial.md) | [上一章：基本骨架](01-skeleton.md) | [下一章：纯 Z2 规范场](03-pure-z2.md)
