# 3. 纯 Z2 规范场的 Heatbath 模拟

Z2（或更一般的 Z_N）是**离散规范群**，每个链接变量只能取 +1 或 -1。CLGLib 对离散群使用专门的离散更新器，而不是 HMC（HMC 要求连续群）。这里我们用 **Heatbath** 更新器。

> 注意：离散规范群目前**不支持 HMC**。热浴（Heatbath）是最自然的更新方式。

## 文件组织

```
Code/Applications/MyProject/
└── MyPureGaugeZ2.cpp

Bin/Debug/
└── MyPureGaugeZ2.yaml
```

> **前置检查**：使用离散规范群之前，请确认 `Code/CLGLib/Core/CLGSetup.h` 中 `_CLG_Z2_GAUGE` 宏已设为 `1`。否则运行时会报 "Unable to create the gauge field!"

## 3.1 C++ 代码（`MyPureGaugeZ2.cpp`）

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../Debug/MyPureGaugeZ2.yaml"), params);

    appSetupLog(params);
    appInitialCLG(params);

    // 热化
    appGetLattice()->m_pUpdator->Update(10, FALSE);

    // 正式演化：Heatbath 每一步都接受，不需要 UpdateUntileAccept
    appClearProfiler();
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, _T("z2_conf"));
    appGetLattice()->m_pUpdator->Update(100, FALSE);

    appDumpProfiler();
    return 0;
}
```

## 3.2 YAML 配置文件（`Bin/Debug/MyPureGaugeZ2.yaml`）

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
    UpdatorType : CHeatbath

Gauge:
    FieldName : CFieldGaugeZ2
    FieldInitialType : EFIT_Random

Action1:
    ActionName : CActionDiscreteZNPlaquette
    Beta : 0.5
    GaugeFields : [1]

Measure1:
    MeasureName : CMeasurePlaqutteEnergy
    GaugeFields : [1]
```

## 3.3 与连续群的关键区别

| 对比项 | SU(3) 连续群 | Z2 离散群 |
|--------|-------------|----------|
| 场类名 | `CFieldGaugeSU3` | `CFieldGaugeZ2` |
| 作用量类名 | `CActionGaugePlaquette` | `CActionDiscreteZNPlaquette` |
| 更新器 | `CHMC` | `CHeatbath` |
| Beta 范围 | 5.5 ~ 7.0 | 0.1 ~ 1.0（小得多） |
| `GaugeFields` 显式指定 | 可选（默认作用于所有 gauge） | **建议显式指定** |

> Z3、Z4、Z5... 只需把 `CFieldGaugeZ2` 改为 `CFieldGaugeZ3` 等即可。D_N（二面体群）系列用 `CFieldGaugeD3`、`CActionDiscreteDNPlaquette`。

> **注意**：`CHeatbath` 更新器的 `Update` 方法内部**不会**调用 `SaveConfiguration`。`SetSaveConfiguration(TRUE, ...)` 这行代码对 Heatbath 实际上不起作用（保留它不会导致错误，只是无效果）。如需保存组态，你需要在每次 `Update` 后手动调用：
> ```cpp
> appGetLattice()->m_pGaugeField[0]->SaveToFile(sFileName);
> ```

---

[< 返回目录](Tutorial.md) | [上一章：纯 SU(3) 规范场](02-pure-su3.md) | [下一章：SU(3) + Wilson Dirac](04-su3-wilson.md)
