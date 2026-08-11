# 6. 对已保存组态进行测量

模拟完成后，你可能已经有了大量保存的组态文件（如 `su3_conf_1.con`、`su3_conf_2.con`...）。

## 文件组织

```
Code/Applications/MyProject/
├── MeasurePolyakov.cpp
├── MeasureWilson.cpp
└── MeasureChiral.cpp

Bin/Debug/
├── MeasurePolyakov.yaml
├── MeasureWilson.yaml
└── MeasureChiral.yaml
```

> **运行时工作目录**：确保在可执行文件所在目录（如 `Bin/Ubuntu/`）运行，这样代码中的 `../Debug/xxx.yaml` 和 `su3_conf_%d.con` 路径才能正确解析。

测量程序与模拟程序的区别：
1. **不配置 Updator**（不需要更新）
2. 规范场的 `FieldInitialType` 设为 **`EFIT_ReadFromFile`**
3. 在 C++ 代码中手动读取每个组态文件，执行测量

## 6.0 测量程序的通用框架

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../Debug/MyMeasurement.yaml"), params);
    appSetupLog(params);

    // 读取测量参数
    INT iStart = 1;
    params.FetchValueINT(_T("StartN"), iStart);
    INT iEnd = 100;
    params.FetchValueINT(_T("EndN"), iEnd);

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    // 获取测量器
    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(
        appGetLattice()->m_pMeasurements->GetMeasureById(1));

    // 遍历所有组态文件
    for (INT i = iStart; i <= iEnd; ++i)
    {
        CCString sFileName;
        sFileName.Format(_T("su3_conf_%d.con"), i);

        // 加载组态
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFileName);

        // 执行测量
        appGetLattice()->m_pMeasurements->OnConfigurationAccepted(
            1, 0,  // 1 个 gauge，0 个 boson
            appGetLattice()->m_pGaugeField,
            NULL,
            appGetLattice()->m_pGaugeField);  // staple 参数，测量时不需要可用 NULL，但这里传入 gauge 自身作为占位

        appGeneral(_T("Measured config %d\n"), i);
    }

    // 输出最终统计结果
    appGetLattice()->m_pMeasurements->Report();

    return 0;
}
```

> 注意：`OnConfigurationAccepted` 的参数较复杂。更简洁的方式是直接调用各个测量器的特定接口，见下节。

---

## 6.1 Polyakov Loop 测量

Polyakov loop 是探测禁闭-退禁闭相变的关键序参量：

$$P(\vec{x}) = \frac{1}{N_c} \text{tr} \prod_{t=0}^{N_t-1} U_t(\vec{x}, t)$$

### YAML 配置

```yaml
Dim : 4
Dir : 4
LatticeLength : [8, 8, 8, 8]
LatticeIndex : CIndexSquare
LatticeBoundary : CBoundaryConditionTorusSquare
ThreadAutoDecompose : 1
RandomType : ER_Schrage
RandomSeed : 1234567

GaugeFieldCount : 1
MeasureListLength : 1
ActionListLength : 0   # 测量不需要作用量

# 注意：测量程序不需要 Updator 配置

Gauge:
    FieldName : CFieldGaugeSU3
    FieldInitialType : EFIT_ReadFromFile
    GaugeFileType : EFFT_CLGBin
    GaugeFileName : ../Debug/su3_conf_1.con

Measure1:
    MeasureName : CMeasurePolyakovXY
    FieldId : 1
    ShowResult : 1
    MeasureDist : 1        # 测量距离分布
```

### C++ 测量代码

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../Debug/MeasurePolyakov.yaml"), params);
    appSetupLog(params);

    INT iStart = 1, iEnd = 50;
    params.FetchValueINT(_T("StartN"), iStart);
    params.FetchValueINT(_T("EndN"), iEnd);

    if (!appInitialCLG(params))
    {
        appCrucial(_T("Initial Failed!\n"));
        return 1;
    }

    CMeasurePolyakovXY* pPL = dynamic_cast<CMeasurePolyakovXY*>(
        appGetLattice()->m_pMeasurements->GetMeasureById(1));

    for (INT i = iStart; i <= iEnd; ++i)
    {
        CCString sFile;
        sFile.Format(_T("su3_conf_%d.con"), i);
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFile);

        // 手动调用 Polyakov 测量
        pPL->OnConfigurationAcceptedSingleField(
            appGetLattice()->m_pGaugeField[0],
            appGetLattice()->m_pGaugeField[0]);  // staple 参数，这里用自身占位

        // 输出当前结果
        if (pPL->m_lstLoop.Num() > 0)
        {
            cuDoubleComplex last = pPL->m_lstLoop[pPL->m_lstLoop.Num() - 1];
            appGeneral(_T("Config %d: Polyakov = %f + %fi\n"),
                       i, last.x, last.y);
        }
    }

    // 最终报告（含统计误差）
    pPL->Report();

    // 导出 CSV（可选）
    pPL->Export(_T("polyakov.csv"), iStart, iEnd, 0, 0);

    return 0;
}
```

### `CMeasurePolyakovXY` 的关键数据成员

| 成员 | 类型 | 含义 |
|------|------|------|
| `m_lstLoop` | `TArray<cuDoubleComplex>` | 每次测量的 Polyakov loop 平均值序列 |
| `m_lstLoopAbs` | `TArray<DOUBLE>` | `|P|` 的序列 |
| `m_lstP` | `TArray<cuDoubleComplex>` | 按距离 `r` 分布的 Polyakov loop |
| `m_lstPAbs` | `TArray<DOUBLE>` | 按距离 `r` 分布的 `|P|` |
| `m_cAverageLoop` | `CLGComplex` | 所有组态的平均值 |

---

## 6.2 Wilson Loop 测量

Wilson loop 用于计算静态夸克势：

$$W(R, T) = \left\langle \text{tr} \prod_{\mathcal{C}} U \right\rangle$$

其中 `C` 是 `R x T` 的矩形路径。

### YAML 配置

```yaml
Dim : 4
Dir : 4
LatticeLength : [8, 8, 8, 8]
LatticeIndex : CIndexSquare
LatticeBoundary : CBoundaryConditionTorusSquare
ThreadAutoDecompose : 1
RandomType : ER_Schrage
RandomSeed : 1234567

GaugeFieldCount : 1
MeasureListLength : 1
ActionListLength : 0

Gauge:
    FieldName : CFieldGaugeSU3
    FieldInitialType : EFIT_ReadFromFile
    GaugeFileType : EFFT_CLGBin
    GaugeFileName : ../Debug/su3_conf_1.con

Measure1:
    MeasureName : CMeasureWilsonLoop
```

### C++ 测量代码

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../Debug/MeasureWilson.yaml"), params);
    appSetupLog(params);

    INT iStart = 1, iEnd = 50;
    params.FetchValueINT(_T("StartN"), iStart);
    params.FetchValueINT(_T("EndN"), iEnd);

    if (!appInitialCLG(params))
        return 1;

    CMeasureWilsonLoop* pWL = dynamic_cast<CMeasureWilsonLoop*>(
        appGetLattice()->m_pMeasurements->GetMeasureById(1));

    for (INT i = iStart; i <= iEnd; ++i)
    {
        CCString sFile;
        sFile.Format(_T("su3_conf_%d.con"), i);
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFile);

        pWL->OnConfigurationAcceptedSingleField(
            appGetLattice()->m_pGaugeField[0],
            appGetLattice()->m_pGaugeField[0]);

        appGeneral(_T("Config %d measured\n"), i);
    }

    pWL->Report();

    // 输出 W(R,T) 到 CSV
    // m_lstAverageC[conf][r] 是对应不同 R 的平均关联函数
    for (INT rIdx = 0; rIdx < pWL->m_lstR.Num(); ++rIdx)
    {
        appGeneral(_T("R = %d:"), pWL->m_lstR[rIdx]);
        for (INT c = 0; c < pWL->m_lstAverageC.Num(); ++c)
        {
            CLGComplex v = pWL->m_lstAverageC[c][rIdx];
            appGeneral(_T(" %f"), v.x);
        }
        appGeneral(_T("\n"));
    }

    return 0;
}
```

### `CMeasureWilsonLoop` 的关键数据成员

| 成员 | 类型 | 含义 |
|------|------|------|
| `m_lstR` | `TArray<UINT>` | 测量的空间距离 `R` 列表 |
| `m_lstAverageC` | `TArray<TArray<CLGComplex>>` | `m_lstAverageC[config_index][r_index]` 为 `W(R,T)` |

> **提示**：Wilson loop 的测量通常需要先对规范场做**规范平滑（smearing）**以获得更好的信噪比。可以在 YAML 中加入 `GaugeSmearing` 配置（参见 `TestSuit_Updator.yaml` 的 `TestWilsonLoop` 示例）。

---

## 6.3 Chiral Condensate 测量

手征凝聚（Chiral Condensate）`\langle\psi\bar{\psi}\rangle` 是探测手征对称性破缺的关键量。CLGLib 使用 **随机源（stochastic source / Z4 噪声源）** 方法计算：

$$\langle \bar{\psi} \psi \rangle = -\frac{1}{V} \text{tr} D^{-1} \approx -\frac{1}{VN_{\text{src}}} \sum_{i=1}^{N_{\text{src}}} \eta_i^\dagger D^{-1} \eta_i$$

### 两种测量器

| 费米子类型 | 测量器类名 |
|-----------|-----------|
| Wilson Dirac | `CMeasureChiralCondensate` |
| Staggered (KS) | `CMeasureChiralCondensateKS` |

### YAML 配置（Wilson Dirac 版）

```yaml
Dim : 4
Dir : 4
LatticeLength : [8, 8, 8, 8]
LatticeIndex : CIndexSquare
LatticeBoundary : CBoundaryConditionTorusSquare
ThreadAutoDecompose : 1
RandomType : ER_Schrage
RandomSeed : 1234567

GaugeFieldCount : 1
FermionFieldCount : 1
MeasureListLength : 1
ActionListLength : 2   # 需要定义费米子作用量，测量代码中 GetActionById(2) 会用到

Gauge:
    FieldName : CFieldGaugeSU3
    FieldInitialType : EFIT_ReadFromFile
    GaugeFileType : EFFT_CLGBin
    GaugeFileName : ../Debug/su3_conf_1.con

# 测量时也需要费米子场定义（但不需要 PoolNumber）
FermionField1:
    FieldName : CFieldFermionWilsonSquareSU3
    FieldInitialType : EFIT_Zero
    Hopping : 0.1575
    FieldId : 2

Solver:
    SolverName : CSLASolverGMRES
    SolverForFieldId : 2
    MaxDim : 20
    Accuracy : 0.0001
    Restart : 15
    AbsoluteAccuracy : 1

Action1:
    ActionName : CActionGaugePlaquette
    Beta : 5.5
    GaugeFields : {1}

Action2:
    ActionName : CActionFermionKS
    FieldId : 2
    GaugeFields : {1}

Measure1:
    MeasureName : CMeasureChiralCondensate
    FieldId : 2          # 测量哪个费米子场
    FieldCount : 4       # Z4 随机源的个数（越多越精确）
    ShowResult : 1
    MeasureDist : 1
```

### C++ 测量代码

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    CYAMLParser::ParseFile(_T("../Debug/MeasureChiral.yaml"), params);
    appSetupLog(params);

    INT iStart = 1, iEnd = 50;
    params.FetchValueINT(_T("StartN"), iStart);
    params.FetchValueINT(_T("EndN"), iEnd);

    // StochasticFieldCount：每个组态使用多少个随机源
    INT iFieldCount = 10;
    params.FetchValueINT(_T("StochasticFieldCount"), iFieldCount);

    if (!appInitialCLG(params))
        return 1;

    CMeasureChiralCondensate* pCC = dynamic_cast<CMeasureChiralCondensate*>(
        appGetLattice()->m_pMeasurements->GetMeasureById(1));

    // 获取费米子场（用于随机源）
    CFieldFermionWilsonSquareSU3* pF1 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(
        appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));
    CFieldFermionWilsonSquareSU3* pF2 = dynamic_cast<CFieldFermionWilsonSquareSU3*>(
        appGetLattice()->GetPooledFieldById(2, _T(__FILE__), __LINE__));

    for (INT i = iStart; i <= iEnd; ++i)
    {
        CCString sFile;
        sFile.Format(_T("su3_conf_%d.con"), i);
        appGetLattice()->m_pGaugeField[0]->InitialFieldWithFile(sFile);

        // ChiralCondensate 使用 Z4 噪声源，需要多次调用
        for (INT src = 0; src < iFieldCount; ++src)
        {
            // 生成 Z4 随机源到 pF1
            pF1->InitialField(EFIT_RandomZ4);

            // 计算 D^{-1} * source -> pF2
            appGetLattice()->GetActionById(2)->PrepareForHMC(
                appGetLattice()->m_pGaugeField[0], pF1, pF2);

            // 调用测量（bStart 和 bEnd 标记首尾）
            pCC->OnConfigurationAcceptedZ4SingleField(
                appGetLattice()->m_pGaugeField[0],
                appGetLattice()->m_pGaugeField[0],
                pF1, pF2,
                0 == src,                    // bStart
                src == iFieldCount - 1);     // bEnd
        }

        appGeneral(_T("Config %d measured with %d sources\n"), i, iFieldCount);
    }

    pCC->Report();

    // 释放池化场
    pF1->Return();
    pF2->Return();

    return 0;
}
```

### `CMeasureChiralCondensate` 的关键数据成员

| 成员 | 类型 | 含义 |
|------|------|------|
| `m_lstCondAll[9]` | `TArray<CLGComplex>[9]` | 9 种不同 gamma 结构的凝聚，All 是所有源的累计 |
| `m_lstCond[9]` | `TArray<CLGComplex>[9]` | 同上，但每次测量后清空（当前组态） |

> **注意**：上述手动调用 `OnConfigurationAcceptedZ4SingleField` 的方式比较复杂。更简单的做法是让框架自动处理随机源——只需在 YAML 中配置好 `FieldCount`，然后在 C++ 中调用 `appGetLattice()->m_pMeasurements->OnConfigurationAccepted(...)` 即可。参见 `Code/Applications/ConstAcc/Measure.cpp` 的完整实现。

---

[< 返回目录](Tutorial.md) | [上一章：SU(3) + Staggered (KS)](05-su3-ks.md) | [下一章：CMake](07-cmake.md)
