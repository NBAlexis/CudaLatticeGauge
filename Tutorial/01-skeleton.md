# 1. CLGLib 程序的基本骨架

CLGLib 的核心设计理念是：**所有物理对象（格点、规范场、费米子场、作用量、测量器、更新器）都通过 YAML 配置文件描述，C++ 代码只负责控制流程。**

## 文件组织

在开始之前，先明确文件放在哪里。CLGLib 的目录结构约定如下：

- **C++ 源文件**：放在 `Code/Applications/你的项目名/` 下（例如 `Code/Applications/MyProject/MyPureGaugeSU3.cpp`）
- **YAML 配置文件**：放在 `Bin/Debug/` 或 `Bin/Ubuntu/` 下（与编译输出的可执行文件同目录或其父目录）
- **`../Debug/xxx.yaml`** 这个相对路径是**相对于可执行文件运行时的工作目录**的。如果你的可执行文件在 `Bin/Ubuntu/`，那么运行时工作目录也应该是 `Bin/Ubuntu/`，这样 `../Debug/xxx.yaml` 才能正确找到 `Bin/Debug/xxx.yaml`

> **提示**：`CLGExample` 的完整参考实现位于 `Code/Applications/CLGExample/CLGExample.cpp`，对应的 YAML 在 `Bin/Debug/CLGExample.yaml`。可以先编译并跑通 CLGExample，确认环境没问题后再写自己的项目。

---

## 最小程序示例

参考 `Code/Applications/CLGExample/CLGExample.cpp`，一个最小的 CLGLib 程序只有不到 40 行：

```cpp
#include "CLGLib.h"

int main(int argc, char * argv[])
{
    CParameters params;
    // 第 1 步：读取 YAML 参数文件
#if _CLG_DEBUG && _CLG_WIN
    CYAMLParser::ParseFile(_T("CLGExample.yaml"), params);
#else
    CYAMLParser::ParseFile(_T("../Debug/CLGExample.yaml"), params);
#endif

    // 第 2 步：初始化日志和 CLG 框架
    appSetupLog(params);
    appInitialCLG(params);

    // 第 3 步：热化（warmup），不做测量，不保存组态
    appGetLattice()->m_pUpdator->SetTestHdiff(TRUE);
    appGetLattice()->m_pUpdator->SetAutoCorrection(FALSE);
    appGetLattice()->m_pUpdator->Update(5, FALSE);

    // 第 4 步：正式演化，每步测量，每接受一次保存组态
    appClearProfiler();
    appGetLattice()->m_pUpdator->SetSaveConfiguration(TRUE, _T("test"));
    appGetLattice()->m_pUpdator->SetAutoCorrection(TRUE);
    appGetLattice()->m_pUpdator->UpdateUntileAccept(20, FALSE);

    appDumpProfiler();
    return 0;
}
```

> **关于代码中的一些符号**：
> - `_T("...")` 是 CLGLib 的跨平台字符串宏，在 Windows 下展开为宽字符字符串 `L"..."`，在 Linux 下为普通字符串 `"..."`。所有涉及字符串的 CLGLib API 都要求用 `_T()` 包裹。
> - `CCString` 是 CLGLib 提供的字符串类（功能类似 `CString`），支持 `Format()` 等格式化方法。它不是 `std::string`，而是框架内部使用的类型。

## 核心 API 速览

| API | 作用 |
|-----|------|
| `CYAMLParser::ParseFile(path, params)` | 读取 YAML 配置文件 |
| `appSetupLog(params)` | 根据参数初始化日志系统 |
| `appInitialCLG(params)` | **核心**：根据参数创建格点、场、作用量、测量器、更新器 |
| `appGetLattice()` | 获取全局格点对象 `CLatticeData*`，所有物理对象都挂在它下面 |
| `appGetLattice()->m_pUpdator` | 获取更新器（HMC / Heatbath 等） |
| `m_pUpdator->Update(n, bMeasure)` | 执行 n 步更新，`bMeasure` 决定是否调用测量器 |
| `m_pUpdator->UpdateUntileAccept(n, bMeasure)` | 执行 n 次**被接受**的更新（跳过被拒绝的） |
| `m_pUpdator->SetSaveConfiguration(bSave, prefix)` | 开启/关闭组态自动保存，文件名为 `prefix_1.con`、`prefix_2.con`... |
| `appDumpProfiler()` | 输出性能分析结果 |

## YAML 文件的通用结构

```yaml
# ============ 格点基本参数 ============
Dim : 4                 # 时空维度
Dir : 4                 # 规范场方向数（通常=Dim）
LatticeLength : [8, 8, 8, 8]   # 各方向格点数
LatticeIndex : CIndexSquare    # 索引类型
LatticeBoundary : CBoundaryConditionTorusSquare  # 边界条件（周期边界）
ThreadAutoDecompose : 1        # 自动分解 CUDA block/thread
RandomType : ER_Schrage        # 随机数生成器
RandomSeed : 1234567           # 随机数种子

# ============ 场、作用量、测量器的数量声明 ============
ActionListLength : 1           # 作用量个数
GaugeFieldCount : 1            # 规范场个数
FermionFieldCount : 0          # 费米子场个数
MeasureListLength : 1          # 测量器个数

# ============ 更新器配置 ============
Updator:
    UpdatorType : CHMC
    Metropolis : 1             # 是否做 Metropolis 检验
    IntegratorType : CIntegratorLeapFrog
    IntegratorStepLength : 1
    IntegratorStep : 60

# ============ 规范场配置 ============
Gauge:
    FieldName : CFieldGaugeSU3
    FieldInitialType : EFIT_Random   # 初始化为随机 SU3 矩阵

# ============ 作用量配置 ============
Action1:
    ActionName : CActionGaugePlaquette
    Beta : 5.5

# ============ 测量器配置 ============
Measure1:
    MeasureName : CMeasurePlaqutteEnergy
```

**关键约定**：
- 场、作用量、测量器都按 `1, 2, 3...` 编号。如果有多个规范场，依次写 `Gauge`、`Gauge2`、`Gauge3`...
- `FieldId` 是场在框架内部的唯一标识，**第一个规范场默认为 1**，后续场（包括费米子场）需要显式指定不冲突的 ID。
- `GaugeFields` / `FieldId` 在 Action 和 Measure 中用于指定该对象作用于哪个场。

## CLGExample.yaml 参考

下面给出 `Bin/Debug/CLGExample.yaml` 的完整内容，供你对比参考。它包含了 SU(3) 规范场 + Wilson Clover 费米子 + plaquette 作用量 + 费米子作用量的完整配置：

```yaml
# VerboseLevel : DETAILED
VerboseLevel : PARANOIAC
VerboseOutput : datetime
Dim : 4
Dir : 4
LatticeLength : [12, 12, 12, 12]
LatticeIndex : CIndexSquare
LatticeBoundary : CBoundaryConditionTorusSquare
ThreadAutoDecompose : 1
RandomType : ER_Schrage
RandomSeed : 1234567
ActionListLength : 2
GaugeFieldCount : 1
FermionFieldCount : 2
MeasureListLength : 1
CacheSolution : 0
MaxThreadPerBlock : 256
ExponentialPrecision : 8
Profiler : 1

Updator:
    UpdatorType : CHMC
    Metropolis : 1
    IntegratorType : CIntegratorNestedForceGradient
    IntegratorStepLength : 1
    IntegratorStep : 10
    NestedStep : 20

Gauge:
    FieldName : CFieldGaugeSU3TreeImproved
    FieldInitialType : EFIT_Random
    RectOverPlaq : -0.05

FermionField1:
    FieldName : CFieldFermionWilsonSquareCloverSU3
    FieldInitialType : EFIT_RandomGaussian
    Hopping : 0.1575
    FieldId : 2
    PoolNumber : 26
    Period : [1, 1, 1, -1]

Solver:
    SolverName : CSLASolverGMRES
    SolverForFieldId : 2
    MaxDim : 50
    Accuracy : 0.000000001
    Restart : 50
    AbsoluteAccuracy : 1

Action1:
    ActionName : CActionGaugePlaquette
    Beta : 7.0

Action2:
    ActionName : CActionFermionKS
    FieldId : 2

Measure1:
    MeasureName : CMeasurePlaqutteEnergy

StapleCache:
    StapleCacheName : CStapleCacheSU3
    FieldId : 1
```

---

[< 返回目录](Tutorial.md) | [下一章：纯 SU(3) 规范场 HMC >](02-pure-su3.md)
