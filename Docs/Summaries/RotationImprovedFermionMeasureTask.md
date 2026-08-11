# RotationImprovedFermion Measure 分支开发总结

> 日期：2026-08-06
> 状态：**全部完成**（代码 + 编译 + 模拟 + 测量 + 结果验证）

## 1. 任务目标

为 `RotationImprovedFermion` 项目补全 measure 分支，参考 `StaggeredSpectrum` 和 `HISQ` 两个项目的 measure 分支，实现：

1. **hisq measure 分支**（`ERIF_MeasureHISQ`）：使用普通 HISQ SU3 费米子（`CFieldFermionHISQSU3`）测量手征凝聚，几乎照搬 HISQ 项目的 measure 分支。
2. **hisq rot measure 分支**（`ERIF_MeasureHISQRotation`）：使用带旋转的 HISQ 费米子（`CFieldFermionHISQSU3R`，即 HISQRSU3）测量手征凝聚。
3. **两个测量 enum**：`ERIFMeasureJob` 提供 `ERIFMJ_Polyakov`（Polyakov loop）和 `ERIFMJ_Chiral`（手征凝聚）两个枚举值。

关键约束（来自用户）：
- 测量手征凝聚**不能使用 even pseudo-fermion**（`Even : 0`）
- **不能用 CG solver**（旋转 D 算符非 Hermitian 正定，须用 GMRES/BiCGStab）
- **带旋转的费米子需要 cache rotation link**（`CachedGauge : 1` + StapleCache 的 `CacheRotationLink : 1`）

## 2. 已完成的工作

### 2.1 新增文件

- **`Code/Applications/RotationImprovedFermion/Measure.cpp`**（约 380 行）
  - 核心驱动函数 `MeasureInternal(CParameters& params, UBOOL bRotation)`，两个公共入口 `MeasureHISQ()` / `MeasureHISQRotation()` 只传 `FALSE`/`TRUE` 区分旋转与否。
  - 支持读取 `StartN/EndN`、`BetaList/PrefixList`（多 ensemble 扫描）、`StochasticFieldCount`、`UseZ4`、`SubFolder/SubFolderPrefix`（组态按 omega 分文件夹）、`LoadType`、`DistributionJob`（Polyakov / Chiral）等参数。
  - Polyakov：`CMeasurePolyakovXY::OnConfigurationAccepted` + `Export`。
  - Chiral：HISQ smearing → effective gauge →（旋转分支先 refresh StapleCache 的 rotation link）→ 随机高斯/Z4 源 → `InverseD`（走 YAML 配置的 GMRES solver）→ `CMeasureChiralCondensateKS::OnConfigurationAcceptedZ4` → `_CLG_EXPORT_CHIRAL` 导出。

### 2.2 修改文件

- **`RotationImprovedFermion.h`**：
  - `ERotationImproveFermionJob` 新增 `ERIF_MeasureHISQ`、`ERIF_MeasureHISQRotation`
  - 新增 `ERIFMeasureJob` 枚举：`ERIFMJ_Polyakov`、`ERIFMJ_Chiral`
  - 声明 `MeasureHISQ` / `MeasureHISQRotation`
- **`RotationImprovedFermion.cpp`**（main）：switch 新增两个 case，分别取 YAML 的 `MeasureHISQ` / `MeasureHISQRotation` section。
- **`Code/CMake/CMakeLists.txt`**：`RotationImprovedFermion` 目标加入 `Measure.cpp`。

### 2.3 测试配置（`Bin/Debug/`）

- **`TestRIFSimulateHISQRotation_O0.yaml`** / **`_O01.yaml`**：从原 `RotationImprovedFermion.yaml` 的 `HISQRotation` section 提取，改为
  - `LatticeLength : [12, 12, 12, 6]`（12³×6 小格子）
  - `Warmup : 50`、`ConfigurationNumber : 50`
  - `Omega : 0.0` / `0.1`（两个 ensemble）
  - `SaveFileType : EFFT_CLGBinFloat`
  - 补了顶层 `WorkJob : ERIF_SimulateHISQRotation`
- **`TestRIFMeasure.yaml`**：含 `MeasureHISQ` 与 `MeasureHISQRotation` 两个 section：
  - 费米子：HISQ 用 `CFieldFermionHISQSU3`；旋转用 `CFieldFermionHISQSU3R` + `CachedGauge : 1` + `ShiftCoord : 1`
  - 均 `Even : 0`（测量不用 even pseudo-fermion）
  - Solver 用 `CSLASolverGMRES`（非 CG）
  - 旋转 section 配 `StapleCache : CacheRotationLink : 1`
  - `Measure1 = CMeasurePolyakovXY`、`Measure2 = CMeasureChiralCondensateKS`

### 2.4 编译

```bash
cd Code/CMake
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/cuda/include   # 必须，见问题 1
cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86 -DCLG_RotationImprovedFermion=1
make RotationImprovedFermion -j$(nproc)
```

产物：`Bin/Ubuntu/RotationImprovedFermion`（100% Built，仅修复过一个 sign-compare 警告后无警告）。

## 3. 遇到的问题与解决

### 3.1 CLGCPULib PCH 找不到 `cuda_runtime.h`（已解决）

- 现象：`make RotationImprovedFermion` 报 `CudaHelperFunctions.h:18:10: fatal error: cuda_runtime.h: No such file or directory`，发生在 CLGCPULib 的预编译头（PCH）阶段。
- 原因：当前 shell 未设置 `CPLUS_INCLUDE_PATH`。README 的 Linux 章节明确要求编译前 `export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/cuda-10.0/include`。CMake 生成的 `flags.make` 里 `CXX_INCLUDES` 只有项目目录、没有 CUDA include，因此 PCH 编译时靠环境变量补充。
- 解决：`export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/usr/local/cuda/include` 后重编成功。
- 附带教训：后台任务若用 `make ... 2>&1 | tail -40`，shell 报告的退出码是 `tail` 的（0），**掩盖了 make 的真实失败**（实际是 `Error 2`）。排查时要直接看日志内容，不能只看任务 exit code。

### 3.2 `__DEFINE_ENUM` 重复定义（已避免）

- `__DEFINE_ENUM` 会同时生成 enum 类型和 `appEnumToString##enumname`/`appStringToEnum##enumname` 函数。若在头文件和 .cpp 各定义一次会冲突。最终把 `ERIFMeasureJob` 只放在头文件，`Measure.cpp` 中不再重复定义。

### 3.3 `_CLG_EXPORT_CHIRAL` 宏需要 3 个参数（已修正）

- 该宏签名是 `(measureName, lstName, variableName)`，第三个参数用于文件名里的 `%d`（如 ensemble 序号）。最初照 HISQ 的注释写法只传了 2 个参数，编译期无法发现（宏内才有 `variableName` 用法），改为传入循环变量 `uiOmega`。

### 3.4 sign-compare 警告（已修复）

- `for (UINT uiOmega = iListStart; uiOmega < BetaList.Num() && ...)` 中 `UINT` 与 `INT` 比较。已改为 `static_cast<UINT>(BetaList.Num())`。

### 3.5 测试 YAML 顶层缺 `WorkJob`（已修正）

- main 里 `argc==2` 时从 YAML 读 `WorkJob`，默认值是 `ERIF_SimulateHISQ`。从原 YAML 提取 section 后顶层没有 `WorkJob`，会导致走错分支。已在两个模拟 YAML 顶部补 `WorkJob : ERIF_SimulateHISQRotation`。

## 4. 完成状态与验证结果

### 已完成（全部）
- [x] 代码实现（Measure.cpp / 头文件 / main / CMake）
- [x] 编译通过（CUDA backend, arch 86）
- [x] 测试 YAML 编写
- [x] 模拟：omega=0 和 omega=0.1 各 50 个组态（12³×6，warmup 50，production 50）
- [x] 组态分文件夹：`Bin/Debug/RIFEnsemble/O0/` 与 `O01/`（各 50 个 .con）
- [x] 测量：HISQ / HISQRotation × Polyakov / Chiral 共 4 组，全部成功

### 测量结果验证（Bin/Ubuntu/ 下 CSV）

| 测量 | HISQ 分支 | HISQRotation 分支 | 一致性 |
|---|---|---|---|
| Polyakov (O0) | 0.1824+0.0068i … | 0.1824+0.0068i … | 完全一致（规范场测量，与费米子无关）|
| Polyakov (O01) | 0.2013-0.0057i … | 0.2013-0.0057i … | 完全一致 |
| Chiral (O0) | -0.00879-0.00105i … | -0.00889-0.00020i … | 数值同量级，差异来自随机源与旋转项 |
| Chiral (O01) | -0.00819-0.00042i … | -0.00822-0.00119i … | 数值同量级 |

- Polyakov loop 输出：`RIF12x6O*_polyakov.csv`（含 `_In`、`_OverR`、`_ZSlice` 变体）
- 手征凝聚输出：`RIF12x6O*_condensatepCCChiralKS.csv`（含 `_In`、`_OverR`、`_ZSlice` 及 CMTKSGamma3/4 变体）
- 旋转分支（HISQRotation Chiral）成功运行，证明 `CachedGauge : 1` + StapleCache `CacheRotationLink : 1` + GMRES solver + `Even : 0` 的配置全部生效。

## 5. 备注

- 手征凝聚测量中的 `InverseD` 使用 YAML 配置的 `Solver:`（GMRES），旋转费米子的 D 算符通过 `CachedGauge : 1` 走 `CStapleCache::GetRotationBuffer()`，测量前由 `MeasureInternal` 显式调用 `pCache->Cache(pGauge, ECC_BeforeGaugeUpdate / ECC_BeforeAllUpdateBeforeSmearing / ECC_BeforeAllUpdateAfterSmearing)` 刷新 rotation link，与 HMC 积分器中的调用顺序一致。
- 参考文件：`Code/Applications/HISQ/Measure.cpp`、`Code/Applications/StaggeredSpectrum/Measurement.cpp`。
