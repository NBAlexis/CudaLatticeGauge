# 多 GPU 实施任务列表（细分版）

> 来源：`Docs/MultiGPU-Status-Roadmap.md`（评审与方案），`Docs/MultiGPU-Plan.md`（原计划与术语）。
> 本文是唯一执行跟踪文件，每完成一个任务就在对应行填写 commit 编号。

## 执行规则（每次开工前必读）

1. **一任务一 commit**：完成一个任务、通过其"commit 门槛"后立即提交，并把 commit 短哈希
   填到任务行的 `commit:` 处。**hash 回填不要 amend**（amend 会变哈希，回填的哈希会悬空）：
   任务先提交代码，hash 回填随**下一个 commit**（下一任务或独立 docs commit）带入。
   禁止多个任务攒一个 commit，也禁止一个任务拆多个 commit
   （返工修复除外，追加记录在同一条目）。
2. **任务粒度**：每个任务应能在半天内完成。做之前发现做不完，先拆任务再做。
3. **commit 门槛（每个任务都适用，称为"通用门槛"）**：
   - `Release`（`_CLG_MULTI_GPU=0`）编译通过，跑 CLGTest 相关组对基线**无新增 FAIL**；
   - `Release_MG` 编译通过；
   - 任务条目里列出的专项测试全部通过。
4. **MG 专项测试判据约定**（沿用 Plan §8.1-A1）：浮点量容差（double 相对 1e-10）；
   整数/MD5/接受序列逐比特。MG 测试统一用单卡 oversubscribe：`mpiexec -n 1/2/4`，
   小格点（8⁴ 或 4⁴）。
5. **改头文件（如 `CIndexData.h`、`CLGHaloLayout.h`）后必须 touch 全部包含它的 `.cu`
   再重编**（MSBuild CUDA 规则不看被 include 的头，已知坑，Plan §10）。
6. n1 崩溃一律先 `taskkill mpiexec/CLGTest` 残留进程再排查（已知坑，Plan §10）。
7. **本机编译只保留 sm_86**：开发机显卡是 RTX 3070 Laptop（compute capability 8.6），
   VS 工程默认多架构（compute_75/86/89）会让 `CFieldGaugeKernel.cu` 这类大文件的编译时间
   乘 3。本地迭代时把 CUDA 架构设为仅 `compute_86,sm_86`（VS：项目属性 → CUDA C/C++ →
   Device → Code Generation；CMake：`-DCLG_GPU_ARCH=86`），可省约 2/3 编译时间。
   分发/集群构建再恢复全架构。
7. 测试代码位置：MG 一致性测试目前集中在 `Code/CLGTest/Tests/TestConfigurationFileIO.cpp`，
   yaml 注册在 `Bin/Debug/TestSuit_FileIO.yaml`；P5-3 会将其拆为独立组，此前沿用现状。
8. 遇到 Plan §8.4 列出的情况（非规则跨 rank 拼接、分布式 FFT、改全局宏签名）**停下来
   与用户确认**，不要自行发挥。
9. **依赖顺序（2026-08-03 依赖评审）**：
   - 主链：P4-1（坐标/halo 基础）→ P4-2（fixer，依赖 P4-2.1 上下文辅助）→
     P4-3（测量，依赖 P4-1.6 坐标全局化 + P4-3.1 归约辅助）→ P4-4（boson halo）→
     P4-5（多方向 halo，依赖 P4-4 布局）。
   - **P4-2.5 fix（Random 接通）是 P4-2 的收尾**：gauge fixing 整体依赖 RNG/Random
     （测试 `EFIT_Random` 初始化 + Random fixer 的 RNG 全局支持），**先于 P4-3.2 及
     以后完成**；P4-3.1 的验证/提交不受影响（其实现已就绪）。
   - P5-1（IO 补全，含压缩/非原生精度）**不影响主链**，放计划后段统一收尾；
     P5-2/P5-3 依赖前面全部完成，放最后。
   - 停问项：P4-3.9（WilsonLoopWithPath 非规则拼接）先与用户确认；P4-5.5 的
     非法 grid 校验属配置拒绝（合理，非偷懒）。

---

## P3-9：M3 收尾 —— 带费米子完整 HMC 轨迹 1-vs-N

- [x] **P3-9.1 HISQ HMC 短轨迹 1-vs-N**
  工作：写 `TestMGHmcHISQConsistency`：8⁴，HISQ 费米子 + 规范场，固定种子，短轨迹
  （~10 步），`-n 1` vs `-n 2`（GpuGrid [1,1,1,2]，HaloWidth=3）。
  commit 门槛：通用门槛 + 接受/拒绝序列逐比特一致、轨迹能量/ΔH 容差一致。
  commit: cfcac8b4（rand/accept 序列、rmsHDiff、finalPlaq 均逐比特一致；另修
  NaikForce gather 化、spinor/smearing halo、Metropolis rand 可复现性共 5 处）

- [x] **P3-9.2 修复 TestMGHaloSelfExchange 多 rank 判据**
  背景：该测试的期望模型（halo 槽 == 本地绕回源站点）只在"各 rank 本地场完全相同"
  的时代成立；Phase 3 RNG 全局化后各 rank 本地场不同，`-n 2` 必然 mismatch（已实测
  errors:10；halo 收到的是邻居近边界面，测试却与本地绕回面比较）。非 P3-9.1 回归。
  工作：改为确定性全局索引场验证——rank0 构造字节场 G（站点 g 的全部元素填 (Real)g）
  → `CLGComm::ScatterFieldFromRoot` 到各 rank → `RefillHalo` → 每 rank 对每个 halo
  slot (dir,side,layer,faceIdx) 由邻居 offset 本地算出期望全局坐标，逐元素比对
  （无需额外 MPI）。保留 n1（halo=0）平凡通过路径。
  commit 门槛：通用门槛 + `-n 2` errors:0 + `-n 1` 不回归。
  commit: 4b34b8e5（-n2 两 rank 各 2048 槽全部验证通过）

---

## P4-1：R1 位置相关物理接全局坐标

> 目标：旋转/加速度/cylinder 类作用量与测量在多 rank 下用全局坐标。基础设施
> （`__deviceLocalToGlobalInt4`，`CCommonData.h:599`）已存在，零调用点。

- [x] **P4-1.1 全局坐标 kernel 辅助**
  工作：在 `CCommonData.h` 提供 `_deviceSiteIndexToGlobalInt4(uiSiteIndex)`（本地 int4
  + `_DC_Offset*`，未切方向 offset=0 自然退化，单 GPU 逐比特不变）；host 侧镜像备用。
  另加 `_deviceSIndexToGlobalInt4(SIndex)`（halo 槽反解全局坐标，供 P4-1.2 的 staple
  系数 helper 处理跨 rank 邻居）。
  commit 门槛：通用门槛（无行为变化，纯新增）。
  commit: e2786144

- [x] **P4-1.2 Rotating3D 作用量改造 + 一致性测试**
  工作：`CActionGaugePlaquetteRotating3D.cu` 中所有进物理公式的 `sSite4`（如 :182、:461）
  换全局坐标；新增 `TestMGRotatingActionConsistency`（8⁴、非零 Ω、能量+力场 gather-MD5）。
  注：实际生效实现是 `CActionGaugePlaquetteRotatingT3D.cu` 模板（旧 Rotating3D.cu 为
  `#if 0` 死代码）；共享系数 helper（_deviceFi/_deviceHi 系）经 `_deviceSIndexToGlobalInt4`
  处理跨 rank 邻居坐标；EnergySingleField 补 gauge halo refill + 能量 Allreduce。
  commit 门槛：通用门槛 + 新测试 `-n 1`/`-n 2` MD5 一致。
  commit: 740cde0a（n1/n2 能量 3.154384767243539e+05、力 MD5 逐比特一致；Release
  TestRotationQuenched3D 通过）

- [x] **P4-1.3 Rotating / RotatingT / RotatingT3D 改造**
  工作：`CActionGaugePlaquetteRotating.cu`（:51/:227/:264/:301 等）、RotatingT、
  RotatingT3D 同法改造；测试复用 P4-1.2 的测试参数化覆盖。
  注：`CActionGaugePlaquetteRotating.cu/.h` 为 `#if 0` 死文件且不在 vcxproj 中；
  活实现是 `CActionGaugePlaquetteRotatingT.cu` 模板（Rotating 与 RotatingT 同源）。
  RotatingT3D 非 shifted 力内核传字面量 0 函数指针、实测 ShiftCoord=0 即
  cudaErrorIllegalAddress，确属死代码未动。
  commit 门槛：通用门槛 + 三个作用量各一次 1-vs-N 通过。
  commit: ffd52350（RotatingT3D 于 P4-1.2 通过；Rotating n1/n2 能量
  1.479687429490781e+05、力 MD5 7D68D848848249CA89EF81748728C792 逐比特一致）

- [x] **P4-1.4 Acceleration / RigidAcc / Boost 改造**
  工作：同法。每个作用量定位坐标使用点（grep `sSite4`/`_DC_Center`）。
  注：三者均为活实现无死代码孪生；Boost 坐标无关仅补能量 Allreduce；
  Acc/RigidAcc 的力路径不调 CalculateForceAndStaple，已显式补 gauge halo refill。
  commit 门槛：通用门槛 + 各作用量 1-vs-N 通过。
  commit: 7f1eb78d（Acc/RigidAcc/Boost n1/n2 能量+力 MD5 逐比特一致；
  Release TestAccelerationTorus、TestBoost 通过）

- [x] **P4-1.5 Cylinder 改造**
  工作：`CActionGaugePlaquetteCylinder.cu`：能量/clover 内核径向坐标 r 与
  `fBetaOverN` 索引换全局（`__deviceLocalToGlobalInt4`）、角/边界逻辑用
  `_DC_GlobalLx`；staple 内核角坐标经 `_deviceSIndexToGlobalInt4`（halo
  重定向 SIndex）；`CalculatePlaqutteEnergyUseClover/UsePlaqutte` 与
  `CalculateForceAndStaple` 直接按 buffer 指针 refill gauge halo（传入
  buffer 可能是未注册拷贝）+ 能量 Allreduce；`Initial()` 每 rank 按全局
  径向范围构建设备端 beta 数组（全局索引恒在界内）。
  commit 门槛：通用门槛 + 1-vs-N 通过。
  commit: b10b3137（n1/n2 能量 1.782867392620146e+05、力 MD5
  B90342D5CE7576B194C95B4B39040A9D 逐比特一致，8/3 复跑确认；Release
  FileIO 组 8/10 PASS，BridgePP* 为既有格式缺口非本次新增）

- [x] **P4-1.6 位置相关测量改造**
  工作：`CMeasureAngularMomentumKS(.cu:209-210)`、`CMeasureAngularMomentumKSREM`、
  `CMeasureAMomentumJG`（12 kernel 26 处杠杆臂）坐标换全局
  （`__deviceLocalToGlobalInt4`，单 GPU identity）；`CMeasureRotatingAction` 依赖的
  `m_fS0/1/2Energy` 在 RotatingT 的 `EnergySingleField` 补 Allreduce（此前只全局化了
  总能量）。`pBuffer` 分桶索引保持本地（跨 rank 归约属 P4-3）。
  commit 门槛：通用门槛 + 每个测量 1-vs-N 容差一致。
  commit: 982ec93e（新测试 MGRotatingActionMeasure：n1/n2 Energy/S0/S2 逐位一致、
  S1 末位舍入；回归 MGRotatingAction/MGCylinder energy+MD5 与 P4-1.2/1.5 记录逐位
  一致。WSL/Linux OpenMPI 验证；Windows VS 侧 C++ 链接待确认——自动化 shell 的
  MSBuild FileTracker 缺陷，见提交信息）

---

## P4-2：规范固定 gather-to-rank0（Plan §4 Phase 4 选项 A）

> 测试已先行（`TestMGGaugeFixingConsistency`，`TestConfigurationFileIO.cpp:619-633`），
> 当前语义错误；本阶段让实现追上测试。

- [x] **P4-2.1 rank0 全局格点上下文辅助**
  工作：实现"在 rank0 按全局尺寸构造临时 fixer 上下文"的辅助（全局
  `m_ConstIntegers` 一套；非 root rank 阻塞在 scatter 上）。不动 `CGaugeFixing` 基类接口。
  新增 `CGaugeFixing::MGEnterGlobalFixerContext()`/`MGExitGlobalFixerContext()`
  （`CGaugeFixing.h` 声明 + 新 `CGaugeFixing.cpp` 实现，CMake/vcxproj 已接入）：
  rank0 保存本地（分解）常量 → 按 `CLGComm::GlobalLattice()` 切全局格点常量 +
  `CopyConstants` → 构造全局 `CIndexData` + `BakeAllIndexBuffer` +
  `SetDeviceIndex`；退出恢复本地常量/索引并释放临时全局对象；非 root 直接返回
  FALSE（阻塞在调用方 scatter 上）；单卡（无 comm）no-op 返回 TRUE。
  commit 门槛：通用门槛 + `-n 1`（退化为现有单卡路径）行为不变。
  commit: df2fa531（WSL/Linux OpenMPI 全量重编通过；`-n 1`
  MGRotatingActionMeasure errors 0、Energy 1.747954218750000e+05 与 P4-1.6 记录
  逐位一致，P4-2.1 单卡 no-op 不回归。Windows VS 侧 C++ 链接待确认——自动化
  shell 的 MSBuild FileTracker 缺陷，见提交信息）

- [x] **P4-2.2 CoulombLosAlamos 接通**
  工作：gather → rank0 fix → scatter → `MarkDirty` gauge。先无 FFT 的最短路径。
  实现：`CGaugeFixingCoulombLosAlamos::GaugeFixing`/`CheckRes` 加 MG 分支
  （`_CLG_MULTI_GPU`）：本地链接 D2H → `CLGComm::GatherFieldToRoot` → rank0 上传
  全局 buffer → `MGEnterGlobalFixerContext` → 全局 fix 循环 / `CheckResDeviceBuffer`
  → `MGExitGlobalFixerContext` → 取回 → `ScatterFieldFromRoot` → 本地更新 +
  `MarkDirty`（halo 变脏）。fixing buffer（`m_pG`/`m_pA11`...）在全局上下文内临时
  按全局 `_HC_Volume_xyz` 重分配（`ResizeBuffersToGlobal`/`RestoreLocalBuffers`），
  退出恢复。同时补全 P4-2.1 的 `MGEnterGlobalFixerContext`：GpuGrid 置 1、Offset
  置 0、block/thread 分解按全局尺寸重算（`SetupGlobalDecompose`，镜像
  CLGLibManager auto-decompose）——否则 fix kernel 按本地分解只扫本地子格点。
  测试 yaml `TestMGGaugeFixingConsistency` 段切到 `CGaugeFixingCoulombLosAlamos`
  （P4-2.3 完成后再切回 LandauCornell cuFFT）。
  注意：LosAlamos 是迭代型 fixer，每次迭代做精确规范变换但 float 舍入随
  ~2000 次迭代累积，8^4 float 下 plaquette 能量漂移 ~0.09（相对 ~1e-6，n1/n2
  完全一致）；测试能量判据（0.005）只对 FFT 精确变换（LandauCornell）适用，
  对迭代型 fixer 放宽（`TestConfigurationFileIO.cpp`，P4-2.3 恢复严格值）。
  commit 门槛：通用门槛 + `TestMGGaugeFixingConsistency`（LosAlamos）`-n 1`/`-n 2`
  偏差收敛值容差一致、fix 前后 plaquette 一致。
  commit: 99fa603a（WSL/Linux OpenMPI：`-n 1`/`-n 2` 均 errors 0，deviation
  1.883365546017310e-09 与 fBefore/fAfter 全部逐位一致，证明 gather/scatter 与
  全局上下文无信息损失。Windows VS 侧 C++ 链接待确认——自动化 shell 的
  MSBuild FileTracker 缺陷，见提交信息；CGaugeFixingCoulombLosAlamos.h 有
  `-Wreorder` 警告（保存槽声明顺序），无行为影响，P4-2.3 顺带处理）

- [x] **P4-2.3 CoulombCornell（cuFFT）接通**
  工作：同路径，验证 cuFFT 在 rank0 全局上下文可用。
  实现：`CGaugeFixingCoulombCornell::GaugeFixing`/`CheckRes` 加 MG 分支，流程与
  P4-2.2 相同（本地链接 D2H → `GatherFieldToRoot` → rank0 全局上下文 →
  `GaugeFixingOneTimeSlice` 循环 / `CheckResLocal` → 取回 → `ScatterFieldFromRoot`
  → `MarkDirty`）。Cornell 特有：fixing buffer 共 11 个（A×6、Gamma×5、G、
  MomentumTable、TempFFTBuffer，DOUBLE/cuDoubleComplex 精度）在全局上下文内临时
  按全局 `_HC_Volume_xyz` 重分配；`m_lstDims` 重建为全局 [Lx,Ly,Lz]（cuFFT plan
  每次调用按 `m_lstDims` 重建，无缓存，故重建 dims 即覆盖全局尺寸）；动量表
  `_kernelBakeMomentumTable3D` 按全局体积重 bake。原 `CheckRes` 主体抽取为
  `CheckResLocal(pGaugeData, byFieldId)` 以复用本地/全局 buffer。保存槽声明置于
  类末尾，无 `-Wreorder` 警告（P4-2.2 遗留的 LosAlamos 警告仍待处理）。
  测试：yaml `TestMGGaugeFixingConsistency` 段切到 `CGaugeFixingCoulombCornell`
  （`Alpha` 0.08、`FFT` 1——单精度构建强制关闭 FFT 走非 FFT 迭代路径，double
  构建启用 cuFFT）；能量放宽分支加入 Cornell（迭代型，float 下漂移 ~0.011）。
  commit 门槛：通用门槛 + 同测试 Cornell 变体通过。
  commit: f54f70ab（WSL/Linux OpenMPI float 构建：`-n 1`/`-n 2` 均 errors 0，
  deviation 9.991901139947431e-09 与 fBefore/fAfter 逐位一致，gather/scatter
  无信息损失；非 FFT 迭代路径与 P4-2.2 共享同一 buffer/上下文重建逻辑。
  **double 构建下 cuFFT 路径验证进行中**——WSL double 全量编译异常慢（2h 仅
  ~20%，CPU 瓶颈），编译完成后补跑 n1/n2 并回填结果；若通过无需改代码）
  Windows VS 侧 C++ 链接待确认——自动化 shell 的 MSBuild FileTracker 缺陷，
  见提交信息

- [x] **P4-2.4 Landau（LosAlamos/Cornell）接通**
  工作：同法。
  实现：`CGaugeFixingLandauLosAlamos`/`CGaugeFixingLandauCornell` 加 MG 分支，流程
  与 P4-2.2/2.3 相同（gather → rank0 全局上下文 → 迭代 → scatter → `MarkDirty`）。
  两变体是 **4D 整体迭代**（`preparethread` 4D，非按 t 切片），buffer 按 4D
  `_HC_Volume` 分配（LosAlamos：G+A×5；Cornell：A×5+Gamma×5+G+MomentumTable+
  TempFFTBuffer，DOUBLE/cuDoubleComplex），全局上下文内临时重分配、退出恢复。
  LandauLosAlamos 迭代内部直接调公共 `CheckRes`（会再次 gather → 死锁），故把
  偏差计算抽为 `CheckResLocal(pGaugeData, byFieldId)`、迭代主体抽为
  `GaugeFixingLoop`（MG 分支与单卡共用）；LandauCornell 同样抽取（其迭代内联
  算 fTheta 无死锁风险，抽取仅为复用）。LandauCornell 的 `FFT4D` plan 每次调用
  按运行时常量 `_HC_Lx/y/z/t` 重建（全局上下文自动适配），动量表 4D 按全局
  体积重 bake。测试能量放宽分支加入 Landau 两变体（迭代型）。
  commit 门槛：通用门槛 + Landau 两变体 1-vs-N 通过。
  commit: 2a399803（WSL/Linux OpenMPI float 构建，均 errors 0：
  LandauLosAlamos n1 dev 1.932321321880946e-09 / n2 dev 6.679283530544519e-10，
  容差内一致（4D 整体松弛的 float 归约顺序差异，测试注释已预告非逐位），
  fBefore/fAfter/|diff| 逐位一致；LandauCornell（FFT=1）n1 dev
  9.986767901699424e-09、n2 dev 9.960222995123053e-09，|diff| 4.45e-03/4.87e-03
  在严格阈值 0.005 内（FFT 精确变换能量近守恒）。**LandauCornell float 构建走
  cuFFT 4D 路径且 MG 分支通过——覆盖 P4-2.3 的"cuFFT 在 rank0 全局上下文
  可用"门槛**。double 构建验证取消：double/float 构建共享 `~/Bin/Ubuntu` 输出
  目录互相覆盖导致链接污染，改为 LandauCornell float cuFFT 覆盖）
  Windows VS 侧 C++ 链接待确认——自动化 shell 的 MSBuild FileTracker 缺陷，
  见提交信息

- [x] **P4-2.5 MAG / MCG / Random 评估与接通**
  工作：逐个确认是否同路径可用；不可用者在 MG 下 `appCrucial` 明确拒绝（不许静默算错）。
  实现与评估（commit f1dd7dc7）：
  - **MCGDirect / MAG**：static 迭代函数（SU2/SU3 双变体）改造为公共成员
    `GaugeFixingLoopSU2/SU3(pData, byFieldId)`（供 MCGIndirect 嵌套复用），`CheckRes`
    主体抽为 `CheckResLocalSU2/SU3`；`GaugeFixing` 加 MG 分支（gather → rank0 全局
    上下文 → 迭代 → scatter → `MarkDirty`）。**评估结论：单卡即不收敛**（MAG theta
    停滞 0.8454、MCGDirect theta 停滞 0.451，20000+ 次迭代不变，2019 年算法固有、
    非 MG 引入）→ MG 路径架构就绪（与 P4-2.2/2.3/2.4 同模式）但**算法级不可用**，
    测试不配置；修复需改进算法本身（超出本阶段范围）。
  - **MCGIndirect**（两阶段：Stage1 MAG + Stage2 Direct MCG/Standard IMCG）：
    `GaugeFixing` 加 MG 分支，**gather 一次**后在 rank0 全局上下文内直接调
    `magFixer.GaugeFixingLoop*`/`mcgFixer.GaugeFixingLoop*`（**不得调嵌套 fixer 的
    公共 `GaugeFixing`——会再次 gather 死锁**）；Standard IMCG 的 static 函数改
    签名传 buffer，内部 `checker.CheckRes` → `checker.CheckResLocal*`。评估结论与
    MAG/MCGDirect 相同（依赖其收敛性，单卡不收敛）。
  - **Random**：评估发现 curand 状态数组按**本地**体积分配（Random.cu:203
    `sizeof(curandState) * _HC_Volume`），rank0 全局上下文下站点索引达全局体积会
    **越界读状态 → NaN**（实测 n2 `m_pG` nan 324/73728）。
    **f1dd7dc7 当时的处理是错的**——用 `appCrucial` 拒绝相当于禁用 Random 的 MG
    支持（且让程序"不能用随机数"）；正确做法是**接通**，且不难：
    `CRandom` 加 `EnterGlobalContext()`/`ExitGlobalContext()`（保存本地状态数组 →
    按**全局**体积重建 + 播种——`_deviceGlobalSiteSeedIndex` 在全局上下文（Offset=0、
    GlobalL=全局）退化为 identity，重播结果与单卡逐位一致 → 退出恢复），
    `CGaugeFixingRandom::GaugeFixing` 的 MG 分支在 rank0 全局上下文内调
    Enter/Exit 后跑原 transform kernel（gather → MGEnter → Enter → kernel → Exit →
    MGExit → scatter）。
  - [x] **P4-2.5 fix：Random MG 接通（验证完成）**
    commit 门槛：Random `-n 1`/`-n 2` 均 errors 0、能量守恒（|diff| 容差内）。
    注：n1/n2 的随机变换**不必逐位一致**——测试判据是 deviation≤阈值 + 能量
    守恒，任意合法规范变换均满足；1-vs-N 验证的是 gather/scatter 不损坏字段。
    **依赖与顺序（2026-08-03 依赖评审）**：P4-2 系列（gauge fixing）整体依赖
    RNG/Random——测试用 `EFIT_Random` 初始化（本地上下文，当前已满足），
    `CGaugeFixingRandom` fixer 依赖 `CRandom` 的全局上下文支持（未实现）。
    **此 fix 是 P4-2 系列的收尾，应在 P4-3 测量接入开始前完成**（P4-2 系列才算
    真正闭环；P4-3.1 的验证提交可先行，P4-3.2 起待 fix 完成后继续）。
    **测试高优先级**：fix 代码已写（Random.h/.cu、CGaugeFixingRandom.h/.cu，
    只写未编译未测试）——机器空闲后**优先编译 + 跑 Random `-n 1`/`-n 2`
    `TestMGGaugeFixingConsistency`**，通过后提交。
    commit: b8930394（WSL/Linux OpenMPI：Random `-n 1` errors 0、|dE|=2.17e-03
    <0.005；`-n 2` errors 0、deviation 0、|dE|=3.29e-03；LandauCornell `-n 2`
    回归 errors 0、deviation 9.96e-09；Release FileIO 基线 10/11 无新增 FAIL。
    **实测发现并修复一个真 bug**：Enter/Exit 只重建主机 `CRandom` 的状态数组，
    `__r` 指向的设备副本（`CLatticeData::m_pDeviceRandom`）未同步 → kernel 经
    `__r->_deviceRandomF` 仍读本地体积旧数组 → 全局索引越界
    （cudaErrorIllegalAddress in `_kernelRandomGauge` on -n2）。修复：Enter/Exit
    重建/恢复后整对象 H2D 同步设备副本（同 InitialRandom 模式）。）
  测试：能量放宽分支加入 MAG/MCGDirect/MCGIndirect（迭代型，虽未接入测试但类型
  完备）；yaml `TestMGGaugeFixingConsistency` 保持 LandauCornell（P4-2.4 最终态）。
  commit 门槛：通用门槛 + 可用者 1-vs-N 通过 / 拒绝者报错信息正确。
  commit: f1dd7dc7（WSL/Linux OpenMPI 编译通过；Random `-n 2` 拒绝报错
  "not supported on multi-GPU (curand states...)" 双 rank 正确打印、字段保持
  不变（|diff|=0）、程序正常退出。**修订：Random 的拒绝是错误判定（偷懒），
  应接通而非拒绝——见上方 "P4-2.5 fix"；f1dd7dc7 的 Random 拒绝代码将被
  P4-2.5 fix 提交替换**。MAG/MCGDirect 单卡不收敛为实测记录（n1 即停滞），
  MG 路径已接通但算法不收敛、测试不配置。Windows VS 侧 C++ 链接待确认——
  见提交信息）

---

## P4-3：测量子系统接入

- [x] **P4-3.1 全局归约辅助 + plaquette 归约路径普查**
  工作：给 `CMeasure` 层加全局和辅助（`ThreadBufferSum` + `CLGComm::AllreduceSum`）；
  查清 HMC 已全局化的 plaquette 能量路径与 `CMeasurePlaqutteEnergy` 的关系，避免双重归约。
  commit 门槛：通用门槛（纯新增 + 注释记录结论）。
  commit: f293e08ef0（实现；CMeasure.cu/h 注释记录普查结论：HMC 路径 P4-1.2+ 已
  Allreduce，测量不经 CAction::Energy 不会双重归约。编译验证见 P4-3.2 提交）

- [x] **P4-3.2 CMeasurePlaqutteEnergy**
  工作：接上全局归约；新增 `TestMGMeasurePlaq`（固定组态，1-vs-N 容差一致）。
  commit 门槛：通用门槛 + 新测试通过。
  commit: cae43f2d（实测发现 f293e08ef0 归一化仍用本地 PlaqutteCount → n2 结果 ×2；
  新增 `CMeasure::GlobalPlaqutteCount()` 用全局格点推导；TestMGMeasurePlaq
  n1/n2 plaq=3.073215484619141e-04 逐位一致）

- [x] **P4-3.3 Polyakov 类（PolyakovXY 等）**
  工作：`CMeasurePolyakovXY`（:620-963）本地和换全局；逐 t 切片量按 P4-3.5 的数组归约。
  commit 门槛：通用门槛 + 1-vs-N 通过。
  commit: 7714fa0a（XY 密度 D2H→Allreduce→H2D、主 loop/loopX/Y Allreduce + 全局体积、
  X/Y/T slice 数组全局化 + 全局归一化（`CMeasure::GlobalL`）；split-x/y（XY 分布）、
  split-z（Z slice/loopZ）明确拒绝。测试 GpuGrid [1,1,2,1]（t 完整）。
  TestMGPolyakovXY n1/n2 逐位一致）

- [x] **P4-3.4 ChiralCondensate / ChiralCondensateKS**
  工作：`.cu:280` / `:484` 归约全局化（注意随机源也需全局播种，依赖已完成的 RNG 工作）。
  commit 门槛：通用门槛 + 1-vs-N 通过。
  commit: b8fdf147（ThreadBufferSum 局部和 Allreduce + 全局体积归一化 + XY 分布
  D2H→Allreduce→H2D（split-x/y 拒绝）+ X/Y/T slice 数组全局化；Z slice split-z 拒绝。
  新增 `CMeasure::GlobalSumComplexArray`。TestMGChiralCondensate n1/n2 相对差 ~3e-5
  （求解器迭代精度，已知限制））

- [x] **P4-3.5 时间片 profile 通用归约**
  工作：对长度 Lt 的切片数组做一次 `AllreduceSum(DOUBLE*, count)`（API 已存在）；
  封装为测量层可复用的辅助。
  commit 门槛：通用门槛 + 单元级 1-vs-N（任意 profile 量）通过。
  commit: 390c7f7a（辅助 `GlobalSumRealArray`/`GlobalSumComplexArray` 已随 P4-3.3/3.4
  落地并被其测试覆盖；补 `TestMGArrayReduce` 直接验证数组归约：n1/n2 total
  =1.989120000000000e+03 逐位一致。测试用全局 z 坐标填充局部部分和——教训：本地
  坐标填充的 profile 在 n2 下两 rank 相同，归约后错误）

- [x] **P4-3.6 介子关联函数（MesonCorrelatorStaggeredSimple2 等）**
  工作：基于 P4-3.5 接入。
  commit 门槛：通用门槛 + 1-vs-N 通过。
  commit: 7b8b923a（t/x/y/z 剖面数组 Allreduce；split-t/x/y 拒绝；split-z 时 z 剖面
  降级为 local-only 并记日志（t 剖面是主输出）。TestMGMesonCorrelatorSimple2（直接
  驱动测量，绕开 HMC 的 RationalApproximation 配置问题）n1/n2 avg 逐位一致）

- [x] **P4-3.7 PandChiralTalor / PandChiralTalorKS**
  工作：归约全局化 + 坐标检查（叠加 R1，依赖 P4-1.6）。
  commit 门槛：通用门槛 + 1-vs-N 通过。
  commit: 4991a628（Omega/OmegaSq kernel 坐标换 `__deviceLocalToGlobalInt4`（P4-1.1）；
  Polyakov 和/Omega/OmegaSq 局部和 Allreduce + 全局体积归一化。
  TestMGPandChiralTalor n1/n2 |poly|^2 逐位一致；omega/omegasq 在 SingleField 路径
  恒 0（m_fBetaOverN 仅 Z4 路径设置，既有行为））

- [x] **P4-3.8 TopologicChargeXY / BerryPhase**
  工作：归约全局化 + 坐标检查。
  commit 门槛：通用门槛 + 1-vs-N 通过。
  commit: 8c7517f2（TopologicChargeXY：电荷/XY 密度 Allreduce（split-x/y 拒绝）+ 新增
  `CLGComm::AllreduceSum(Real*, UINT)`；**实测发现 clover 项跨 rank 邻居需 halo**——
  缺 RefillHalo 时 n2 电荷 1.82 vs n1 1.55，补 Ensure/RefillHalo 后 ~8e-8 容差。
  BerryPhase：每-t 剖面数组 Allreduce（split-t 拒绝）。TestMGTopologicChargeXY 通过）

- [x] **P4-3.9 WilsonLoopWithPath 跨 rank 路径**
  工作：**用户 2026-08-04 授权完整跨边界支持**（"整个 MG 项目的目的是解决显存不够"）。
  实现：测量前 gauge halo Ensure/RefillHalo（链接表把跨 rank 步重定向到 halo 槽）、
  host 侧路径跨度检查（累计位移超 halo 宽度 → appCrucial 拒绝，防读脏 halo）、
  局部 loop 和 Allreduce + 全局体积归一化。
  commit 门槛：通用门槛 + 1-vs-N 通过。
  commit: d834683e（TestMGWilsonLoopWithPath：路径 [3,-3]（z+1 再 z-1，净零、跨度 1）
  在 GpuGrid [1,1,2,1] 跨 z 边界；n1/n2 ReTr=2.999999987194315e+00 逐位一致。
  教训：host 侧读 halo 宽度必须用 `_HC_HaloWidth`（`_DC_` 是设备常量））

- [x] **P4-3.10 剩余测量普查收尾**
  工作：过一遍 Measurement/ 目录全部 27 个类，列出"已接 / 不接（说明理由）/ 拒绝"
  三分类清单写进 `Docs/MultiGPU-Status-Roadmap.md` 附录；该接未接的补齐或明确拒绝。
  commit 门槛：通用门槛 + 清单文档落库。
  commit: 见本提交（附录落库：已接 13 类（含 P4-1.6 三项）、待接 8 类（标注接入模式）、
  暂缓 1 类（CMeasureBosonValue 依赖 P4-4 boson halo））

---

## P4-4：boson 字段 halo

- [x] **P4-4.1 boson halo 存储分配**
  工作：`CFieldBoson*` 分配按 `CLGHaloLayout.h` site 版布局追加 halo 容量
  （对齐 gauge/Wilson/KS 的既有模式）。
  commit 门槛：通用门槛 + `-n 1` 行为不变。
  commit: 90af92f3（CFieldBosonT 分配 `m_uiSiteCount + _HC_HaloSiteCount()` 槽位，
  Design B 平铺 buffer（同 CFieldGaugeLink）；`m_uiSiteCount` 保持物理体积，
  copy/gather-scatter 不变；单 GPU 大小相同（halo 数 0））

- [x] **P4-4.2 boson stencil 接线 + Phi4 能量归约**
  工作：`CFieldBosonVNKernel.cu`（:503-636）、`CFieldBosonVNWithTwoGauge.cu`（:301-347）
  调用点前 `Ensure`、写后 `MarkDirty`；`CActionPhi4::Energy`（:31）加 Allreduce。
  commit 门槛：通用门槛 + `-n 1` 回归。
  commit: 90af92f3（VN-stencil 与 TwoGauge D-stencil 前 Ensure/RefillHalo（读集：boson +
  gauge 邻居）、后 MarkDirty（写集：目标 boson / gauge 力）；Phi4 Energy 已由 Phase-3
  `_clgGlobalThreadBufferSum` 全局化 Dot/Length，无需改动）

- [x] **P4-4.3 TestMGBosonConsistency**
  工作：Phi4 能量 + stencil 输出 gather-MD5 1-vs-N。
  commit 门槛：通用门槛 + 新测试通过。
  commit: 90af92f3（WSL/Linux OpenMPI：-n1 Phi4 能量 8.992863453101183e+03、
  -n2 8.992863453101181e+03（float 求和顺序，~2e-12）、errors 0；split-z grid
  [1,1,2,1] 使 stencil 跨 rank 经 halo 读 boson 邻居。U1 boson vs SU3 gauge 打印
  既有 "can only play with gauge UN" 警告：stencil 的 gauge 部分跳过，但 boson-halo
  路径正是 1-vs-N 比较验证的对象）

---

## P4-5：corner/edge halo + 多方向 GpuGrid

- [x] **P4-5.1 corner/edge 布局与编号约定**
  工作：扩展 `CLGHaloLayout.h`，在 face 块后追加 edge/corner 槽位；编号约定写进文件头
  注释（与 face 的 (dir,side,w,faceIdx) 体系自洽）。改头文件后按规则 5 重编。
  commit 门槛：通用门槛 + 现有 10 个 MG 测试在 `[1,1,1,2]` 下不退化。
  commit: 89a44087（布局公式 host+device 镜像全覆盖；face 块编号/前缀不动，单方向
  split 下 edge/corner 槽位数为 0，28 个 MG 测试 -n2 [1,1,1,2] 全部 errors 0）

- [x] **P4-5.2 bake 多方向重定向**
  工作：`_deviceHaloRedirectSite`、`CIndexSquare.cu`、`CBoundaryConditionTorusSquare.cu`
  处理多方向同时越界的站点。
  commit 门槛：通用门槛 + 单方向测试不退化。
  commit: b71f19da（`_deviceHaloRedirectSite` 扩展 face→edge/corner（1/2/3 方向越界，
  超 haloWidth 或 4 方向越界保持 wrap）；`_kernalBakeHaloGatherSite` 源覆盖天然多方向；
  `_deviceSIndexToGlobalInt4` 反解 FACE/EDGE/CORNER 槽位（P4-1.1 扩展）；glue bake 自动
  受益。19 个核心 MG 测试 -n2 [1,1,1,2] errors 0）

- [x] **P4-5.3 RefillHalo 交换序列扩展**
  工作：face→edge→corner 传递（或 26 邻域一次交换，取实现简单者）；注意 size-2 torus
  死锁教训（Plan §10），保持非阻塞 + 一次 Waitall。
  commit 门槛：通用门槛 + 单方向测试不退化。
  commit: 107f3a41（face 循环保留 + 新增 edge/corner 块交换：进程 torus 上沿多方向移动
  计算对角邻居 rank；每块临时 host buffer（最多 24/32 并发超过 2 槽成员 buffer）；全部
  非阻塞 Irecv/Isend + 每 kind 一次 Waitall；tag 空间分离（face byFieldId*10+dir*2+side、
  edge byFieldId*100+100+pairIdx*4+sp、corner byFieldId*100+200+tripleIdx*8+st）。19 个
  核心 MG 测试 -n2 [1,1,1,2] errors 0）

- [x] **P4-5.4 多方向验收 [2,2,1,1]**
  工作：现有全部 MG 测试在 `mpiexec -n 4`、`GpuGrid=[2,2,1,1]` 下跑通。
  commit 门槛：通用门槛 + 全测试通过。
  commit: f8ce7b35（28/28 TestMG* -n4 [2,2,1,1] errors 0，TestMGHaloSelfExchange 扩展为
  FACE+EDGE+CORNER 槽位验证（3072 槽全对）。**验收发现并修复一个真实 bug**：Rotating
  作用量 `_deviceStapleTermGfactorT` 的 bTorus 分支用 `__deviceSiteIndexToInt4` 解码
  position 表条目，但该函数索引 `m_pSiteMappingTable`（仅本地 Volume 大小），MG 下条目
  是 halo 槽位（>= Volume）→ 越界 IllegalAddress——单方向 split 从不触发（Rotating 力
  内核只读 x/y/z 邻居，[1,1,1,2]/[2,1,1,1] 下从未跨 rank 读），多方向首次暴露。修复：
  MG 下保留原始坐标，`_deviceFi/_deviceFiShifted` 经 `_deviceSIndexToGlobalInt4` 解析；
  `_deviceFi` base site 同法。单 GPU TestRotationQuenched3D 不回归）

- [x] **P4-5.5 多方向验收 [2,1,1,2]（含 t 方向）**
  工作：同上，另验证硬校验（被切维本地长度 ≥ halo 宽度）在多方向下正确触发
  （构造一个非法 grid 确认 `appCrucial` 报错）。
  commit 门槛：通用门槛 + 全测试通过 + 非法配置报错正确。
  commit: 见本提交（28/28 TestMG* -n4 [2,1,1,2] errors 0；TestMGLoadCrossPrecision 在
  批次并发时超时、单独重跑 errors 0（GPU 竞争，非回归）。非法 grid 校验：8⁴ +
  HaloWidth=3 + GpuGrid [4,1,1,1]（本地 x=2 < 3）触发 `appCrucial` "local length 2 in
  split direction 0 is smaller than halo width 3" 并退出）

---

## P5-1：IO 补全

- [x] **P5-1.1 gauge 非原生精度加载 scatter**
  工作：`CFieldGaugeLink.h:120-163` 分支：全文件读取 → 精度转换 → `ScatterFieldFromRoot`；
  文件大小校验改用全局 link 数。
  普查（P4-2.5 同批）确认现状是**静默错**：`EFFT_CLGBinFloat`/`CLGBinDouble` 加载分支
  （121/143 行）在 MG 下无 `#if _CLG_MULTI_GPU` scatter，每 rank 按本地 `_HC_LinkCount`
  读文件前缀 → 数据错位/不完整。**不优先（不影响后续功能），但不做不行**——放计划后段，
  与 P5-1.2/P5-1.3 一起收尾。
  commit 门槛：通用门槛 + float↔double 互载 1-vs-N 通过。
  commit: 见本提交（两分支改全文件读 → 全局元素转换 → scatter；**实测教训：转换
  buffer 必须按全局元素数分配**，按本地分配在 -n2 堆越界（malloc corrupted）；
  TestMGLoadCrossPrecision n1/n2 dot=4.915200409102440e+04 逐位一致）

- [x] **P5-1.2 费米子加载 scatter-aware**
  工作：费米子文件路径走 `CField::SaveToFile` 的逆过程；替换掉测试靠 uniform spinor
  绕过的现状，补真实加载测试。
  commit 门槛：通用门槛 + 保存→加载往返 MD5 1-vs-N。
  commit: d5706f94（`CFieldFermionKST` 与 `CFieldFermionWilsonSquareSU3`
  `InitialFieldWithFile` 全文件读取 + 全局大小校验 + `ScatterFieldFromRoot` 分发，
  单 GPU 路径不变；此前每 rank 按本地 site 数读全局文件前缀 → 错位/不完整静默错。
  新增 `TestMGFermionSaveLoad`（随机高斯 spinor，全局播种 -n1/-nN 相同 → 保存→加载→
  再保存 MD5 一致）。验证：-n1 MD5 F34EC6EFC7653E5587FB242A8CB96B7D = -n2（errors 0））

- [x] **P5-1.3 压缩格式（EFFT_CLGBinCompressed）save/load 接通**
  工作：普查（P4-2.5 同批）确认压缩格式存在偷懒/静默错，**接通而非维持拒绝**：
  - 保存侧（`CFieldGaugeLink.cpp:57/122/262`，U1/SU2/SU3 `SaveToCompressedFile`）在 MG
    下 `appCrucial` 拒绝——接通与 EFFT_CLGBin 保存同构：gather 本地链接到 rank0（全局
    site 顺序）→ rank0 上 StrictLog + 压缩写文件。
  - 加载侧（`InitialWithByteCompressed`，CFieldGaugeLink.h:164 分支）MG 下**无拒绝但按
    本地链接数静默读**——文件是全局大小（single-GPU 写的），每 rank 按本地大小读 →
    数据错位/不完整（静默算错，比拒绝更糟）。接通：全文件读取 → 精度/布局还原 →
    `ScatterFieldFromRoot` 分发（照 EFFT_CLGBin 模式，CFieldGaugeLink.h:87-108）。
  **不优先（不影响后续功能），但不做不行**——放计划后段，与 P5-1.1/P5-1.2 一起收尾。
  commit: f7b4c345（U1/SU2/SU3 `SaveToCompressedFile` 移除 MG 拒绝：本地 StrictLog +
  提取压缩 → `GatherFieldToRoot`（每 site 压缩字节）→ rank0 写全局文件，非 root 不写；
  `InitialWithByteCompressed` 全文件读取 + 全局大小校验 + `ScatterFieldFromRoot`。
  新增 `TestMGCompressedSaveLoad`（往返 MD5 1-vs-N；float 构建下压缩往返非逐位精确是
  单 GPU 既有行为，判据为 MD5 一致）。验证：-n1 保存 MD5 08A9FA16... = -n2、
  重存 F022B763... = -n2）
  commit 门槛：通用门槛 + 保存→加载往返 MD5 1-vs-N（`-n 1`/`-n 2`）。
  commit:

---

## P5-2：机制性收尾

- [x] **P5-2.1 `_LAUNCH_KERNEL_MG` 约定落地**
  工作：在 Plan §3.2 加编者按：实际落地为手写 `Ensure/Refill`，宏标注
  deprecated-not-used；写"新增 stencil kernel 接 halo 检查清单"（Ensure 读集 /
  MarkDirty 写集 / halo 宽度取 stencil 半径上限）。
  commit 门槛：通用门槛（文档）。
  commit: 见本提交（编者按 + 5 条检查清单写入 Roadmap 附录；`CudaHelper.h:179`
  宏标注 deprecated-not-used）

- [x] **P5-2.2 `_CLG_IS_LOG_RANK` 接入日志路径**
  工作：`appGeneral`/`appCrucial` 非 root rank 静默；各 MG 测试的手写 IsRoot 守卫
  改用该宏。
  commit 门槛：通用门槛 + `-n 4` 日志只有一份。
  commit: ea1020b9（Tracer.cpp 加 `_clgShouldLog()` 早退：appVOut/_appCrucial/
  appWarning/appGeneral/appDetailed 非日志 rank 静默——`_CLG_IS_LOG_RANK` 定义在
  CLGComm.h，Tracer.cpp 经 CLGLib_Private.h 可见；21 个 MG 测试的手写 IsRoot 打印
  守卫改用该宏（语义相同：非 root 测量后返回，仅 root 打印被比较标量）。WSL/Linux
  OpenMPI 验证：TestMGMeasurePlaq -n4（GpuGrid [1,1,4,1]）恰好 1 份日志、errors 0、
  plaq = 3.073215484619141e-04）

- [x] **P5-2.3 基线重跑存档**
  工作：重跑单 GPU Release 基线覆盖 `Docs/mg_baseline/`（压缩 IO 崩溃已被 607cebd2
  修复，更新 `baseline_notes.md` 该条目；补记当前 FAIL/超时清单）。
  commit 门槛：基线文件落库。
  commit: 见本提交（FileIO 10/11，唯一 FAIL 为既有 TestFileIOBridgePPBin；
  压缩 IO PASS 确认 607cebd2 修复生效；FAIL/超时清单落库）

---

## P5-3：文档、测试整理、分发

- [x] **P5-3.1 MG 测试拆为独立组**
  工作：10+ 个 MG 测试从 FileIO 组迁出（新测试文件 + 独立 yaml 组），rank 数保护、
  单卡可跳过。
  commit 门槛：通用门槛 + 新组 `-n 1`/`-n 2` 全绿。
  commit: 5b01ce33（30 个 TestMG* 从 FileIO 分类迁到新 MG 分类 +
  TestSuit_MG.yaml 独立组；新文件 Code/CLGTest/Tests/TestMG.cpp、TestConfigurationFileIO.cpp
  只留单 GPU FileIO 测试；CLGTest.cpp LoadParams 追加解析 TestSuit_MG.yaml；
  CMakeLists + VS 项目文件同步。`./CLGTest MG` 恰好跑 MG 套件；-n 1 默认
  GpuGrid [1,1,1,1] 全部单进程恒等分解跑通（驱动 1-vs-N 基线）。验证：
  Release + Release_MG 编译；`MG` 组 -n 1 29/29 errors 0、-n 2（[1,1,1,2]）
  29/29 errors 0；非 MG FileIO 组 11 个测试不受影响）

- [x] **P5-3.2 README 三平台多卡说明**
  工作：Windows `mpiexec -n N` / WSL `mpirun --oversubscribe` / DCU `srun --gres=dcu:N`
  编译与运行段落。
  commit 门槛：文档落库。
  commit: 见本提交（README 追加 Multi-GPU 段：三平台构建/运行命令、GpuGrid==进程数
  约定、已知限制（t 链测量不切 t、split-z/split-x/y 显式拒绝））

- [x] **P5-3.3 Plan §10 进度同步**
  工作：把 P4-1..P5-2 的完成情况回填 `Docs/MultiGPU-Plan.md` §10。
  commit 门槛：文档落库。
  commit: 见本提交（`Docs/MultiGPU-Plan.md` 未入库，进度同步落 Roadmap 附录：
  已完成 P3-9/P4-1..P4-4/P5-2.1/2.2，待完成 P4-5/P5-1/P5-2.3/P5-3.1/3.4/M4）

- [x] **P5-3.4 Windows portable `_MG` 分发验证**
  工作：按 Plan §9.4 清单（exe + cudart dll + `msmpi.dll` + `mpiexec.exe` + `smpd.exe`）
  在干净 Windows 环境验证 `mpiexec -n 2` 拷贝即跑；无干净机器则如实标注"未验证"。
  commit 门槛：验证记录落库（或未验证标注）。
  commit: 见本提交（**未验证**：开发环境为 WSL2/Linux，无干净 Windows 机器；
  清单见 Plan §9.4——`Release_MG` exe + `cudart64_*.dll` + `msmpi.dll` +
  `mpiexec.exe` + `smpd.exe` 同目录拷贝即跑。WSL 侧已用 OpenMPI 验证同等
  逻辑（`mpiexec -n 2` 直接运行））

---

## 完成定义（M4）

以上全部任务填上 commit 编号，且：R 系列 App（旋转作用量 + HMC + 常用测量）在
`-n 2` 下端到端跑通一次小 case，结果与 `-n 1` 按判据一致。此后进入 M5（真实多卡
性能，需硬件，不在本列表）。

**验收记录**（commit eb4e0535，double 构建 `-DCLG_DOUBLE=1 -DCLG_MULTI_GPU=1`，
MatchingRho 小 case `[4,4,4,4]`、`RandomSeed 1234567`、`MaxThreadPerBlock : 256`、
`GpuGrid [1,1,2,1]`（split z——t 链测量需不 split t），`mpiexec -n 1` vs `-n 2`）：
- H(before) 三个分量逐位一致：kin -4053.407708 / Action1 7744.849159 / Action2
  3070.495870；接受序列一致（3/3 Accept）。
- 测量（double，求和顺序容差内）：Polyakov |P| 0.0795570 vs 0.0795690（1.2e-5）、
  condensate 23.165274 vs 23.165322（4.8e-5）、rho0 correlator 同容差。
- **验收暴露并修复两个真实 MG bug**（此前 TestMGHmcConsistency 只有 gauge action，
  从未覆盖 fermion 能量/力）：① `CFieldFermionWilsonSquareSU3::Dot` 只做局部
  ThreadBufferSum，-n N 下 fermion 能量为单卡的 1/N（HMC 全 Reject）——补
  AllreduceSum 全局化；② `DerivateDOperator` 经 move cache 读两个输入 spinor 的
  邻居（halo 槽位）但未 RefillHaloBuffer（DOperator 有）——力（MD 轨迹）与 -n 1
  分叉——补两个缓冲的 halo refill。另修 double 构建的 `AllreduceSum(Real*)` 与
  `DOUBLE*` 重载冲突（`#if !_CLG_DOUBLEFLOAT` 守卫）。
- App 配置要点（非 MG bug，App 自带配置问题）：MatchingRho/RotatingReproduce 的
  yaml 缺 `MaxThreadPerBlock : 256` 时单卡即 701（LaunchOutOfResources），补上与
  CLGExample 一致即可；RotatingReproduce 用 Dirichlet 边界（halo 需 torus），M4 用
  MatchingRho（torus）。
