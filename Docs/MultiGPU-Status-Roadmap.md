# 多 GPU 支持：进度总结、可行性评估与后续实施计划

> 日期：2026-07-30。分支 `multi-GPU`（领先 `gitee/multi-GPU` 12 个提交）。
> 本文是 `Docs/MultiGPU-Plan.md`（下称"计划"）的配套评审文档：总结已落地的部分，
> 评估继续实施的可行性，列出需要改进/补齐的方案，并给出剩余工作的详细实施计划。
> 计划文档的术语（rank、halo、本地/全局坐标、"逐比特 vs 数值一致"等）沿用 §8.3。

---

## 1. 项目进度总结

### 1.1 总体状态

| 里程碑（计划 §7） | 状态 | 证据 |
|---|---|---|
| M1：Phase 0+1 分解/自交换/gather-scatter | ✅ 达成 | commit `0057b6ca`–`bff39871`；`TestMGGatherScatterRoundTrip`、`TestMGHaloSelfExchange` |
| M2：Phase 2 纯规范 stencil 多卡一致 | ✅ 达成 | commit `d1fe7066`；`TestMGGaugeForceConsistency` MD5 `9E8DEA66...` 1-vs-N 逐比特 |
| M3：Phase 3 带费米子完整 HMC 一致 | 🔶 核心达成 | 费米子 D（Wilson/KS/HISQ）与纯规范 HMC 已验证；**带费米子的完整 HMC 轨迹尚未做 1-vs-N 验收** |
| M4：Phase 4+5 规范固定/测量/IO/文档 | ❌ 未开始（IO 部分有基础） | 见 §3 缺口清单 |
| M5：Phase 6 性能优化 | ❌ 未开始 | 需真实多卡硬件 |

**一句话结论**：通信与分解的"地基+主体"已完工并经 1-vs-N 数值验收；剩余是
"按子系统逐个接线"的扫尾（测量、规范固定、boson、位置相关物理）与多方向切分、IO 补全。

### 1.2 各 Phase 落地情况

**Phase 0 — 脚手架**（完成）：`_CLG_MULTI_GPU` 开关；`CLGComm`（拓扑/Allreduce/
Gather/Scatter/Barrier）；`CHaloManager`；VS `Debug_MG`/`Release_MG` 配置（14 个
vcxproj 脚本化生成，MS-MPI 引用条件化）；CMake `-DCLG_MULTI_GPU=1`；MS-MPI 本机验证。
开关关闭时单 GPU 行为零变化（基线文档 `Docs/mg_baseline/baseline_notes.md` 据此建立）。

**Phase 1 — 分解/索引/字段**（完成）：每 rank 子格点常量就地改写；extended flat
buffer halo 布局（Design B，`Core/Distributed/CLGHaloLayout.h`）；越界邻居 `SIndex`
重定向到 halo 槽（kernel 透明）；gather∘scatter == 恒等（MD5 1-vs-N 逐比特）；
单卡 self-exchange 单元测试。全局坐标基础设施（`_DC_Offset*`、
`__deviceLocalToGlobalInt4`）已备——但**只有 RNG 和 halo bake 在用**（见 §3-缺口 1）。

**Phase 2 — halo 交换接入 stencil**（完成）：gauge staple/plaquette/force、Wilson D、
plain staggered D、HISQ/Naik D（HaloWidth=3）、Fat357Lepage/Naik smearing，全部
1-vs-N gather-MD5 逐比特一致（`TestMGGaugeForceConsistency`、`TestMGWilsonDConsistency`、
`TestMGStaggeredDConsistency`、`TestMGHISQDConsistency`）。
关键设计兑现：fermion halo refill 挂在不可覆盖的 `CFieldFermionKS::DOperator`
chokepoint，一处覆盖所有 KS 子类；smearing 走 baked mapping table 零 kernel 改动。

**Phase 3 — 全局归约/求解器/RNG/HMC**（代码完成，验收部分完成）：
- `CCommonKernelField::Dot/LengthSq/Sum` 经 `_clgGlobalThreadBufferSum` 全局化
  （`CFieldCommonKernel.cu`）→ **SparseLinearAlgebra 全部求解器自动继承**（该目录
  无任何独立归约路径，已核实）。
- HMC 能量与 Metropolis 全局化（`CHMC.cpp`、`CActionGaugePlaquette.cpp`）；
  `TestMGHmcConsistency`：纯规范短轨迹**接收数精确一致**，rmsHDiff/plaquette 容差一致。
- RNG 全局播种（`Random.cu` 按全局 site index），`TestMGRngGlobalSeed` 逐比特分解不变。
- integrator-copy 的 gauge halo 修复（commit `49a68b58`，HMC energy+force 1-vs-N 验证）。
- **未做的验收**：带费米子（HISQ）的完整 HMC 轨迹 1-vs-N——M3 的最后一公里。

**额外收获**：基线文档记录的 `TestFileIOCLGCompressed`（32⁴）既有崩溃，已在
commit `607cebd2` 修复（三个先于多 GPU 存在的 OOB/未初始化 bug：`m_byFieldId` 未初始化、
`StrictExp` 越界读 `m_me[9]`、`InitialWithByteCompressed` 堆溢出）。`Docs/mg_baseline/
baseline_notes.md` 中该条分类已过时，基线应重跑一次存档（见 §4-Task P5-3）。

### 1.3 验证基座

单卡 oversubscribe（`mpiexec -n N` 全指 device 0）方法学运转良好，10 个 MG 测试
全部注册在 FileIO 组（`TestConfigurationFileIO.cpp`）。已知操作陷阱（MS-MPI stdout
缓冲造成的"假死锁"、改头文件后需 touch 相关 `.cu` 重编）已记录在计划 §10。

---

## 2. 继续实施的可行性评估

**结论：高度可行。** 理由：

1. **架构赌注已全部兑现**。计划 §2.6 的三个关键判断——邻居查表使 kernel 零改动、
   算子入口收敛（chokepoint）、extended flat buffer 使 halo 对 kernel 透明——在
   gauge/Wilson/KS/HISQ/smearing 五个子系统上逐一验证成功。剩余子系统（测量、boson、
   规范固定）的耦合模式与已完成的同类，没有新的技术风险。
2. **通信层足够简单且够用**。同步 `cudaMemcpy` + host-staging + 非阻塞 MPI 的实现
   已被全部 10 个测试覆盖；CUDA-aware/NCCL 只是性能选项，不挡正确性。
3. **回归纪律可执行**。`_CLG_MULTI_GPU=0` 下所有 MG 代码编译出局，单 GPU 路径零改动；
   每个 Task 的验收都可以在本机单卡 oversubscribe 下完成，无需等硬件。

**真正的风险不在通信，而在"漏接"**：MG 接入目前是**手动的**（`_LAUNCH_KERNEL_MG`
宏定义了但零调用点，实际是各调用点手写 `Ensure/RefillHaloBuffer`）。凡是"写代码的人
没想到要接"的子系统就静默算错——这正是 §3 缺口清单的来源。因此后续计划的核心原则是
**按子系统枚举、逐个接线、每个配 1-vs-N 测试**，而不是再改架构。

**排序考量**：
- 缺口 1（R1 位置相关物理）是 R 系列 App（本项目主力）多卡化的**前提**，优先级最高。
- 缺口 3（规范固定）测试已先行（TDD），实现路径明确（计划 §4 Phase 4 选项 A）。
- 缺口 5（corner/edge halo）只影响多方向 `GpuGrid`；单方向 `[N,1,1,1]` 已正确。
  多方向是扩展到更多卡时降低通信/计算比的关键（计划 §9.2），但不阻塞 2–4 卡使用。
- Phase 6（stream 重叠）必须等真实多卡硬件，本机无法验证，不在本计划内排期。

---

## 3. 需要改进/补齐的方案（缺口清单）

> 按严重度排序。每条附代码证据。

### 缺口 1（严重）：R1 位置相关物理未接全局坐标
计划 §1.4-R1 要求 Phase 1 完成，实际只落了基础设施：
`__deviceLocalToGlobalInt4()`（`CCommonData.h:599`）全仓库**零调用点**。
受影响的文件（均 0 处 `_CLG_MULTI_GPU`/`_DC_Offset`）：
- 作用量：`CActionGaugePlaquetteRotating/Rotating3D/RotatingT/RotatingT3D/Acceleration/
  RigidAcc/Boost/Cylinder`（如 `CActionGaugePlaquetteRotating.cu:51` 直接用本地
  `sSite4.x - _DC_Centerx`）。
- 测量：`CMeasureRotatingAction`、`CMeasureAngularMomentum*`、`CMeasurePandChiralTalor(KS)`。
后果：多 rank 下旋转/加速度类物理**静默算错**且 plaquette 类测试暴露不出。

### 缺口 2（严重）：测量子系统整体未接入
`Code/CLGLib/Measurement/` 全部 27 个 `CMeasure*` 类 0 处 MG 代码：
- 直接调 `ThreadBufferSum/ReduceReal` 的（`CMeasureChiralCondensate`、
  `CMeasureWilsonLoopWithPath`、`CMeasurePolyakovXY` 等）：`CFieldCommonKernel.cu:26`
  注释明确**故意不**全局化 `ThreadBufferSum` → 整格和约为真值的 1/N。
- 时间片 profile 类（介子关联函数等）：需要按 t 切片跨 rank 归并，完全未处理。
- `CMeasurePlaqutteEnergy` 调用的 `CalculatePlaqutteEnergy` 只做本地和（测试里的
  Allreduce 是手动加的），HMC 外的 plaquette 测量会错。
- 位置相关测量叠加缺口 1，双重错误。

### 缺口 3（高）：Phase 4 规范固定实现缺失（测试先行）
`TestMGGaugeFixingConsistency`（`TestConfigurationFileIO.cpp:619-633`）的注释描述了
"gather 到 rank0 → 临时全局格点上下文 → fix → scatter 回来"的实现，但该实现**不存在**；
当前测试实际在各 rank 本地子格点上跑单卡 fixer，多 rank 语义错误。
cuFFT 用户：`CGaugeFixingCoulombCornell.cu`、`CGaugeFixingLandauCornell.cu`；
其余（LosAlamos/MAG/MCG*/Random）不用 FFT，但同样需要全局格点。

### 缺口 4（高）：boson 字段无 halo 存储
`CFieldBoson*` 的 kernel 是 stencil 型（读 `m_pMoveCache`），但 halo 尾部扩容只做了
gauge/Wilson/KS 三类字段；`CActionPhi4::Energy` 也无 Allreduce。多 rank 下 Phi4 类
物理静默算错。

### 缺口 5（中）：corner/edge halo 未填，多方向 GpuGrid 不可用
`_deviceHaloRedirectSite` 只处理 face 单元，corner/edge 保持周期 wrap
（`CLGHaloLayout.h:92-101`）。`[1,1,1,N]` 正确；`[2,2,1,1]` 等多方向切分会错。
计划 §9.2 已把多方向列为"跑通后启用"，现状与此一致，但扩展性叙事（§9.2 通信量
提示）依赖此项完成。

### 缺口 6（中）：IO 路径不完整
- gauge **非原生精度**加载分支（`CFieldGaugeLink.h:120-163`）无 scatter，且按本地
  `_HC_LinkCount` 校验文件大小 → MG 下加载全局文件 size mismatch。
- 压缩格式：`SaveToCompressedFile` MG 下直接拒绝（`appCrucial`），
  `InitialWithByteCompressed` 无 scatter。
- 费米子文件加载不 scatter-aware（测试靠 uniform spinor 绕过）。
- 保存路径已正确（`CField::SaveToFile` gather 到 rank0 写文件）。

### 缺口 7（低）：机制性欠账
- `_LAUNCH_KERNEL_MG` 宏零调用点：实际靠手写 `Ensure/Refill`。要么按 §3.2 接入
  （防"新 stencil kernel 忘记接 halo"），要么在文档中明确"手写 Ensure"是既定约定
  并把宏降级为文档。二选一，当前状态最糟（宏在、没人用、新人会误以为有自动保障）。
- `_CLG_IS_LOG_RANK` 定义后零调用点：非 root rank 仍在刷日志，N 大时日志不可读。
- `Docs/mg_baseline/baseline_notes.md` 过时（缺口已修复），基线需重跑存档。

---

## 4. 详细实施计划

> 每个 Task 遵循计划 §8.2 的强制回归：`Release`（MG=0）构建 + CLGTest 对基线无新增
> FAIL → `Release_MG` `-n 1` → `-n 2`（小格点 oversubscribe）。全部验收在本机单卡完成。
> 浮点判据按 §8.1-A1：浮点容差（double 相对 1e-10）、整数/MD5 逐比特。

### Task P4-1：R1 全局坐标接入（预计 2–3 天）

**目标**：位置相关作用量/测量在多 rank 下用全局坐标，`GpuGrid=[1,1,1,2]` 与 `-n 1`
数值一致。

**步骤**：
1. 在 `CCommonData.h` 的 `__deviceLocalToGlobalInt4` 基础上，为 kernel 内使用提供
   简洁形式（如 `_deviceSiteIndexToGlobalInt4(uiSiteIndex)`），未切方向直接加 offset，
   周期方向对全局 `Lμ` 取模。**不改 `_DC_Center*` 语义**——它本来就是全局量。
2. 逐个改造作用量 kernel（`CActionGaugePlaquetteRotating*.cu`、`Acceleration/RigidAcc/
   Boost/Cylinder`）：把进物理公式的 `sSite4` 换成全局坐标。用 `#if _CLG_MULTI_GPU`
   保护或在单 rank 下 offset=0 自然退化（offset=0 时加 0，单 GPU 逐比特不变，
   优先选无条件写法减少分叉）。
3. 同步改造 `CMeasureAngularMomentum*`、`CMeasureRotatingAction` 的坐标来源。
4. **验收**：新增 `TestMGRotatingActionConsistency`：8⁴、`CActionGaugePlaquetteRotating3D`
   非零 Ω，能量与力场 gather-MD5 1-vs-N 一致。这是计划 §4 Phase 1 遗留的 R1 验收。

**依赖**：无。**这是所有 R 系列 App 多卡运行的前置条件。**

### Task P4-2：规范固定 gather-to-rank0（预计 2–3 天）

**目标**：实现计划 §4 Phase 4 选项 A，使已存在的 `TestMGGaugeFixingConsistency` 真正通过。

**步骤**：
1. 在 `CGaugeFixing` 基类（或 `CLGLibManager` 层）加 MG 分支：
   `GatherFieldToRoot` → rank0 用**全局格点几何**构造临时 fixer 上下文（需要一套
   全局尺寸的 `m_ConstIntegers`，最简做法：rank0 单独 new 一个按全局格点初始化的
   子管理器/辅助对象，fix 完销毁，显存峰值只在该步）→ `ScatterFieldFromRoot`。
2. 非 root rank 在 fix 期间阻塞等待 scatter（`Barrier` 或直接靠 scatter 的通信同步）。
3. fix 后 `MarkDirty` gauge 场（halo 过期）。
4. **验收**：`TestMGGaugeFixingConsistency` `-n 1`/`-n 2`：Coulomb 偏差收敛值容差一致、
   fix 前后 plaquette（规范不变量）一致。先 CoulombLosAlamos（无 FFT，路径短），
   再 CoulombCornell（cuFFT，验证 FFT 在全局上下文下可用）。

**注意**：分布式 FFT（选项 B）按计划 §8.4 需显式授权，本 Task 不做。

### Task P4-3：测量子系统接入（预计 3–5 天）

**目标**：常用测量在多 rank 下给出正确的全局值。

**步骤**：
1. **整格标量和**（工作量小、收益大，先做）：给 `CMeasure` 层提供一个全局化归约辅助
   （如 `CMeasure::GlobalSum(...)`，内部 `ThreadBufferSum` + `CLGComm::AllreduceSum`），
   逐个替换 `CMeasurePlaqutteEnergy`、`CMeasureChiralCondensate(KS)`、
   `CMeasureTopologicChargeXY`、`CMeasurePolyakovXY`、`CMeasureBosonValue` 等的本地和。
   `CMeasurePlaqutteEnergy` 所依赖的 `CalculatePlaqutteEnergy` 走 Task 已在 HMC 验证的
   全局化路径或在此补 Allreduce（二选一，避免双重归约——先查清 HMC 里怎么做的再定）。
2. **时间片 profile**（介子关联函数、`PandChiralTalor`）：本地按 t 切片求和 →
   对长度为 Lt 的数组做一次 `AllreduceSum(DOUBLE*, count)`（API 已存在）。
3. **Wilson loop / 非规则路径**（`CMeasureWilsonLoopWithPath`）：路径可能跨 rank 边界——
   按计划 §8.4，遇到非规则跨 rank 拼接先停下来与用户确认范围，首版可只支持
   "路径在 rank 本地或可分解为整格归约"的子集并明确报错。
4. **验收**：每个接入的测量加 1-vs-N 对照（小格点、固定组态，容差一致）。
   接入顺序按用户实际使用频率排：plaq energy → Polyakov → chiral condensate → 介子。

### Task P4-4：boson 字段 halo（预计 1–2 天）

**步骤**：
1. 按 gauge/Wilson/KS 的既有模式给 `CFieldBoson*` 的分配追加 halo 容量
   （`CLGHaloLayout.h` 的 site 版布局直接可用）。
2. boson stencil kernel（`CFieldBosonVNKernel.cu` 等）的调用点前加
   `Ensure(bosonFieldId, 1)` + 写后 `MarkDirty`。
3. `CActionPhi4::Energy` 加 Allreduce。
4. **验收**：`TestMGBosonConsistency`（Phi4 能量 + D 类 stencil 输出 gather-MD5 1-vs-N）。

### Task P4-5：corner/edge halo + 多方向 GpuGrid（预计 3–4 天）

**步骤**：
1. 扩展 `CLGHaloLayout.h`：在 face 块之后追加 edge/corner 槽位，编号约定补进文档
   （与 face 块同一套 (dir,side,w,faceIdx) 体系的延伸，bake 与 pack/unpack 逐字一致）。
2. `_deviceHaloRedirectSite` 处理多方向同时越界的站点；bake（`CIndexSquare.cu`、
   `CBoundaryConditionTorusSquare.cu`）生成 corner/edge 重定向。
3. `RefillHalo` 的交换序列扩展为 4 维网格的 face→edge→corner 传递（或直接按
   26 邻域一次性交换——维度低时两者等价，取实现简单的）。
4. **验收**：现有 10 个 MG 测试在 `GpuGrid=[2,2,1,1]`（`-n 4`）下全部通过；
   再抽 `[2,1,1,2]` 验证含 t 方向的多方向组合。计划 §9.2 的硬校验（被切维本地长度
   ≥ halo 宽度）已存在，补测试确认它在多方向下正确触发。

### Task P5-1：IO 补全（预计 2–3 天）

1. gauge 非原生精度加载：先读全文件到各 rank → 转换 → `ScatterFieldFromRoot`
   （校验改为全局 link 数）。
2. 费米子加载 scatter-aware（`CFieldFermion` 的文件路径走 `CField::SaveToFile` 的
   逆过程）。
3. 压缩格式：首版维持"MG 下拒绝"但把 `appCrucial` 信息写清楚（指引用非压缩格式）；
   或按 `607cebd2` 修复后的路径补 scatter。取工作量小的。
4. **验收**：`TestSaveLoad*` 的 MG 变体（保存→加载往返 MD5，1-vs-N）。

### Task P5-2：机制性收尾（预计 1 天）

1. `_LAUNCH_KERNEL_MG` 二选一：推荐**文档化手写 Ensure 约定**（在计划 §3.2 加编者按，
   说明实际落地形式与"新增 stencil kernel 必须手动接 halo"的检查清单），把宏标注为
   deprecated-not-used；真改成自动宏收益小、触碰面大。
2. `_CLG_IS_LOG_RANK` 接到 `appGeneral`/`appCrucial` 日志路径，非 root rank 静默。
3. 重跑基线并存档（覆盖 `baseline_notes.md` 中已修复的压缩 IO 条目；补记
   `TestFileIOCLGCompressed` 现已 PASS）。

### Task P5-3：文档、分发、CI（预计 2 天，部分可与上面并行）

1. 更新 `Docs/MultiGPU-Plan.md` §10（本轮已补 Phase 3 记录）。
2. README 补三平台多卡编译/运行说明（`mpiexec -np` / `srun --gres=dcu:N`）。
3. Windows portable `_MG` 分发清单验证（计划 §9.4；需一台干净 Windows 机，
   无则标注"未验证"）。
4. MG 回归测试整理进 CLGTest 的独立组（现在挂在 FileIO 组），rank 数保护、单卡跳过。

### 不在本期排期

- **Phase 6 性能优化**（stream 重叠、interior/boundary 拆分）：需真实多卡硬件验证。
- **分布式 FFT**（Phase 4 选项 B）：需显式授权（计划 §8.4）。
- **多机跨节点**：Linux 超算侧按 §9.4 现场编译即可，Windows 跨节点不做。

### 总排期估算

| 阶段 | 内容 | 预估 |
|---|---|---|
| P4-1 | R1 全局坐标 | 2–3 天 |
| P4-2 | 规范固定 | 2–3 天 |
| P4-3 | 测量接入 | 3–5 天 |
| P4-4 | boson halo | 1–2 天 |
| P4-5 | corner/edge + 多方向 | 3–4 天 |
| P5-1/2/3 | IO/机制/文档 | 4–6 天（部分可并行） |

P4-1 → P4-2 → P4-3 是向"R 系列 App 多卡生产"的最短路径，建议按此顺序；
P4-4/P4-5 按用户实际物理需求插排。全部完成后 M4 达成，进入等硬件的 M5。

---

## 附录：Measurement 目录测量类 MG 普查（P4-3.10，2026-08-05）

普查范围：`Code/CLGLib/Measurement/` 全部测量类。三分类：已接（MG 下结果全局正确，
有 1-vs-N 验证）/ 待接（实现就绪但尚未做 MG 归约，需按标注的同类模式补齐）/
不接（明确理由）。沿 t 的链乘积类（Polyakov、WilsonLoop 的 t 方向）要求**t 不切分**
的 grid；沿切分方向的链/剖面需跨 rank 拼接，由 halo 或显式拒绝处理。

### 已接（有 1-vs-N 测试）
| 类 | 处理 | 验证 |
|---|---|---|
| CMeasureRotatingAction | P4-1.6 坐标全局化 + 能量/分量 Allreduce | MGRotatingActionMeasure |
| CMeasureAngularMomentumKS / KSREM | P4-1.6 杠杆臂坐标全局化（12 kernel 26 处） | P4-1.6 记录 |
| CMeasureAMomentumJG | P4-1.6 同法 | P4-1.6 记录 |
| CMeasurePlaqutteEnergy | P4-3.2 GlobalSumReal + 全局 PlaqutteCount | MGMeasurePlaq |
| CMeasurePolyakovXY | P4-3.3 XY 密度/loop/slice 全局化 + 全局归一化 | MGPolyakovXY |
| CMeasureChiralCondensate / KS | P4-3.4 求和全局化 + 全局体积 + XY 分布/切片 | MGChiralCondensate |
| CMeasureMesonCorrelatorStaggeredSimple2 | P4-3.6 t/x/y/z 剖面数组全局化 | MGMesonCorrelatorSimple2 |
| CMeasurePandChiralTalor / KS | P4-3.7 全局坐标 + Polyakov/Omega 归约 | MGPandChiralTalor |
| CMeasureTopologicChargeXY | P4-3.8 电荷/密度全局化 + clover halo refill | MGTopologicChargeXY |
| CMeasureBerryPhase | P4-3.8 每-t 剖面数组全局化 | （与 TopologicChargeXY 同批） |
| CMeasureWilsonLoopWithPath | P4-3.9 halo refill + 路径跨度检查 + 全局化 | MGWilsonLoopWithPath |

### 待接（未做 MG 归约，按标注模式补齐）
| 类 | 模式 |
|---|---|
| CMeasureMesonCorrelator | 同 Simple2（P4-3.6）：时间/空间剖面数组全局化 |
| CMeasureMesonCorrelatorStaggered | 同 Simple2（P4-3.6） |
| CMeasurePolyakovXY3D | 同 PolyakovXY（P4-3.3） |
| CMeasureConnectedChiralSusceptibilityKS | 同 ChiralCondensateKS（P4-3.4） |
| CMeasureWilsonLoop / WilsonLoopXY / WilsonNality | t 方向链乘积：t 不切分时本地归约全局化（P4-3.3 模式）；矩形路径跨 rank 拼接同 WilsonLoopWithPath（P4-3.9） |
| CMeasureAMomentumStochastic | 随机源测量，依赖 Z4 + 求解器；全局化模式同 ChiralCondensate（P4-3.4），但验证需求解器精度容差 |
| CMeasureChargeAndCurrents | 归约点普查未完成，按测量输出类型套用 GlobalSumReal/ComplexArray |
| CMeasureBosonValue | **依赖 P4-4（boson halo）**：boson 场 MG 下 halo 未接通前不接 |

### 不接
无（所有类均可按上述模式接通；暂缓项仅为依赖未完成基础设施的 CMeasureBosonValue）。

---

## 附录：`_LAUNCH_KERNEL_MG` 约定（P5-2.1，2026-08-05）

**编者按**：原 Plan §3.2 设想的 `_LAUNCH_KERNEL_MG` 宏（launch 时自动 halo 交换）
**实际未落地**——所有 MG stencil 路径都是**手写 `Ensure`/`RefillHalo`（读前）+
`MarkDirty`（写后）**。`CudaHelper.h:179` 的 `_LAUNCH_KERNEL_MG` 宏保留但标注
deprecated-not-used；不要在新代码中使用。

**新增 stencil kernel 接 halo 检查清单**（照此接线，否则 MG 下跨 rank 邻居读脏数据）：
1. 确认 stencil 的**读集**：kernel 直接索引 `m_pDeviceData[localIdx ± offset]` 或经
   `__idx->m_pDeviceIndexLinkToSIndex/PositionToSIndex` 重定向到 halo 槽的所有字段；
   读集内的每个 field 都要 `Ensure(fieldId, width)` + `RefillHalo(fieldId)`。
2. 确认 stencil 的**写集**：被 kernel 改写的字段在 launch 后 `MarkDirty(fieldId)`；
   写集字段的 halo 尾段随后续 Refill 自动刷新（除非先被读）。
3. halo 宽度取 stencil 半径上限：单格点 stencil 用 `_HC_HaloWidth`；Naik/HISQ
   类 3-link stencil 必须 `Ensure(..., 3)` 且 yaml 配 `HaloWidth: 3`（见
   TestMGHISQDConsistency）。
4. host 侧读 halo 相关常量一律用 `_HC_` 前缀宏（`_DC_` 是设备常量，host 侧恒 0，
   P4-3.9 实测教训）。
5. 完成接线后必须跑一次该字段的 1-vs-N 测试（split 方向覆盖跨 rank 邻居路径）。

---

## 附录：执行进度同步（P5-3.3，2026-08-05）

> 原计划 `Docs/MultiGPU-Plan.md` 未入库（仓库只有本 Roadmap），§10 进度同步
> 落在此处。跟踪文件 `multi-gpu.md` 是唯一执行清单，逐任务带 commit 回填。

已完成并验证（commit 见 multi-gpu.md）：
- P3-9（M3 收尾）、P4-1 全部（坐标/halo 基础）、P4-2 全部（gauge fixing，
  含 P4-2.5 fix Random 接通）、P4-3 全部（测量子系统）、P4-4 全部（boson halo）、
  P5-2.1（`_LAUNCH_KERNEL_MG` 约定）。

待完成：
- P4-5（corner/edge halo + 多方向 GpuGrid，5 项）
- P5-1（IO 补全，3 项）
- P5-2.2（`_CLG_IS_LOG_RANK` 接入）、P5-2.3（基线重跑存档）
- P5-3.1（MG 测试拆组）、P5-3.2（README，见本提交）、P5-3.4（Windows 分发验证）
- M4 完成定义（R 系列 App -n 2 端到端）
