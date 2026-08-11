# CUDA Kernel 宏参考

CLGLib 封装了所有 CUDA kernel 启动和线程分解操作。禁止直接使用 `<<< >>>` 语法。

## 文件位置

| 宏类别 | 定义文件 |
|--------|---------|
| Kernel 启动 | `Core/CudaHelper.h` |
| 线程分解 + Kernel 入口 | `CLGLib_Private.h` |
| 模板参数辅助 | `Core/CLGDefine.h` |

---

## Kernel 启动宏

### `_LAUNCH_KERNEL(func, block, threads, ...)`

标准 kernel 启动。所有参数通过 host 栈打包后传入设备。

```cpp
// 非模板 kernel
_LAUNCH_KERNEL(_kernelMyFunc, block, threads, arg1, arg2);

// 模板 kernel — 必须用 TMPARG 包裹模板参数
_LAUNCH_KERNEL(_kernelMyFunc TMPARG(deviceSU3), block, threads, arg1, arg2);
_LAUNCH_KERNEL(_kernelMyFunc TMPARG(deviceVector, deviceGauge), block, threads, arg1);
```

**规则**：
- `func` — `__global__` kernel 函数指针
- `block` — `dim3` 或 `UINT`，block 数量
- `threads` — `dim3` 或 `UINT`，每 block 线程数
- 变参 — kernel 的参数列表，宏内部自动计数和打包

**底层**：使用 `launchKernel()` 函数 + 4KB host 栈分配参数数组，调用后自动 `FreeStack()`。

### `_LAUNCH_KERNEL0(func, block, threads)`

无参数 kernel 启动。

```cpp
_LAUNCH_KERNEL0(_kernelInit, block, threads);
```

### `_LAUNCH_KERNELS(func, block, threads, sharemem, ...)`

带动态共享内存的 kernel 启动。

```cpp
_LAUNCH_KERNELS(_kernelReduce, block, threads, sizeof(Real) * threads, buffer, count);
```

### `TMPARG(...)`

模板参数包裹宏，定义在 `CLGDefine.h`：

```cpp
#define TMPARG(...) <__VA_ARGS__>
```

**用途**：让 `_LAUNCH_KERNEL` 宏正确解析模板参数列表。

---

## Multi-GPU：唯一推荐写法与自动 guard

多 GPU（`_CLG_MULTI_GPU`）构建下，普通 kernel **保持本节上述写法不变**——
kernel body 不写任何 MPI / rank / GpuGrid / halo 分支，host 侧也不需要 MG 标注。
halo 一致性由 launch 自动 guard 保证：

- `_LAUNCH_KERNEL*` 展开为带 guard 的 `<<<>>>` 启动：registry 按浅层 pointer 参数
  发现已注册 buffer handle，launch 前统一 `Ensure` 读集 halo、`Commit` 时 invalidate
  写集；CUDA 编译器看到的仍是类型完整的 `<<<>>>` 调用，类型检查不丢失；
- 参数表达式只做单次求值（single-evaluation probe），guard 不改变参数语义；
- 未注册的 raw pointer（scratch、常量表等）按 `LocalOnly` 处理，不触发交换。

### `_LAUNCH_KERNEL_RANK_LOCAL(func, block, threads, ...)`

rank 非对称 launch 的**唯一**合法形式（例外分类见 [multi-gpu.md](multi-gpu.md) §8）：与 `_LAUNCH_KERNEL` 走同一 backend
raw 路径、保留 CUDA 类型检查、对写 buffer 做 invalidation，但**绝不**执行 MPI halo
refill。它只接受 registry 标记为 `LocalOnly` / `RankLocalScratch` 的 handle；参数一旦
命中 `HaloCapable` handle 立即 fail-fast。分布式 field 的 rank 非对称更新必须改造成
collective / scatter 语义，不能靠该宏放行。

### 禁止事项

- 业务代码禁止直接使用 `<<<>>>`；raw `_LAUN_KERNEL`（注意拼写）是
  `CudaHelper.cu` 基础设施私有，不对外开放；
- `_LAUNCH_KERNEL_MG`、`_MGEXCH`、field-ID 版 `_R/_W` 标注宏是迁移期产物，
  已从生产调用路径删除（I10），新代码不得使用；
- 嵌套 pointer / 隐藏 pointer / 动态 stencil 等例外场景的完整清单与处置见
  [multi-gpu.md](multi-gpu.md) 的「halo guard 例外清单」。

---

## 线程分解宏（Host 端）

所有宏自动根据格点尺寸和 `_HC_` 常量计算 `block` 和 `threads`。

### 标准分解

| 宏 | 用途 | 生成的变量 |
|----|------|-----------|
| `preparethread` | 全格点分解 | `block = _HC_DecompBlock`, `threads = _HC_DecompThread` |
| `preparethreadHalf` | even/odd 半格点 | `block = _HC_DecompBlockHalf`, `threads = _HC_DecompThreadHalf` |
| `preparethreadDir` | per-link 分解 | `block = _HC_DecompBlockDir`, `threads = _HC_DecompThreadDir` |

```cpp
preparethread;
_LAUNCH_KERNEL(_kernelFunc, block, threads, pGauge);
```

### 单 block 覆盖（`preparethreadE` 系列）

确保所有 sites/links 落在一个 block 内，同时每个 site/link 处理多个元素。

| 宏 | 用途 |
|----|------|
| `preparethreadE(element_count)` | 所有 sites 在一个 block，每 site `element_count` 个元素 |
| `preparethreadEDir(element_count)` | 所有 links 在一个 block，每 link `element_count` 个元素 |
| `preparethreadEVar(element_count, blockvar, threadvar)` | 同上，但自定义变量名 |

```cpp
// element_count = 4 (e.g. 4 个方向各一个 staple)
preparethreadE(4);
_LAUNCH_KERNEL(_kernelMultiElement, block, threads, pData);
```

### 3D 线程分解

用于固定某一维度、在其余三维上并行。

| 宏 | 固定维度 | 并行维度 |
|----|---------|---------|
| `preparethread_S` | `t` | x, y, z |
| `preparethread_Sxyt` | `z` | x, y, t |
| `preparethread_Sxzt` | `y` | x, z, t |
| `preparethread_Syzt` | `x` | y, z, t |

```cpp
preparethread_S;
_LAUNCH_KERNEL(_kernelSlice, block3d, threads3d, pGauge, uiT);
```

---

## Kernel 入口宏（Device 端）

所有宏自动计算 `uiSiteIndex` 并做边界检查（`>= _DC_Volume` 则 `return`）。

### 基础入口

| 宏 | 计算内容 | 额外变量 |
|----|---------|---------|
| `intokernal` | `uiSiteIndex` | 无 |
| `intokernalInt4` | `uiSiteIndex` + `sSite4`（4D 坐标） | `SSmallInt4 sSite4` |
| `intokernalInt4EO` | `uiSiteIndex` + `sSite4` + even/odd 检查 | `sSite4`，需 `bEven` 参数 |

```cpp
__global__ void _kernelFunc(...)
{
    intokernalInt4;
    // 现在可用 uiSiteIndex 和 sSite4
    deviceSU3 link = pData[uiSiteIndex];
}
```

### Even-Odd 半格点

| 宏 | 说明 |
|----|------|
| `intokernalEOHalf` | 半格点（volume/2 线程），自动处理奇偶偏移，提供 `eta` 和 `mask` |
| `intokernalEO` | 全格点线程，但跳过不符合奇偶要求的 site |

`intokernalEOHalf` 关键逻辑：
- 线程覆盖 `volume/2` 个 sites
- 根据 `bEven` 参数自动选择 even 或 odd 子格点
- 若邻居 site 的奇偶性不匹配，自动偏移到正确 site
- 提供 `eta = pEtaTable[uiSiteIndex]`（staggered 相位）

### Per-Link 入口

| 宏 | 计算内容 |
|----|---------|
| `intokernalDir` | `uiSiteIndex` + `dir`（0.._DC_Dir-1） |
| `intokernalDirInt4` | 同上 + `sSite4` |
| `intokernalDirEO` | 同上 + even/odd 检查 |
| `intokernalDirInt4EO` | 同上 + `sSite4` + even/odd 检查 |

```cpp
preparethreadDir;
_LAUNCH_KERNEL(_kernelPerLink, block, threads, pGauge);

__global__ void _kernelPerLink(...)
{
    intokernalDir;
    // 可用 uiSiteIndex 和 dir
}
```

### 多元素入口（`intokernalE` 系列）

与 `preparethreadE` 配套使用。

| 宏 | 计算内容 |
|----|---------|
| `intokernalE(element_count)` | `uiSiteIndex` + `elementIdx`（0..element_count-1） |
| `intokernalEDir(element_count)` | `uiSiteIndex` + `dir` + `elementIdx` |
| `intokernalEDir_NoDir(element_count)` | `uiSiteIndex` + `elementIdx`（不含 dir） |

### 3D Block 入口

与 `preparethread_S*` 配套使用。

| 宏 | 参数 | 计算内容 |
|----|------|---------|
| `intokernalInt4_S(uiT)` | 固定时间 `t` | `sSite4(x,y,z,uiT)` + `uiSiteIndex` + `uiSiteIndex3D` |
| `intokernalInt4_S_Only3D(uiT)` | 固定时间 `t` | 仅 `sSite4` + `uiSiteIndex3D`（无 `uiSiteIndex`） |
| `intokernalInt4_Sxyt(uiZ)` | 固定 z | `sSite4(x,y,uiZ,t)` + `uiSiteIndex` + `uiSiteIndex3DXYT` |
| `intokernalInt4_Sxzt(uiY)` | 固定 y | `sSite4(x,uiY,z,t)` + `uiSiteIndex` + `uiSiteIndex3DXZT` |
| `intokernalInt4_Syzt(uiX)` | 固定 x | `sSite4(uiX,y,z,t)` + `uiSiteIndex` + `uiSiteIndex3DYZT` |

---

## 废弃的宏

以下宏已标记废弃，新代码不应使用：

| 宏 | 状态 | 替代方案 |
|----|------|---------|
| `intokernaldir` | 废弃 | `intokernalDir` |
| `intokernalInt4dirC` | 废弃 | `intokernalDirInt4` |
| `intokernaldirEO` | 废弃 | `intokernalDirEO` |

---

## 完整示例

### 示例 1：标准 site-wise kernel

```cpp
// Host
void MyClass::MyMethod(CFieldGaugeSU3* pGauge)
{
    preparethread;
    _LAUNCH_KERNEL(_kernelSiteWise TMPARG(deviceSU3), block, threads,
                   pGauge->m_pDeviceData);
}

// Device
__global__ void _kernelSiteWise(deviceSU3* pData)
{
    intokernalInt4;
    // sSite4 可用，uiSiteIndex 已检查边界
    pData[uiSiteIndex] = makeSU3Id();
}
```

### 示例 2：Even-odd 半格点 kernel

```cpp
// Host
preparethreadHalf;
_LAUNCH_KERNEL(_kernelEOEven TMPARG(deviceSU3), block, threads,
               pGauge->m_pDeviceData, pEtaTable, TRUE);  // bEven = TRUE
_LAUNCH_KERNEL(_kernelEOOdd TMPARG(deviceSU3), block, threads,
               pGauge->m_pDeviceData, pEtaTable, FALSE); // bEven = FALSE

// Device
__global__ void _kernelEOEven(deviceSU3* pData, const BYTE* pEtaTable, UBOOL bEven)
{
    intokernalEOHalf;
    // uiSiteIndex 指向正确的 even site
    // eta 包含 staggered 相位信息
    pData[uiSiteIndex] = makeSU3Id();
}
```

### 示例 3：Per-link kernel

```cpp
// Host
preparethreadDir;
_LAUNCH_KERNEL(_kernelPerLink TMPARG(deviceSU3), block, threads,
               pGauge->m_pDeviceData);

// Device
__global__ void _kernelPerLink(deviceSU3* pData)
{
    intokernalDirInt4;
    // dir: 0=x, 1=y, 2=z, 3=t
    const UINT uiLinkIndex = uiSiteIndex * _DC_Dir + dir;
    pData[uiLinkIndex] = makeSU3Id();
}
```

### 示例 4：3D slice kernel

```cpp
// Host
preparethread_S;
for (UINT t = 0; t < _HC_Lt; ++t)
{
    _LAUNCH_KERNEL(_kernelSlice TMPARG(deviceSU3), block3d, threads3d,
                   pGauge->m_pDeviceData, t);
}

// Device
__global__ void _kernelSlice(deviceSU3* pData, UINT uiT)
{
    intokernalInt4_S(uiT);
    // sSite4.x, .y, .z 来自线程坐标，.w = uiT
    pData[uiSiteIndex] = makeSU3Id();
}
```

---

## 相关常量

线程分解宏依赖的 host 端常量（定义在 `Data/CCommonData.h`）：

| 常量 | 说明 |
|------|------|
| `_HC_Volume` | 格点总体积 |
| `_HC_VolumeHalf` | 半格点体积（even/odd） |
| `_HC_Dir` | 维度数（通常 4） |
| `_HC_DecompBlock` / `_HC_DecompThread` | 标准分解参数 |
| `_HC_DecompBlockHalf` / `_HC_DecompThreadHalf` | 半格点分解 |
| `_HC_DecompBlockDir` / `_HC_DecompThreadDir` | per-link 分解 |

设备端常量（定义在 `CudaHelper.h`，通过 `__constant__` 内存传递）：

| 常量 | 说明 |
|------|------|
| `_DC_Volume` | 设备端格点体积 |
| `_DC_Dir` | 设备端维度数 |
| `_DC_MultX/Y/Z` | 索引乘法因子 |
| `_DC_Lx/Ly/Lz/Lt` | 各方向尺寸 |
