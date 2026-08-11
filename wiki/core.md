# Core 模块

CUDA 基础设施、内存管理、kernel 启动宏、全局管理器、RTTI 系统。

## 文件清单

| 文件 | 路径 | 说明 |
|------|------|------|
| `CLGDefine.h` | `Core/CLGDefine.h` | 基础宏（namespace、安全删除、断言） |
| `CLGFloat.h` | `Core/CLGFloat.h` | 精度抽象层（单/双精度切换） |
| `CLGSetup.h` | `Core/CLGSetup.h` | 编译期配置开关 |
| `CBase.h` | `Core/CBase.h` | RTTI 基类系统 |
| `CCudaBuffer.h` | `Core/CCudaBuffer.h` | GPU 内存池管理器 |
| `CLGLibManager.h` | `Core/CLGLibManager.h` | 全局管理器单例 |
| `CudaHelper.h` | `Core/CudaHelper.h` | CUDA 设备管理、常量内存、kernel 启动 |
| `CudaHelperFunctions.h` | `Core/CudaHelperFunctions.h` | CUDA 运行时包装、错误字符串 |
| `CLGLib_Private.h` | `CLGLib_Private.h` | 预编译头：线程分解宏、kernel 入口宏、字段头文件集合 |

## 核心概念

### 精度抽象

`CLGFloat.h` 定义了 `Real` 和 `CLGComplex`，在编译时映射为 `float`/`cuComplex` 或 `double`/`cuDoubleComplex`：

```cpp
#if _CLG_DOUBLEFLOAT
    typedef double Real;
    typedef cuDoubleComplex CLGComplex;
#else
    typedef float Real;
    typedef cuComplex CLGComplex;
#endif
```

所有数学函数通过宏或模板统一处理两种精度。物理常量（`PI`、`PI2`、`PI_s2` 等）同样按精度定义。

### RTTI 与对象工厂

`CBase.h` 提供轻量级自定义 RTTI，无需依赖 C++ 的 `typeid`/`dynamic_cast`：

```cpp
// 头文件
__CLG_REGISTER_HELPER_HEADER(CMyClass)

class CLGAPI CMyClass : public CBase
{
    __CLGDECLARE_CLASS(CMyClass)
    // ...
};

// 实现文件
__CLGIMPLEMENT_CLASS(CMyClass)
```

- `CClass` — 元类描述符（类名、大小、工厂函数）
- `CBase` — 所有对象的抽象基类，提供 `GetClass()` 和 `GetInfos()`
- `CClassGather` — 全局类注册表（单链表 `m_pClasses`），支持按名称查找和创建对象
- `__CLG_REGISTER_HELPER_HEADER` — 在头文件中生成辅助结构体，强制 Linux 静态链接时保留注册构造函数

**用途**：参数文件反序列化时，通过类名字符串创建对应实例：
```cpp
CBase* pObj = appCreate(_T("CActionGaugePlaquette"));
// 实际由 CLGLibManager 内部调用
```

**注意**：静态库在 Linux 下链接时，若某个 `.cpp` 中的对象未被直接引用，其全局构造函数可能被优化掉，导致类未注册。`__CLG_REGISTER_HELPER_HEADER` 生成的辅助对象被显式引用，从而保证注册代码被链接。

### GPU 内存池

`CCudaBuffer` 预分配一块大设备内存，小分配（<4KB）从中子分配，减少 `cudaMalloc` 碎片。大分配直接 fallback 到 `cudaMalloc`。

```cpp
__cudaMalloc(ptr, size);   // 带追踪的分配（记录文件名/行号）
__cudaFree(ptr);           // 释放
appSafeCudaFree(ptr);      // 安全释放 + 置 NULL
```

### 全局管理器

`CLGLibManager` 是单例，整个模拟的生命周期由它管理：

1. 读取 YAML 参数 → 创建 `CIndex` → 烘焙索引表
2. 按参数创建场（gauge/fermion/boson）
3. 按参数创建作用量、更新器、测量器、求解器
4. 运行时通过全局函数访问：
   - `appGetLattice()` — 格点数据容器
   - `appGetCudaHelper()` — CUDA 辅助
   - `appGetRandom()` — 随机数生成器

## CUDA Kernel 启动基础设施

### `_LAUNCH_KERNEL` 宏（强制使用）

**禁止裸 `<<< >>>`**。所有 kernel 启动必须通过：

```cpp
// 非模板 kernel
_LAUNCH_KERNEL(_kernelFunc, block, threads, arg1, arg2, ...);

// 模板 kernel — 使用 TMPARG 宏
_LAUNCH_KERNEL(_kernelFunc TMPARG(deviceSU3), block, threads, arg1, arg2);
_LAUNCH_KERNEL(_kernelFunc TMPARG(deviceVector, deviceGauge), block, threads, ...);
```

`TMPARG` 定义在 `CLGDefine.h`：`#define TMPARG(...) <__VA_ARGS__>`

底层通过 `launchKernel()` 函数实现，支持跨平台（CUDA/DCU）。`_LAUNCH_KERNEL` 使用栈分配参数数组，调用后自动 `FreeStack()`。

### 线程分解宏

在 `CLGLib_Private.h` 中定义：

| 宏 | 用途 | 变量 |
|----|------|------|
| `preparethread` | 标准格点分解 | `block = _HC_DecompBlock`, `threads = _HC_DecompThread` |
| `preparethreadHalf` | even/odd 半格点 | `block = _HC_DecompBlockHalf`, `threads = _HC_DecompThreadHalf` |
| `preparethreadDir` | per-link 分解 | `block = _HC_DecompBlockDir`, `threads = _HC_DecompThreadDir` |
| `preparethreadE(n)` | 确保每个 block 覆盖所有 sites，每 site n 个元素 | 自定义 block/threads |
| `preparethreadEDir(n)` | 同上，per-link | 自定义 block/threads |

Host 端常量（`_HC_` 前缀）在 `CLGLibManager::InitialWithParameter` 中根据格点大小自动计算。Device 端常量（`_DC_` 前缀）通过 `__constant__` 内存传递。

### Kernel 入口宏

| 宏 | 计算内容 | 边界检查 |
|----|---------|---------|
| `intokernal` | `uiSiteIndex` | `uiSiteIndex >= _DC_Volume` 则 return |
| `intokernalInt4` | `uiSiteIndex` + `sSite4`（4D 坐标） | 同上 |
| `intokernalInt4dirC` | `uiSiteIndex` + `sSite4` + `uiDir = _DC_Dir` | 同上（废弃中） |
| `intokernalEOHalf` | even/odd half-lattice + `eta` 表查值 | 自动处理奇偶偏移 |
| `intokernalDir` | `uiLinkIndex` → `uiSiteIndex` + `dir` | `uiSiteIndex >= _DC_Volume` |
| `intokernalDirInt4` | 同上 + `sSite4` | 同上 |

**3D 线程分解**（用于特定平面/体积切片）：
- `preparethread_S` / `intokernalInt4_S(uiT)` — x,y,z 3D block，固定 t
- `preparethread_Sxyt` / `intokernalInt4_Sxyt(uiZ)` — x,y,t 3D block，固定 z
- 类似的有 `_Sxzt`, `_Syzt` 变体

### 设备端常量内存

`CudaHelper.h` 声明了以下 `__device__ __constant__` 数组：

| 常量 | 类型 | 说明 |
|------|------|------|
| `_constIntegers[]` | `UINT` | 格点尺寸、分解参数、场数量等（`EConstIntId` 索引） |
| `_constSignedIntegers[]` | `INT` | 格点中心坐标 |
| `_constFloats[]` | `Real` | 物理常数（如规范动量因子） |
| `__fieldPointers[]` | `CField*` | 设备端场指针数组（按 field ID） |
| `__boundaryFieldPointers[]` | `CFieldBoundaryParent*` | 边界场指针 |
| `__r` | `CRandom*` | 设备端随机数生成器 |
| `__idx` | `CIndexData*` | 设备端索引数据缓存 |
| `__chiralGamma[]` | `gammaMatrix` | Gamma 矩阵（chiral 基） |
| `__SU3Generators[9]` | `deviceSU3` | SU(3) 生成元 |
| `_plaq_idx[6][2]` | `SCHAR` | 6 个 plaquette 方向对 |

Host 端通过 `CCudaHelper::CopyConstants()` 一次性拷贝到设备。

### 并行归约

`CCudaHelper` 提供：
- `ReduceReal(DOUBLE* buffer, UINT count)` — 实数数组求和
- `ReduceComplex(cuDoubleComplex* buffer, UINT count)` — 复数数组求和
- `ThreadBufferSum(...)` / `ThreadBufferZero(...)` — 线程局部 buffer 操作

## 关键设计模式

- **预编译头聚合**：`CLGLib_Private.h` 是项目的"超级头文件"，按模块顺序包含所有基础设施头文件。`.cu` 文件通常只包含它。
- **常量内存查表**：格点尺寸、分解参数、索引指针等全部通过 `__constant__` 内存传递，kernel 中零开销访问。
- **参数栈分配**：`_LAUNCH_KERNEL` 使用固定大小的 host 栈（4KB）打包 kernel 参数，避免动态分配。
- **工厂 + RTTI**：参数驱动的对象创建（`CCLGLibManager::CreateAction` 等）依赖 `CClassGather` 的全局注册表。

## 常用宏速查

| 宏 | 位置 | 说明 |
|----|------|------|
| `__BEGIN_NAMESPACE` / `__END_NAMESPACE` | `CLGDefine.h` | `namespace CLGLib { }` |
| `__CLGDECLARE_CLASS(name)` | `CBase.h` | 声明 RTTI |
| `__CLGIMPLEMENT_CLASS(name)` | `CBase.h` | 实现 RTTI |
| `_LAUNCH_KERNEL(func, block, threads, args)` | `CudaHelper.h` | 跨平台 kernel 启动 |
| `TMPARG(...)` | `CLGDefine.h` | 模板参数宏 `#define TMPARG(...) <__VA_ARGS__>` |
| `preparethread` | `CLGLib_Private.h` | 标准分解 |
| `intokernal` | `CLGLib_Private.h` | 计算索引 + 边界检查 |
| `_HC_Volume` | `Data/CCommonData.h` | host 端格点体积 |
| `_DC_Volume` | `Data/CCommonData.h` | device 端格点体积（常量内存） |
| `F(x)` | `CLGDefine.h` | 字面量宏，单精度 `x##f`，双精度 `x` |
| `appSafeFree` / `appSafeDelete` | `CLGDefine.h` | 安全释放 + 置 NULL |
