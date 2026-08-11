# Platform Tools 模块

平台无关工具：动态数组、哈希表、字符串、YAML 解析、性能分析、计时器、位标志。

## 文件清单

### 数据结构

| 文件 | 路径 | 说明 |
|------|------|------|
| `TArray.h` | `Tools/Data/TArray.h` | 动态数组（类似 MFC CArray） |
| `THashMap.h` | `Tools/Data/THashMap.h` | 哈希表 |
| `CCString.h` | `Tools/Data/CCString.h` | 字符串类 |
| `CBitFlag.h` | `Tools/Data/CBitFlag.h` | 64 位位标志 |
| `CUUID.h` | `Tools/Data/CUUID.h` | UUID 生成 |
| `CLinkedList.h` | `Tools/Data/CLinkedList.h` | 链表 |
| `MemStack.h` | `Tools/Data/MemStack.h` | 内存栈分配器 |
| `TemplateFunctions.h` | `Tools/Data/TemplateFunctions.h` | 通用模板辅助函数 |
| `STDStringFunctions.h` | `Tools/Data/STDStringFunctions.h` | 字符串转换函数 |

### 配置与解析

| 文件 | 路径 | 说明 |
|------|------|------|
| `CYAMLParser.h` | `Tools/CYAMLParser.h` | YAML 参数文件解析 |
| `EnumGather.h` | `Tools/EnumGather.h` | 枚举收集辅助 |

### 性能与调试

| 文件 | 路径 | 说明 |
|------|------|------|
| `CSimpleProfiler.h` | `Tools/CSimpleProfiler.h` | 简单函数级性能分析 |
| `Timer.h` | `Tools/Timer.h` | 高精度计时器 |
| `Tracer.h` | `Tools/Tracer.h` | 代码跟踪/日志 |

## TArray — 动态数组

**文件**: `Tools/Data/TArray.h`

模板动态数组，接口类似 MFC `CArray` + 部分 STL 风格：

| 方法 | 说明 |
|------|------|
| `SetSize(n, growBy)` | 设置大小，自动分配/释放 |
| `AddItem(element)` | 追加元素 |
| `AddSize(n)` / `AddZeroed(n)` | 扩展并清零 |
| `RemoveAt(idx, count)` / `RemoveItem(item)` | 删除 |
| `InsertAt(idx, element, count)` | 插入 |
| `GetAt(idx)` / `operator[]` | 访问 |
| `GetData()` | 原始指针访问 |
| `Num()` / `GetSize()` / `IsEmpty()` | 大小查询 |
| `ReserveSpace(n)` | 预留容量 |
| `FreeExtra()` / `Shrink()` | 收缩到实际大小 |
| `RemoveAll(bFreeMemory)` / `Reset()` | 清空 |
| `FindItemIndex(item)` | 查找 |
| `Copy(other)` / `operator=` | 复制 |
| `Append(src)` | 追加数组 |
| `PushBack(element)` / `Pop()` | 栈操作 |

**内存管理**：使用 `new BYTE[]` 分配原始内存 + placement new 构造，支持带析构的清理。

## THashMap — 哈希表

**文件**: `Tools/Data/THashMap.h`

基于字符串键的哈希表，广泛用于：
- 参数查找（`CParameters`）
- 类名到工厂函数的映射（RTTI）
- 字段 ID 到字段对象的映射

| 方法 | 说明 |
|------|------|
| `operator[key]` | 插入/访问 |
| `Exist(key)` | 存在检查 |
| `Lookup(key, value)` | 安全查找 |
| `RemoveAll()` | 清空 |

## CCString — 字符串

**文件**: `Tools/Data/CCString.h`

项目内字符串类，封装 `std::basic_string<TCHAR>`。

**常用全局函数**：
- `appToString(val)` — 任意类型转字符串
- `appStrToINT/Real/FLOAT/DOUBLE(s)` — 字符串转数值
- `appGetPath(buf, len)` — 获取可执行路径

## CBitFlag — 位标志

**文件**: `Tools/Data/CBitFlag.h`

64 位（`QWORD`）位标志封装：

| 方法 | 说明 |
|------|------|
| `SetFlag(flag)` / `ClearFlag(flag)` | 设置/清除 |
| `HasFlag(flag)` | 检查 |
| `ToggleFlagBy(flag, bSet)` | 条件切换 |

## CParameters / CYAMLParser — YAML 参数系统

**文件**: `Tools/CYAMLParser.h`

`CParameters`：分层参数存储，支持标量、数组、嵌套参数表。

```cpp
CParameters param;
CYAMLParser::ParseFile("config.txt", param);

// 获取标量
INT steps;
param.FetchValueINT("StepCount", steps);

// 获取数组
TArray<Real> betas;
param.FetchValueArrayReal("BetaList", betas);

// 获取子表
CParameters actionParam = param.GetParameter("Action");
```

**获取方法**：
- `FetchValueINT/Real/FLOAT/DOUBLE(key, value)` — 标量
- `FetchStringValue(key, value)` — 字符串
- `FetchStringVectorValue(key, value)` — 字符串数组
- `FetchParameterValue(key, value)` — 嵌套参数表
- `Exist(key)` — 存在检查
- `GetParameter(key)` — 获取子表（不存在则报错）

**宏辅助**：`_FetchFunction(typen)` 和 `_FetchFunctionArray(typen)` 自动生成类型化的获取方法。

## CSimpleProfiler — 性能分析

**文件**: `Tools/CSimpleProfiler.h`

函数级性能分析器，记录调用次数和总耗时（纳秒）。

**使用方式**：
```cpp
void MyFunction()
{
    _RECORD(MyFunction);  // 自动记录从构造到析构的时间
    // ... 函数体
}
```

**全局函数**：
- `appDumpProfiler()` — 打印所有记录
- `appClearProfiler()` — 清空记录

**开关**：由 `_CLG_PROFILER` 宏控制，关闭时 `_RECORD` 展开为空。**默认关闭**（`CLGSetup.h` 中 `_CLG_PROFILER 0`），可用 CMake 选项 `-DCLG_PROFILER=1` 开启。

**注意**：开启后 `CSimpleProfileRecord` 构造时会 `cudaDeviceSynchronize()` 以准确测量 kernel 耗时，即每个 `_RECORD` 点一次全设备同步，开销很大，仅用于性能分析。

**相关开关 `_CLG_CHECKSYNCHRONIZE`**（`Core/CudaHelperFunctions.h`）：控制 `_CHECKCUDA`（`launchKernel` 及若干检查点中的 `appSynchronize()`），默认关闭（`CLGSetup.h` 中 `_CLG_CHECKSYNCHRONIZE 0`），CMake 选项 `-DCLG_CHECKSYNCHRONIZE=1` 开启。开启后每次 kernel launch 都同步设备，用于尽早暴露 kernel 错误，会串行化所有 launch，仅供调试。

## 关键设计模式

- **TArray 为默认容器**：几乎所有列表/数组都使用 `TArray` 而非 `std::vector`。
- **THashMap 为默认字典**：类工厂、参数表、字段映射都使用 `THashMap`。
- **YAML 为唯一配置格式**：所有模拟参数通过 `CYAMLParser` 从 `.txt` 文件读取。
- **Profiler 为编译期开关**：通过 `_CLG_PROFILER` 宏完全消除性能分析开销。
