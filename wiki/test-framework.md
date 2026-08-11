# 测试框架

CLGTest：单元测试注册、运行和报告系统。

## 文件清单

| 文件 | 路径 | 说明 |
|------|------|------|
| `CLGTest.h` | `CLGTest/CLGTest.h` | 测试注册宏、测试基类、断言宏 |
| `CLGTest.cpp` | `CLGTest/CLGTest.cpp` | 测试运行器、主函数 |

## 测试注册宏

**文件**: `CLGTest/CLGTest.h`

```cpp
__REGIST_TEST(functionname, category, paramName, showName)
```

将测试函数注册到全局测试表中。参数：
- `functionname` — 测试函数名（`UINT func(CParameters&)` 签名）
- `category` — 测试类别（如 `GaugeFixing`, `Action`, `Fermion`）
- `paramName` — YAML 参数文件中的测试配置名（对应 YAML 中的顶级 key）
- `showName` — 显示名称

**断言宏**：
- `test_EQUAL(a, b)` — 相等断言
- `test_NEAR(a, b, eps)` — 近似相等断言
- `test_TRUE(cond)` / `test_FALSE(cond)` — 布尔断言

## 测试运行器

**文件**: `CLGTest/CLGTest.cpp`

### 编译

```bash
cd Code/CMake
cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86
make -j$(nproc)
```

> **注意**：完整 clean build 耗时较长（大型 CUDA 项目），如果终端有超时限制，建议把编译放到后台：
> ```bash
> nohup make -j$(nproc) > build.log 2>&1 &
> tail -f build.log
> ```
> 或者减少并行度（如 `make -j4`）以降低峰值耗时。

编译产物在 `Bin/Ubuntu/CLGTest`（Release）或 `Bin/UbuntuDebug/CLGTest`（Debug）。

### 运行测试

从 `Bin/Ubuntu/`（或 `Bin/UbuntuDebug/`）目录运行：

```bash
./CLGTest <testname>           # 按名称运行单个测试（匹配 ParamName，大小写不敏感）
./CLGTest <category>           # 运行某类别所有测试（类别名小写，如 updator）
./CLGTest all                  # 运行所有测试
```

`<testname>` 对应 `__REGIST_TEST` 的第三个参数（`paramName`），也对应 YAML 文件的顶级 key。

不带参数运行时进入交互模式，支持以下命令：
- `<number>` — 按序号运行测试
- `<testname>` — 按名称运行测试
- `<category>` — 按类别名运行（小写）
- `r` — 运行所有测试
- `l` — 列出所有测试
- `p` — 重新加载参数文件
- `q` — 退出

### 执行流程

1. 初始化 CLG（格点、场、作用量、测量器等）
2. 查找并调用注册的测试函数
3. 收集断言结果
4. 报告通过/失败
## 参数文件

测试所需的格点参数、场配置、作用量参数等通过 YAML 参数文件传入。参数文件放置在 `Bin/Debug/`（或 `Bin/Release/`）目录下，扩展名为 `.yaml`。

**典型参数文件结构**：
```yaml
Lattice:
    Name: CIndexSquare
    LatticeLength: [8, 8, 8, 8]

Gauge:
    FieldName: CFieldGaugeSU3
    ...

TestMyFeature:
    TestConfig:
        ProjectionType: MAG
```

## 添加新测试

### Step 1: 编写测试函数

测试函数签名必须是 `UINT func(CParameters& params)`，返回错误计数（0 = 通过）。

在 `Code/CLGTest/Tests/` 下创建 `.cpp` 文件（CPU 代码，不允许依赖 CUDA），参考：

```cpp
#include "CLGTest.h"

UINT TestMyFeature(CParameters& params)
{
    UINT uiErrors = 0;

    // 通过 framework 获取场对象
    CFieldGauge* pGauge = dynamic_cast<CFieldGauge*>(
        appGetLattice()->m_pGaugeField[0]);

    // 通过 framework 获取临时场（pooled field，用后需 Return）
    CFieldGauge* pCopy = dynamic_cast<CFieldGauge*>(
        appGetLattice()->GetPooledFieldById(1, _T(__FILE__), __LINE__));

    // 执行操作（framework 方法，非 CUDA API）
    pGauge->CopyTo(pCopy);
    pCopy->AxpyMinus(pGauge);

    // 验证
    DOUBLE diff = pCopy->GetLength();
    if (diff > F(1e-8))
    {
        ++uiErrors;
        appCrucial(_T("FAILED: diff = %e\n"), diff);
    }

    pCopy->Return();
    return uiErrors;
}
```

### Step 2: 注册测试

在同一个 `.cpp` 文件底部调用注册宏：

```cpp
// __REGIST_TEST(函数名, 类别, YAML参数名, 显示名)
__REGIST_TEST(TestMyFeature, GaugeFixing, TestMyFeature, "My Feature Test")

// ___REGIST_TEST 变体支持额外 tag：
//   _TEST_NOCHECK — 默认跳过（需手动触发）
//   _TEST_BOUND   — 仅当 _CLG_USE_LAUNCH_BOUND 时运行
//   _TEST_DOUBLE  — 仅当 _CLG_DOUBLEFLOAT 时运行
___REGIST_TEST(TestSlowFeature, Tools, TestSlowFeature, "Slow Test", _TEST_NOCHECK)
```

### Step 3: 创建 YAML 参数文件

在 `Bin/Debug/` 目录下创建 YAML 文件，顶级 key 必须与 `__REGIST_TEST` 的 `paramName` 一致：

```yaml
TestMyFeature:
    Lattice:
        Name: CIndexSquare
        LatticeLength: [4, 4, 4, 4]
    GaugeFieldInit:
        - FieldType: EFT_GaugeSU3
          FieldName: CFieldGaugeSU3
    Beta: 6.0
```

### Step 4: 注册 YAML 文件

在 `CLGTest.cpp` 的 `LoadParams()` 函数中添加相应文件：

```cpp
CYAMLParser::ParseFile(_T("../Debug/TestSuit_MyFeature.yaml"), params);
```

### Step 5: 更新 CMakeLists.txt

在 `Code/CMake/CMakeLists.txt` 的 `CLGTest` 源文件列表中添加新文件：

```cmake
add_executable(CLGTest 
    ...
    ${PROJECT_SOURCE_DIR}/CLGTest/Tests/TestMyFeature.cpp
)
```

### 关键约束

- **CLGTest 不允许依赖 CUDA**：测试函数只能通过 CLGLib framework（`appGetLattice()`、`CField::CopyTo()`、`CAction::Energy()` 等）操作数据，禁止直接调用 `cudaMalloc`、`cudaMemcpy`、`<<<>>>` 等 CUDA API。
- 测试函数返回错误计数（非零 = 失败），不使用断言宏（gtest 风格）。
- 通过 `appGetLattice()->GetPooledFieldById(...)` 获取的场必须调用 `->Return()` 归还。
- 参数通过 YAML 文件传入，框架自动解析和初始化。

## 常用测试模式

- **能量守恒测试**：HMC 轨迹前后能量差应在 Metropolis 接受范围内
- **规范不变性测试**：gauge fixing 或 gauge transformation 后可观测量不变
- **厄米性测试**：D 算子满足 `⟨x|D|y⟩ = ⟨y|D|x⟩*`（或 `γ5`-厄米）
- **反厄米性测试**：staggered D 算子应反厄米
- **解析解比对**：小格点上与已知解析结果比对
