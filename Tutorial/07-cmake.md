# 7. 如何把你的项目加入 CMake 编译

假设你的项目文件放在 `Code/Applications/MyProject/` 目录下：

```
Code/Applications/MyProject/
├── MyPureGaugeSU3.cpp
├── MySU3Wilson.cpp
└── MeasurePolyakov.cpp
```

每个 `.cpp` 文件都有自己的 `main()` 函数，因此需要**分别编译成独立可执行文件**。

## 7.1 在 CMakeLists.txt 中添加

打开 `Code/CMake/CMakeLists.txt`，在文件末尾（`CLGTest` 之后）添加：

```cmake
# ==================== 
# MyProject 
# =================

if (CLG_MyProject)
    include_directories(${PROJECT_SOURCE_DIR}/Applications/MyProject)

    # 纯 SU(3) HMC
    add_executable(MyPureGaugeSU3 
        ${PROJECT_SOURCE_DIR}/Applications/MyProject/MyPureGaugeSU3.cpp)
    target_link_libraries(MyPureGaugeSU3 CLGLib CLGCPULib)

    # SU(3) + Wilson Dirac
    add_executable(MySU3Wilson 
        ${PROJECT_SOURCE_DIR}/Applications/MyProject/MySU3Wilson.cpp)
    target_link_libraries(MySU3Wilson CLGLib CLGCPULib)

    # 测量程序
    add_executable(MeasurePolyakov 
        ${PROJECT_SOURCE_DIR}/Applications/MyProject/MeasurePolyakov.cpp)
    target_link_libraries(MeasurePolyakov CLGLib CLGCPULib)
else()
    message("project MyProject not built, to enable it using -DCLG_MyProject=1")
endif()
```

## 7.2 编译

```bash
cd Code/CMake
cmake CMakeLists.txt -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86 -DCLG_MyProject=1
make -j4
```

编译完成后，可执行文件输出到 `Bin/Ubuntu/`（Linux）或 `Bin/Debug/`（Windows）。

## 7.3 运行

**关键：运行时的工作目录必须是可执行文件所在目录**，否则相对路径无法正确解析。

```bash
cd Bin/Ubuntu

# 运行模拟
./MyPureGaugeSU3

# 运行测量
./MeasurePolyakov
```

运行前确保：
1. YAML 配置文件已放在 `Bin/Debug/` 目录下
2. 如果是测量程序，组态文件（`.con`）也在可执行文件同目录下

---

[< 返回目录](Tutorial.md) | [上一章：测量](06-measurements.md) | [下一章：FAQ](08-faq.md)
