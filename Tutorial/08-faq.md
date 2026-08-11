# 8. 常见问题速查

### Q1: 程序启动时报 "Unable to create the gauge field!"

**原因**：`FieldName` 对应的类未在编译时启用。CLGLib 通过 `CLGSetup.h` 中的宏控制哪些场类被编译。检查 `_CLG_SU3_GAUGE`、`_CLG_Z2_GAUGE` 等宏是否为 1。

### Q2: 费米子模拟时求解器不收敛

**排查步骤**：
1. 检查 `Accuracy` 是否设得太小（如 `1e-10`），可以先放宽到 `1e-4` 测试
2. 检查 `PoolNumber` 是否足够（Wilson + GMRES 至少需要 `6 + MaxDim`）
3. 检查 `Period` 时间方向是否为 `-1`（费米子需要反周期边界）
4. 检查 `Hopping` / `Mass` 参数是否合理

### Q3: 如何只读取一个组态文件做测量，不配置 Updator？

YAML 中**不配置 `Updator` 节**即可。`appInitialCLG` 会跳过更新器初始化。但务必把规范场的 `FieldInitialType` 设为 `EFIT_ReadFromFile`，并提供 `GaugeFileName`。

### Q4: 保存的组态文件格式

| 格式 | 说明 | 适用场景 |
|------|------|----------|
| `EFFT_CLGBin`（默认） | 二进制格式，含 MD5 校验 | 标准保存 |
| `EFFT_CLGBinCompressed` | 压缩二进制 | 大格点节省空间 |
| `EFFT_CLGBinDouble` | 双精度二进制 | 需要高精度 |
| `EFFT_CLGBinFloat` | 单精度二进制 | 节省空间 |
| `EFFT_BridgePPTXT` | Bridge++ 文本格式 | 与其他代码交换 |
| `EFFT_BridgePPBin` | Bridge++ 二进制格式 | 与其他代码交换 |

在 YAML 中指定：
```yaml
Gauge:
    FieldName : CFieldGaugeSU3
    FieldInitialType : EFIT_ReadFromFile
    GaugeFileType : EFFT_CLGBin
    GaugeFileName : ../Debug/my_conf.con
```

### Q5: `Update` 和 `UpdateUntileAccept` 的区别

- `Update(n, bMeasure)` —— 执行 `n` 次 HMC trajectory，**不管是否被 Metropolis 接受**
- `UpdateUntileAccept(n, bMeasure)` —— 执行直到产生 `n` 次**被接受**的组态。对于纯规范场 HMC 接受率通常 60% ~ 90%，所以实际执行的步数约为 `n / accept_rate`

### Q6: 离散规范群（Z_N, D_N）能否做 HMC？

**不能**。HMC 要求作用量对场变量的导数（力）有良好定义，离散群不存在连续导数。离散群必须使用 `CHeatbath` 更新器。

### Q7: 为什么程序运行后控制台没有输出？

`appSetupLog(params)` 会将所有输出重定向到以时间戳命名的 `.log` 文件（如 `21-05-2026 02-13-48.log`），控制台不再显示。如需同时输出到控制台，可在 `appSetupLog` 前加 `printf` 调试，或检查生成的日志文件。

### Q8: 程序报 "unable to read parameter file"

检查以下几点：
1. **运行时工作目录**：`../Debug/xxx.yaml` 是相对于**工作目录**的路径，不是相对于可执行文件的路径。确保运行时 `pwd` 是可执行文件所在目录（如 `Bin/Ubuntu/`）
2. **YAML 文件位置**：确保 `Bin/Debug/xxx.yaml` 确实存在
3. **大小写敏感**：Linux 文件名区分大小写

### Q9: 我把公共参数放在 YAML 顶层，子节点里为什么读取不到？

`CParameters::GetParameter("NodeName")` 返回子树的参数引用，但它**不会继承**父节点的参数。例如：

```yaml
Dim : 4
LatticeLength : [8, 8, 8, 8]

JobSU3:
    ActionListLength : 1
```

如果你用 `params.GetParameter(_T("JobSU3"))` 获取子树，子树中只有 `ActionListLength`，**没有 `Dim` 和 `LatticeLength`**。每个 YAML 块必须是自包含的。

### Q10: 组态文件编号从几开始？

框架自动保存的组态文件编号从 **1** 开始，扩展名为 `.con`。例如 `su3_conf_1.con`、`su3_conf_2.con`。注意不是从 0 开始，扩展名也不是 `.con_`。

---

[< 返回目录](Tutorial.md) | [上一章：CMake](07-cmake.md)
