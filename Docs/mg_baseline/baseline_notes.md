# MG 单 GPU Release 基线（P5-2.3，2026-08-05）

构建：`cmake -DCLG_BACKEND=CUDA -DCLG_GPU_ARCH=86`（无 `-DCLG_MULTI_GPU`），
commit 基线：multi-GPU 分支 P5-2.3 提交。

## FileIO 组结果（Release，单进程）

- 10 success / 11 total（1 skipped）/ 2 errors
- 唯一 FAIL：**TestFileIOBridgePPBin**（errors:2）——既有格式缺口
  （Bridge++ 二进制格式未实现，非 MG 相关，P4-1.5 起持续存在）。
- 压缩 IO（TestFileIOCLGCompressed）PASS：`607cebd2` 修复后无崩溃
  （原 P5-2.3 记录的"压缩 IO 崩溃"已消除）。

## FAIL / 超时清单

| 测试 | 状态 | 说明 |
|---|---|---|
| TestFileIOBridgePPBin | FAIL (2 errors) | Bridge++ 格式缺口，既有 |
| 其余 FileIO 测试 | PASS | 含全部 TestMG* 单卡退化路径（无 comm 时 no-op） |

## 备注

- 本基线在 WSL2/Linux + RTX 3090（sm_86）+ CUDA 12.9 + OpenMPI 4.1.6 上录制。
- Windows VS 构建基线未录制（无 Windows 环境，见 P5-3.4）。
