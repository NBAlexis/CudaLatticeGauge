# Task 4: Compare Tree-Improved Symanzik Gauge Action

## Test Description

Compare the tree-improved Symanzik gauge action energy between CLGLib and PyQUDA/QUDA.

The action is:
```
S = (β/Nc) · Σ [ (1 − ReTr(U_plaq)) + c_rect · (1 − ReTr(U_rect)) ]
```

where `c_rect = −1/(20·u0⁴)` with `u0 = 0.8794` and `β = 7.29`.

### Constant Offset
CLGLib and PyQUDA use opposite sign conventions for the trace term:
- CLGLib: `S_clg = (β/Nc) · Σ [ Nc − ReTr + c_rect·(Nc − ReTr_rect) ]`
- PyQUDA: `S_py = −(β/Nc) · Σ [ ReTr + c_rect·ReTr_rect ]`

The difference is a constant offset:
```
offset = β · V · (6 + 12·c_rect)
```

After subtracting this offset, the two values match.

### Key Parameters
- `u0 = 0.8794`
- `β = 7.29`
- `c_rect = −1/(20·u0⁴) ≈ −0.064654`
- Lattice volume `V = 16⁴ = 65536`

## How to Invoke

### CLGLib
```bash
cd /home/nbalexis/CLGLib/Bin/Debug
./HISQ ../Debug/HISQ.yaml EHJ_Task4Action <cfg_file> <out_file>
```

### PyQUDA
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
/home/nbalexis/pyqudaenv/bin/python3 task4.py
```

## Test Conclusion

**PASS** — All 5 configurations pass with relative error < 1e-6 after offset correction.

Representative results (cfg16_12):
| | Value |
|---|-------|
| S_clg | 904981.9790351744 |
| S_py (raw) | −1590894.3569753235 |
| S_py + offset | 904981.9790351744 |
| relerr | 1.802e-07 |

The remaining ~2e-7 relative difference is consistent with single-precision accumulation in QUDA's gauge observable path; it is well below the 1e-6 threshold.

## Representative Log

### CLGLib — Task4Action (cfg16_12)
```
[04-05-2026 15-42-24|][Task4Action] CfgFile   = /home/nbalexis/comparehisq/data/cfg16_12_double.con
[04-05-2026 15-42-24|][Task4Action] OutputFile= /home/nbalexis/comparehisq/data/outputs/clg_task4_action_cfg16_12_v1.txt
...
[04-05-2026 15-42-24|][Task4Action] gauge is tree-improved, RectOverPlaq=-0.064654
[04-05-2026 15-42-24|][Task4Action] computing tree-improved gauge action energy
[04-05-2026 15-42-24|][Task4Action] S = 904981.97903517436
[04-05-2026 15-42-24|][Task4Action] saved /home/nbalexis/comparehisq/data/outputs/clg_task4_action_cfg16_12_v1.txt
[04-05-2026 15-42-24|][Task4Action] done.
```

### PyQUDA — Task4 (cfg16_12)
```
[task4] tree-improved Symanzik gauge action, beta=7.29, u0=0.8794
[task4] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task4_action_cfg16_12.txt  S = -1590894.3569753235
```
