# Task 1: Compare HISQ Dirac Operator M·v

## Test Description

Compare the HISQ staggered Dirac operator application `M·v` between CLGLib and PyQUDA/QUDA on identical gauge configurations and random fermion vectors.

Two sub-modes are tested:
- **Task1Full**: Full lattice application `pFerm <- (2am + D)·v` via `DWithMass`
- **Task1EO**: Even-even preconditioned application `pFerm <- (M†M)_ee·v_e` via `DDdagger` with `m_bEvenPseudofermion=TRUE`

### Key Parameters
- Lattice: `16^4`
- `2am = 0.22` (CLGLib), PyQUDA mass = `-0.11` (see Mass Sign Convention below)
- `naik_epsilon = -0.11`, Naik = `-(1+epsilon)/24 ≈ -0.037083`
- `t_boundary = -1`, anisotropy = `1.0`
- MILC staggered phase convention

### Mass Sign Convention Fix
QUDA's staggered operator is `M_QUDA(m)·v = 2m·v − H(v)`, while CLGLib's is `M_clg·v = 2am·v + H(v)`. To make them match:
- PyQUDA uses `mass = -0.11`
- Output is negated: `−M_QUDA(−m)·v = M_clg(m)·v`
- `M†M` is sign-invariant, so Task1EO is unaffected.

## How to Invoke

### CLGLib
```bash
cd /home/nbalexis/CLGLib/Bin/Debug
./HISQ ../Debug/HISQ.yaml EHJ_Task1Full <cfg_file> <v_file> <out_file>
./HISQ ../Debug/HISQ.yaml EHJ_Task1EO   <cfg_file> <v_file> <out_file>
```

### PyQUDA
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
/home/nbalexis/pyqudaenv/bin/python3 task1.py
```

### Comparison
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
python3 compare.py
```

## Test Conclusion

**PASS** — All 5 configurations × 3 vectors × 2 modes pass with relative error < 1e-6.

Representative results (cfg16_12, v1):
| Mode | |Δ| | ratio | maxabs |
|------|------|-------|--------|
| Full | 3.280e-12 | 1.635e-15 | 3.553e-14 |
| EO | 1.014e-11 | 2.043e-15 | 1.365e-13 |

Both modes achieve machine-precision agreement (~1e-15 relative).

## Representative Log

### CLGLib — Task1Full (cfg16_12, v1)
```
[04-05-2026 15-42-14|][Task1Full] CfgFile   = /home/nbalexis/comparehisq/data/cfg16_12_double.con
[04-05-2026 15-42-14|][Task1Full] VFile     = /home/nbalexis/comparehisq/data/v1_double.con
[04-05-2026 15-42-14|][Task1Full] OutputFile= /home/nbalexis/comparehisq/data/outputs/clg_task1_full_cfg16_12_v1.con
...
[04-05-2026 15-42-15|][Task1Full] before run: m_f2am=0.220000 m_bEvenPseudofermion=0
[04-05-2026 15-42-15|][Task1Full] applying (2am + D)*v with 2am=0.220000
[04-05-2026 15-42-15|][Task1Full] saved /home/nbalexis/comparehisq/data/outputs/clg_task1_full_cfg16_12_v1.con md5=5569D105C2DEB8C307344A89A5532C9A
[04-05-2026 15-42-15|][Task1Full] done.
```

### CLGLib — Task1EO (cfg16_12, v1)
```
[04-05-2026 15-42-15|][Task1EO] CfgFile   = /home/nbalexis/comparehisq/data/cfg16_12_double.con
[04-05-2026 15-42-15|][Task1EO] VFile     = /home/nbalexis/comparehisq/data/v1_double.con
[04-05-2026 15-42-15|][Task1EO] OutputFile= /home/nbalexis/comparehisq/data/outputs/clg_task1_eo_cfg16_12_v1.con
...
[04-05-2026 15-42-16|][Task1EO] before run: m_f2am=0.220000 m_bEvenPseudofermion=1
[04-05-2026 15-42-16|][Task1EO] applying (M^dM)_ee*v_e with 2am=0.220000
[04-05-2026 15-42-16|][Task1EO] saved /home/nbalexis/comparehisq/data/outputs/clg_task1_eo_cfg16_12_v1.con md5=8A844E95C118D075D8A0FED1F9F726A4
[04-05-2026 15-42-16|][Task1EO] done.
```

### PyQUDA — Task1 (cfg16_12)
```
========== cfg16_12 ==========
[gauge] reading /home/nbalexis/comparehisq/data/cfg16_12.npy
[v1] norm² = 3.927244e+05
[v2] norm² = 3.936435e+05
[v3] norm² = 3.926523e+05
[task1] dirac with mass=-0.11 (so −M_QUDA(−m) = M_clg)
[task1/full v1] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task1_full_cfg16_12_v1.npy  norm²=4.027979e+06
[task1/eo   v1] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task1_eo_cfg16_12_v1.npy  norm²(even)=2.464137e+07
[task1/full v2] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task1_full_cfg16_12_v2.npy  norm²=4.033028e+06
[task1/eo   v2] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task1_eo_cfg16_12_v2.npy  norm²(even)=2.477091e+07
[task1/full v3] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task1_full_cfg16_12_v3.npy  norm²=4.024546e+06
[task1/eo   v3] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task1_eo_cfg16_12_v3.npy  norm²(even)=2.467219e+07
```
