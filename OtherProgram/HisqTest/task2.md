# Task 2: Compare Rational Approximation (M†M)^{1/4}·v

## Test Description

Compare the rational approximation application of `(M†M)^{+1/4}` on the HISQ staggered operator between CLGLib and PyQUDA/QUDA.

Two sub-modes are tested:
- **Task2Full**: Full lattice `pFerm <- (M†M)^{1/4}·v` via `D_MC` with `m_bEvenPseudofermion=FALSE`
- **Task2EO**: Even-even preconditioned `pFerm <- (M†M)_ee^{1/4}·v_e` via `D_MC` with `m_bEvenPseudofermion=TRUE`

### Rational Approximation Setup
The Remez approximation is generated for `f(x) = (x + 4m²)^{+1/4}` with `m = 0.11` (so `4m² = 0.0484`), on the interval `[1e-3, 50]` with degree 12. The coefficients are stored in `rational_pos14.json`.

Both sides set **mass = 0** in the Dirac operator; the mass shift is entirely carried by the rational approximation offsets. This ensures:
- Full: `f(D†D)·v = (D†D + 4m²)^{1/4}·v = (M†M)^{1/4}·v`
- EO: `f(D_eo·D_oe)·v_e = (D_eo·D_oe + 4m²)^{1/4}·v_e = (M†M)_ee^{1/4}·v_e`

### Key Parameters
- Same lattice, gauge, and `v` as Task 1
- `Mass = 0` in both CLGLib YAML and PyQUDA `HISQDirac`
- Naik, epsilon, boundary conditions identical to Task 1

## How to Invoke

### CLGLib
```bash
cd /home/nbalexis/CLGLib/Bin/Debug
./HISQ ../Debug/HISQ.yaml EHJ_Task2Full <cfg_file> <v_file> <out_file>
./HISQ ../Debug/HISQ.yaml EHJ_Task2EO   <cfg_file> <v_file> <out_file>
```

### PyQUDA
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
/home/nbalexis/pyqudaenv/bin/python3 task2.py
```

## Test Conclusion

**PASS** — All 5 configurations × 3 vectors × 2 modes pass with relative error < 1e-6.

Representative results (cfg16_12, v1):
| Mode | |Δ| | ratio | maxabs |
|------|------|-------|--------|
| Full | 3.193e-10 | 2.890e-13 | 2.634e-12 |
| EO | 1.559e-12 | 1.996e-15 | 2.522e-14 |

The EO mode achieves machine-precision agreement (~2e-15). The Full mode is slightly looser (~3e-13) because PyQUDA's full-lattice path uses a custom Python multi-shift CG instead of QUDA's native multi-shift solver, but still well within tolerance.

## Representative Log

### CLGLib — Task2Full (cfg16_12, v1)
```
[04-05-2026 15-42-16|][Task2Full] CfgFile   = /home/nbalexis/comparehisq/data/cfg16_12_double.con
[04-05-2026 15-42-16|][Task2Full] VFile     = /home/nbalexis/comparehisq/data/v1_double.con
[04-05-2026 15-42-16|][Task2Full] OutputFile= /home/nbalexis/comparehisq/data/outputs/clg_task2_full_cfg16_12_v1.con
...
[04-05-2026 15-42-17|][Task2Full] before run: m_f2am=0.000000 m_bEvenPseudofermion=0
[04-05-2026 15-42-17|][Task2Full] applying (M^dM)^{1/4}*v via D_MC with 2am=0.000000
[04-05-2026 15-42-19|][Task2Full] saved /home/nbalexis/comparehisq/data/outputs/clg_task2_full_cfg16_12_v1.con md5=00A826973408BB48D559547F1AB0B0DC
[04-05-2026 15-42-19|][Task2Full] done.
```

### CLGLib — Task2EO (cfg16_12, v1)
```
[04-05-2026 15-42-19|][Task2EO] CfgFile   = /home/nbalexis/comparehisq/data/cfg16_12_double.con
[04-05-2026 15-42-19|][Task2EO] VFile     = /home/nbalexis/comparehisq/data/v1_double.con
[04-05-2026 15-42-19|][Task2EO] OutputFile= /home/nbalexis/comparehisq/data/outputs/clg_task2_eo_cfg16_12_v1.con
...
[04-05-2026 15-42-20|][Task2EO] before run: m_f2am=0.000000 m_bEvenPseudofermion=1
[04-05-2026 15-42-20|][Task2EO] applying (M^dM)_ee^{1/4}*v_e via D_MC with 2am=0.000000
[04-05-2026 15-42-21|][Task2EO] saved /home/nbalexis/comparehisq/data/outputs/clg_task2_eo_cfg16_12_v1.con md5=6C331F73FA0DF4C4777D1A74D88BA3F6
[04-05-2026 15-42-21|][Task2EO] done.
```

### PyQUDA — Task2Full multi-shift CG (cfg16_12, v1)
```
[task2/3] dirac with mass=0 (rational carries the +4m² shift)
[task2/full v1] (M†M)^(1/4)·v via custom multi-shift CG (12 shifts)
    shift[ 0] sigma=5.445258e-02 residue=-8.827397e-04  CG iters=  186 rel_resid=8.84e-13  (3.29s)
    shift[ 1] sigma=7.965177e-02 residue=-3.676156e-03  CG iters=  167 rel_resid=8.73e-13  (2.76s)
    shift[ 2] sigma=1.401896e-01 residue=-1.093447e-02  CG iters=  139 rel_resid=8.67e-13  (2.43s)
    shift[ 3] sigma=2.745223e-01 residue=-3.008906e-02  CG iters=  106 rel_resid=8.94e-13  (1.82s)
    shift[ 4] sigma=5.685767e-01 residue=-8.105869e-02  CG iters=   78 rel_resid=8.94e-13  (1.40s)
    shift[ 5] sigma=1.213337e+00 residue=-2.180630e-01  CG iters=   56 rel_resid=6.22e-13  (1.10s)
    shift[ 6] sigma=2.641536e+00 residue=-5.943394e-01  CG iters=   39 rel_resid=6.28e-13  (0.70s)
    shift[ 7] sigma=5.880862e+00 residue=-1.677783e+00  CG iters=   27 rel_resid=6.02e-13  (0.44s)
    shift[ 8] sigma=1.362508e+01 residue=-5.134075e+00  CG iters=   19 rel_resid=3.92e-13  (0.43s)
    shift[ 9] sigma=3.449013e+01 residue=-1.890993e+01  CG iters=   13 rel_resid=5.83e-13  (0.27s)
    shift[10] sigma=1.097714e+02 residue=-1.104013e+02  CG iters=    9 rel_resid=3.93e-13  (0.18s)
    shift[11] sigma=7.996778e+02 residue=-3.002811e+03  CG iters=    6 rel_resid=5.17e-14  (0.11s)
[task2/full v1] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task2_full_cfg16_12_v1.npy  norm²=1.220496e+06
```
