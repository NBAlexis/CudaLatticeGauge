# Task 5: Compare Average Plaquette

## Test Description

Compare the average plaquette value computed from identical gauge configurations.

- **CLGLib**: Uses `CFieldGaugeSU3::CalculatePlaqutteEnergyOriginal(1.0/3.0)` to get the pure plaquette energy, then computes `avgPlaq = 1.0 − energy/(6·V)`.
- **PyQUDA**: Uses `gauge_action.dirac.plaquette()` which returns `(total, spatial, temporal)` plaquette values; the `total` component is used for comparison.

### Key Parameters
- Same gauge configurations as Tasks 1-4
- No fermion fields involved

## How to Invoke

### CLGLib
```bash
cd /home/nbalexis/CLGLib/Bin/Debug
./HISQ ../Debug/HISQ.yaml EHJ_Task5Plaq <cfg_file> <out_file>
```

### PyQUDA
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
/home/nbalexis/pyqudaenv/bin/python3 task5.py
```

## Test Conclusion

**PASS** — All 5 configurations pass with relative error effectively zero.

Representative results (cfg16_12):
| | Value |
|---|-------|
| avgPlaq_clg | 0.6034114169406844 |
| avgPlaq_py | 0.6034114169406844 |
| relerr | 0.000e+00 |

The plaquette agrees to machine precision across all configurations.

## Representative Log

### CLGLib — Task5Plaq (cfg16_12)
```
[04-05-2026 15-42-24|][Task5Plaq] CfgFile   = /home/nbalexis/comparehisq/data/cfg16_12_double.con
[04-05-2026 15-42-24|][Task5Plaq] OutputFile= /home/nbalexis/comparehisq/data/outputs/clg_task5_plaq_cfg16_12_v1.txt
...
[04-05-2026 15-42-25|][Task5Plaq] gauge is tree-improved, RectOverPlaq=-0.064654
[04-05-2026 15-42-25|][Task5Plaq] computing average plaquette
[04-05-2026 15-42-25|][Task5Plaq] energy=155944.97627625184 avgPlaq=0.60341141694068445
[04-05-2026 15-42-25|][Task5Plaq] saved /home/nbalexis/comparehisq/data/outputs/clg_task5_plaq_cfg16_12_v1.txt
[04-05-2026 15-42-25|][Task5Plaq] done.
```

### PyQUDA — Task5 (cfg16_12)
```
[task5] plaquette
[task5] saved /home/nbalexis/comparehisq/data/outputs/pyquda_task5_plaq_cfg16_12.txt
  total=0.6034114169  spatial=0.6029718708  temporal=0.6038509631
```
