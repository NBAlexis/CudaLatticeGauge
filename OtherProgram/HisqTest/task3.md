# Task 3: Compare Fermion Action S = v_e†·(M†M)_ee^{-1/4}·v_e

## Test Description

Compare the even-even preconditioned fermion action scalar value:
```
S = v_e† · (M†M)_ee^{-1/4} · v_e
```

This is computed differently on each side:
- **CLGLib**: Uses `CActionFermionHISQCombined::Energy()` which internally calls `D_MD` (the `MD` rational approximation, i.e. `(M†M)^{-1/4}`) and returns `phi†·D_MD(phi)`.
- **PyQUDA**: Uses `HMC.fermionAction()` which calls `loadGaugeMomSmeared()` → `monomial.action()` → `updateFatLong()` → `invertMultiShift("force")` with `compute_action=1`. The returned value omits the `norm * <b, b>` term, matching CLGLib's `bAction=TRUE` behavior.

The test is split into two independent sub-tests:

### Test A: Manually injected neg14 rational coefficients
All rational slots (force/sample/action) are injected with `rational_neg14.json` coefficients. This makes PyQUDA's rational identical to CLGLib's MD rational.

### Test B: Native `staggeredFermionRationalParam` force rational
Uses PyQUDA's `staggeredFermionRationalParam(mass=[0.11], num_flavor=[1], ...)` to generate its own force rational. The CLGLib side uses matching native MD coefficients in `Task3TestB.yaml`.

### Key Parameters
- Same lattice, gauge, and `v` as Tasks 1-2
- `Mass = 0` in both sides; mass shift carried by rational offsets
- `m_bEvenPseudofermion = TRUE` / odd sites zeroed before computation
- Rational: `f(x) = (x + 4m²)^{-1/4}`

## File Layout (Independent Tests)

Each sub-test is completely independent:

| File | Purpose |
|------|---------|
| `task3_testA.py` | PyQUDA driver for Test A |
| `task3_testB.py` | PyQUDA driver for Test B |
| `Task3TestA.yaml` | CLGLib YAML config for Test A (neg14 MD coeffs) |
| `Task3TestB.yaml` | CLGLib YAML config for Test B (native MD coeffs) |
| `run_all_clglib.py` | Runs all CLGLib standalone executables |

## How to Invoke

### CLGLib — Test A
```bash
cd /home/nbalexis/CLGLib/Bin/Ubuntu
./HISQ /home/nbalexis/CLGLib/OtherProgram/HisqTest/Task3TestA.yaml EHJ_Task3EO \
    <cfg_file> <v_file> <out_file>
```

### CLGLib — Test B
```bash
cd /home/nbalexis/CLGLib/Bin/Ubuntu
./HISQ /home/nbalexis/CLGLib/OtherProgram/HisqTest/Task3TestB.yaml EHJ_Task3EO \
    <cfg_file> <v_file> <out_file>
```

Or run all CLGLib tests at once:
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
python3 run_all_clglib.py
```

### PyQUDA — Test A
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
/home/nbalexis/pyqudaenv/bin/python3 task3_testA.py
```

### PyQUDA — Test B
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
/home/nbalexis/pyqudaenv/bin/python3 task3_testB.py
```

## Test Conclusion

**PASS** — All 5 configurations × 3 vectors × 2 sub-tests pass with ratio = 1.000000.

Representative results (Test A, cfg16_12):
| v | S_clg | S_py | ratio |
|---|-------|------|-------|
| v1 | 114898.161740 | 114898.161740 | 1.000000 |
| v2 | 115720.162140 | 115720.162140 | 1.000000 |
| v3 | 115105.660896 | 115105.660896 | 1.000000 |

Representative results (Test B, cfg16_12):
| v | S_clg | S_py | ratio |
|---|-------|------|-------|
| v1 | 114898.161740 | 114898.161740 | 1.000000 |
| v2 | 115720.162140 | 115720.162140 | 1.000000 |
| v3 | 115105.660896 | 115105.660896 | 1.000000 |

Both sub-tests achieve machine-precision agreement (~1e-16).

## Representative Log

### CLGLib — Task3TestA (cfg16_12, v1)
```
[Task3EO] CfgFile   = /home/nbalexis/comparehisq/data/cfg16_12_double.con
[Task3EO] VFile     = /home/nbalexis/comparehisq/data/v1_double.con
[Task3EO] OutputFile= /home/nbalexis/comparehisq/data/outputs/clg_task3a_eo_cfg16_12_v1.txt
...
[Task3EO] S = 114898.16173995406
[Task3EO] saved .../clg_task3a_eo_cfg16_12_v1.txt
[Task3EO] done.
```

### PyQUDA — Test A (cfg16_12, v1)
```
========== cfg16_12 ==========
MultiShiftCG: Converged after 138 iterations
...
  v1: clg=114898.161740  py=114898.161740  ratio=1.000000
```
