# Task 6: Compare Combined Multi-Flavor Fermion Action

## Test Description

Compare the combined fermion action for **3 independent staggered fermion fields** using `CActionFermionHISQCombined` in CLGLib against the equivalent PyQUDA computation.

The action is:
```
S = Σ_i v_i† · (M†M)_ee^{-1/4} · v_i   (i = 1, 2, 3)
```

- **CLGLib**: `CActionFermionHISQCombined` with `FieldIds: [2, 3, 4]` calls `EnergyS` on each field internally, summing the individual `phi†·D_MD(phi)` contributions.
- **PyQUDA**: Uses `HMC.fermionAction()` with 3 `HISQAction` monomials. Each monomial's `action()` calls `updateFatLong()` → `invertMultiShift("force")` with `compute_action=1`. The returned value omits the `norm * <b, b>` term, matching CLGLib's `bAction=TRUE` behavior.

### Key Parameters
- Same lattice and gauge as Tasks 1-5
- 3 independent random vectors `v1`, `v2`, `v3`
- `Mass = 0` on all fields; mass shift carried by rational offsets
- `m_bEvenPseudofermion = TRUE` on all fields; odd sites zeroed via `ZeroOnEvenOdd(FALSE)`
- Each field id (2, 3, 4) must have its own `CMultiShiftBiCGStab` solver registered in YAML

### Critical Fix: Per-Field Multi-Shift Solver
Initially only field id 2 had a solver registered. Fields 3 and 4 silently failed their `D_MD` calls, returning incorrect results (~508k vs expected ~346k). The fix was adding `MSSolver3` and `MSSolver4` blocks to the YAML configuration.

## File Layout (Independent Test)

| File | Purpose |
|------|---------|
| `task6.py` | PyQUDA driver |
| `Task6.yaml` | CLGLib YAML config for 3-field combined action |
| `run_all_clglib.py` | Runs all CLGLib standalone executables |

## How to Invoke

### CLGLib
```bash
cd /home/nbalexis/CLGLib/Bin/Ubuntu
./HISQ /home/nbalexis/CLGLib/OtherProgram/HisqTest/Task6.yaml EHJ_Task6FermionAction \
    <cfg_file> <v1_file> <v2_file> <v3_file> <out_file>
```

Or via the unified runner:
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
python3 run_all_clglib.py
```

### PyQUDA
```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
/home/nbalexis/pyqudaenv/bin/python3 task6.py
```

## Test Conclusion

**PASS** — All 5 configurations pass with ratio = 1.000000.

Representative results (cfg16_12):
| | Value |
|---|-------|
| S_clg | 345723.984776 |
| S_py | 345723.984776 |
| ratio | 1.000000 |

Per-field breakdown (PyQUDA, cfg16_12):
| Field | S_i |
|-------|-----|
| v1 | 114898.161740 |
| v2 | 115720.162140 |
| v3 | 115105.660896 |
| **Total** | **345723.984776** |

The combined action achieves machine-precision agreement (~1e-16).

## Representative Log

### CLGLib — Task6 (cfg16_12)
```
[Task6FermionAction] CfgFile=/home/nbalexis/comparehisq/data/cfg16_12_double.con
[Task6FermionAction] V1File=/home/nbalexis/comparehisq/data/v1_double.con
[Task6FermionAction] V2File=/home/nbalexis/comparehisq/data/v2_double.con
[Task6FermionAction] V3File=/home/nbalexis/comparehisq/data/v3_double.con
[Task6FermionAction] OutputFile=/home/nbalexis/comparehisq/data/outputs/clg_task6_fermion_cfg16_12_v1.txt
...
[Task6FermionAction] computing combined fermion action energy (3 fields)
[Task6FermionAction] S = 345723.98477570503
[Task6FermionAction] saved .../clg_task6_fermion_cfg16_12_v1.txt
[Task6FermionAction] done.
```

### PyQUDA — Task6 (cfg16_12)
```
========== cfg16_12 ==========
MultiShiftCG: Converged after 138 iterations
...
  v1: S_i = 114898.161740
  v2: S_i = 115720.162140
  v3: S_i = 115105.660896
  total: clg=345723.984776  py=345723.984776  ratio=1.000000
```
