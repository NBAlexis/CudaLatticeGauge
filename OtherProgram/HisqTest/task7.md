# Task 7: Gauge Momentum Convention Verification

## Objective

Verify that setting `GaugeMomentumFactor = 2.0` in CLGLib makes its gauge momentum conventions fully consistent with QUDA/PyQUDA.

## Background

Both CLGLib and QUDA use the same HMC update equations (`P += dt*F`, `U = exp(dt*P)*U`), but their momentum definitions differ by a factor:

- **QUDA** momentum action per link: `0.5*tr(P†P) - 4.0`
- **CLGLib** `LengthSq` per link: `tr(P†P) - 8.0`

With `GaugeMomentumFactor = 2.0`:
- CLGLib random momentum generation multiplies by `sqrt(2.0)`
- CLGLib `CalculateKinematicEnergy` divides `LengthSq` by `2.0`

This makes:
```
CLGLib CalculateKinematicEnergy = LengthSq / 2.0
                                = sum[tr(P†P) - 8.0] / 2.0
                                = sum[0.5*tr(P†P) - 4.0]
                                = QUDA momAction
```

## Two Verifications

### Verification 1: Random momentum norm comparison

Generate random Gaussian momentum independently in both libraries and compare the per-link average of `tr(P†P)`.

- **QUDA**: `gaussMomQuda(seed, sigma=1.0)` (Box-Muller with `radius = sqrt(-log(u))`)
- **CLGLib**: `InitialField(EFIT_RandomGenerator)` (Box-Muller with `radius = sqrt(-0.5*log(u))`, then scaled by `sqrt(2.0)`)

Both should produce statistically equivalent norms (same distribution, not exact match due to different RNGs).

### Verification 2: Same-field energy comparison

Use the **same** momentum field in both libraries and verify the energy/action ratio is exactly 1.0.

1. QUDA generates random momentum
2. Python reconstructs the full 3x3 matrices from QUDA's reconstruct=10 format
3. The field is saved in CLGLib binary format via `transformToCon`
4. CLGLib loads it and computes `CalculateKinematicEnergy`
5. Compare with QUDA `momActionQuda` (which includes the -4.0 per-link offset)

## Key Parameters

- Lattice: `4x4x4x12`
- `t_boundary = -1`, anisotropy = `1.0`
- Seed: `12345` (QUDA), `1234567` (CLGLib YAML)
- `GaugeMomentumFactor = 2.0`

## How to Run

```bash
cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
/home/nbalexis/pyqudaenv/bin/python3 task7.py
```

The script automatically:
1. Injects parameters into `Task7Mom.yaml`
2. Runs QUDA/PyQUDA for both verifications
3. Calls CLGLib `HISQ` binary for both verifications
4. Compares results and prints PASS/FAIL

## Results

### Verification 1: Random momentum norm

| Metric | QUDA | CLGLib |
|--------|------|--------|
| Per-link avg `tr(P†P)` | 7.9558 | 7.9380 |
| Relative difference | 0.22% | |
| Verdict | **PASS** (statistically consistent) | |

### Verification 2: Same-field energy/action

| Metric | Value |
|--------|-------|
| QUDA `momAction` (with -4 offset) | 37.6967204211 |
| CLGLib `CalculateKinematicEnergy` | 37.6967204211 |
| Ratio | **1.0000000000** |
| Difference from 1.0 | 0.000e+00 |
| Verdict | **PASS** (exact match to machine precision) |

## Implementation Details

### CLGLib side

- `Code/Applications/HISQ/Task7Mom.cpp` — Task implementation with two modes:
  - `GenerateRandom=1`: generate random Gaussian momentum, save to `OutMomFile`
  - `GenerateRandom=0`: load momentum from `MomFile`, compute `CalculateKinematicEnergy`
- `Code/Applications/HISQ/HISQ.h` — added `EHJ_Task7Mom` enum + declaration
- `Code/Applications/HISQ/HISQ.cpp` — added switch case for `EHJ_Task7Mom`
- `Code/CMake/CMakeLists.txt` — added `Task7Mom.cpp` to `add_executable(HISQ ...)`

### Python side

- `OtherProgram/HisqTest/task7.py` — end-to-end verification script
- `OtherProgram/HisqTest/Task7Mom.yaml` — standalone YAML config with `GaugeMomentumFactor: 2.0`

## Bug Fixed

**Issue**: `Task7Mom.cpp` originally used `EFIT_RandomGaussian`, but `_kernelInitialLink` (the gauge field initialization kernel) does **not** handle this case — it only has `EFIT_RandomGenerator`. As a result, the field was not reinitialized with random momentum; it retained the dummy gauge configuration loaded by `SetupGaugeOnly`, giving `tr(P†P) = 6.0` per link instead of the expected ~8.0.

**Fix**: Changed `EFIT_RandomGaussian` to `EFIT_RandomGenerator` in `Task7Mom.cpp`.
