"""Generate shared inputs for the CLGLib vs PyQUDA HISQ comparison.

Outputs (under data/):
- v{1,2,3}.npy                random staggered fermion lexico [lt, lz, ly, lx, 3] complex128
- v{1,2,3}_double.con         same v in CLGLib double-precision binary layout
- {CFG_PREFIX}_NN_double.con  for NN in CFG_LIST
- rational_pos14.json         (norm, residue, offset) for f(x) = (x + 4m^2)^(+1/4)
- rational_neg14.json         (norm, residue, offset) for f(x) = (x + 4m^2)^(-1/4)
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, "/home/nbalexis/CLGLib/OtherProgram/NpyConfiguration")
from Converter.NpyConfigurationToCon import (  # noqa: E402
    transformStaggeredFermionFromCon,
    transformStaggeredFermionToCon,
    transformToCon,
)

from pyquda_utils.alg_remez import AlgRemez  # noqa: E402
from pyquda_utils.hmc_param import _StaggeredMD  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
# 16^4 cfgs are produced locally by test_hmc_hisq_16.py and saved under DATA.
PYQUDA_DATA = DATA

LX, LY, LZ, LT = 16, 16, 16, 16
MASS = 0.11  # PyQUDA convention; 2am = 0.22

# Post-warmup HMC trajectories saved by test_hmc_hisq_16.py: cfg16_{12,14,16,18,20}.npy
CFG_LIST = [12, 14, 16, 18, 20]
CFG_PREFIX = "cfg16"

# Three independent random fermion fields v1,v2,v3 with distinct seeds.
V_LIST = [1, 2, 3]
V_SEEDS = {1: 1234, 2: 5678, 3: 9012}


def gen_v():
    for v_idx in V_LIST:
        seed = V_SEEDS[v_idx]
        rng = np.random.default_rng(seed=seed)
        v_lex = (rng.standard_normal((LT, LZ, LY, LX, 3))
                 + 1j * rng.standard_normal((LT, LZ, LY, LX, 3))).astype(np.complex128)

        v_npy = os.path.join(DATA, f"v{v_idx}.npy")
        np.save(v_npy, v_lex)
        print(f"[v{v_idx}] saved {v_npy} shape={v_lex.shape} dtype={v_lex.dtype} "
              f"seed={seed} norm={np.linalg.norm(v_lex):.6e}")

        v_con = os.path.join(DATA, f"v{v_idx}_double.con")
        transformStaggeredFermionToCon(v_con, v_lex, double=True)
        print(f"[v{v_idx}] saved {v_con}")

        # Round-trip self-check
        v_back = transformStaggeredFermionFromCon(v_con, LX, LY, LZ, LT, double=True)
        diff = np.abs(v_back.astype(np.complex128) - v_lex).max()
        print(f"[v{v_idx}] round-trip max abs diff = {diff:.3e}")
        assert diff < 1e-6, f"v{v_idx} round-trip mismatch"


def gen_cfg():
    for n in CFG_LIST:
        src = os.path.join(PYQUDA_DATA, f"{CFG_PREFIX}_{n}.npy")
        dst = os.path.join(DATA, f"{CFG_PREFIX}_{n}_double.con")
        a = np.load(src)
        transformToCon(dst, a, double=True)
        print(f"[cfg] {src} -> {dst}, shape={a.shape} dtype={a.dtype}")


def gen_remez():
    """Generate Remez approximations for (x+4m^2)^{nf/4}, nf in {+1, -1}.

    Bounds chosen for HISQ on small lattices: (M^dagger M) eigenvalues lie roughly
    between (2am)^2 ~ 0.0484 (low) and ~ a few (high). We pad both sides for safety.
    """
    LOW = 1e-3
    HIGH = 50.0
    PRECISION_BITS = 40
    DEGREE = 12

    # f(x) = (x + 4m^2)^(+1/4) -> PFE
    remez = AlgRemez(LOW, HIGH, PRECISION_BITS)
    remez.generateApprox(DEGREE, DEGREE, _StaggeredMD([MASS], [+1]))
    norm_p, res_p, off_p = remez.getPFE()
    print(f"\n[remez +1/4] norm={norm_p}")
    for r, b in zip(res_p, off_p):
        print(f"  r={r:.10e}  b={b:.10e}")
    err_p = _check_remez_error(norm_p, res_p, off_p, MASS, +1, LOW, HIGH)
    print(f"[remez +1/4] max relative err on [{LOW},{HIGH}]: {err_p:.3e}")

    with open(os.path.join(DATA, "rational_pos14.json"), "w") as fp:
        json.dump({
            "exponent": "+1/4",
            "function": "(x + 4*m^2)^{+1/4}",
            "mass": MASS,
            "low": LOW, "high": HIGH, "degree": DEGREE,
            "norm": norm_p, "residue": res_p, "offset": off_p,
            "clglib_flat": [norm_p] + res_p + off_p,
            "max_relative_error": err_p,
        }, fp, indent=2)

    # f(x) = (x + 4m^2)^(-1/4) -> approximated directly (NOT via getIPFE on +1/4 Remez)
    remez = AlgRemez(LOW, HIGH, PRECISION_BITS)
    remez.generateApprox(DEGREE, DEGREE, _StaggeredMD([MASS], [-1]))
    norm_n, res_n, off_n = remez.getPFE()
    print(f"\n[remez -1/4] norm={norm_n}")
    for r, b in zip(res_n, off_n):
        print(f"  r={r:.10e}  b={b:.10e}")
    err_n = _check_remez_error(norm_n, res_n, off_n, MASS, -1, LOW, HIGH)
    print(f"[remez -1/4] max relative err on [{LOW},{HIGH}]: {err_n:.3e}")

    with open(os.path.join(DATA, "rational_neg14.json"), "w") as fp:
        json.dump({
            "exponent": "-1/4",
            "function": "(x + 4*m^2)^{-1/4}",
            "mass": MASS,
            "low": LOW, "high": HIGH, "degree": DEGREE,
            "norm": norm_n, "residue": res_n, "offset": off_n,
            "clglib_flat": [norm_n] + res_n + off_n,
            "max_relative_error": err_n,
        }, fp, indent=2)


def _check_remez_error(norm, residue, offset, m, nf, low, high, samples=200):
    """Sample relative error of norm + sum r_i / (x + b_i) vs (x+4m^2)^{nf/4}."""
    xs = np.geomspace(low, high, samples)
    target = (xs + 4 * m ** 2) ** (nf / 4)
    approx = norm + np.zeros_like(xs)
    for r, b in zip(residue, offset):
        approx = approx + r / (xs + b)
    return float(np.max(np.abs(approx - target) / np.abs(target)))


def main():
    os.makedirs(DATA, exist_ok=True)
    gen_v()
    gen_cfg()
    gen_remez()
    print("\nDone.")


if __name__ == "__main__":
    main()
