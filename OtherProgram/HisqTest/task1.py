#!/usr/bin/env python3
"""Task 1: Compare HISQ Dirac operator M·v.

Task1Full: pFerm <- (2am + D)·v  (full lattice)
Task1EO:   pFerm <- (M†M)_ee · v_e  (even-even preconditioned)

QUDA mass sign convention:
  M_QUDA(m)·v = 2m·v − H(v)
  M_clg·v     = 2am·v + H(v)
So PyQUDA uses mass=-0.11 and negates the output for Full mode.
M†M is sign-invariant, so EO mode is unaffected.
"""

import os

import numpy as np

from pyquda import init
from pyquda.dirac.hisq import HISQDirac
from pyquda.enum_quda import QudaMassNormalization, QudaMatPCType, QudaSolutionType
from pyquda.field import LatticeInfo, LatticeStaggeredFermion
from pyquda_utils.io import readNPYGauge

LATT_SIZE = [16, 16, 16, 16]
T_BOUNDARY = -1
ANISOTROPY = 1.0
MASS = 0.11
NAIK_EPS = -0.11
TOL = 1e-12
MAXITER = 10000
CFG_LIST = [12, 14, 16, 18, 20]
CFG_PREFIX = "cfg16"
V_LIST = [1, 2, 3]

HERE = os.path.dirname(os.path.abspath(__file__))
# Data and outputs live in the comparehisq workspace
COMPAREHISQ = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), "comparehisq")
DATA = os.path.join(COMPAREHISQ, "data")
OUT = os.path.join(DATA, "outputs")


def make_v(latt_info, v_idx):
    v_lex = np.load(os.path.join(DATA, f"v{v_idx}.npy")).astype(np.complex128)
    return LatticeStaggeredFermion(latt_info, latt_info.evenodd(v_lex, False))


def lexico_eo_only(field_lex, latt_info):
    out = field_lex.copy()
    Lt, Lz, Ly, Lx = latt_info.size[3], latt_info.size[2], latt_info.size[1], latt_info.size[0]
    for t in range(Lt):
        for z in range(Lz):
            for y in range(Ly):
                for x in range(Lx):
                    if (x + y + z + t) % 2 == 1:
                        out[t, z, y, x, :] = 0
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    init([1, 1, 1, 1], backend="numpy", resource_path=os.path.join(COMPAREHISQ, ".cache"))

    for cfg_idx in CFG_LIST:
        print(f"\n========== {CFG_PREFIX}_{cfg_idx} ==========")
        latt_info = LatticeInfo(LATT_SIZE, t_boundary=T_BOUNDARY, anisotropy=ANISOTROPY)
        gauge = readNPYGauge(os.path.join(DATA, f"{CFG_PREFIX}_{cfg_idx}.npy"))

        v_fields = {v_idx: make_v(latt_info, v_idx) for v_idx in V_LIST}

        print("[task1] dirac with mass=-0.11 (so -M_QUDA(-m) = M_clg)")
        d1 = HISQDirac(latt_info, -MASS, TOL, MAXITER, NAIK_EPS)
        d1.invert_param.mass_normalization = QudaMassNormalization.QUDA_MASS_NORMALIZATION
        d1.loadGauge(gauge)

        for v_idx in V_LIST:
            v_field = v_fields[v_idx]
            # Full
            d1.invert_param.solution_type = QudaSolutionType.QUDA_MAT_SOLUTION
            out = d1.mat(v_field)
            out_lex = -out.lexico()
            out_path = os.path.join(OUT, f"pyquda_task1_full_{CFG_PREFIX}_{cfg_idx}_v{v_idx}.npy")
            np.save(out_path, out_lex)
            print(f"  [full v{v_idx}] saved {out_path}  norm²={float(np.vdot(out_lex, out_lex).real):.6e}")

            # EO
            d1.invert_param.solution_type = QudaSolutionType.QUDA_MATPC_SOLUTION
            d1.invert_param.matpc_type = QudaMatPCType.QUDA_MATPC_EVEN_EVEN
            out = d1.mat(v_field)
            out_lex = lexico_eo_only(out.lexico(), latt_info)
            out_path = os.path.join(OUT, f"pyquda_task1_eo_{CFG_PREFIX}_{cfg_idx}_v{v_idx}.npy")
            np.save(out_path, out_lex)
            print(f"  [eo   v{v_idx}] saved {out_path}  norm²(even)={float(np.vdot(out_lex, out_lex).real):.6e}")

        d1.freeGauge()

    print("\nDone.")


if __name__ == "__main__":
    main()
