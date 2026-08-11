#!/usr/bin/env python3
"""Task 2: Compare rational approximation (M†M)^{1/4}·v.

Task2Full: f(D†D)·v  on full lattice (custom Python multi-shift CG)
Task2EO:   f(D†_eo·D_oe)·v_e  via HISQAction.invertMultiShift("action")

Rational: f(x) = (x + 4m²)^{+1/4} from rational_pos14.json.
"""

import json
import os
import time

import numpy as np

from pyquda import init
from pyquda.action import HISQAction
from pyquda.action.abstract import RationalParam
from pyquda.dirac.hisq import HISQDirac
from pyquda.enum_quda import QudaMassNormalization, QudaSolutionType
from pyquda.field import LatticeInfo, LatticeStaggeredFermion, MultiLatticeStaggeredFermion
from pyquda_utils.io import readNPYGauge

LATT_SIZE = [16, 16, 16, 16]
T_BOUNDARY = -1
ANISOTROPY = 1.0
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


def load_rational(name):
    with open(os.path.join(DATA, name)) as fp:
        d = json.load(fp)
    return float(d["norm"]), [float(r) for r in d["residue"]], [float(b) for b in d["offset"]]


def cg_shifted(matvec, b_vec, sigma, tol, maxiter):
    x = np.zeros_like(b_vec)
    r = b_vec.copy()
    p = b_vec.copy()
    rsq = float(np.vdot(r, r).real)
    rsq_init = rsq
    if rsq_init == 0.0:
        return x, 0, 0.0
    for it in range(maxiter):
        Ap = matvec(p) + sigma * p
        pAp = float(np.vdot(p, Ap).real)
        alpha = rsq / pAp
        x += alpha * p
        r -= alpha * Ap
        rsq_new = float(np.vdot(r, r).real)
        if np.sqrt(rsq_new / rsq_init) < tol:
            return x, it + 1, np.sqrt(rsq_new / rsq_init)
        beta = rsq_new / rsq
        p = r + beta * p
        rsq = rsq_new
    return x, maxiter, np.sqrt(rsq / rsq_init)


def rational_apply_full(dirac, v_full, norm, residue, offset, tol, maxiter):
    dirac.invert_param.solution_type = QudaSolutionType.QUDA_MAT_SOLUTION

    def matvec(x_lex):
        latt_info = v_full.latt_info
        x_field = LatticeStaggeredFermion(latt_info, latt_info.evenodd(x_lex.reshape(*v_full.lexico().shape), False))
        out = dirac.matDagMat(x_field)
        return out.lexico()

    v_lex = v_full.lexico()
    result = norm * v_lex.copy()
    for i, (r_i, b_i) in enumerate(zip(residue, offset)):
        t0 = time.perf_counter()
        x_i, iters, rrel = cg_shifted(matvec, v_lex, b_i, tol, maxiter)
        dt = time.perf_counter() - t0
        result += r_i * x_i
        print(f"    shift[{i:2d}] sigma={b_i:.6e} residue={r_i:+.6e}  CG iters={iters:5d} rel_resid={rrel:.2e}  ({dt:.2f}s)")
    return result


def make_action(latt_info, norm, residue, offset):
    rp = RationalParam(
        norm_force=norm, residue_force=residue, offset_force=offset,
        norm_sample=norm, residue_sample=residue, offset_sample=offset,
        norm_action=norm, residue_action=residue, offset_action=offset,
    )
    act = HISQAction(latt_info, rp, TOL, MAXITER, naik_epsilon=NAIK_EPS)
    act.quark = MultiLatticeStaggeredFermion(latt_info, act.max_num_offset)
    return act


def apply_rational_eo(act, v):
    act.phi.data[:] = 0
    act.eta.data[:] = 0
    act.phi.even = v.even
    act.invert_param.compute_action = 0
    act.invertMultiShift("action")
    return act.eta


def main():
    os.makedirs(OUT, exist_ok=True)
    init([1, 1, 1, 1], backend="numpy", resource_path=os.path.join(COMPAREHISQ, ".cache"))

    norm_p, res_p, off_p = load_rational("rational_pos14.json")

    for cfg_idx in CFG_LIST:
        print(f"\n========== {CFG_PREFIX}_{cfg_idx} ==========")
        latt_info = LatticeInfo(LATT_SIZE, t_boundary=T_BOUNDARY, anisotropy=ANISOTROPY)
        gauge = readNPYGauge(os.path.join(DATA, f"{CFG_PREFIX}_{cfg_idx}.npy"))

        v_fields = {v_idx: make_v(latt_info, v_idx) for v_idx in V_LIST}

        # Task 2 Full
        print("[task2] dirac with mass=0 (rational carries the +4m² shift)")
        d0 = HISQDirac(latt_info, 0.0, TOL, MAXITER, NAIK_EPS)
        d0.invert_param.mass_normalization = QudaMassNormalization.QUDA_MASS_NORMALIZATION
        d0.loadGauge(gauge)

        for v_idx in V_LIST:
            v_field = v_fields[v_idx]
            print(f"  [full v{v_idx}] (M†M)^(1/4)·v via custom multi-shift CG ({len(off_p)} shifts)")
            out_lex = rational_apply_full(d0, v_field, norm_p, res_p, off_p, TOL, MAXITER)
            out_path = os.path.join(OUT, f"pyquda_task2_full_{CFG_PREFIX}_{cfg_idx}_v{v_idx}.npy")
            np.save(out_path, out_lex)
            print(f"  [full v{v_idx}] saved {out_path}  norm²={float(np.vdot(out_lex, out_lex).real):.6e}")
        d0.freeGauge()

        # Task 2 EO
        print("  [eo  ] (M†M)_ee^(1/4)·v_e via HISQAction.invertMultiShift(action)")
        act_pos = make_action(latt_info, norm_p, res_p, off_p)
        act_pos.dirac.loadGauge(gauge)
        for v_idx in V_LIST:
            v_field = v_fields[v_idx]
            out = apply_rational_eo(act_pos, v_field)
            out_lex = lexico_eo_only(out.lexico(), latt_info)
            out_path = os.path.join(OUT, f"pyquda_task2_eo_{CFG_PREFIX}_{cfg_idx}_v{v_idx}.npy")
            np.save(out_path, out_lex)
            print(f"  [eo   v{v_idx}] saved {out_path}  norm²(even)={float(np.vdot(out_lex, out_lex).real):.6e}")
        act_pos.dirac.freeGauge()

    print("\nDone.")


if __name__ == "__main__":
    main()
