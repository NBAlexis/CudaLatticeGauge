#!/usr/bin/env python3
"""Task 6 via HMC flow (combined 3-field action)."""

import json
import os

import numpy as np

from pyquda import init
from pyquda.action.abstract import RationalParam
from pyquda.action.hisq import HISQAction
from pyquda.field import LatticeInfo, LatticeStaggeredFermion, MultiLatticeStaggeredFermion
from pyquda.hmc import HMC, Integrator
from pyquda_utils.io import readNPYGauge

LATT_SIZE = [16, 16, 16, 16]
T_BOUNDARY = -1
ANISOTROPY = 1.0
NAIK_EPS = -0.11
CFG_PREFIX = "cfg16"
CFG_LIST = [12, 14, 16, 18, 20]
V_LIST = [1, 2, 3]

HERE = os.path.dirname(os.path.abspath(__file__))
COMPAREHISQ = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), "comparehisq")
DATA_DIR = os.path.join(COMPAREHISQ, "data")
OUT_DIR = os.path.join(DATA_DIR, "outputs")


class DummyIntegrator(Integrator):
    def integrate(self, updateGauge, updateMom, t):
        pass


def load_rational(name):
    path = os.path.join(DATA_DIR, name)
    with open(path) as f:
        d = json.load(f)
    return d["norm"], d["residue"], d["offset"]


def read_clg(cfg_idx):
    path = os.path.join(OUT_DIR, f"clg_task6_fermion_{CFG_PREFIX}_{cfg_idx}_v1.txt")
    with open(path) as f:
        return float(f.read().strip())


def main():
    init([1, 1, 1, 1], backend="numpy", resource_path=os.path.join(COMPAREHISQ, ".cache"))

    norm_n, res_n, off_n = load_rational("rational_neg14.json")

    for cfg_idx in CFG_LIST:
        print(f"\n========== {CFG_PREFIX}_{cfg_idx} ==========")
        latt_info = LatticeInfo(LATT_SIZE, t_boundary=T_BOUNDARY, anisotropy=ANISOTROPY)
        gauge = readNPYGauge(os.path.join(DATA_DIR, f"{CFG_PREFIX}_{cfg_idx}.npy"))

        v_fields = {}
        for v_idx in V_LIST:
            v_lex = np.load(os.path.join(DATA_DIR, f"v{v_idx}.npy"))
            v_fields[v_idx] = LatticeStaggeredFermion(latt_info, latt_info.evenodd(v_lex, False))

        rp = RationalParam(
            norm_force=norm_n, residue_force=res_n, offset_force=off_n,
            norm_sample=norm_n, residue_sample=res_n, offset_sample=off_n,
            norm_action=norm_n, residue_action=res_n, offset_action=off_n,
        )

        # Create 3 HISQActions, one per field
        acts = []
        for _ in V_LIST:
            act = HISQAction(latt_info, rp, 1e-12, 10000, naik_epsilon=NAIK_EPS)
            act.quark = MultiLatticeStaggeredFermion(latt_info, act.max_num_offset)
            acts.append(act)

        hmc = HMC(latt_info=latt_info, monomials=acts, integrator=DummyIntegrator(1))
        gauge.toDevice()
        hmc.initialize(1234, gauge)

        for i, v_idx in enumerate(V_LIST):
            acts[i].phi.data[:] = 0
            acts[i].phi.even = v_fields[v_idx].even
            acts[i].quark = MultiLatticeStaggeredFermion(latt_info, acts[i].max_num_offset)

        total_py = hmc.fermionAction()
        S_clg = read_clg(cfg_idx)
        print(f"  total: clg={S_clg:.6f}  py={total_py:.6f}  ratio={total_py/S_clg:.6f}")

    print("\nDone.")


if __name__ == "__main__":
    main()
