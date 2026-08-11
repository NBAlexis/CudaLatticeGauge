#!/usr/bin/env python3
"""Task 5: Compare average plaquette.

PyQUDA: gauge_action.dirac.plaquette() returns (total, spatial, temporal).
CLGLib: CalculatePlaqutteEnergyOriginal then avgPlaq = 1 - energy/(6*V).
"""

import os

import numpy as np

from pyquda import init
from pyquda.action.gauge import GaugeAction
from pyquda.field import LatticeInfo
from pyquda_utils.hmc_param import symanzikTreeGaugeLoopParam
from pyquda_utils.io import readNPYGauge

LATT_SIZE = [16, 16, 16, 16]
T_BOUNDARY = -1
ANISOTROPY = 1.0
NAIK_EPS = -0.11
CFG_LIST = [12, 14, 16, 18, 20]
CFG_PREFIX = "cfg16"

U0 = 0.8794
BETA_TASK4 = 7.29

HERE = os.path.dirname(os.path.abspath(__file__))
# Data and outputs live in the comparehisq workspace
COMPAREHISQ = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), "comparehisq")
DATA = os.path.join(COMPAREHISQ, "data")
OUT = os.path.join(DATA, "outputs")


def main():
    os.makedirs(OUT, exist_ok=True)
    init([1, 1, 1, 1], backend="numpy", resource_path=os.path.join(COMPAREHISQ, ".cache"))

    for cfg_idx in CFG_LIST:
        print(f"\n========== {CFG_PREFIX}_{cfg_idx} ==========")
        latt_info = LatticeInfo(LATT_SIZE, t_boundary=T_BOUNDARY, anisotropy=ANISOTROPY)
        gauge = readNPYGauge(os.path.join(DATA, f"{CFG_PREFIX}_{cfg_idx}.npy"))

        print("[task5] plaquette")
        loop_param = symanzikTreeGaugeLoopParam(U0)
        gauge_action = GaugeAction(latt_info, loop_param, BETA_TASK4)
        gauge_action.dirac.loadGauge(gauge)
        _ = gauge_action.dirac.plaquette()
        plaq = gauge_action.dirac.plaquette()
        total_plaq = float(plaq[0])
        out_path = os.path.join(OUT, f"pyquda_task5_plaq_{CFG_PREFIX}_{cfg_idx}.txt")
        with open(out_path, "w") as fp:
            fp.write(f"{total_plaq:.17g}\n")
        print(f"  saved {out_path}  total={total_plaq:.10f}  spatial={float(plaq[1]):.10f}  temporal={float(plaq[2]):.10f}")
        gauge_action.dirac.freeGauge()

    print("\nDone.")


if __name__ == "__main__":
    main()
