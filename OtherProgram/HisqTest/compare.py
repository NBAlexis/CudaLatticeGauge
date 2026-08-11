"""Compare CLGLib vs PyQUDA HISQ outputs across 5 cfgs x 3 v's x 5 sub-tasks
plus 2 gauge-only tasks (Task4 action, Task5 plaquette).

Reads:
- data/outputs/clg_<task>_<mode>_<CFG_PREFIX>_<NN>_v<i>.con      (CLG vector outputs, double precision)
- data/outputs/clg_<task>_<mode>_<CFG_PREFIX>_<NN>_v<i>.txt      (CLG scalar action, Task3 only)
- data/outputs/pyquda_<task>_<mode>_<CFG_PREFIX>_<NN>_v<i>.npy   (PyQUDA vector outputs, lexico)
- data/outputs/pyquda_<task>_<mode>_<CFG_PREFIX>_<NN>_v<i>.txt   (PyQUDA scalar action, Task3 only)
- data/outputs/clg_task4_action_<CFG_PREFIX>_<NN>_v<i>.txt       (CLG tree-improved action)
- data/outputs/pyquda_task4_action_<CFG_PREFIX>_<NN>.txt         (PyQUDA tree-improved action)
- data/outputs/clg_task5_plaq_<CFG_PREFIX>_<NN>_v<i>.txt         (CLG plaquette)
- data/outputs/pyquda_task5_plaq_<CFG_PREFIX>_<NN>.txt           (PyQUDA plaquette)

Reports |CLG - PyQUDA| / |PyQUDA| for each, and the largest 5 |Delta| sites.
For Task4, subtracts the known constant offset before comparison.
"""
import os
import sys

import numpy as np

sys.path.insert(0, "/home/nbalexis/CLGLib/OtherProgram/NpyConfiguration")
from Converter.NpyConfigurationToCon import transformStaggeredFermionFromCon  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
# Data lives in the comparehisq workspace (sibling of CLGLib/OtherProgram)
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), "comparehisq", "data")
OUT = os.path.join(DATA, "outputs")

LX, LY, LZ, LT = 16, 16, 16, 16
CFG_LIST = [12, 14, 16, 18, 20]
CFG_PREFIX = "cfg16"
V_LIST = [1, 2, 3]

# Task 4 constant offset: CLGLib uses (beta/Nc)*[sum(Nc-ReTr) + Crect*sum(Nc-ReTr_rect)]
# while PyQUDA uses -(beta/Nc)*[sum(ReTr) + c_rect*sum(ReTr_rect)].
# The difference is beta * V * (6 + 12*Crect).
U0 = 0.8794
BETA_TASK4 = 7.29
VOLUME = LX * LY * LZ * LT
CRECT = -1.0 / (20.0 * U0 * U0)
TASK4_OFFSET = BETA_TASK4 * VOLUME * (6.0 + 12.0 * CRECT)


def cmp_vec(label: str, clg_path: str, py_path: str) -> None:
    if not os.path.isfile(clg_path):
        print(f"[{label}] MISSING CLG: {clg_path}")
        return
    if not os.path.isfile(py_path):
        print(f"[{label}] MISSING PY : {py_path}")
        return
    clg = transformStaggeredFermionFromCon(clg_path, LX, LY, LZ, LT, double=True)
    py = np.load(py_path)
    if clg.shape != py.shape:
        print(f"[{label}] SHAPE mismatch clg={clg.shape} py={py.shape}")
        return
    diff = clg - py
    py_norm = float(np.linalg.norm(py))
    diff_norm = float(np.linalg.norm(diff))
    max_abs = float(np.abs(diff).max())
    rel = diff_norm / py_norm if py_norm > 0 else float("nan")
    flag = "OK " if rel < 1e-6 else "BAD"
    print(f"[{flag}] {label:44s} |Δ|={diff_norm:.3e} ratio={rel:.3e} maxabs={max_abs:.3e} |py|={py_norm:.3e}")
    if rel >= 1e-6 or max_abs > 1e-4:
        flat = np.abs(diff).reshape(-1)
        top = np.argsort(flat)[::-1][:5]
        for idx in top:
            t, z, y, x, c = np.unravel_index(idx, diff.shape)
            print(f"    site=(x={x},y={y},z={z},t={t},c={c}) "
                  f"clg={clg[t, z, y, x, c]} py={py[t, z, y, x, c]} Δ={diff[t, z, y, x, c]}")


def cmp_scalar(label: str, clg_path: str, py_path: str, offset: float = 0.0) -> None:
    if not os.path.isfile(clg_path):
        print(f"[{label}] MISSING CLG: {clg_path}")
        return
    if not os.path.isfile(py_path):
        print(f"[{label}] MISSING PY : {py_path}")
        return
    s_clg = float(open(clg_path).read().strip())
    s_py = float(open(py_path).read().strip())
    s_clg_adj = s_clg - offset
    rel = abs(s_clg_adj - s_py) / abs(s_py) if abs(s_py) > 0 else float("nan")
    abs_err = abs(s_clg_adj - s_py)
    flag = "OK " if rel < 1e-6 else "BAD"
    print(f"[{flag}] {label:44s} S_clg={s_clg:.10f} S_py={s_py:.10f} relerr={rel:.3e} abserr={abs_err:.3e}")


def main():
    print("=" * 88)
    print(f"{'task/mode/cfg/v':44s} {'metric':>40s}")
    print("=" * 88)
    for cfg in CFG_LIST:
        for v in V_LIST:
            for mode in ["full", "eo"]:
                for task_n in [1, 2]:
                    label = f"task{task_n}/{mode}/{CFG_PREFIX}_{cfg}/v{v}"
                    clg = os.path.join(OUT, f"clg_task{task_n}_{mode}_{CFG_PREFIX}_{cfg}_v{v}.con")
                    py = os.path.join(OUT, f"pyquda_task{task_n}_{mode}_{CFG_PREFIX}_{cfg}_v{v}.npy")
                    cmp_vec(label, clg, py)
            # Task 3 is EO scalar only
            label = f"task3/eo/{CFG_PREFIX}_{cfg}/v{v}"
            clg = os.path.join(OUT, f"clg_task3_eo_{CFG_PREFIX}_{cfg}_v{v}.txt")
            py = os.path.join(OUT, f"pyquda_task3_eo_{CFG_PREFIX}_{cfg}_v{v}.txt")
            cmp_scalar(label, clg, py)

    print("=" * 88)
    for cfg in CFG_LIST:
        # Task 4 (gauge-only, run once per cfg; v=1 in CLG filenames)
        label = f"task4/action/{CFG_PREFIX}_{cfg}"
        clg = os.path.join(OUT, f"clg_task4_action_{CFG_PREFIX}_{cfg}_v1.txt")
        py = os.path.join(OUT, f"pyquda_task4_action_{CFG_PREFIX}_{cfg}.txt")
        cmp_scalar(label, clg, py, offset=TASK4_OFFSET)

        # Task 5 (gauge-only, run once per cfg; v=1 in CLG filenames)
        label = f"task5/plaq/{CFG_PREFIX}_{cfg}"
        clg = os.path.join(OUT, f"clg_task5_plaq_{CFG_PREFIX}_{cfg}_v1.txt")
        py = os.path.join(OUT, f"pyquda_task5_plaq_{CFG_PREFIX}_{cfg}.txt")
        cmp_scalar(label, clg, py)

        # Task 6 (combined fermion action, 3 flavors, run once per cfg)
        label = f"task6/fermion/{CFG_PREFIX}_{cfg}"
        clg = os.path.join(OUT, f"clg_task6_fermion_{CFG_PREFIX}_{cfg}_v1.txt")
        py = os.path.join(OUT, f"pyquda_task6_fermion_{CFG_PREFIX}_{cfg}.txt")
        cmp_scalar(label, clg, py)


if __name__ == "__main__":
    main()
