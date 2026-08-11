#!/usr/bin/env python3
"""Run all CLGLib HISQ comparison tasks via the unified HISQ executable."""

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
# Data and outputs live in the comparehisq workspace
COMPAREHISQ = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), "comparehisq")
DATA = os.path.join(COMPAREHISQ, "data")
OUT = os.path.join(DATA, "outputs")
LOGS = os.path.join(HERE, "logs")

CLG_BIN = "/home/nbalexis/CLGLib/Bin/Ubuntu/HISQ"

CFG_LIST = [12, 14, 16, 18, 20]
CFG_PREFIX = "cfg16"
V_LIST = [1, 2, 3]

os.makedirs(OUT, exist_ok=True)
os.makedirs(LOGS, exist_ok=True)


def run_task(task_name, yaml_name, job_enum, extra_args_fn):
    yaml = os.path.join(HERE, yaml_name)
    print(f"\n========== {task_name} (CLGLib) ==========")
    failed = []
    for cfg_idx in CFG_LIST:
        for v_idx in V_LIST:
            stem, args = extra_args_fn(cfg_idx, v_idx)
            out_file = args[-1]
            log_path = os.path.join(LOGS, f"clg_{stem}.log")
            cmd = [CLG_BIN, yaml, job_enum] + args
            print(f"  [{stem}] running...", end="", flush=True)
            t0 = time.perf_counter()
            with open(log_path, "w") as logfp:
                proc = subprocess.run(cmd, cwd=os.path.dirname(CLG_BIN), stdout=logfp, stderr=subprocess.STDOUT)
            dt = time.perf_counter() - t0
            ok = (proc.returncode == 0) and os.path.isfile(out_file)
            status = "OK" if ok else f"FAIL(rc={proc.returncode})"
            print(f" {status} ({dt:.2f}s)")
            if not ok:
                failed.append(stem)
    return failed


def run_task1():
    failed = []
    for mode, enum, yaml in [("full", "EHJ_Task1Full", "Task1Full.yaml"),
                              ("eo", "EHJ_Task1EO", "Task1EO.yaml")]:
        def make_args(cfg, v):
            stem = f"task1_{mode}_{CFG_PREFIX}_{cfg}_v{v}"
            cfg_file = os.path.join(DATA, f"{CFG_PREFIX}_{cfg}_double.con")
            v_file = os.path.join(DATA, f"v{v}_double.con")
            out_file = os.path.join(OUT, f"clg_{stem}.con")
            return stem, [cfg_file, v_file, out_file]
        failed.extend(run_task(f"Task 1 {mode.upper()}", yaml, enum, make_args))
    return failed


def run_task2():
    failed = []
    for mode, enum, yaml in [("full", "EHJ_Task2Full", "Task2Full.yaml"),
                              ("eo", "EHJ_Task2EO", "Task2EO.yaml")]:
        def make_args(cfg, v):
            stem = f"task2_{mode}_{CFG_PREFIX}_{cfg}_v{v}"
            cfg_file = os.path.join(DATA, f"{CFG_PREFIX}_{cfg}_double.con")
            v_file = os.path.join(DATA, f"v{v}_double.con")
            out_file = os.path.join(OUT, f"clg_{stem}.con")
            return stem, [cfg_file, v_file, out_file]
        failed.extend(run_task(f"Task 2 {mode.upper()}", yaml, enum, make_args))
    return failed


def run_task3a():
    def make_args(cfg, v):
        stem = f"task3a_eo_{CFG_PREFIX}_{cfg}_v{v}"
        cfg_file = os.path.join(DATA, f"{CFG_PREFIX}_{cfg}_double.con")
        v_file = os.path.join(DATA, f"v{v}_double.con")
        out_file = os.path.join(OUT, f"clg_{stem}.txt")
        return stem, [cfg_file, v_file, out_file]
    return run_task("Task 3 Test A", "Task3TestA.yaml", "EHJ_Task3EO", make_args)


def run_task3b():
    def make_args(cfg, v):
        stem = f"task3b_eo_{CFG_PREFIX}_{cfg}_v{v}"
        cfg_file = os.path.join(DATA, f"{CFG_PREFIX}_{cfg}_double.con")
        v_file = os.path.join(DATA, f"v{v}_double.con")
        out_file = os.path.join(OUT, f"clg_{stem}.txt")
        return stem, [cfg_file, v_file, out_file]
    return run_task("Task 3 Test B", "Task3TestB.yaml", "EHJ_Task3EO", make_args)


def run_task4():
    def make_args(cfg, v):
        stem = f"task4_action_{CFG_PREFIX}_{cfg}_v{v}"
        cfg_file = os.path.join(DATA, f"{CFG_PREFIX}_{cfg}_double.con")
        out_file = os.path.join(OUT, f"clg_{stem}.txt")
        return stem, [cfg_file, out_file]
    return run_task("Task 4", "Task4.yaml", "EHJ_Task4Action", make_args)


def run_task5():
    def make_args(cfg, v):
        stem = f"task5_plaq_{CFG_PREFIX}_{cfg}_v{v}"
        cfg_file = os.path.join(DATA, f"{CFG_PREFIX}_{cfg}_double.con")
        out_file = os.path.join(OUT, f"clg_{stem}.txt")
        return stem, [cfg_file, out_file]
    return run_task("Task 5", "Task5.yaml", "EHJ_Task5Plaq", make_args)


def run_task6():
    def make_args(cfg, v):
        stem = f"task6_fermion_{CFG_PREFIX}_{cfg}_v1"
        cfg_file = os.path.join(DATA, f"{CFG_PREFIX}_{cfg}_double.con")
        v1_file = os.path.join(DATA, "v1_double.con")
        v2_file = os.path.join(DATA, "v2_double.con")
        v3_file = os.path.join(DATA, "v3_double.con")
        out_file = os.path.join(OUT, f"clg_{stem}.txt")
        return stem, [cfg_file, v1_file, v2_file, v3_file, out_file]
    return run_task("Task 6", "Task6.yaml", "EHJ_Task6FermionAction", make_args)


def main():
    t0 = time.perf_counter()
    failed = []
    failed.extend(run_task1())
    failed.extend(run_task2())
    failed.extend(run_task3a())
    failed.extend(run_task3b())
    failed.extend(run_task4())
    failed.extend(run_task5())
    failed.extend(run_task6())
    print(f"\nAll CLGLib runs done in {time.perf_counter() - t0:.1f}s.")
    if failed:
        print(f"FAILED: {failed}")
        sys.exit(1)
    print("All CLGLib runs passed.")


if __name__ == "__main__":
    main()
