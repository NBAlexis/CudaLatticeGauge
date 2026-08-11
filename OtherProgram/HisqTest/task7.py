#!/usr/bin/env python3
"""Task 7: Verify gauge momentum conventions with GaugeMomentumFactor=2.0.

Two verifications:
1. Random momentum: CLGLib (GaugeMomentumFactor=2.0) and QUDA gaussMomQuda
   generate momenta with similar norms (statistical equivalence).
2. Same field: CLGLib CalculateKinematicEnergy with GaugeMomentumFactor=2.0
   equals QUDA momAction raw value (0.5*tr(P^dagger P)) exactly,
   i.e. ratio = 1.0, no scaling needed.

Usage:
    cd /home/nbalexis/CLGLib/OtherProgram/HisqTest
    /home/nbalexis/pyqudaenv/bin/python3 task7.py
"""

import os
import sys
import numpy as np

os.environ["QUDA_RESOURCE_PATH"] = "/home/nbalexis/comparehisq/.cache"

from pyquda import init
from pyquda.field import LatticeInfo, LatticeMom, LatticeGauge
from pyquda.quda import momResidentQuda, gaussMomQuda, momActionQuda
from pyquda.dirac.gauge import GaugeDirac

sys.path.insert(0, "/home/nbalexis/CLGLib/OtherProgram/NpyConfiguration")
from Converter.NpyConfigurationToCon import transformToCon, transformFromCon

LATT_SIZE = [4, 4, 4, 12]
T_BOUNDARY = -1
ANISOTROPY = 1.0
SEED = 12345
GAUGE_MOM_FACTOR = 2.0

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "..", "..", "..", "comparehisq", "data", "outputs")
os.makedirs(OUT_DIR, exist_ok=True)

CLG_BIN = "/home/nbalexis/CLGLib/Bin/Ubuntu/HISQ"
CLG_YAML = os.path.join(HERE, "Task7Mom.yaml")


def reconstruct_mom_matrix(data10):
    """Reconstruct 3x3 anti-Hermitian traceless matrix from QUDA reconstruct=10."""
    b01_re, b01_im = data10[0], data10[1]
    b02_re, b02_im = data10[2], data10[3]
    b12_re, b12_im = data10[4], data10[5]
    a0, a1, a2 = data10[6], data10[7], data10[8]

    P = np.zeros((3, 3), dtype=np.complex128)
    P[0, 0] = 1j * a0
    P[1, 1] = 1j * a1
    P[2, 2] = 1j * a2
    P[0, 1] = b01_re + 1j * b01_im
    P[1, 0] = -b01_re + 1j * b01_im
    P[0, 2] = b02_re + 1j * b02_im
    P[2, 0] = -b02_re + 1j * b02_im
    P[1, 2] = b12_re + 1j * b12_im
    P[2, 1] = -b12_re + 1j * b12_im
    return P


def length_sq(P):
    """tr(P^dagger P) for a single 3x3 matrix."""
    return np.trace(P.conj().T @ P).real


def run_clglib_task7(generate_random, mom_file, out_file, out_mom_file=None):
    """Run CLGLib Task7Mom via subprocess."""
    import subprocess

    yaml_lines = []
    with open(CLG_YAML) as f:
        for line in f:
            yaml_lines.append(line)

    # Inject job parameters
    DUMMY_CFG = "/home/nbalexis/comparehisq/data/cfg16_12_double.con"
    job_lines = [
        f"    GenerateRandom : {1 if generate_random else 0}\n",
        f"    CfgFile : {DUMMY_CFG}\n",
        f"    MomFile : {mom_file}\n",
        f"    OutputFile : {out_file}\n",
    ]
    if out_mom_file:
        job_lines.append(f"    OutMomFile : {out_mom_file}\n")

    # Inject job parameters into the existing JobTask7Mom block
    new_lines = []
    in_job = False
    for line in yaml_lines:
        new_lines.append(line)
        if line.strip().startswith("JobTask7Mom"):
            in_job = True
        elif in_job and line.strip() == "":
            # After JobTask7Mom header line, insert parameters before first blank
            # Actually, insert right after JobTask7Mom line
            pass

    # Simpler approach: insert job params right after "JobTask7Mom:" line
    new_lines = []
    inserted = False
    for line in yaml_lines:
        new_lines.append(line)
        if not inserted and line.strip().startswith("JobTask7Mom"):
            new_lines.extend(job_lines)
            inserted = True

    # Append WorkJob at the end if not present
    has_workjob = any(line.strip().startswith("WorkJob") for line in yaml_lines)
    if not has_workjob:
        new_lines.append(f"\nWorkJob : EHJ_Task7Mom\n")
    temp_yaml = os.path.join(OUT_DIR, "task7_temp.yaml")
    with open(temp_yaml, "w") as f:
        f.writelines(new_lines)

    cmd = [CLG_BIN, temp_yaml, "EHJ_Task7Mom"]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.dirname(CLG_BIN))
    if result.returncode != 0:
        print(f"CLGLib stderr: {result.stderr}")
        raise RuntimeError(f"CLGLib failed with code {result.returncode}")

    with open(out_file) as f:
        return float(f.read().strip())


def main():
    init([1, 1, 1, 1], backend="numpy", resource_path=os.path.join(HERE, ".cache"))

    latt_info = LatticeInfo(LATT_SIZE, t_boundary=T_BOUNDARY, anisotropy=ANISOTROPY)
    vol = np.prod(LATT_SIZE)
    nlinks = vol * 4

    print("=" * 70)
    print("Task 7: Verify GaugeMomentumFactor=2.0 makes CLGLib/QUDA consistent")
    print("=" * 70)
    print(f"Lattice: {LATT_SIZE}, volume = {vol}, links = {nlinks}")
    print(f"GaugeMomentumFactor = {GAUGE_MOM_FACTOR}")

    # ========================================================================
    # Verification 1: Random momentum norm comparison
    # ========================================================================
    print("\n" + "=" * 70)
    print("Verification 1: Random momentum norm comparison")
    print("=" * 70)

    # 1a. Generate QUDA random momentum and compute LengthSq
    mom_quda = LatticeMom(latt_info)
    gauge = LatticeGauge(latt_info)
    gauge.setAntiPeriodicT()
    gd = GaugeDirac(latt_info)
    gd.loadGauge(gauge)

    momResidentQuda(mom_quda.data_ptrs, gd.gauge_param)
    gaussMomQuda(SEED, 1.0)

    # Read back
    param = gd.gauge_param
    param.make_resident_mom = 0
    param.return_result_mom = 1
    momResidentQuda(mom_quda.data_ptrs, param)

    lex = mom_quda.lexico()
    ndir, lt, lz, ly, lx, ncomp = lex.shape

    P_quda = np.zeros((ndir, lt, lz, ly, lx, 3, 3), dtype=np.complex128)
    for mu in range(ndir):
        for t in range(lt):
            for z in range(lz):
                for y in range(ly):
                    for x in range(lx):
                        P_quda[mu, t, z, y, x] = reconstruct_mom_matrix(lex[mu, t, z, y, x, :])

    quda_length_sq = 0.0
    for mu in range(ndir):
        for t in range(lt):
            for z in range(lz):
                for y in range(ly):
                    for x in range(lx):
                        quda_length_sq += length_sq(P_quda[mu, t, z, y, x])

    print(f"\n[QUDA random momentum]")
    print(f"  Total LengthSq (tr(P^dagger P)) = {quda_length_sq:.10f}")
    print(f"  Per-link average                = {quda_length_sq / nlinks:.10f}")

    # 1b. Generate CLGLib random momentum and compute LengthSq
    clg_out_mom = os.path.join(OUT_DIR, "task7_clg_random_mom.bin")
    clg_out_result = os.path.join(OUT_DIR, "task7_clg_random_result.txt")

    # Clean up old files
    for f in [clg_out_mom, clg_out_result]:
        if os.path.exists(f):
            os.remove(f)

    clg_kin = run_clglib_task7(
        generate_random=True,
        mom_file="",
        out_file=clg_out_result,
        out_mom_file=clg_out_mom,
    )
    # CLGLib: kin = LengthSq / factor, where LengthSq = sum_links(tr(P†P) - 8)
    # So: avg_tr(P†P) = 8 + (kin * factor) / nlinks
    clg_length_sq = clg_kin * GAUGE_MOM_FACTOR
    clg_avg_tr = 8.0 + clg_length_sq / nlinks

    print(f"\n[CLGLib random momentum (GaugeMomentumFactor={GAUGE_MOM_FACTOR})]")
    print(f"  CalculateKinematicEnergy        = {clg_kin:.10f}")
    print(f"  LengthSq (kin * factor)         = {clg_length_sq:.10f}")
    print(f"  Per-link avg tr(P†P)            = {clg_avg_tr:.10f}")

    # Compare per-link tr(P†P) averages
    # QUDA: quda_length_sq is already sum(tr(P†P))
    quda_avg_tr = quda_length_sq / nlinks
    avg_diff = abs(clg_avg_tr - quda_avg_tr) / quda_avg_tr if quda_avg_tr != 0 else 0.0
    print(f"\n[Comparison]")
    print(f"  QUDA per-link avg tr(P†P)       = {quda_avg_tr:.10f}")
    print(f"  CLGLib per-link avg tr(P†P)     = {clg_avg_tr:.10f}")
    print(f"  Relative difference             = {avg_diff:.3e}")

    # The two random number generators differ, so we expect statistical
    # agreement (same distribution), not exact numerical match.
    if avg_diff < 0.1:  # 10% tolerance for single sample
        print("  PASS: Norms are statistically consistent")
    else:
        print("  WARN: Norms differ significantly")

    # ========================================================================
    # Verification 2: Same field -> kinetic energy / momAction ratio = 1.0
    # ========================================================================
    print("\n" + "=" * 70)
    print("Verification 2: Same field energy comparison (ratio should be 1.0)")
    print("=" * 70)

    # 2a. QUDA: compute momAction on GPU and raw action in Python
    # Re-generate with same seed for consistency
    mom_quda2 = LatticeMom(latt_info)
    gd.gauge_param.make_resident_mom = 1
    gd.gauge_param.return_result_mom = 0
    momResidentQuda(mom_quda2.data_ptrs, gd.gauge_param)
    gaussMomQuda(SEED + 1, 1.0)  # Different seed to avoid correlation

    # Compute momAction BEFORE releasing resident momentum
    nullptr = np.empty((0, 0), "<c16")
    action_quda_gpu = momActionQuda(nullptr, gd.gauge_param)

    # Now read back the momentum
    gd.gauge_param.make_resident_mom = 0
    gd.gauge_param.return_result_mom = 1
    momResidentQuda(mom_quda2.data_ptrs, gd.gauge_param)
    gd.gauge_param.make_resident_mom = 1
    gd.gauge_param.return_result_mom = 0

    lex2 = mom_quda2.lexico()
    P_quda2 = np.zeros((ndir, lt, lz, ly, lx, 3, 3), dtype=np.complex128)
    for mu in range(ndir):
        for t in range(lt):
            for z in range(lz):
                for y in range(ly):
                    for x in range(lx):
                        P_quda2[mu, t, z, y, x] = reconstruct_mom_matrix(lex2[mu, t, z, y, x, :])

    quda_raw_action = 0.0
    for mu in range(ndir):
        for t in range(lt):
            for z in range(lz):
                for y in range(ly):
                    for x in range(lx):
                        P = P_quda2[mu, t, z, y, x]
                        # QUDA formula per link: 0.5*|P_ii.im|^2 + sum|P_ij|^2 - 4.0
                        s = 0.0
                        s += 0.5 * P[0, 0].imag ** 2
                        s += 0.5 * P[1, 1].imag ** 2
                        s += 0.5 * P[2, 2].imag ** 2
                        s += P[0, 1].real ** 2 + P[0, 1].imag ** 2
                        s += P[0, 2].real ** 2 + P[0, 2].imag ** 2
                        s += P[1, 2].real ** 2 + P[1, 2].imag ** 2
                        quda_raw_action += s  # No -4.0 offset

    print(f"\n[QUDA same-field computation]")
    print(f"  momAction (GPU, with -4 offset) = {action_quda_gpu:.10f}")
    print(f"  Raw action (no -4 offset)       = {quda_raw_action:.10f}")
    print(f"  Offset contribution (-4*nlinks) = {-4.0 * nlinks:.10f}")
    print(f"  GPU matches CPU?                = {np.isclose(action_quda_gpu, quda_raw_action - 4.0 * nlinks)}")

    # 2b. Save QUDA momentum in CLGLib format
    clg_mom_file = os.path.join(OUT_DIR, "task7_quda_mom_for_clg.bin")
    transformToCon(clg_mom_file, P_quda2, double=True)

    # 2c. Run CLGLib to compute kinematic energy
    clg_energy_file = os.path.join(OUT_DIR, "task7_clg_energy.txt")
    if os.path.exists(clg_energy_file):
        os.remove(clg_energy_file)

    clg_kin_same = run_clglib_task7(
        generate_random=False,
        mom_file=clg_mom_file,
        out_file=clg_energy_file,
    )

    print(f"\n[CLGLib same-field computation (GaugeMomentumFactor={GAUGE_MOM_FACTOR})]")
    print(f"  CalculateKinematicEnergy        = {clg_kin_same:.10f}")
    print(f"  Equivalent to LengthSq / factor = {clg_kin_same * GAUGE_MOM_FACTOR:.10f}")

    # 2d. Compare
    # CLGLib CalculateKinematicEnergy = (tr(P†P) - 8) / GaugeMomentumFactor
    # QUDA momAction = sum(0.5*tr(P†P) - 4.0) = sum(tr(P†P) - 8) / 2
    # With GaugeMomentumFactor=2.0, they are mathematically equivalent.
    ratio = clg_kin_same / action_quda_gpu if action_quda_gpu != 0 else 0.0
    print(f"\n[Comparison]")
    print(f"  CLGLib kin / QUDA momAction     = {ratio:.10f}")
    print(f"  Expected                        = 1.0000000000")
    print(f"  Difference from 1.0             = {abs(ratio - 1.0):.3e}")

    if np.isclose(ratio, 1.0, rtol=1e-12, atol=1e-12):
        print("  PASS: ratio is exactly 1.0 to machine precision")
    else:
        print(f"  FAIL: ratio deviates from 1.0 by {abs(ratio - 1.0):.3e}")

    # ========================================================================
    # Summary
    # ========================================================================
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"GaugeMomentumFactor = {GAUGE_MOM_FACTOR}")
    print(f"")
    print(f"Verification 1 (random momentum norm):")
    print(f"  QUDA per-link tr(P†P)    = {quda_avg_tr:.6f}")
    print(f"  CLGLib per-link tr(P†P)  = {clg_avg_tr:.6f}")
    print(f"  Relative diff            = {avg_diff:.3e}")
    print(f"")
    print(f"Verification 2 (same field energy/action):")
    print(f"  QUDA momAction   = {action_quda_gpu:.6f}")
    print(f"  CLGLib kinEnergy = {clg_kin_same:.6f}")
    print(f"  Ratio            = {ratio:.10f}")
    print(f"  (Expected: 1.0, diff = {abs(ratio - 1.0):.3e})")
    print("=" * 70)


if __name__ == "__main__":
    main()
