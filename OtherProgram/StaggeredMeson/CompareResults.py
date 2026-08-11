"""
CompareResults.py
=================
Compare correlators from propav2.py (data/testc24p31/) against reference
(data/correlators/c24p31/).

Checks:
  1. Combined files:  {basestr}_{channel}.npy  vs  reference
  2. Per-operator files:  correlationp2p_{type}_{j}.npy  vs reference combined
     (recombined from per-operator to verify consistency)
"""

import numpy as np
import os

basestr = "c24p31"
test_dir = f"data/testc24p31"
ref_dir = f"data/correlators/{basestr}"

# ── 1. Compare combined files ──
print("=" * 60)
print("Comparing combined files:  test vs reference")
print("=" * 60)

all_pass = True
for channel in range(20):
    test_path = f"{test_dir}/{basestr}_{channel}.npy"
    ref_path = f"{ref_dir}/{basestr}_{channel}.npy"

    if not os.path.exists(test_path):
        print(f"  [MISSING] {test_path}")
        all_pass = False
        continue
    if not os.path.exists(ref_path):
        print(f"  [MISSING] ref {ref_path}")
        all_pass = False
        continue

    test = np.load(test_path)
    ref = np.load(ref_path)

    if test.shape != ref.shape:
        print(f"  channel {channel}: SHAPE MISMATCH  test={test.shape}  ref={ref.shape}")
        all_pass = False
        continue

    max_diff = np.max(np.abs(test - ref))
    if max_diff == 0:
        print(f"  channel {channel}: EXACT MATCH  shape={test.shape}")
    else:
        rel_diff = np.max(np.abs(test - ref) / (np.abs(ref) + 1e-30))
        print(f"  channel {channel}: max_abs_diff={max_diff:.6e}  max_rel_diff={rel_diff:.6e}  shape={test.shape}")
        if max_diff > 1e-10:
            all_pass = False

if all_pass:
    print("\n>>> All combined files match. <<<")
else:
    print("\n>>> Some combined files differ. <<<")

# ── 2. Recombine per-operator files and compare ──
print("\n" + "=" * 60)
print("Recombining per-operator files and comparing to reference")
print("=" * 60)

# Need all_signs to know how many operators per channel
import sys
sys.path.insert(0, ".")
from MesonStructures import all_signs, all_deltas, delta_to_deltaidx
all_signs_lst = all_signs()

all_pass2 = True
for channel in range(20):
    ref_path = f"{ref_dir}/{basestr}_{channel}.npy"
    if not os.path.exists(ref_path):
        continue
    ref = np.load(ref_path)

    n_ops = len(all_signs_lst[channel])
    combined = None
    for j in range(n_ops):
        path = f"{test_dir}/correlationp2p_{channel}_{j}.npy"
        if not os.path.exists(path):
            print(f"  [MISSING] {path}")
            all_pass2 = False
            continue
        data = np.load(path)
        if combined is None:
            combined = data.copy()
        else:
            combined = combined + data

    if combined is None:
        continue

    combined_real = np.real(combined)
    if combined_real.shape != ref.shape:
        print(f"  channel {channel}: SHAPE MISMATCH  recombined={combined_real.shape}  ref={ref.shape}")
        all_pass2 = False
        continue

    max_diff = np.max(np.abs(combined_real - ref))
    if max_diff == 0:
        print(f"  channel {channel}: EXACT MATCH  ({n_ops} operators)")
    else:
        rel_diff = np.max(np.abs(combined_real - ref) / (np.abs(ref) + 1e-30))
        print(f"  channel {channel}: max_abs_diff={max_diff:.6e}  max_rel_diff={rel_diff:.6e}  ({n_ops} operators)")
        if max_diff > 1e-10:
            all_pass2 = False

if all_pass2:
    print("\n>>> Recombined per-operator files match reference. <<<")
else:
    print("\n>>> Some recombined files differ. <<<")
