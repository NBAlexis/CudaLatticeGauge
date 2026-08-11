#!/bin/bash
# runtest.sh TestName GpuCount -- multi-GPU-improve1.md 3.10 test driver.
#
# Runs exactly one _TEST_MULTIGPU-registered test under mpiexec:
#   1. rankCount and GpuGrid come from `CLGTest --mg-config TestName`
#      (this script keeps no test-name table of its own);
#   2. GpuCount (devices per node, >= 1, <= rankCount) is validated here;
#      the worker re-verifies the visible device count with the CUDA API;
#   3. `mpiexec -n rankCount CLGTest TestName --mg-worker --gpu-grid ...
#      --device-per-node GpuCount` runs the test; its exit code is
#      propagated unchanged (no log polling, no YAML edits, no kills).
#
# Environment overrides: CLGTEST_BIN (default Bin/Ubuntu/CLGTest),
# MPIEXEC (default mpiexec). Paths containing spaces are supported.
set -u

if [ $# -ne 2 ]; then
    echo "usage: $0 TestName GpuCount" >&2
    exit 2
fi
TEST_NAME=$1
GPU_COUNT=$2

case "$GPU_COUNT" in
    ''|*[!0-9]*)
        echo "runtest: GpuCount must be a positive integer, got '$GPU_COUNT'" >&2
        exit 2
        ;;
esac
if [ "$GPU_COUNT" -lt 1 ]; then
    echo "runtest: GpuCount must be >= 1" >&2
    exit 2
fi

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
CLGTEST_BIN=${CLGTEST_BIN:-"$SCRIPT_DIR/Bin/Ubuntu/CLGTest"}
MPIEXEC=${MPIEXEC:-mpiexec}
# Open MPI refuses to launch rankCount > slot count (cores) unless told to
# oversubscribe; probe for Open MPI/OpenRTE before adding the flag (MPICH and
# other MPIEXEC overrides do not know it).
MPIEXEC_FLAGS=
if "$MPIEXEC" --version 2>/dev/null | grep -qiE 'open[ -]?mpi|openrte'; then
    MPIEXEC_FLAGS=--oversubscribe
fi
# Resolve a relative CLGTEST_BIN against the caller's cwd before we cd.
case "$CLGTEST_BIN" in
    /*) ;;
    *) CLGTEST_BIN="$(cd "$(dirname "$CLGTEST_BIN")" && pwd)/$(basename "$CLGTEST_BIN")" ;;
esac
if [ ! -x "$CLGTEST_BIN" ]; then
    echo "runtest: CLGTest binary not found or not executable: $CLGTEST_BIN" >&2
    exit 2
fi

# CLGTest loads the test YAML blocks via paths relative to Bin/Ubuntu.
cd "$SCRIPT_DIR/Bin/Ubuntu" || { echo "runtest: cannot cd $SCRIPT_DIR/Bin/Ubuntu" >&2; exit 2; }

CONFIG=$("$CLGTEST_BIN" --mg-config "$TEST_NAME")
if [ $? -ne 0 ]; then
    echo "runtest: --mg-config failed for $TEST_NAME (unknown test, not _TEST_MULTIGPU, or invalid metadata)" >&2
    exit 1
fi
# Expect exactly: rankCount gx gy gz gt
set -- $CONFIG
if [ $# -ne 5 ]; then
    echo "runtest: malformed --mg-config output for $TEST_NAME: '$CONFIG'" >&2
    exit 1
fi
RANK_COUNT=$1
GRID="$2,$3,$4,$5"
case "$RANK_COUNT" in
    ''|*[!0-9]*)
        echo "runtest: malformed rankCount in --mg-config output: '$CONFIG'" >&2
        exit 1
        ;;
esac
if [ "$GPU_COUNT" -gt "$RANK_COUNT" ]; then
    echo "runtest: GpuCount $GPU_COUNT > rankCount $RANK_COUNT for $TEST_NAME" >&2
    exit 2
fi

echo "runtest: $TEST_NAME ranks=$RANK_COUNT grid=[$GRID] devicePerNode=$GPU_COUNT" >&2
"$MPIEXEC" $MPIEXEC_FLAGS -n "$RANK_COUNT" "$CLGTEST_BIN" "$TEST_NAME" --mg-worker --gpu-grid "$GRID" --device-per-node "$GPU_COUNT"
exit $?
