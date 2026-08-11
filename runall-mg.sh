#!/bin/bash
# Run every _TEST_MULTIGPU test in Bin/Debug/TestSuit_MG.yaml on real GPUs
# via runtest.sh (mpiexec -n rankCount, GpuCount devices per node).
# usage: runall-mg.sh [GpuCount]   (default 2)
set -u
GPU_COUNT=${1:-2}
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
cd "$SCRIPT_DIR"

# WSL2: /usr/lib/wsl/lib libcuda shim is stale; use the live driver-store copy.
WSL_DRIVER_STORE=$(dirname "$(readlink -f "$SCRIPT_DIR/Libs/wsl-cuda/libcuda.so.1")")
export PATH="$SCRIPT_DIR/Libs/openmpi/bin:/usr/local/cuda/bin:$PATH"
export LD_LIBRARY_PATH="$SCRIPT_DIR/Libs/wsl-cuda:$WSL_DRIVER_STORE:$SCRIPT_DIR/Libs/openmpi/lib:/usr/local/cuda/lib64:${LD_LIBRARY_PATH:-}"
export MPIEXEC="$SCRIPT_DIR/Libs/openmpi/bin/mpiexec"

LOGDIR="$SCRIPT_DIR/Bin/Ubuntu/mg_logs_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$LOGDIR"
echo "logs: $LOGDIR  GpuCount=$GPU_COUNT"

# Extract test names and rank counts from the yaml (test name line starts at
# column 0; rank count follows inside the block).
mapfile -t TESTS < <(awk '
    /^[ \t]*#/ { next }
    /^TestMG[A-Za-z0-9_]*:/ { name=$1; sub(":$","",name) }
    /MultiGPUTestRankCount[ \t]*:[ \t]*[0-9]+/ { if ($3 ~ /^[0-9]+$/) print name, $3 }
' Bin/Debug/TestSuit_MG.yaml)

PASS=0; FAIL=0; FAILED=""
for entry in "${TESTS[@]}"; do
    set -- $entry
    NAME=$1; RANKS=$2
    if [ "$GPU_COUNT" -gt "$RANKS" ]; then GC=$RANKS; else GC=$GPU_COUNT; fi
    echo "=== $NAME (ranks=$RANKS, gpuPerNode=$GC) ==="
    if ./runtest.sh "$NAME" "$GC" > "$LOGDIR/$NAME.log" 2>&1; then
        ERR=$(grep -oE "errors:[ ]*[0-9]+" "$LOGDIR/$NAME.log" | tail -1)
        echo "PASS $NAME  ($ERR)"
        PASS=$((PASS+1))
    else
        echo "FAIL $NAME  (exit=$?)"
        FAIL=$((FAIL+1)); FAILED="$FAILED $NAME"
    fi
done

echo
echo "===== SUMMARY: PASS=$PASS FAIL=$FAIL ====="
[ -n "$FAILED" ] && echo "failed:$FAILED"
echo "logs: $LOGDIR"
