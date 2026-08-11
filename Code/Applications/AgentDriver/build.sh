#!/bin/bash
set -e

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
CLG_ROOT="$SCRIPT_DIR/../../.."
CLG_CMAKE_DIR="$SCRIPT_DIR/../../CMake"
AGENT_BUILD_DIR="$SCRIPT_DIR/build"
BIN_DIR="$CLG_ROOT/Bin/Ubuntu"

BACKEND="${1:-CUDA}"
GPU_ARCH="${2:-86}"
PRECISION_RAW="${3:-single}"
PRECISION=$(printf '%s' "$PRECISION_RAW" | tr '[:upper:]' '[:lower:]')

case "$PRECISION" in
    single|float|fp32)
        PRECISION_CMAKE_ARGS=(-UCLG_DOUBLE)
        PRECISION_LABEL="single"
        ;;
    double|fp64)
        PRECISION_CMAKE_ARGS=(-DCLG_DOUBLE=1)
        PRECISION_LABEL="double"
        ;;
    *)
        echo "Usage: $0 [BACKEND] [GPU_ARCH] [single|double]" >&2
        echo "Unknown precision: $PRECISION_RAW" >&2
        exit 2
        ;;
esac

echo "Building CLGLib (backend=$BACKEND, arch=$GPU_ARCH, precision=$PRECISION_LABEL)..."
mkdir -p "$CLG_CMAKE_DIR/build"
cd "$CLG_CMAKE_DIR/build"
cmake .. -DCLG_BACKEND="$BACKEND" -DCLG_GPU_ARCH="$GPU_ARCH" "${PRECISION_CMAKE_ARGS[@]}"
make -j$(nproc)

echo "Building AgentDriver (precision=$PRECISION_LABEL)..."
mkdir -p "$AGENT_BUILD_DIR"
cd "$AGENT_BUILD_DIR"
cmake .. -DCLG_GPU_ARCH="$GPU_ARCH" "${PRECISION_CMAKE_ARGS[@]}"
make -j$(nproc)

echo "Done. Binary: $BIN_DIR/AgentDriver"
