#!/bin/bash
#=============================================================================
# check_raw_launch_typefail.sh
#
# Improve-1 I5 negative compile check (multi-GPU-improve1.md 3.3 gate), run
# by the opt-in CheckRawLaunchTypeFail cmake target -- NOT part of ALL.
#
# Compiles CLGTest/Tests/TestImprove5BadLaunch.cu (an intentionally
# type-mismatched _LAUNCH_KERNEL call) with the same defines/includes as the
# configured build, and EXPECTS the compile to fail: the raw <<<>>> backend
# must keep compile-time type checking of kernel arguments.
#
#   compile FAILS  -> OK (raw backend type-checks; this is a raw build)
#   compile PASSES -> either not a _CLG_LAUNCH_KERNEL=0 build, or the raw
#                     backend lost its type checking (a bug). Exit 1.
#
# Usage: check_raw_launch_typefail.sh <build-dir> <source-dir>
#=============================================================================
set -u

BUILD_DIR="$1"
SRC_DIR="$2"
FLAGS="$BUILD_DIR/CMakeFiles/CLGTest.dir/flags.make"
OUT="$BUILD_DIR/TestImprove5BadLaunch"

if [ ! -f "$FLAGS" ]; then
    echo "flags.make not found at $FLAGS -- configure the build first."
    exit 2
fi

DEFINES=$(grep '^CXX_DEFINES' "$FLAGS" | cut -d= -f2-)
INCLUDES=$(grep '^CXX_INCLUDES' "$FLAGS" | cut -d= -f2-)
#GPU arch flags come from the CUDA target (CLGLib); strip embedded quotes so
#plain word-splitting passes them to nvcc the way the Makefile does.
CUDA_FLAGS=$(grep '^CUDA_FLAGS' "$BUILD_DIR/CMakeFiles/CLGLib.dir/flags.make" | cut -d= -f2- | tr -d '"')
NVCC=${CUDACXX:-nvcc}

# shellcheck disable=SC2086
# -rdc=true matches the real build (and is what makes the extern __constant__
# arrays in CudaHelper.h legal); -forward-unknown-to-host-compiler likewise.
if "$NVCC" -forward-unknown-to-host-compiler -std=c++17 -rdc=true -x cu \
        -c "$SRC_DIR/CLGTest/Tests/TestImprove5BadLaunch.cu" \
        $DEFINES $INCLUDES $CUDA_FLAGS -o "$OUT.o" 2> "$OUT.err"; then
    echo "UNEXPECTED: the bad launch COMPILED. On a _CLG_LAUNCH_KERNEL=0 (raw)"
    echo "build this means the <<<>>> backend lost its type checking -- a bug."
    echo "(On a _CLG_LAUNCH_KERNEL=1 build type erasure makes this compile; the"
    echo "check is only meaningful on a raw build.)"
    exit 1
fi

echo "OK: bad launch rejected at compile time. First diagnostics from $OUT.err:"
head -5 "$OUT.err"
exit 0
