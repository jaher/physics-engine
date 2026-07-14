#!/usr/bin/env bash
# Build and run every test with g++ directly (no CMake needed) and tally results.
set -u
cd "$(dirname "$0")/.."
CXX=${CXX:-g++}
FLAGS="-std=c++17 -O2 -pthread -Iinclude -Itests"
fail=0
for t in core particle rigidbody collision resolution stacking features harmonics sph leaf convex articulation contacts2 softbody constraint2 query robotics advanced ros compound em explosion; do
    if ! $CXX $FLAGS "tests/test_$t.cpp" -o "/tmp/pe_test_$t" 2>/tmp/pe_build_$t.log; then
        echo "BUILD FAILED: $t"; cat /tmp/pe_build_$t.log; fail=1; continue
    fi
    "/tmp/pe_test_$t" || fail=1
done
# GPU (CUDA) suite — only if nvcc is present
if command -v nvcc >/dev/null 2>&1; then
    if nvcc -O3 -std=c++17 -gencode arch=compute_90,code=compute_90 -Iinclude -Itests tests/test_gpu.cu -o /tmp/pe_test_gpu 2>/tmp/pe_build_gpu.log; then
        /tmp/pe_test_gpu || fail=1
    else echo "BUILD FAILED: gpu"; cat /tmp/pe_build_gpu.log; fail=1; fi
else echo "gpu                           (skipped — no nvcc/CUDA)"; fi
if [ $fail -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; fi
exit $fail
