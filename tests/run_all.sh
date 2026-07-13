#!/usr/bin/env bash
# Build and run every test with g++ directly (no CMake needed) and tally results.
set -u
cd "$(dirname "$0")/.."
CXX=${CXX:-g++}
FLAGS="-std=c++17 -O2 -pthread -Iinclude -Itests"
fail=0
for t in core particle rigidbody collision resolution stacking features harmonics sph leaf convex articulation contacts2 softbody constraint2 query robotics advanced; do
    if ! $CXX $FLAGS "tests/test_$t.cpp" -o "/tmp/pe_test_$t" 2>/tmp/pe_build_$t.log; then
        echo "BUILD FAILED: $t"; cat /tmp/pe_build_$t.log; fail=1; continue
    fi
    "/tmp/pe_test_$t" || fail=1
done
if [ $fail -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; fi
exit $fail
