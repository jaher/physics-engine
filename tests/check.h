// Minimal test harness: CHECK(cond) and CHECK_NEAR(a,b,tol) with a pass/fail tally.
#pragma once
#include <cstdio>
#include <cmath>
#include <string>

namespace test {
inline int& passes() { static int p = 0; return p; }
inline int& fails()  { static int f = 0; return f; }
inline void ok(bool cond, const char* expr, const char* file, int line) {
    if (cond) passes()++;
    else { fails()++; std::printf("  FAIL  %s:%d   %s\n", file, line, expr); }
}
inline int report(const char* name) {
    std::printf("%-28s  %d passed, %d failed\n", name, passes(), fails());
    return fails() == 0 ? 0 : 1;
}
inline bool near(double a, double b, double tol = 1e-6) { return std::fabs(a - b) <= tol; }
}
#define CHECK(cond)            test::ok((cond), #cond, __FILE__, __LINE__)
#define CHECK_NEAR(a,b,tol)    test::ok(test::near((a),(b),(tol)), #a " ~= " #b, __FILE__, __LINE__)
