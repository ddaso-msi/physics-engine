#pragma once
// Tiny dependency-free test runner.
#include <cmath>
#include <cstdio>
#include <functional>
#include <vector>

namespace test {

struct Case { const char* name; std::function<void()> fn; };
inline std::vector<Case>& registry() { static std::vector<Case> r; return r; }
inline int& failures() { static int f = 0; return f; }

struct Register {
    Register(const char* name, std::function<void()> fn) { registry().push_back({name, std::move(fn)}); }
};

inline void check(bool ok, const char* expr, const char* file, int line) {
    if (!ok) { std::printf("    FAIL %s:%d: %s\n", file, line, expr); ++failures(); }
}
inline void check_near(double a, double b, double eps, const char* expr, const char* file, int line) {
    if (!(std::fabs(a - b) <= eps)) {
        std::printf("    FAIL %s:%d: %s (%.9g vs %.9g)\n", file, line, expr, a, b);
        ++failures();
    }
}

}  // namespace test

#define TEST(name) \
    static void test_##name(); \
    static test::Register reg_##name(#name, test_##name); \
    static void test_##name()

#define CHECK(expr) test::check((expr), #expr, __FILE__, __LINE__)
#define CHECK_NEAR(a, b, eps) test::check_near((a), (b), (eps), #a " ~ " #b, __FILE__, __LINE__)
