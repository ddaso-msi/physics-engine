#include "test.hpp"

int main() {
    for (auto& c : test::registry()) {
        int before = test::failures();
        c.fn();
        std::printf("[%s] %s\n", test::failures() == before ? " ok " : "FAIL", c.name);
        std::fflush(stdout);  // so a crash inside a later test does not swallow this result
    }
    std::printf("%zu tests, %d failed checks\n", test::registry().size(), test::failures());
    return test::failures() == 0 ? 0 : 1;
}
