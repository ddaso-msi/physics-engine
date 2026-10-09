#include "levels/levels.hpp"

#include "levels/knock_down.hpp"

namespace puzzle {

namespace {

template <class T>
std::unique_ptr<Level> make() {
    return std::make_unique<T>();
}

struct Registered {
    std::unique_ptr<Level> (*make)();
    LevelEntry entry;
};

const std::vector<Registered>& registry() {
    static const std::vector<Registered> levels = [] {
        std::vector<Registered> r;
        for (auto make_fn : {&make<KnockDown>}) {
            const std::unique_ptr<Level> sample = make_fn();
            r.push_back({make_fn, {sample->name(), sample->objective()}});
        }
        return r;
    }();
    return levels;
}

}  // namespace

int level_count() { return static_cast<int>(registry().size()); }
const LevelEntry& level_entry(int index) { return registry()[static_cast<size_t>(index)].entry; }

std::unique_ptr<Level> make_level(int index) {
    if (index < 0 || index >= level_count()) return nullptr;
    return registry()[static_cast<size_t>(index)].make();
}

}  // namespace puzzle
