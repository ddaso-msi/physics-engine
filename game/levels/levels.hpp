#pragma once
// The list of levels, in menu order.
#include "level.hpp"

#include <memory>

namespace puzzle {

struct LevelEntry {
    const char* name;
    const char* objective;
};

int level_count();
const LevelEntry& level_entry(int index);
// A fresh instance of level `index` (0-based), ready to play. Null if the index is out of range.
std::unique_ptr<Level> make_level(int index);

}  // namespace puzzle
