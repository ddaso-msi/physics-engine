#include "level.hpp"

#include <scenes.hpp>

namespace puzzle {

void add_room(World& w) {
    using phys::scenes::add_static_box;
    const Real t = kFloorTop * 0.5f;
    add_static_box(w, {kWorldW * 0.5f, t}, {kWorldW * 0.5f, t});
    add_static_box(w, {t, kWorldH * 0.5f}, {t, kWorldH * 0.5f});
    add_static_box(w, {t, kWorldH * 0.5f}, {kWorldW - t, kWorldH * 0.5f});
    add_static_box(w, {kWorldW * 0.5f, t}, {kWorldW * 0.5f, kWorldH + t});
}

}  // namespace puzzle
