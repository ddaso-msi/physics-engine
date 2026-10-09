// Stage 5 demo: collision response sandbox.
//
// Click to drop a shape at the cursor. Things to try:
//   - drop a rubber ball and a wooden one: bounce height scales with restitution (e^2 per bounce)
//   - drop a crate on a ramp, then switch to ice: friction decides if it slides or holds
//   - stack crates: the solver's repeated sweeps are what keep the tower up
//
//   1 box   2 circle   3 hexagon   4 triangle     M cycle material (rubber / wood / ice)
//   C show contacts    X clear dynamic bodies     R reset            Esc quit
#include <SDL3/SDL.h>
#include <phys/timestep.hpp>
#include <phys/world.hpp>

#include <cmath>
#include <vector>

using namespace phys;

namespace {

constexpr int kWidth = 960, kHeight = 640;
constexpr float kPixelsPerMeter = 60.0f;
constexpr float kWorldW = static_cast<float>(kWidth) / kPixelsPerMeter;
constexpr float kWorldH = static_cast<float>(kHeight) / kPixelsPerMeter;
constexpr size_t kMaxBodies = 160;  // the narrow phase is O(n^2) until stage 7

SDL_FPoint to_screen(Vec2 w) { return {w.x * kPixelsPerMeter, (kWorldH - w.y) * kPixelsPerMeter}; }
Vec2 to_world(float sx, float sy) { return {sx / kPixelsPerMeter, kWorldH - sy / kPixelsPerMeter}; }

struct Material {
    const char* name;
    Real restitution, friction;
};
constexpr Material kMaterials[] = {{"rubber", 0.8f, 0.6f}, {"wood", 0.2f, 0.5f}, {"ice", 0.0f, 0.03f}};

enum class Kind { Box, Circle, Hexagon, Triangle };
const char* kind_name(Kind k) {
    switch (k) {
        case Kind::Box: return "box";
        case Kind::Circle: return "circle";
        case Kind::Hexagon: return "hexagon";
        default: return "triangle";
    }
}

Body make_body(Kind kind, Vec2 pos, Real angle, const Material& m) {
    Shape shape;
    switch (kind) {
        case Kind::Box: shape = Shape::make_polygon(Polygon::box(0.5f, 0.5f)); break;
        case Kind::Circle: shape = Shape::make_circle(0.5f); break;
        case Kind::Hexagon: {
            std::vector<Vec2> pts;
            for (int k = 0; k < 6; ++k) {
                float a = 2.0f * kPi * static_cast<float>(k) / 6.0f;
                pts.push_back({0.6f * std::cos(a), 0.6f * std::sin(a)});
            }
            shape = Shape::make_polygon(*Polygon::from_points(pts));
            break;
        }
        case Kind::Triangle: {
            const Vec2 tri[3] = {{0, 0}, {1.3f, 0}, {0.4f, 1.1f}};
            shape = Shape::make_polygon(*Polygon::from_points(tri));
            break;
        }
    }
    Body b(shape, pos, angle);
    b.restitution = m.restitution;
    b.friction = m.friction;
    return b;
}

Body make_static_box(Vec2 half, Vec2 pos, Real angle) {
    Body b(Shape::make_polygon(Polygon::box(half.x, half.y)), pos, angle, BodyType::Static);
    b.restitution = 0.0f;
    b.friction = 0.5f;
    return b;
}

void build_level(World& world, size_t& static_count) {
    world.bodies.clear();
    world.add(make_static_box({kWorldW * 0.5f, 0.25f}, {kWorldW * 0.5f, 0.25f}, 0));            // floor
    world.add(make_static_box({0.25f, kWorldH * 0.5f}, {0.25f, kWorldH * 0.5f}, 0));            // left wall
    world.add(make_static_box({0.25f, kWorldH * 0.5f}, {kWorldW - 0.25f, kWorldH * 0.5f}, 0));  // right wall
    world.add(make_static_box({2.8f, 0.15f}, {3.6f, 6.2f}, -0.3f));                             // ramp, slopes right
    world.add(make_static_box({2.8f, 0.15f}, {11.6f, 4.4f}, 0.3f));                             // ramp, slopes left
    static_count = world.bodies.size();

    for (int i = 0; i < 4; ++i)  // a tower of crates
        world.add(make_body(Kind::Box, {8.0f, 1.0f + static_cast<float>(i) * 1.001f}, 0, kMaterials[1]));
    world.add(make_body(Kind::Circle, {2.6f, 8.5f}, 0, kMaterials[0]));
}

void draw_body(SDL_Renderer* ren, const Body& b, SDL_Color color) {
    SDL_SetRenderDrawColor(ren, color.r, color.g, color.b, 255);
    if (b.shape.type == Shape::Type::Circle) {
        constexpr int kSegments = 32;
        SDL_FPoint pts[kSegments + 1];
        const float r = b.shape.circle.radius;
        for (int i = 0; i <= kSegments; ++i) {
            float a = 2.0f * kPi * static_cast<float>(i) / kSegments;
            pts[i] = to_screen(b.pos + Vec2{r * std::cos(a), r * std::sin(a)});
        }
        SDL_RenderLines(ren, pts, kSegments + 1);
        SDL_FPoint c = to_screen(b.pos), tip = to_screen(b.pos + b.q.x_axis() * r);  // spoke shows spin
        SDL_RenderLine(ren, c.x, c.y, tip.x, tip.y);
    } else {
        const Polygon& p = b.shape.polygon;
        SDL_FPoint pts[Polygon::kMaxVertices + 1];
        for (int i = 0; i <= p.count; ++i) pts[i] = to_screen(apply(b.transform(), p.vertices[i % p.count]));
        SDL_RenderLines(ren, pts, p.count + 1);
    }
}

}  // namespace

int main() {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    if (!SDL_CreateWindowAndRenderer("physics-engine: collision response", kWidth, kHeight, 0, &window,
                                     &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    World world;
    size_t static_count = 0;
    build_level(world, static_count);

    Kind kind = Kind::Box;
    size_t material = 1;
    bool show_contacts = true;
    int spawned = 0;

    FixedTimestep timestep(1.0f / 120.0f);
    Uint64 last = SDL_GetTicksNS();
    bool running = true;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = false;
            if (e.type == SDL_EVENT_KEY_DOWN) {
                switch (e.key.key) {
                    case SDLK_ESCAPE: running = false; break;
                    case SDLK_1: kind = Kind::Box; break;
                    case SDLK_2: kind = Kind::Circle; break;
                    case SDLK_3: kind = Kind::Hexagon; break;
                    case SDLK_4: kind = Kind::Triangle; break;
                    case SDLK_M: material = (material + 1) % (sizeof(kMaterials) / sizeof(kMaterials[0])); break;
                    case SDLK_C: show_contacts = !show_contacts; break;
                    case SDLK_X: world.bodies.resize(static_count); break;
                    case SDLK_R: build_level(world, static_count); break;
                    default: break;
                }
            }
            if (e.type == SDL_EVENT_MOUSE_BUTTON_DOWN && e.button.button == SDL_BUTTON_LEFT &&
                world.bodies.size() < kMaxBodies) {
                // A small angle offset per spawn so stacked drops don't land perfectly aligned.
                const Real tilt = 0.17f * static_cast<Real>(spawned++ % 5) - 0.3f;
                world.add(make_body(kind, to_world(e.button.x, e.button.y), tilt, kMaterials[material]));
            }
        }

        Uint64 now = SDL_GetTicksNS();
        Real frame = static_cast<Real>(static_cast<double>(now - last) * 1e-9);
        last = now;
        timestep.advance(frame, [&](Real dt) { world.step(dt); });

        SDL_SetRenderDrawColor(renderer, 20, 22, 28, 255);
        SDL_RenderClear(renderer);
        for (const Body& b : world.bodies)
            draw_body(renderer, b, b.type == BodyType::Static ? SDL_Color{110, 114, 130, 255} : SDL_Color{120, 200, 255, 255});

        if (show_contacts) {
            for (const ContactPair& c : world.contacts()) {
                for (int k = 0; k < c.manifold.count; ++k) {
                    const SDL_FPoint p = to_screen(c.manifold.points[k].point);
                    const SDL_FPoint tip = to_screen(c.manifold.points[k].point + c.manifold.normal * 0.4f);
                    SDL_SetRenderDrawColor(renderer, 110, 220, 140, 255);
                    SDL_RenderLine(renderer, p.x, p.y, tip.x, tip.y);
                    SDL_SetRenderDrawColor(renderer, 255, 210, 90, 255);
                    SDL_FRect dot{p.x - 2.5f, p.y - 2.5f, 5.0f, 5.0f};
                    SDL_RenderFillRect(renderer, &dot);
                }
            }
        }

        SDL_SetRenderDrawColor(renderer, 150, 154, 170, 255);
        SDL_RenderDebugText(renderer, 70.0f, 12.0f, "click: drop shape | 1-4 shape | M material | C contacts | X clear | R reset");
        const Material& m = kMaterials[material];
        SDL_RenderDebugTextFormat(renderer, 70.0f, 28.0f, "%s, %s (e=%.2f, mu=%.2f) | bodies %d | contacts %d",
                                  kind_name(kind), m.name, static_cast<double>(m.restitution),
                                  static_cast<double>(m.friction), static_cast<int>(world.bodies.size() - static_count),
                                  static_cast<int>(world.contacts().size()));
        SDL_RenderPresent(renderer);
    }

    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
