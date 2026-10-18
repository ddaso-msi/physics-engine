// Physics Puzzle Lab 3D: the puzzle game on the 3D engine.
//
// Every level is a physical problem. Some run all the time (aim and fire); the others have a setup
// phase, where the world is frozen and you arrange things, and a run phase, where physics decides.
//
//   left mouse  act (aim and fire, drag a piece, pull a weight)    right-drag  turn the camera   wheel  zoom
//   Space  run / back to setup    R  reset the level    Tab  pause    V  reset the camera
//   Q / E  turn a held piece      F  stand a held piece on end      Esc  cancel, then back to the menu
//   N  next level (after solving)  H  help   C  show contacts
//
// Checking a build without a display: arguments are run in order as a script, then the program exits.
//   --level N          open level N
//   --mouse X Y        put the cursor at a window pixel
//   --press X Y / --drag X Y / --release X Y / --click X Y   the left button, at window pixels
//   --orbit DYAW DPITCH   turn the camera (radians)       --zoom F   multiply the camera distance
//   --steps N          advance the simulation N fixed steps
//   --key NAME         run, reset, pause, help, contacts, next, escape, left, right, tip, view
//   --capture FILE     draw the current frame and save it as a BMP
//   --record FILE      from here on, --steps and --glide also write video (raw RGB24 at 60 frames a second)
//                      to FILE and a line per audible event to FILE.events
//   --glide X Y N      move the cursor to a window pixel over N steps
//   --at X Y Z         move the cursor onto a world point (as seen by the camera now)
//   --glide-at X Y Z N the same, over N steps
//   --down / --up      press / release the left button wherever the cursor is
//   --sweep DYAW DPITCH N   turn the camera that far over N steps
#include "hud.hpp"
#include "scene.hpp"

#include <SDL3/SDL.h>
#include <levels/levels.hpp>
#include <phys/timestep.hpp>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_set>

using namespace puzzle3d;
using namespace puzzle3d::app;

namespace {

constexpr float kBarHeight = 88.0f;
constexpr float kFooterHeight = 36.0f;
constexpr float kFooterTop = static_cast<float>(kHeight) - kFooterHeight;

enum class Action { None, Run, Reset, Pause, Help, Menu, Next, Contacts };

struct App {
    std::unique_ptr<Level> level;  // null on the menu
    int index = 0;
    demo3d::Camera camera;
    bool paused = false;
    bool show_help = false;
    bool show_contacts = false;
    bool running = true;
    bool pressed = false;   // left button went down over the scene and has not come up yet
    bool orbiting = false;  // right button held: turning the camera
    float mouse_x = 0, mouse_y = 0;
    phys::FixedTimestep timestep{kTimeStep};
};

// Video capture for the scripted run. Off unless --record was given.
struct Recorder {
    std::FILE* video = nullptr;
    std::FILE* events = nullptr;
    int tick = 0;
    std::unordered_set<std::uint64_t> touching;
    const Level* watched = nullptr;
    Status status = Status::Playing;
};
Recorder g_recorder;

void note(const char* what, double strength = 1.0) {
    if (g_recorder.events) std::fprintf(g_recorder.events, "%.4f %s %.4f\n", g_recorder.tick / 120.0, what, strength);
}

struct BarButton {
    Button button;
    Action action;
};
constexpr BarButton kButtons[] = {
    {{{452.0f, 10.0f, 176.0f, 36.0f}, "Space Run"}, Action::Run},
    {{{636.0f, 10.0f, 136.0f, 36.0f}, "R Reset"}, Action::Reset},
    {{{780.0f, 10.0f, 176.0f, 36.0f}, "Tab Pause"}, Action::Pause},
    {{{964.0f, 10.0f, 120.0f, 36.0f}, "H Help"}, Action::Help},
    {{{1092.0f, 10.0f, 168.0f, 36.0f}, "Esc Menu"}, Action::Menu},
};

SDL_FRect card_rect(int i) { return {140.0f, 250.0f + 78.0f * static_cast<float>(i), 1000.0f, 64.0f}; }
bool inside(const SDL_FRect& r, float x, float y) { return x >= r.x && x < r.x + r.w && y >= r.y && y < r.y + r.h; }
bool over_scene(const App& app) { return app.mouse_y > kBarHeight && app.mouse_y < kFooterTop; }

// The ray from the camera through a window pixel.
Ray pointer_ray(const demo3d::Camera& cam, float x, float y) {
    const Real f = cam.focal();
    const Real sx = (x - static_cast<float>(cam.width) * 0.5f) / f, sy = (static_cast<float>(cam.height) * 0.5f - y) / f;
    return {cam.eye(), (cam.forward() + cam.right() * sx + cam.up() * sy).normalized()};
}

void reset_view(App& app) {
    const View v = app.level->view();
    app.camera.target = v.target;
    app.camera.yaw = v.yaw;
    app.camera.pitch = v.pitch;
    app.camera.distance = v.distance;
}

void start_level(App& app, int index) {
    std::unique_ptr<Level> level = make_level(index);
    if (!level) return;
    app.level = std::move(level);
    app.index = index;
    app.paused = app.show_help = app.pressed = false;
    app.timestep = phys::FixedTimestep(kTimeStep);
    reset_view(app);
    note("select");
}

void leave_level(App& app) {
    app.level.reset();
    app.pressed = false;
}

void apply(App& app, Action action) {
    Level& level = *app.level;
    switch (action) {
        case Action::Run:
            if (!level.has_setup()) break;
            if (level.phase == Phase::Setup) level.start();
            else level.stop();
            note(level.phase == Phase::Running ? "run" : "stop");
            app.paused = app.pressed = false;
            app.timestep = phys::FixedTimestep(kTimeStep);
            break;
        case Action::Reset:
            level.reset();
            note("stop");
            app.paused = app.pressed = false;
            app.timestep = phys::FixedTimestep(kTimeStep);
            break;
        case Action::Pause: app.paused = !app.paused; break;
        case Action::Help: app.show_help = !app.show_help; break;
        case Action::Contacts: app.show_contacts = !app.show_contacts; break;
        case Action::Next:
            if (level.status != Status::Won) break;
            if (app.index + 1 < level_count()) start_level(app, app.index + 1);
            else leave_level(app);
            break;
        case Action::Menu:  // Esc backs out one thing at a time: the overlay, a held piece, then the level
            if (app.show_help) app.show_help = false;
            else if (level.cancel()) app.pressed = false;
            else leave_level(app);
            break;
        case Action::None: break;
    }
}

void press(App& app, float x, float y) {
    app.mouse_x = x;
    app.mouse_y = y;
    if (!app.level) {
        for (int i = 0; i < level_count(); ++i)
            if (inside(card_rect(i), x, y)) start_level(app, i);
        return;
    }
    if (app.show_help) {
        app.show_help = false;
        return;
    }
    for (const BarButton& b : kButtons)
        if (b.button.hit(x, y)) {
            apply(app, b.action);
            return;
        }
    if (over_scene(app) && !app.paused) {
        app.pressed = true;
        const size_t bodies = app.level->world.bodies.size();
        app.level->press(pointer_ray(app.camera, x, y));
        if (app.level->phase == Phase::Running && app.level->world.bodies.size() > bodies) note("fire");
        else if (app.level->phase == Phase::Setup) note("pick", 0.5);
    }
}

void move(App& app, float x, float y) {
    const float dx = x - app.mouse_x, dy = y - app.mouse_y;
    app.mouse_x = x;
    app.mouse_y = y;
    if (app.orbiting) app.camera.orbit(-dx * 0.008f, dy * 0.006f);
    if (app.level && app.pressed) app.level->drag(pointer_ray(app.camera, x, y));
}

void release(App& app, float x, float y) {
    app.mouse_x = x;
    app.mouse_y = y;
    if (app.level && app.pressed) {
        app.level->release(pointer_ray(app.camera, x, y));
        if (app.level->phase == Phase::Setup) note("place");
    }
    app.pressed = false;
}

void key(App& app, SDL_Keycode k) {
    if (!app.level) {
        if (k == SDLK_Q) app.running = false;
        if (k == SDLK_RETURN) start_level(app, 0);
        if (k >= SDLK_1 && k <= SDLK_9) start_level(app, static_cast<int>(k - SDLK_1));
        return;
    }
    switch (k) {
        case SDLK_SPACE: apply(app, Action::Run); break;
        case SDLK_R: apply(app, Action::Reset); break;
        case SDLK_TAB: apply(app, Action::Pause); break;
        case SDLK_H: apply(app, Action::Help); break;
        case SDLK_C: apply(app, Action::Contacts); break;
        case SDLK_N: apply(app, Action::Next); break;
        case SDLK_ESCAPE: apply(app, Action::Menu); break;
        case SDLK_V: reset_view(app); break;
        case SDLK_Q: app.level->rotate(1); break;
        case SDLK_E: app.level->rotate(-1); break;
        case SDLK_F: app.level->tip(); break;
        default: break;
    }
    if (app.pressed) app.level->drag(pointer_ray(app.camera, app.mouse_x, app.mouse_y));  // show the turned piece at once
}

void handle_event(App& app, const SDL_Event& e) {
    switch (e.type) {
        case SDL_EVENT_QUIT: app.running = false; break;
        case SDL_EVENT_KEY_DOWN:
            if (!e.key.repeat) key(app, e.key.key);
            break;
        case SDL_EVENT_MOUSE_MOTION: move(app, e.motion.x, e.motion.y); break;
        case SDL_EVENT_MOUSE_BUTTON_DOWN:  // one event per press: holding the button does not repeat
            if (e.button.button == SDL_BUTTON_LEFT) press(app, e.button.x, e.button.y);
            if (e.button.button == SDL_BUTTON_RIGHT) {
                app.mouse_x = e.button.x;
                app.mouse_y = e.button.y;
                app.orbiting = true;
            }
            break;
        case SDL_EVENT_MOUSE_BUTTON_UP:
            if (e.button.button == SDL_BUTTON_LEFT) release(app, e.button.x, e.button.y);
            if (e.button.button == SDL_BUTTON_RIGHT) app.orbiting = false;
            break;
        case SDL_EVENT_MOUSE_WHEEL:
            if (app.level) app.camera.scroll(e.wheel.y);
            break;
        default: break;
    }
}

void update(App& app, Real frame_seconds) {
    if (!app.level || app.paused || app.show_help) return;
    app.timestep.advance(frame_seconds, [&](Real dt) { app.level->step(dt); });
}

// ---- drawing ----

void draw_hud(SDL_Renderer* ren, const App& app) {
    const Level& level = *app.level;
    const float width = static_cast<float>(kWidth), cx = width * 0.5f;
    fill_rect(ren, {0, 0, width, kBarHeight}, color::kPanel);
    set_color(ren, color::kPanelEdge);
    SDL_RenderLine(ren, 0, kBarHeight, width, kBarHeight);

    text(ren, 20.0f, 20.0f, 2.0f, color::kText, "%d %s", app.index + 1, level.name());
    const bool setup = level.phase == Phase::Setup;
    const char* mode = app.paused ? "PAUSED" : level.status == Status::Won ? "SOLVED" : level.status == Status::Failed ? "FAILED" : setup ? "SETUP" : "RUNNING";
    const SDL_Color mode_color = app.paused ? color::kTarget : level.status == Status::Won ? color::kGood : level.status == Status::Failed ? color::kBad : setup ? color::kToolEdge : color::kText;
    const SDL_FRect pill{320.0f, 12.0f, 124.0f, 32.0f};
    outline_rect(ren, pill, mode_color);
    text_centered(ren, pill.x + pill.w * 0.5f, 20.0f, 2.0f, mode_color, mode);

    for (const BarButton& b : kButtons) {
        Button button = b.button;
        const bool enabled = b.action != Action::Run || level.has_setup();
        if (b.action == Action::Run && !setup && level.has_setup()) button.label = "Space Stop";
        if (b.action == Action::Pause && app.paused) button.label = "Tab Resume";
        draw_button(ren, button, button.hit(app.mouse_x, app.mouse_y), enabled);
    }

    text(ren, 20.0f, 60.0f, 2.0f, color::kDim, "%s", level.objective());
    const std::string stats = level.stats();
    text(ren, width - 20.0f - text_width(stats.c_str(), 2.0f), 60.0f, 2.0f, color::kText, "%s", stats.c_str());

    fill_rect(ren, {0, kFooterTop, width, kFooterHeight}, color::kPanel);
    text(ren, 20.0f, kFooterTop + 10.0f, 2.0f, color::kText, "%s", level.hint());
    char counters[64];
    SDL_snprintf(counters, sizeof(counters), "Try %d   %5.1f s", level.attempt, level.sim_time);
    text(ren, width - 20.0f - text_width(counters, 2.0f), kFooterTop + 10.0f, 2.0f, color::kDim, "%s", counters);

    if (app.show_help) {
        const std::span<const char* const> lines = level.help();
        const char* keys[] = {
            "Right-drag  turn the camera   Wheel  zoom   V  reset it",
            "Space  run, or go back to setup",
            "R      reset the level        Tab  pause or resume",
            "Q / E  turn a held piece      F    stand it on end",
            "N      next level, once solved",
            "C      show contact points",
            "Esc    cancel, close this, or back to the menu",
        };
        const float h = 110.0f + 30.0f * static_cast<float>(lines.size() + std::size(keys));
        const SDL_FRect r{cx - 480.0f, 110.0f, 960.0f, h};
        draw_panel(ren, r, color::kToolEdge);
        text_centered(ren, cx, r.y + 20.0f, 3.0f, color::kText, "HOW TO PLAY");
        float y = r.y + 66.0f;
        for (const char* l : lines) {
            text(ren, r.x + 40.0f, y, 2.0f, color::kText, "%s", l);
            y += 30.0f;
        }
        y += 20.0f;
        for (const char* l : keys) {
            text(ren, r.x + 40.0f, y, 2.0f, color::kDim, "%s", l);
            y += 30.0f;
        }
    } else if (level.status == Status::Won) {
        const SDL_FRect r{cx - 360.0f, 104.0f, 720.0f, 124.0f};
        draw_panel(ren, r, color::kGood);
        text_centered(ren, cx, r.y + 14.0f, 4.0f, color::kGood, "SOLVED!");
        const std::string line = level.result() + ", try " + std::to_string(level.attempt);
        text_centered(ren, cx, r.y + 60.0f, 2.0f, color::kText, line.c_str());
        const bool last = app.index + 1 >= level_count();
        text_centered(ren, cx, r.y + 92.0f, 2.0f, color::kDim,
                      level.has_setup() ? (last ? "N menu   Space back to setup   R reset" : "N next level   Space back to setup   R reset")
                                        : (last ? "N menu   R play again" : "N next level   R play again"));
    } else if (level.status == Status::Failed) {
        const SDL_FRect r{cx - 360.0f, 104.0f, 720.0f, 112.0f};
        draw_panel(ren, r, color::kBad);
        text_centered(ren, cx, r.y + 14.0f, 3.0f, color::kBad, "NOT SOLVED");
        text_centered(ren, cx, r.y + 50.0f, 2.0f, color::kText, level.hint());
        text_centered(ren, cx, r.y + 80.0f, 2.0f, color::kDim, level.has_setup() ? "Space back to setup   R reset" : "R try again");
    }
}

void draw_menu(SDL_Renderer* ren, const App& app) {
    const float cx = static_cast<float>(kWidth) * 0.5f;
    text_centered(ren, cx, 64.0f, 6.0f, color::kText, "PHYSICS PUZZLE LAB");
    text_centered(ren, cx, 124.0f, 4.0f, color::kToolEdge, "3D");
    text_centered(ren, cx, 172.0f, 2.0f, color::kDim, "Real physics, in three dimensions. Experiment until it works.");
    text_centered(ren, cx, 210.0f, 2.0f, color::kToolEdge, "Pick a level");
    for (int i = 0; i < level_count(); ++i) {
        const SDL_FRect r = card_rect(i);
        const bool hovered = inside(r, app.mouse_x, app.mouse_y);
        fill_rect(ren, r, hovered ? SDL_Color{48, 56, 78, 255} : SDL_Color{33, 38, 52, 255});
        outline_rect(ren, r, hovered ? color::kToolEdge : color::kPanelEdge);
        text(ren, r.x + 20.0f, r.y + 20.0f, 3.0f, color::kTarget, "%d", i + 1);
        text(ren, r.x + 70.0f, r.y + 12.0f, 2.0f, color::kText, "%s", level_entry(i).name);
        text(ren, r.x + 70.0f, r.y + 38.0f, 2.0f, color::kDim, "%s", level_entry(i).objective);
        text(ren, r.x + r.w - 90.0f, r.y + 24.0f, 2.0f, hovered ? color::kGood : color::kDim, "PLAY");
    }
    text_centered(ren, cx, 660.0f, 2.0f, color::kDim, "Click a level or press its number.   Q quit");
}

void draw(SDL_Renderer* ren, const App& app) {
    set_color(ren, color::kSky);
    SDL_RenderClear(ren);
    if (!app.level) {
        draw_menu(ren, app);
        return;
    }
    const bool aiming = over_scene(app) && !app.paused && !app.show_help && !app.orbiting;
    const Ray ray = pointer_ray(app.camera, app.mouse_x, app.mouse_y);
    draw_level(ren, app.camera, *app.level, aiming ? &ray : nullptr, app.show_contacts);
    draw_hud(ren, app);
}

// ---- recording ----

void draw_cursor(SDL_Renderer* ren, const App& app) {
    const SDL_FPoint p{app.mouse_x, app.mouse_y};
    fill_disc(ren, p, app.pressed ? 11.0f : 8.0f, {10, 12, 18, 255});
    fill_disc(ren, p, app.pressed ? 8.0f : 6.0f, app.pressed ? color::kTarget : color::kText);
}

bool tick(App& app, SDL_Renderer* ren) {
    update(app, kTimeStep * 1.0001f);
    Recorder& rec = g_recorder;
    if (!rec.video) return true;
    const Level* level = app.level.get();
    if (level != rec.watched) {
        rec.watched = level;
        rec.touching.clear();
        rec.status = Status::Playing;
    }
    if (level && level->phase == Phase::Running && !app.paused) {
        std::unordered_set<std::uint64_t> now;
        for (const phys3d::ContactPair& c : level->world.contacts()) {
            if (c.dormant || c.manifold.count == 0) continue;
            const std::uint64_t pair = (static_cast<std::uint64_t>(static_cast<std::uint32_t>(c.a)) << 32) | static_cast<std::uint32_t>(c.b);
            now.insert(pair);
            if (rec.touching.count(pair)) continue;
            double impulse = 0;
            for (int k = 0; k < c.manifold.count; ++k) impulse += static_cast<double>(c.manifold.points[k].normal_impulse);
            if (impulse > 0.02 && level->sim_time > 0.1) note("impact", impulse);
        }
        rec.touching = std::move(now);
    } else {
        rec.touching.clear();
    }
    if (level && level->status != rec.status) {
        rec.status = level->status;
        if (rec.status == Status::Won) note("win");
        if (rec.status == Status::Failed) note("fail");
    }
    if (++rec.tick % 2 != 0) return true;
    draw(ren, app);
    draw_cursor(ren, app);
    SDL_Surface* shot = SDL_RenderReadPixels(ren, nullptr);
    SDL_Surface* rgb = shot ? SDL_ConvertSurface(shot, SDL_PIXELFORMAT_RGB24) : nullptr;
    bool ok = rgb != nullptr;
    for (int y = 0; ok && y < rgb->h; ++y)
        ok = std::fwrite(static_cast<const Uint8*>(rgb->pixels) + static_cast<size_t>(y) * static_cast<size_t>(rgb->pitch), 3,
                         static_cast<size_t>(rgb->w), rec.video) == static_cast<size_t>(rgb->w);
    SDL_DestroySurface(rgb);
    SDL_DestroySurface(shot);
    return ok;
}

bool run_script(App& app, SDL_Renderer* ren, int argc, char** argv) {
    auto number = [&](int i) { return i < argc ? static_cast<float>(std::atof(argv[i])) : 0.0f; };
    auto is = [](const char* a, const char* b) { return std::strcmp(a, b) == 0; };
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (is(a, "--level")) {
            start_level(app, static_cast<int>(number(++i)) - 1);
        } else if (is(a, "--mouse") || is(a, "--press") || is(a, "--drag") || is(a, "--release") || is(a, "--click")) {
            const float x = number(i + 1), y = number(i + 2);
            i += 2;
            if (is(a, "--press") || is(a, "--click")) press(app, x, y);
            if (is(a, "--mouse") || is(a, "--drag")) move(app, x, y);
            if (is(a, "--release") || is(a, "--click")) release(app, x, y);
        } else if (is(a, "--orbit")) {
            app.camera.orbit(number(i + 1), number(i + 2));
            i += 2;
        } else if (is(a, "--zoom")) {
            app.camera.zoom(number(++i));
        } else if (is(a, "--steps")) {
            const int n = static_cast<int>(number(++i));
            for (int s = 0; s < n; ++s)
                if (!tick(app, ren)) return false;
        } else if (is(a, "--glide")) {
            const float tx = number(i + 1), ty = number(i + 2);
            const int n = std::max(1, static_cast<int>(number(i + 3)));
            i += 3;
            const float x0 = app.mouse_x, y0 = app.mouse_y;
            for (int s = 1; s <= n; ++s) {
                float t = static_cast<float>(s) / static_cast<float>(n);
                t = t * t * (3.0f - 2.0f * t);
                move(app, x0 + (tx - x0) * t, y0 + (ty - y0) * t);
                if (!tick(app, ren)) return false;
            }
        } else if (is(a, "--at") || is(a, "--glide-at")) {
            const demo3d::ScreenPoint to = app.camera.to_screen(app.camera.to_view({number(i + 1), number(i + 2), number(i + 3)}));
            const int n = is(a, "--at") ? 1 : std::max(1, static_cast<int>(number(i + 4)));
            i += is(a, "--at") ? 3 : 4;
            const float x0 = app.mouse_x, y0 = app.mouse_y;
            for (int s = 1; s <= n; ++s) {
                float t = static_cast<float>(s) / static_cast<float>(n);
                t = t * t * (3.0f - 2.0f * t);
                move(app, x0 + (to.x - x0) * t, y0 + (to.y - y0) * t);
                if (n > 1 && !tick(app, ren)) return false;
            }
        } else if (is(a, "--sweep")) {
            const Real dyaw = number(i + 1), dpitch = number(i + 2);
            const int n = std::max(1, static_cast<int>(number(i + 3)));
            i += 3;
            for (int s = 0; s < n; ++s) {
                app.camera.orbit(dyaw / static_cast<Real>(n), dpitch / static_cast<Real>(n));
                if (app.pressed) app.level->drag(pointer_ray(app.camera, app.mouse_x, app.mouse_y));
                if (!tick(app, ren)) return false;
            }
        } else if (is(a, "--down")) {
            press(app, app.mouse_x, app.mouse_y);
        } else if (is(a, "--up")) {
            release(app, app.mouse_x, app.mouse_y);
        } else if (is(a, "--record") && i + 1 < argc) {
            const std::string path = argv[++i];
            g_recorder.video = std::fopen(path.c_str(), "wb");
            g_recorder.events = std::fopen((path + ".events").c_str(), "w");
            if (!g_recorder.video || !g_recorder.events) {
                SDL_Log("cannot record to %s", path.c_str());
                return false;
            }
        } else if (is(a, "--key") && i + 1 < argc) {
            const char* k = argv[++i];
            key(app, is(k, "run") ? SDLK_SPACE : is(k, "reset") ? SDLK_R : is(k, "pause") ? SDLK_TAB : is(k, "help") ? SDLK_H
                     : is(k, "contacts") ? SDLK_C : is(k, "next") ? SDLK_N : is(k, "left") ? SDLK_Q : is(k, "right") ? SDLK_E
                     : is(k, "tip") ? SDLK_F : is(k, "view") ? SDLK_V : SDLK_ESCAPE);
        } else if (is(a, "--capture") && i + 1 < argc) {
            draw(ren, app);
            SDL_Surface* shot = SDL_RenderReadPixels(ren, nullptr);
            const bool ok = shot && SDL_SaveBMP(shot, argv[++i]);
            SDL_DestroySurface(shot);
            if (!ok) {
                SDL_Log("capture failed: %s", SDL_GetError());
                return false;
            }
        } else {
            SDL_Log("unknown argument: %s", a);
            return false;
        }
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        SDL_Log("SDL_Init failed: %s", SDL_GetError());
        return 1;
    }
    SDL_Window* window = nullptr;
    SDL_Renderer* renderer = nullptr;
    if (!SDL_CreateWindowAndRenderer("Physics Puzzle Lab 3D", kWidth, kHeight, 0, &window, &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    App app;
    app.camera.width = kWidth;
    app.camera.height = kHeight;
    int result = 0;
    if (argc > 1) {
        result = run_script(app, renderer, argc, argv) ? 0 : 1;
    } else {
        Uint64 last = SDL_GetTicksNS();
        while (app.running) {
            SDL_Event e;
            while (SDL_PollEvent(&e)) handle_event(app, e);
            const Uint64 now = SDL_GetTicksNS();
            update(app, static_cast<Real>(static_cast<double>(now - last) * 1e-9));
            last = now;
            draw(renderer, app);
            SDL_RenderPresent(renderer);
        }
    }

    app.level.reset();
    if (g_recorder.video) std::fclose(g_recorder.video);
    if (g_recorder.events) std::fclose(g_recorder.events);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return result;
}
