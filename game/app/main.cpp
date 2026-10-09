// Physics Puzzle Lab: a puzzle game on top of the 2D engine.
//
// Every level is a physical problem. Some run all the time (aim and fire); the others have a setup
// phase, where the world is frozen and you arrange things, and a run phase, where physics decides.
//
//   Space  run / back to setup      R  reset the level       Tab  pause
//   Q / E or mouse wheel  rotate the piece you are holding   right click  put it back
//   N  next level (after solving)   H  help    C  show contacts    Esc  cancel, then back to the menu
//
// Checking a build without a display: arguments are run in order as a script, then the program exits.
//   --level N          open level N
//   --mouse X Y        put the cursor at a world position (metres)
//   --press X Y        press the left button there     --drag X Y   move with it held
//   --release X Y      let go there                    --click X Y  press and release
//   --steps N          advance the simulation N fixed steps
//   --key NAME         run, reset, pause, help, contacts, next, escape, left, right
//   --capture FILE     draw the current frame and save it as a BMP
#include "render.hpp"
#include "scene.hpp"

#include <SDL3/SDL.h>
#include <levels/levels.hpp>
#include <phys/timestep.hpp>

#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

using namespace puzzle;
using namespace puzzle::app;
using phys::Vec2;

namespace {

constexpr float kBarHeight = 88.0f;
constexpr float kFooterTop = static_cast<float>(kHeight) - kFloorTop * kPixelsPerMeter;

enum class Action { None, Run, Reset, Pause, Help, Menu, Next, Contacts };

struct App {
    std::unique_ptr<Level> level;  // null on the menu
    int index = 0;                 // which level `level` is
    bool paused = false;
    bool show_help = false;
    bool show_contacts = false;
    bool running = true;
    bool pressed = false;          // left button went down over the world and has not come up yet
    float mouse_x = 0, mouse_y = 0;
    phys::FixedTimestep timestep{kTimeStep};
};

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
bool over_world(const App& app) { return app.mouse_y > kBarHeight && app.mouse_y < kFooterTop; }
Vec2 pointer(const App& app) { return to_world(app.mouse_x, app.mouse_y); }

void start_level(App& app, int index) {
    std::unique_ptr<Level> level = make_level(index);
    if (!level) return;
    app.level = std::move(level);
    app.index = index;
    app.paused = app.show_help = app.pressed = false;
    app.timestep = phys::FixedTimestep(kTimeStep);
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
            app.paused = app.pressed = false;
            app.timestep = phys::FixedTimestep(kTimeStep);  // no leftover time from the previous run
            break;
        case Action::Reset:
            level.reset();
            app.paused = app.pressed = false;
            app.timestep = phys::FixedTimestep(kTimeStep);
            break;
        case Action::Pause: app.paused = !app.paused; break;
        case Action::Help: app.show_help = !app.show_help; break;
        case Action::Contacts: app.show_contacts = !app.show_contacts; break;
        case Action::Next:  // only once solved; after the last level, back to the menu
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
    if (app.show_help) {  // any click dismisses the overlay and does nothing else
        app.show_help = false;
        return;
    }
    for (const BarButton& b : kButtons)
        if (b.button.hit(x, y)) {
            apply(app, b.action);
            return;
        }
    if (over_world(app) && !app.paused) {
        app.pressed = true;
        app.level->press(to_world(x, y));
    }
}

void move(App& app, float x, float y) {
    app.mouse_x = x;
    app.mouse_y = y;
    if (app.level && app.pressed) app.level->drag(to_world(x, y));
}

void release(App& app, float x, float y) {
    app.mouse_x = x;
    app.mouse_y = y;
    if (app.level && app.pressed) app.level->release(to_world(x, y));
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
        case SDLK_Q: app.level->rotate(1); break;
        case SDLK_E: app.level->rotate(-1); break;
        default: break;
    }
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
            if (e.button.button == SDL_BUTTON_RIGHT && app.level && app.level->cancel()) app.pressed = false;
            break;
        case SDL_EVENT_MOUSE_BUTTON_UP:
            if (e.button.button == SDL_BUTTON_LEFT) release(app, e.button.x, e.button.y);
            break;
        case SDL_EVENT_MOUSE_WHEEL:
            if (app.level && e.wheel.y != 0) app.level->rotate(e.wheel.y > 0 ? 1 : -1);
            break;
        default: break;
    }
}

void update(App& app, phys::Real frame_seconds) {
    // Paused, on the menu or reading the help: bank no time, so play resumes exactly where it stopped.
    if (!app.level || app.paused || app.show_help) return;
    app.timestep.advance(frame_seconds, [&](phys::Real dt) { app.level->step(dt); });
}

// ---- drawing ----

void draw_panel(SDL_Renderer* ren, const SDL_FRect& r, SDL_Color edge) {
    fill_rect(ren, r, color::kPanel);
    outline_rect(ren, r, edge);
    outline_rect(ren, {r.x + 1, r.y + 1, r.w - 2, r.h - 2}, edge);
}

void draw_hud(SDL_Renderer* ren, const App& app) {
    const Level& level = *app.level;
    const float width = static_cast<float>(kWidth), cx = width * 0.5f;
    fill_rect(ren, {0, 0, width, kBarHeight}, color::kPanel);
    set_color(ren, color::kPanelEdge);
    SDL_RenderLine(ren, 0, kBarHeight, width, kBarHeight);

    text(ren, 20.0f, 20.0f, 2.0f, color::kText, "%d %s", app.index + 1, level.name());
    // What the game is doing right now, so setup, running, paused and finished never look alike.
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

    // The floor strip doubles as a footer: what to do now on the left, the run's counters on the right.
    text(ren, 48.0f, kFooterTop + 8.0f, 2.0f, color::kText, "%s", level.hint());
    char counters[64];
    SDL_snprintf(counters, sizeof(counters), "Try %d   %5.1f s", level.attempt, level.sim_time);
    text(ren, width - 48.0f - text_width(counters, 2.0f), kFooterTop + 8.0f, 2.0f, color::kDim, "%s", counters);

    if (app.show_help) {
        const std::span<const char* const> lines = level.help();
        const char* keys[] = {
            "Space  run, or go back to setup",
            "R      reset the level",
            "Tab    pause or resume",
            "Q / E  rotate a held piece (or the wheel)",
            "RMB    put a held piece back (right mouse button)",
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
    draw_grid(ren);
    const float cx = static_cast<float>(kWidth) * 0.5f;
    text_centered(ren, cx, 70.0f, 6.0f, color::kText, "PHYSICS PUZZLE LAB");
    text_centered(ren, cx, 140.0f, 2.0f, color::kDim, "Real physics. Experiment until it works.");
    text_centered(ren, cx, 200.0f, 2.0f, color::kToolEdge, "Pick a level");
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
    set_color(ren, color::kBackground);
    SDL_RenderClear(ren);
    if (!app.level) {
        draw_menu(ren, app);
        return;
    }
    const bool aiming = over_world(app) && !app.paused && !app.show_help;
    draw_level(ren, *app.level, pointer(app), aiming, app.show_contacts);
    draw_hud(ren, app);
}

// ---- scripted run (see the top of the file) ----

bool run_script(App& app, SDL_Renderer* ren, int argc, char** argv) {
    auto number = [&](int i) { return i < argc ? static_cast<float>(std::atof(argv[i])) : 0.0f; };
    auto is = [](const char* a, const char* b) { return std::strcmp(a, b) == 0; };
    for (int i = 1; i < argc; ++i) {
        const char* a = argv[i];
        if (is(a, "--level")) {
            start_level(app, static_cast<int>(number(++i)) - 1);
        } else if (is(a, "--mouse") || is(a, "--press") || is(a, "--drag") || is(a, "--release") || is(a, "--click")) {
            const SDL_FPoint p = to_screen({number(i + 1), number(i + 2)});
            i += 2;
            if (is(a, "--press") || is(a, "--click")) press(app, p.x, p.y);
            if (is(a, "--mouse") || is(a, "--drag")) move(app, p.x, p.y);
            if (is(a, "--release") || is(a, "--click")) release(app, p.x, p.y);
        } else if (is(a, "--steps")) {
            const int n = static_cast<int>(number(++i));
            for (int s = 0; s < n; ++s) update(app, kTimeStep * 1.0001f);
        } else if (is(a, "--key") && i + 1 < argc) {
            const char* k = argv[++i];
            key(app, is(k, "run") ? SDLK_SPACE : is(k, "reset") ? SDLK_R : is(k, "pause") ? SDLK_TAB : is(k, "help") ? SDLK_H
                     : is(k, "contacts") ? SDLK_C : is(k, "next") ? SDLK_N : is(k, "left") ? SDLK_Q : is(k, "right") ? SDLK_E : SDLK_ESCAPE);
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
    if (!SDL_CreateWindowAndRenderer("Physics Puzzle Lab", kWidth, kHeight, 0, &window, &renderer)) {
        SDL_Log("window creation failed: %s", SDL_GetError());
        SDL_Quit();
        return 1;
    }
    SDL_SetRenderDrawBlendMode(renderer, SDL_BLENDMODE_BLEND);

    App app;
    int result = 0;
    if (argc > 1) {
        result = run_script(app, renderer, argc, argv) ? 0 : 1;
    } else {
        Uint64 last = SDL_GetTicksNS();
        while (app.running) {
            SDL_Event e;
            while (SDL_PollEvent(&e)) handle_event(app, e);
            const Uint64 now = SDL_GetTicksNS();
            update(app, static_cast<phys::Real>(static_cast<double>(now - last) * 1e-9));
            last = now;
            draw(renderer, app);
            SDL_RenderPresent(renderer);
        }
    }

    app.level.reset();
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return result;
}
