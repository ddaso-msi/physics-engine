#pragma once
// An orbit camera and perspective projection for the wireframe demo. No SDL in here, so it can be tested.
//
// The camera sits on a sphere around `target`: `yaw` turns it about the vertical, `pitch` raises it, and
// `distance` is the sphere's radius. It looks at the target with world +y up.
#include <phys3d/math.hpp>

#include <algorithm>
#include <cmath>

namespace demo3d {

using phys3d::Real;
using phys3d::Vec3;

struct ScreenPoint {
    float x = 0, y = 0;
};

struct Camera {
    Vec3 target{0, 1.5f, 0};
    Real yaw = 0.6f;       // radians about +y; 0 looks along -z from +z
    Real pitch = 0.35f;    // radians above the horizon
    Real distance = 14;
    Real fov_y = 0.9f;     // vertical field of view, radians
    Real near_plane = 0.1f;
    int width = 960, height = 640;

    Vec3 eye() const {
        return target + Vec3{std::cos(pitch) * std::sin(yaw), std::sin(pitch), std::cos(pitch) * std::cos(yaw)} * distance;
    }
    // Unit vectors of the camera's own frame, in world space.
    Vec3 forward() const { return (target - eye()).normalized(); }
    Vec3 right() const { return cross(forward(), Vec3{0, 1, 0}).normalized(); }
    Vec3 up() const { return cross(right(), forward()); }

    // World point -> camera space: x to the right, y up, z = distance in front of the camera.
    Vec3 to_view(Vec3 world) const {
        const Vec3 d = world - eye();
        return {dot(d, right()), dot(d, up()), dot(d, forward())};
    }

    // Pixels per unit of (x/z): chosen so the vertical field of view fills the window height.
    Real focal() const { return static_cast<Real>(height) * 0.5f / std::tan(fov_y * 0.5f); }

    // Camera-space point (z > 0) -> pixel. Screen y grows downward.
    ScreenPoint to_screen(Vec3 view) const {
        const Real f = focal();
        return {static_cast<float>(width) * 0.5f + view.x / view.z * f, static_cast<float>(height) * 0.5f - view.y / view.z * f};
    }

    // Projects the world-space segment a-b, cutting off whatever lies behind the near plane. Returns false
    // if nothing of it is in front of the camera.
    bool project_segment(Vec3 a, Vec3 b, ScreenPoint& out_a, ScreenPoint& out_b, Real* depth = nullptr) const {
        Vec3 va = to_view(a), vb = to_view(b);
        if (va.z < near_plane && vb.z < near_plane) return false;
        if (va.z < near_plane) va = va + (vb - va) * ((near_plane - va.z) / (vb.z - va.z));
        if (vb.z < near_plane) vb = vb + (va - vb) * ((near_plane - vb.z) / (va.z - vb.z));
        out_a = to_screen(va);
        out_b = to_screen(vb);
        if (depth) *depth = (va.z + vb.z) * 0.5f;
        return true;
    }

    void orbit(Real d_yaw, Real d_pitch) {
        yaw += d_yaw;
        pitch = std::clamp(pitch + d_pitch, static_cast<Real>(-1.5), static_cast<Real>(1.5));  // stop short of the poles
    }
    void zoom(Real factor) { distance = std::clamp(distance * factor, static_cast<Real>(2), static_cast<Real>(80)); }
};

}  // namespace demo3d
