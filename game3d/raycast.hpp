#pragma once
// Rays against bodies, for pointing at things with the mouse. The engine has no ray queries (it never
// needed one), so the game has its own, for each shape the engine knows.
#include <phys3d/world.hpp>

#include <optional>

namespace puzzle3d {

struct Ray {
    phys3d::Vec3 origin, dir;  // dir is a unit vector
    phys3d::Vec3 at(phys3d::Real t) const { return origin + dir * t; }
};

struct RayHit {
    int body = -1;
    phys3d::Real t = 0;   // distance along the ray
    phys3d::Vec3 point, normal;
};

// Where the ray first enters the body's shape (t >= 0), if it does. A ray starting inside reports t = 0.
std::optional<RayHit> raycast(const phys3d::Body& body, const Ray& ray);
// The nearest body the ray enters, skipping body `skip` (or none, if -1).
std::optional<RayHit> raycast(const phys3d::World& world, const Ray& ray, int skip = -1);
// Where the ray crosses the horizontal plane at height y, going forward.
std::optional<phys3d::Vec3> hit_plane_y(const Ray& ray, phys3d::Real y);

}  // namespace puzzle3d
