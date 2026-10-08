# physics-engine

A 2D rigid-body physics engine written from scratch in C++20, built to understand every line.
No third-party physics or math libraries. The engine is a pure library; a separate SDL3 demo app
draws shapes, contacts and AABBs as a visual debugger.

Planned stages: math, integration, rigid bodies, narrow-phase collision (SAT + clipping), impulse
response with friction, sequential-impulse solver with warm starting, broadphase, joints, sleeping
and CCD. A 3D version follows once the 2D engine is done.

Status: project skeleton only.

## Layout

```
engine/   static library "phys" (no graphics dependency)
demo/     SDL3 debug renderer and scenes
tests/    dependency-free test runner
```
