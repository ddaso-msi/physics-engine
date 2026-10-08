# physics-engine

A 2D rigid-body physics engine written from scratch in C++20, built to understand every line.
No third-party physics or math libraries. The engine is a pure library; a separate SDL3 demo app
draws shapes, contacts and AABBs as a visual debugger.

Planned stages: math, integration, rigid bodies, narrow-phase collision (SAT + clipping), impulse
response with friction, sequential-impulse solver with warm starting, broadphase, joints, sleeping
and CCD. A 3D version follows once the 2D engine is done.

Status: stages 0-2 done (build setup, 2D math, integrators + fixed timestep).

## Layout

```
engine/   static library "phys" (no graphics dependency)
demo/     SDL3 debug renderer and scenes
tests/    dependency-free test runner
```

## Build

```bash
cmake --preset debug          # Ninja, ASan + UBSan
cmake --build --preset debug
ctest --preset debug
./build/debug/demo/demo       # needs SDL3 (brew install sdl3)
```

## Stage 2: integration

`engine/include/phys/integrate.hpp`, `timestep.hpp`. Three ways to advance `x'' = a`, differing only
in update order:

| scheme | update | behaviour |
|---|---|---|
| explicit Euler | `x += v dt; v += a dt` | energy grows every step (x535 over 100 oscillator periods at dt=0.01) |
| semi-implicit Euler | `v += a dt; x += v dt` | symplectic: energy stays within +/-0.5%; first order |
| velocity Verlet | `x += v dt + a dt^2/2; v += (a0+a1) dt/2` | symplectic, second order, exact for constant `a` |

Semi-implicit Euler is what the rest of the engine will use (cheap, one force evaluation, stable).
`FixedTimestep` banks real frame time and spends it in constant `dt` chunks, with a cap against the
spiral of death.

Demo: three identical springs, one per scheme. Watch the red one fly off while the others stay put.
