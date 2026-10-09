# physics-engine

A 2D rigid-body physics engine written from scratch in C++20, built to understand every line.
No third-party physics or math libraries. The engine is a pure library; a separate SDL3 demo app
draws shapes, contacts and AABBs as a visual debugger.

Planned stages: math, integration, rigid bodies, narrow-phase collision (SAT + clipping), impulse
response with friction, sequential-impulse solver with warm starting, broadphase, joints, sleeping
and CCD. A 3D version follows once the 2D engine is done.

Status: stages 0-3 done (build setup, 2D math, integrators + fixed timestep, rigid body state + shapes).

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
./build/debug/demo/demo       # rigid bodies; needs SDL3 (brew install sdl3)
./build/debug/demo/demo_integrators   # stage 2 spring comparison
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

## Stage 3: rigid body state

`engine/include/phys/shapes.hpp`, `body.hpp`. A `Body` is a centre-of-mass position, an angle, linear
and angular velocity, force/torque accumulators and inverse mass/inertia. Inverse values of 0 mean
"infinite", which is how static bodies are represented: no special cases in the maths.

- **Shapes** live in body space with the centre of mass at the origin. `Polygon::from_points` accepts
  either winding, recentres on the centroid, computes outward normals, and rejects concave,
  collinear, degenerate or >8-vertex input.
- **Mass and inertia** from uniform density. Circle: `m = rho pi r^2`, `I = m r^2 / 2`. Polygon: fan
  of triangles from the centroid, each adding `rho * D/12 * (e1.e1 + e1.e2 + e2.e2)` with `D = e1 x e2`.
  Tests check box, regular n-gons (`I = m R^2/6 (1 + 2cos^2(pi/n))`) and an off-origin box.
- **Force at a point** adds torque `r x F`; **impulse at a point** changes velocity by `j/m` and
  angular velocity by `(r x j)/I`. A test confirms the impulse changes angular momentum about the
  world origin by exactly `p x j`.
- Gravity is an acceleration (mass-independent). Integration is semi-implicit Euler for both
  linear and angular state.

Demo: drag from a point on a body and release. Hit it through the centre and it only translates; hit
it off-centre and it spins. The long thin bar resists spinning (large `I`), the triangle is the easiest.
