# physics-engine

A 2D rigid-body physics engine written from scratch in C++20, built to understand every line.
No third-party physics or math libraries. The engine is a pure library; a separate SDL3 demo app
draws shapes, contacts and AABBs as a visual debugger.

Planned stages: math, integration, rigid bodies, narrow-phase collision (SAT + clipping), impulse
response with friction, sequential-impulse solver with warm starting, broadphase, joints, sleeping
and CCD. A 3D version follows once the 2D engine is done.

Status: stages 0-5 done (build setup, 2D math, integrators + fixed timestep, rigid body state + shapes, narrow-phase collision, contact response).

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
./build/debug/demo/demo               # latest stage (physics sandbox); needs SDL3 (brew install sdl3)
./build/debug/demo/demo_collision     # stage 4 contact manifold debug view
./build/debug/demo/demo_bodies        # stage 3 impulses at a point
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

## Stage 4: collision detection

`engine/include/phys/collision.hpp`. `collide(a, b, manifold)` handles every circle/polygon pairing.
A `Manifold` is a unit normal **pointing from A to B**, the minimum translation `depth`, and one or
two contact points, each with its own depth and placed halfway through the overlap.

- **Circle vs circle**: compare centre distance to the radius sum. Concentric circles pick an
  arbitrary direction instead of dividing by zero.
- **Polygon vs circle**: work in the polygon's frame. Find the face the centre is furthest in front of;
  if that exceeds the radius there is a gap. Otherwise the Voronoi region of that face decides the
  feature: a vertex (normal = vertex to centre) or the face itself (normal = face normal). A centre
  inside the polygon exits through the nearest face.
- **Polygon vs polygon**: separating axis test over both polygons' face normals (enough in 2D). The
  face with the *least* penetration becomes the reference face. The incident polygon's most
  anti-parallel edge is clipped to the reference edge's side planes (Sutherland-Hodgman), and the
  points that end up behind the reference face are the contacts. A small bias keeps A's face as the
  reference when two are nearly tied, so contacts do not flicker between frames.

Tests: closed-form cases (face-to-face boxes give 2 points, a 45-degree diamond gives 1, clipped
overlap, every circle region), plus three property tests over thousands of random shape pairs:
the reported (normal, depth) is exactly the minimum translation that separates the pair; any point
sampled inside both bodies implies a reported collision; and results are invariant under rigid
motion of the whole scene.

Demo: drag bodies together and watch contact points, normals and depths. Wheel or Q/E rotates the
body under the cursor.

## Stage 5: collision response

`engine/include/phys/solver.hpp`, `world.hpp`. `World::step(dt)` runs: find contacts, apply gravity to
velocities, solve contacts, move bodies. Only the velocity change is solved; positions follow.

Each contact point is a velocity constraint solved by **sequential impulses**:

- **Normal**: the speed along the normal must reach a target `bias >= 0`. The target is the bounce,
  `-e * approach_speed` (computed once from the arrival speed, and only above a 1 m/s threshold so
  resting bodies don't chatter), or a push-out speed `0.2/dt * (depth - slop)` (Baumgarte
  stabilisation) when the bodies overlap.
- **Friction**: removes tangential speed, clamped to the Coulomb cone `|jt| <= mu * jn`.
- Each constraint is fixed in isolation by `lambda = -(v - bias) * effective_mass`, where the
  effective mass `1 / (1/mA + 1/mB + (r x n)^2/IA + (r x n)^2/IB)` accounts for spin. Fixing one
  contact disturbs the others, so the solver sweeps all of them (10 times by default).
- The **accumulated** impulse is clamped (push, never pull; friction inside the cone), not each
  sweep's increment, so a later sweep can undo an earlier overshoot.
- Materials mix per contact: restitution = max of the two, friction = sqrt of the product.

Tests use closed forms: elastic equal masses swap velocities; momentum is conserved and separation
speed is `e` times closing speed; a dropped ball leaves at `e` times its arrival speed and rises to
`e^2 h`; a block on a 20-degree slope accelerates at `g (sin t - mu cos t)` and holds when
`mu > tan t`; a sliding disc settles into rolling at exactly 2/3 of its initial speed with `w = -v/r`.
A tower of four crates stands, and a single solver sweep demonstrably does not hold it.

Known limitations (addressed later): resting stacks creep sideways by a millimetre or two per second
(Stage 6 warm starting, Stage 9 sleeping); contacts are rebuilt from scratch every step, so there is no
warm start yet; the broad phase is still all pairs.

Demo: click to drop shapes; `M` switches rubber / wood / ice; drop things on the ramps.
