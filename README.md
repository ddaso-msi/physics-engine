# physics-engine

A 2D rigid-body physics engine written from scratch in C++20, built to understand every line.
No third-party physics or math libraries. The engine is a pure library; a separate SDL3 demo app
draws shapes, contacts and AABBs as a visual debugger.

Planned stages: math, integration, rigid bodies, narrow-phase collision (SAT + clipping), impulse
response with friction, sequential-impulse solver with warm starting, broadphase, joints, sleeping
and CCD. A 3D version follows once the 2D engine is done.

Status: stages 0-9 done (build setup, 2D math, integrators + fixed timestep, rigid body state + shapes, narrow-phase collision, contact response, warm starting + contact persistence, broad phase, joints, sleeping + continuous collision + scene gallery).

## Layout

```
engine/   static library "phys" (no graphics dependency)
demo/     SDL3 debug renderer and scenes
tests/    dependency-free test runner
scenes/   library of gallery scenes (Newton's cradle, rag doll, car) shared by the demo and the tests
bench/    broad phase / World::step timings (build the release preset)
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

Recording controls in `demo`: `H` hides or shows all on-screen text, `Backspace` reloads the current scene
from its initial state (solver, sleep, CCD and broad phase toggles are left as set), `Tab` pauses and resumes.

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

Demo: click to drop shapes; `M` switches rubber / wood / ice; drop things on the ramps.

## Stage 6: warm starting and contact persistence

Contacts are rebuilt from scratch every step, so by default each step's solver starts from zero
impulse and has to rediscover how hard every contact is pushing. A resting contact needs nearly the
same impulse every frame, so we carry it over instead.

- **Feature ids** (`ContactPoint::id`, `collision.cpp`): each contact point is tagged with the pair of
  geometric features that made it: the reference face, the incident edge, and a tag for whether the
  point is an original vertex or one created by clipping (and which plane clipped it). Circle contacts
  use the polygon's face or vertex index. The id is unchanged while the same features touch, so a
  crate sliding along a floor keeps its two ids, and a box tipping onto a corner keeps that corner's id.
- **Persistence** (`transfer_impulses`, `World::step`): after finding this step's contacts, each point
  looks up last step's point with the same body pair and the same id and inherits its accumulated
  normal and tangent impulse. Points with no match start from zero.
- **Warm starting** (`solve_contacts`): the inherited impulses are applied to the velocities before
  the sweeps (as a separate pass after all constraints are built, so no contact's bounce is measured
  from velocities another contact's old impulse already changed). The sweeps then only have to
  correct a small change, and the final impulses are written back for the next step.

Measured effect (crates, dt = 1/120, wood):

| scene | warm starting off | warm starting on |
|---|---|---|
| tower of 4, 1 sweep per step | collapses (top at 0.5 m, should be 3.5) | stands (3.49 m) |
| tower of 10, 4 sweeps | collapses | stands (9.46 m, about 1 cm wobble) |
| pyramid of 21, 4 sweeps | 13 cm of slump | 2 cm |
| pyramid of 55, 10 sweeps | n/a | stands, max displacement 4 cm, speed 0 |
| resting tower, sideways creep | about 1.5 mm/s, never stops | frozen after about 3 s |

Both settings are in `SolverSettings` (`warm_starting`, `iterations`). In the demo press `P` for a
28-crate pyramid, `W` to switch warm starting off, and `-`/`=` to change the number of sweeps.

Not done yet (Stage 9): bodies never sleep, so a settled pile is still simulated every step. Overlap
is still corrected by Baumgarte velocity bias, which adds a little energy; a separate position solve
could remove that.

## Stage 7: broad phase

The narrow phase is exact but costs real work per pair, and testing every pair is O(n^2). The broad
phase compares cheap AABBs (`aabb.hpp`) first and hands the narrow phase only the survivors. There are
three interchangeable implementations in `broadphase.hpp`, selected by `World::broadphase`:

- **Brute force**: every pair. The reference the others are tested against.
- **Sweep and prune**: bodies kept sorted by the left edge of their box on x. Body i can only overlap
  the following bodies whose left edge is no further right than i's right edge, so the scan stops
  early. The order persists between frames and insertion sort fixes it up in near-linear time (it falls
  back to a full sort after a teleport or a burst of new bodies).
- **Dynamic AABB tree** (`dynamic_tree.hpp`, the default): a balanced binary tree whose leaves are
  bodies' boxes padded by 10 cm and whose interior nodes enclose their children. Inserting walks down
  choosing the child that grows the total perimeter least; AVL-style rotations keep the height near
  log2(n) (1000 boxes inserted in sorted order give height 10, not 999). A body only touches the tree
  when its real box escapes its padding, so most steps cost one containment test per body. Each moving
  body then queries the tree for overlaps.

Two things hold for all three, and are tested: **pairs come out sorted by (a, b)**, so the narrow phase
sees the same pairs in the same order whichever runs, and a stepped scene is bit-for-bit identical
across broad phases; and **pairs of two immovable bodies are never reported**.

**Collision filtering**: each body has a `category` and a `mask`; two bodies collide only if each one's
mask admits the other's category (`should_collide`). The filter runs after the broad phase, before the
narrow phase. In the demo, ghost mode (`G`) makes shapes pass through each other but still hit the floor.

Benchmark (`phys_bench`, release build, Apple clang, ms per frame; boxes jitter a little every frame):

| bodies | brute force | sweep and prune | tree |
|---|---|---|---|
| 1000 scattered | 2.03 | 0.12 | 0.25 |
| 5000 scattered | 23.6 | 0.60 | 3.30 |
| 10000 scattered | (skipped) | 1.55 | 10.2 |
| 5000 in a tall column | 8.24 | 24.4 | 0.42 |
| 10000 in a tall column | (skipped) | 121.8 | 0.94 |

What this says, honestly: both beat brute force by a wide margin. Brute force quadruples its work
when bodies double; sweep and prune and the tree need about 2.3x-2.8x as many comparisons per doubling
(work growing like n^1.2 to n^1.5). On evenly scattered bodies sweep and prune is the fastest, 2x-7x
quicker than the tree depending on size. But sweep and prune sorts on
one axis, so when many boxes share an x range (a tall stack, a wide floor) it compares nearly every
pair and ends up slower than brute force: 130x slower than the tree on a 10000-box column. The tree
has no such pathological case, which is why it is the default. In a dense pile it is not faster than
the others: a settled pile of 2000 crates in a closed room (the whole `World::step`, release build) costs
3.24 ms with brute force, 2.22 ms with sweep and prune and 3.18 ms with the tree, so below a few thousand
dense bodies the broad phase choice hardly matters and sweep and prune is the cheapest. The tree's query
cost per body also grows faster than log n (about 60 node visits per body at 2000 bodies, 160 at 10000),
so there is room to improve its quality (for example a periodic bulk rebuild); not done.

*Correction:* an earlier version of this section quoted whole-step times (brute force 3.7 ms, sweep and
prune 0.33 ms, tree 1.5 ms at 2000 crates) from a benchmark scene with a floor but no walls. Crates knocked
sideways fell out of the world and were never colliding again, so most of what it timed was free fall.
The scene now has walls and starts the crates on a non-overlapping grid, and the numbers above replace
those. The broad-phase-only tables are unaffected: they time bare boxes.

Demo: `F` rains 100 shapes in from the top, one every 7 steps into a free slot so nothing starts out
overlapping (it stops once the pile reaches the top: about 115 shapes in an empty room), `B` switches broad phase, `T` draws the tree's boxes
(or each body's tight box in the other modes), and the readout shows box tests, candidate pairs and
confirmed contacts for the last step.

## Stage 8: joints

`engine/include/phys/joints.hpp`. A joint takes away some of two bodies' relative freedom. Every
joint is turned into a few **scalar velocity rows** in one shared form,

    J v = lin_a . v_a + ang_a w_a + lin_b . v_b + ang_b w_b          (the rate of change of the constrained quantity)
    lambda = -(J v + bias) / (J M^-1 J^T),      bodies receive the impulse  M^-1 J^T lambda

and the solver neither knows nor cares which joint a row came from. Rows are solved in the same sweeps
as contacts (joints first), with accumulated-impulse clamping and warm starting. A row can be an
equality (rigid), one-sided (a limit: accumulated impulse >= 0), capped (a motor: +-max force x dt) or
soft (a spring-damper in the implicit-Euler form, stable however stiff). Index -1 means "the fixed
world", so a joint can pin a body to a point.

| joint | rows | extras |
|---|---|---|
| Distance | 1 along the line between anchors | rigid, or a spring (`frequency_hz`, `damping_ratio`) |
| Revolute | 2 (anchors coincide: one per world axis) | angle limit, motor with a torque cap |
| Prismatic | 2 (stay on the rail, no relative rotation) | translation limit, motor with a force cap |
| Mouse | 2 (soft, to a target point) | force cap; used for dragging |

Jointed bodies do not collide with each other unless `collide_connected` is set. `World::truncate`
drops joints that lose a body.

**Position error is not fed back as velocity.** The first version added a Baumgarte term to each row
(like contacts do) and a hinge chain with a weight hanging on it tore itself apart: even a 16:1 mass
ratio gave a 1.1 m gap at rest. A warm-started impulse carries the "push" that closed last frame's error
into the next frame, after the error is gone, so each frame overshoots a little more. The fix is
Box2D's: rows carry no position bias, and `solve_joint_positions` runs after the bodies move and
shifts positions and angles directly (never velocities, so it cannot add energy; at most 20 cm or
8 degrees per correction). That pass must be **Gauss-Seidel**: joints are corrected one at a time, each
from the poses the previous correction left. My first version built every joint's rows once and applied
all the fixes, and in a chain the shared link received both fixes and overshot: a settled rope bridge
that should hold to a millimetre had a 25 cm gap. Contacts still use Baumgarte for overlap; joints do not.

*Correction:* "one joint at a time" was not enough; it has to be one **row** at a time. The rows of a
single joint share its two bodies, and I was still correcting all of them from one measurement. With a
hinge at its limit that goes wrong: closing the pin turns the arm off the limit, turning it back opens
the pin, and the two overshoot against each other without settling. An arm dropped onto its lower limit
came to hang a few degrees above the stop with the pin pulled apart, and because it was "not at the
limit" nothing stopped its speed building up. I found it while writing the same joint for the 3D engine;
a sweep of 108 arms (three lengths, thicknesses, limits and timesteps) had 56 failing in 2D. The limit
test here had only checked how far the arm swung, never where it came to rest. Each row is now measured
again after the one before it; all 108 rest on the stop, pinned, and that sweep is a test.

Tests: every row's Jacobian is checked against a numerical derivative of the constraint it represents
(including a prismatic joint with the bodies slid off the anchor, where one term only shows up), plus
closed forms: a pinned disc swings with the physical-pendulum period (`T0 (1 + theta^2/16 + ...)`, within
1%), a spring oscillates at its frequency, damping ratio sets the decay and critical damping does not
overshoot, a rail box slides at `g sin(theta)` and stops at its limits, motors reach their speed or
accelerate at `tau / I`, a mouse joint's force cap holds. Position correction is tested not to add energy,
to ease large errors back rather than teleport, and not to make a chain worse.

Measured (dt = 1/120, 10 sweeps, 4 position passes, 8-link chain of 0.05 kg links, worst hinge gap):

| weight on the end vs a link | hanging at rest | swung at 3 m/s |
|---|---|---|
| 16:1 | 2 mm | 3 mm |
| 79:1 | 10 mm | 17 mm |
| 314:1 | 40 mm | 64 mm |
| 1571:1 | 320 mm | torn apart |

Warm starting is essential for hanging loads: at 79:1 with it off the same chain stretches 7 m. A rope
bridge pinned at both ends and dragged 2.2 m down by a mouse joint keeps every hinge within 5 mm.

Limits of this solver (not fixed): point constraints are two independent scalar rows, not a 2x2 block;
sequential impulses move information one joint per sweep, so a very heavy swinging end needs more sweeps
(Box2D's own advice is to keep mass ratios modest); angle limits are only meaningful inside (-pi, pi]; there is
no joint breaking.

Demo: `J` builds a joint level (rope bridge, rope pendulum, spring crate, motorised paddle, tilted rail
with end stops). Right-drag grabs any dynamic body with a mouse joint.

## Stage 9: sleeping, continuous collision, gallery

**Sleeping** (`World::allow_sleep`, `time_to_sleep` 0.5 s, thresholds 0.01 m/s and 0.035 rad/s). Each step
builds *islands*, the connected groups of dynamic bodies linked by contacts or joints, with union-find;
static bodies link nothing, so a floor of separate piles is not one island. A body that is slower than
both thresholds accumulates sleep time; an island goes to sleep when its *most restless* member has been
still long enough, so a pile sleeps all at once or not at all. Sleeping bodies are not integrated or
solved, and a pair where neither body is active is not even collided. Their contacts are kept as *dormant*
(with their impulses), so waking a sleeping pyramid with a tap warm starts it instead of letting it
settle again from nothing: the crates below the poked one move less than 2 cm.

Waking: an awake body touching a sleeping island joins it and wakes it (this falls out of building the
islands from contacts that include the sleepers); setting a body's velocity or applying a force wakes
it; a mouse joint keeps its body awake; and a motor told to turn (nonzero speed and torque) wakes its
bodies. That last rule exists because commanding a motor touches no velocity and no force: my first parked
car simply ignored its throttle until I made commands wake the mechanism. Moving a body by hand needs
`World::wake`, since the world cannot see a teleport.

**Continuous collision** (`World::continuous`), against static bodies only. A body that moves, *or
turns*, more than half its inscribed radius in a step is swept along its path: the pose is tested at
samples no further apart than the inscribed radius, which a convex body at least that wide cannot step
over however thin the obstacle, and the first hit is refined by bisection. A body that was clear is
placed just inside the wall so the next step's contact solver responds to a real contact. A fast body
that starts *already in contact* may sink no deeper than a quarter of its inscribed radius per step.

Two things went wrong before it worked, and are now regression tests. First, a spinning box hit a thin
wall corner-first: the solver stopped the corner but turned most of the momentum into a 205 rad/s swing,
and the box's centre ended up on the far side of a 4 cm wall within one step; the sweep had skipped it
because it "started in contact". Hence the depth cap, and rotation counting as motion (a pinned rod
spinning at 40 rad/s sweeps 33 cm of tip per step without its centre moving). Second, the sweep must ask
the narrow phase's exact question, argument order included: at razor-thin overlaps `collide(a, b)` and
`collide(b, a)` can disagree through rounding. Stress test: 3000 random shots (discs, boxes, triangles, up
to 200 m/s, walls 1 to 10 cm thick, spin up to 30 rad/s): 64% tunnel without the pass, none with it.

*Correction:* that stress test checked only that nothing got through, and so missed a bug that stopped
bodies too well. A harsher run (20000 shots, plates down to 1:10, spin up to 80 rad/s) that also checked
energy and whether the body ever left the wall found one shot in ten hanging against the wall for ever,
still moving into it at full speed. The sweep took the body's turn during the step from its start and end
angles, and angles are stored wrapped into (-pi, pi]. An off-centre impact can leave a box spinning at
400 rad/s, over half a revolution per step, and then that difference is the short way round: the sweep
tested poses the box never passed through and put it back in the same pose every step. The turn is now
`w * dt`. After the fix: 100000 shots at each of 60, 120 and 240 Hz, none through, none trapped, none
gaining energy. I first blamed the depth cap and rewrote it (measuring depth by the deepest point's travel,
as the 3D engine does, and adding an impulse at the moment of impact); with the angle fixed, removing each
of those changes made no difference to any run, so they are not in the code.

**Gallery** (`scenes/`, demo key `L`): Newton's cradle (rigid ropes, elastic frictionless balls; lifting one
ball sends one out the far side to 99.6% of its height with the middle three staying put, and two balls
send two), a rag doll (ten bodies, nine limited hinges, no self-collision), a car (a prismatic slide with
end stops plus a soft spring per wheel, motors in the axles; 5.4 m/s at 12 rad/s, brakes to a stop in
0.6 s), and a shooting range for the continuous pass. The car's suspension carrier needs real mass: a
spring-damper joint's stiffness is frequency times the effective mass of the two bodies it joins, so a
feather-light carrier gave a feather-weak spring and the suspension bottomed out.

Sleeping, measured (release build, ms per step, a pile of crates dropped into a closed room and settled
for 20 s, tree broad phase): 250 crates 0.120 ms awake, 0.068 ms asleep (64 still awake); at 500 to 2000
crates it saved nothing, because one restless crate keeps its whole island awake and large piles keep
shifting. Whether a given big pile fully sleeps within 20 s also varied between builds, so do not count on
it: it is a saving for scenes that really come to rest.

Not done: continuous collision between two moving bodies (a fast bullet can still pass through a thin
*dynamic* plank); splitting a big island so its quiet half could sleep; angular-velocity clamping.
