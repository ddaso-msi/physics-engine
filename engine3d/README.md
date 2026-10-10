# engine3d: the 3D engine (work in progress, branch `3d`)

A separate library (`phys3d`, namespace `phys3d`) next to the finished 2D engine (tag `2d-complete`).
The 2D code is not modified. `phys3d` borrows only what has no dimension: `Real` from `phys/math.hpp`
today, and later `timestep.hpp` and `solver_settings.hpp`.

Why separate rather than one engine templated over the dimension: in 3D the pieces change in kind, not just
in size. A rotation is a quaternion with an axis instead of an angle; inertia is a 3x3 tensor, not a number,
and angular momentum is not parallel to angular velocity; contacts need edge-edge cases and face clipping on
polygons in 3D. Keeping the two codebases apart keeps each one readable and keeps the 2D tests meaningful.

Build and test (same presets as the 2D engine; `ctest` runs both suites):

```bash
cmake --preset debug && cmake --build --preset debug && ctest --preset debug
./build/debug/tests3d/phys3d_tests
```

## Step 1: math (`math.hpp`, `inertia.hpp`)

- **Conventions**: right-handed, y up; a positive rotation about an axis is counter-clockwise looking from the
  axis tip toward the origin (x -> y -> z -> x). `Mat3` is stored by columns and multiplies column vectors.
  `Quat` is (x, y, z, w) with w the scalar part; `q` and `-q` are the same rotation.
- **`Vec3`**: arithmetic, `dot`, `cross`, `length`, a zero-safe `normalized`, indexing by axis.
- **`Mat3`**: product, transpose, determinant, inverse (each row of the inverse is a cross product of two
  columns over the determinant; singular gives zeros, not NaN), `diagonal`, `outer`, and `skew(a)`, the matrix
  that does a cross product (`skew(a) * b == cross(a, b)`).
- **`Quat`**: Hamilton product (`a * b` = "b first, then a"), `from_axis_angle`, `rotate` (no matrix needed),
  `to_mat3`, `from_mat3` (branches on the largest component so it never divides by something near zero),
  `rotation_angle`, and `integrate_orientation`: one Euler step of `dq/dt = 1/2 (omega, 0) q` followed by a
  renormalise, which keeps a spinning body's orientation a valid rotation indefinitely.
- **Inertia**: `sphere_inertia`, `box_inertia`, `parallel_axis` and `to_world_frame` (`R I R^T`, used for the
  inertia tensor and for its inverse, which is what the solver stores).

Tests (`tests3d/`, 23): the cross product against the right-hand rule and against five identities over 300 random
vector triples (Lagrange, BAC-CAB, Jacobi...); a hand-checked 3x3 inverse; matrix algebra over random matrices;
`skew`; 90-degree rotations about each axis; quaternion composition, inverse, double cover, and agreement with
the matrix, with the matrix round trip forced through all four branches; integration of constant angular velocity
against the closed form, and 20000 steps staying a valid rotation; and the inertia formulas checked against an
independent brute-force integral (box, sphere, and a box about its corner, which also tests the parallel axis
theorem), a cube's inertia being the same in every orientation, and angular momentum not being parallel to
angular velocity except along a principal axis. 17 deliberate breakages of the maths (a flipped cross-product
component, wrong quaternion product terms, a missing factor of 2, a transposed inverse, a missing
`transposed()` in the world-frame rotation...) are each caught.

## Step 2: rigid body state (`body.hpp`, `src/body.cpp`)

A `Body` is the 2D one extended: position of the centre of mass, a `Quat` orientation, linear velocity, an
angular velocity **vector** `w` (world frame), force and torque accumulators, inverse mass, and the inertia
tensor and its inverse stored in the body frame and rotated into the world frame on demand. Zero inverse mass
and inertia mean immovable. `apply_force_at` adds the torque `r x F`; `apply_impulse_at` changes `v` by `j/m`
and `w` by `I^-1 (r x j)`.

What has no 2D counterpart is the **gyroscopic term**. A free body conserves angular momentum `L = I_world w`,
and `I_world` turns with the body, so `w` has to change even with no torque: `dw/dt = I^-1 (torque - w x (I w))`.
`integrate_velocity` offers three treatments (`Gyroscopic`), measured on a tumbling brick over 10 s:

| mode | energy | angular momentum direction |
|---|---|---|
| `Off` (w constant) | exactly conserved | wanders by 0.98 rad: wrong physics |
| `Explicit` (term from the start of the step) | grows: x1.59 at dt 1/60, x1.11 at 1/240; a fast thin plate reaches 10^6 in 5 s | held to 0.006 rad |
| `Implicit` (one Newton step in the body frame; the default) | never rises; loses 29% at 1/60, 9% at 1/240; the fast plate keeps 99.98% | held to 0.011 rad |

So the default is stable but dissipative: a tumbling body slowly spins down, faster at coarse steps. That is
the usual trade (it is Catto's formulation) and is far better than gaining energy, but it is not conservation.

*Correction (found in step 7):* "never rises" was only true at moderate spin. The implicit step is a single
Newton iteration, and when a body turns several radians per step (hundreds of rad/s, as after a hard corner
impact) that iteration can overshoot badly: a stress test saw one step hand back 780 times the energy it was
given. A gyroscopic torque does no work, so `integrate_velocity` now compares the rotational energy before and
after and scales the result back if it rose. With that guard the statement holds at any spin (tested up to
900 rad/s).

Tests (16, 39 in total): mass and inertia of the solid shapes; a fixed body ignores everything, including a
velocity written into it; free fall; force at a point; an impulse changes linear momentum by `j` and angular
momentum about the origin by `hit x j` in any orientation; a torque on a turned body meets the world-frame
inertia; spin about a principal axis is steady in every mode; a sphere needs no gyroscopic term; the table
above; a symmetric top's body-frame `w` circles its axis at the closed-form rate `(I3 - I1)/I1 * w3`
(within 2% at dt 1/1000); and the tennis-racket effect: spin about the intermediate axis with a 1% wobble flips
right over (5 to -4.9 rad/s) while spin about the other two axes stays at 5.000. 14 deliberate breakages of the
body code are each caught; three survived at first (explicit term's sign, torque using the body-frame inertia,
fixed bodies moving) and got their own tests.

## Step 3: shapes and contact generation (`shapes.hpp`, `collision.hpp`, `src/collision.cpp`)

A body now has a `Shape` (sphere or box). `collide(a, b, manifold)` keeps the 2D contract: the normal points
from A to B, `depth` is how far B must move along it to separate, and each point sits halfway through the
overlap. A 3D manifold holds up to **four** points, because two faces meet in a polygon and three or more
points are needed to stop a box rocking.

- **Sphere vs sphere**: centre distance against the radius sum.
- **Box vs sphere**: in the box's frame, clamp the sphere's centre to the box to get the nearest point. How
  many coordinates were clamped says whether that is a face (1), an edge (2) or a corner (3). A centre inside
  the box leaves through the nearest face.
- **Box vs box**: the separating axis test over **15 axes**: the 3 face normals of each box and the 9 cross
  products of one edge direction from each. The 9 are not optional: two cubes turned so an edge of each faces
  the other are separated by no face normal at all for centre distances between 2.83 and 3.83, only by the
  axis perpendicular to both edges. The axis of least overlap decides the contact:
  - a **face** axis: the other box's most opposed face is clipped (Sutherland-Hodgman) against the four side
    planes of the reference face; the points left behind the reference face are the patch. Up to 8 points come
    out, reduced to 4: the deepest, the farthest from it, the one making the biggest triangle with those, and
    the farthest on the other side.
  - an **edge** axis: one point, midway between the closest points of the two edges.

Two tuning choices, each pinned by a test. A face of A is kept as the reference unless another axis is
*clearly* shallower (5% and 1 cm): equal cubes stacked have A's face and B's face tied to within rounding, and
the strict minimum would flip the reference from frame to frame. And the side planes are pushed out 2 mm before
clipping: with equal boxes the corners of one face lie exactly on the side planes of the other, and without the
margin each corner flickers between "this corner" and "a point cut by plane k", which are different features
with different ids. (I checked the 2D engine for the same flicker: in settled stacks 0 of 12000 contact pairs
changed ids between frames, so it is left alone.)

Tests (20, 59 in total): closed forms for every feature (sphere on a face, edge, corner, inside; a 4-point patch
whose corners and area are known; an overhang clipped at the edge; separation along each face axis; the edge-edge
case above, both separated and touching; a corner pressed into a face giving one point; an octagonal overlap
reduced to four points covering more than 2.0 of its 3.31 area) and properties over thousands of random pairs:
the reported normal and depth really separate the pair and nothing much smaller does; every contact point lies in
the overlap; any point sampled inside both bodies implies a collision; the result turns with the scene;
swapping the bodies flips the normal. Of 24 deliberate breakages, 23 are caught (no edge axes, an unnormalised
edge axis, each clipping mistake, the wrong incident face, a point off the midpoint, a naive reduction...).
Two survived at first and showed real gaps, which is how the reference-face and clip-margin tests came about.
One still survives: dropping only the 1 cm absolute part of the face preference, leaving the 5% relative part.
The relative part alone settles the tie that is tested; the absolute part is there for overlaps so shallow
that 5% of them is below rounding noise, and no test exercises that yet.

## Step 4: the solver and the world (`solver.hpp`, `world.hpp`)

`World::step` finds contacts (every pair, after a bounding-sphere rejection), matches each point to last
step's by feature id, applies gravity, solves the contacts and moves the bodies. The solver is the 2D one
carried over: sequential impulses, the accumulated impulse clamped (push, never pull), a bounce fixed up front
from the arrival speed, Baumgarte push-out of overlap, warm starting. What is 3D-specific:

- **Effective mass with a tensor**: `1/m_eff = 1/mA + 1/mB + dir.((IA^-1 (ra x dir)) x ra) + (same for B)`,
  with each body's inverse inertia rotated into the world frame once per solve.
- **Friction in a plane**: two tangent rows per point, and the *length* of the combined friction impulse is
  limited to `mu * jn` (a circle). Clamping the rows separately makes a square, 41% stronger along diagonals.
- **A tangent basis that depends only on the normal** (`tangent_basis`), so a steady contact keeps the same
  tangents and its stored friction impulses can be reused.

Measured against closed forms:

| scene | expected | got |
|---|---|---|
| crate on a 20-degree slope, mu 0.1 | `g (sin - mu cos)` = 2.442 m/s^2 | 2.442 |
| slab sliding along a diagonal, mu 0.4, after 0.5 s | 4.038 m/s, no curving | 4.038, direction kept to 1 part in 10^5 |
| sliding solid sphere settling into rolling | 5/7 of its speed, `w = v/r` | 5.000 of 7, 10.000 rad/s |
| sphere rolling down a slope | `5/7 g sin(theta)` = 2.071 m/s^2 | 2.070 |
| resting crate, sum of its four normal impulses | `m g dt` = 0.08175 | 0.08175 |

Stacking: a tower of 4 and of 10 crates and a 4-3-2-1 pyramid of 30 stand still. Warm starting matters exactly
as in 2D: a tower of 4 holds on one sweep per step when warm and collapses when cold; a tower of 10 on four
sweeps likewise.

Tests (26, 85 in total): the table above; elastic exchange along arbitrary lines; momentum and restitution;
bounce to `e^2 h`; a crate dropped on a corner ends flat and still; one solve conserves total momentum and
angular momentum for 600 random pairs of tilted bricks (which is what proves the world-frame inertia is used);
the stored friction impulses on a held crate add up to `m g sin(theta) dt` uphill, for both tangent components;
friction mixing is symmetric; a body starting inside the floor is pushed out; friction is warm started; a
tumbling body in a World uses the gyroscopic setting; determinism. 22 deliberate breakages are each caught.
Six survived at first and each exposed a missing test, which is where the last six of those came from. One
more was caught only by crashing the test binary: my own tests indexed the contact list after a failed size
check. That is guarded now, and while chasing it I made the 8-to-4 point reduction safe against NaN positions.

Not here yet: joints, sleeping, continuous collision, a broad phase (it tests all pairs), and any way to look
at it.

## Step 5: a wireframe demo (`demo3d/`)

`./build/debug/demo3d/demo3d` draws the 3D world as lines through an orbit camera: boxes as their 12 edges,
spheres as three circles fixed in the sphere's own frame (so spin is visible), contacts as a yellow dot with a
green normal, and a ground grid. Further lines are dimmer; nothing is hidden. `demo3d/camera.hpp` holds the
camera and projection with no SDL in it, so it is unit tested (frame orthonormal and right-handed, target at the
screen centre, perspective halving with doubled distance, the field of view filling the window, segments
clipped at the near plane).

Scenes: `1` a tower of eight crates, `2` the 4-3-2-1 pyramid, `3` forty mixed boxes and spheres raining down,
`4` a slope with a crate that holds, a crate that slides and a ball that rolls, `5` three bricks tumbling in
zero gravity, one spun about each of its axes (the middle one keeps flipping; `G` cycles the gyroscopic mode,
and "explicit" makes them slowly gain speed, "off" stops the middle one flipping at all). Left-drag orbits, the wheel zooms, `Space` shoots a ball
from the camera, `F` drops ten more bodies, `W` and `-`/`=` change warm starting and the sweep count,
`Tab` pauses, `H` hides the text, `R` or `Backspace` reloads.

For checking without a window: `SDL_VIDEODRIVER=dummy demo3d --scene 2 --steps 240 --capture out.bmp`.

## Step 6: broad phase (`aabb.hpp`, `dynamic_tree.hpp`, `broadphase.hpp`)

The world no longer runs the narrow phase on every pair. Each body gets an AABB (for a turned box the extent
along world axis i is `sum_j |R_ij| h_j`), and `World::broadphase` picks how candidate pairs are found:

- **BruteForce**: every pair of AABBs. The reference.
- **DynamicTree** (the default): the 2D engine's balanced AABB tree, ported unchanged apart from the box type
  and its cost (surface area in place of perimeter). Boxes are stored padded by 10 cm and re-inserted only
  when a body escapes its padding; each movable body queries the tree.

As in 2D, both return pairs sorted by `(a, b)` and never a pair of two immovable bodies, so a scene steps
bit-for-bit identically with either (tested over 360 steps of a pyramid with balls). There is no sweep and
prune in 3D.

Measured (`-O2`, ms):

| scene | bodies | brute force | tree |
|---|---|---|---|
| broad phase alone, boxes scattered at constant density | 2000 | 6.10 | 0.83 |
| | 5000 | 28.97 | 3.38 |
| | 10000 | (skipped) | 8.01 |
| whole `World::step`, crates settled in a closed room | 1000 | 1.78 | 1.94 |
| | 2000 | 3.94 | 6.43 |

So the tree is 7 to 9 times faster where bodies are spread out, and that gap widens with size. In a dense
settled pile it is no faster up to 1000 bodies and slower at 2000, even though it does a third as many
comparisons: a brute-force test is six comparisons on contiguous memory, a tree visit chases a pointer, and
in a pile every padded box overlaps its neighbours so queries cannot prune much (374 node visits per body at
2000). The same thing showed up in 2D. The tree is the default because it is the one that scales; for a single
dense pile of a thousand bodies brute force is as good. Improving the tree's quality in piles (a bulk
rebuild, smaller padding for resting bodies) is not done.

Tests (10, 99 in total): AABBs of spheres and turned boxes, tight around all eight corners for 300 random
orientations; overlap, contains and merge checked along each axis separately; the tree stays valid and
balanced through 4000 random inserts, moves and removals and answers queries exactly like a linear scan; the
tree broad phase finds every real overlap while bodies jitter, teleport, vanish and appear; its work grows
sub-quadratically; the simulation is identical with either broad phase. 13 deliberate breakages are each caught.

The demo's body cap is now 1200, `B` switches the broad phase, and the readout shows box tests, candidate
pairs and confirmed contacts.

## Step 7: sleeping and continuous collision (ported from the 2D engine)

**Sleeping** works as in 2D: islands (bodies linked by contacts) are found with union-find, the floor links
nothing, and an island sleeps when its most restless body has been slower than 0.01 m/s and 0.035 rad/s for
half a second. Sleeping bodies are not integrated, solved or collided with each other; their contacts are
kept dormant with their impulses, so a sleeping pyramid that is poked wakes warm started (crates below the
poked one move under 2 cm). An awake body touching a sleeping island wakes it; so does setting a velocity,
force or torque. Moving a sleeper by hand needs `World::wake`.

**Continuous collision** against static bodies also follows 2D: a body that moves or turns more than half
its inscribed radius in a step is swept at poses no further apart than that radius; a body that was clear is
placed just inside the first thing it touches so the solver responds next step.

A 3000-shot stress run (spheres and boxes 8 to 40 cm, 30 to 200 m/s, spin up to 30 rad/s per axis, walls 1 to
10 cm thick) found three things the 150-shot test had not:

1. **8 shots got through with the sweep on.** For a body already in contact I had ported the 2D rule "may
   sink no deeper than a quarter radius", measured by the overlap depth `collide()` reports. For a plate
   thinner than the wall that depth stops growing as the plate slices in, then measures the way out the far
   side. The rule now follows the body's farthest point along the direction it entered by, and stops it when
   that point has advanced a quarter radius.
2. **Some shots left a collision with far more energy than they arrived with** (one, x780). That was not the
   sweep or the solver but the implicit gyroscopic step overshooting at extreme spin; see the correction
   under step 2.
3. **One shot appeared to deadlock**, pinned at the wall spinning at 500 rad/s. I first added an impulse to
   break it; measuring afterwards showed the impulse changed nothing (the stall had been the gyroscopic
   overshoot again), so it was removed.

After the fixes: 0 of 3000 get through (59% do without the sweep), none is left pinned at the wall, none ever
exceeds its starting energy, and every shot has turned back within 23 steps. The ten worst shots are kept as a
regression test, along with the first 400.

Tests (17, 116 in total): the 2D sleeping tests in 3D (asleep means bit-for-bit still; an island sleeps and
wakes as one; thresholds on speed and on spin; velocity, force and torque wake; a teleport needs `wake`;
sleeping does not change where a pyramid ends up); a fast ball and 150 random shots against a thin wall with
the sweep on and off; a fast fall onto a thin floor; a fast spinner stopped at a post it only crosses between
steps; slow bodies left bit-for-bit alone; and the regressions above. Of 20 deliberate breakages 19 are
caught; the one that survives (placing a capped body at the other end of a bracket 1/65536 of a step wide) is
equivalent. Two pieces of code that mutants showed to be doing nothing were deleted rather than kept: the
impulse above, and a "take the short way round" quaternion negation that was wrong for turns over half a
revolution.

In the demo `S` toggles sleeping (sleepers are drawn dim), `K` toggles the sweep, and the readout shows awake
bodies, islands, and how many bodies were swept and stopped.

## Step 8: joints (`joints.hpp`, `src/joints.cpp`)

Three joints, on the 2D engine's method: each is a few scalar velocity rows `J v = lin_a.v_a + ang_a.w_a +
lin_b.v_b + ang_b.w_b`, solved by impulses in the same sweeps as the contacts (joints first), with the
accumulated impulse clamped and kept to warm start the next step, and position error removed by a separate
pass that moves poses, never velocities.

- **Distance**: one row along the line between two anchors. A rod, or with `frequency_hz > 0` a spring-damper.
- **Ball and socket**: the two anchors coincide. A point has three coordinates, so three rows, one per world
  axis; all three rotations stay free.
- **Hinge**: the ball's three rows, plus two that stop relative turning about the two directions
  perpendicular to the axis. Each body carries its own copy of the axis (`a1`, `b1`); for a small
  misalignment `a1 x b1` is the rotation that has pulled them apart, which gives both rows their error. The
  angle about the axis is measured between a reference direction carried by each body, with `atan2`, so it
  lies in (-pi, pi]. On it sit an optional motor (target speed, torque limit) and optional limits (one-sided
  rows).

What changes from 2D: the angular part of a row is a vector, and the effective mass goes through the inertia
tensor, `J M^-1 J^T = |lin_a|^2/m_a + ang_a.(I_a^-1 ang_a) + ...`. Each row stores `I^-1 ang` for both bodies
once, which is also exactly the spin a unit impulse gives them.

In the world: jointed bodies do not collide unless the joint says so; a joint links its bodies into one island,
so they sleep and wake together; a motor told to turn wakes its mechanism; `truncate` drops joints whose
bodies went.

**A bug found here, which the 2D engine had too.** The position pass corrected joints one at a time but all
the rows of a joint from one measurement. For a hinge at its limit that fails: closing the pin turns the arm
off the limit, turning it back opens the pin, and the two overshoot against each other. An arm dropped onto
its lower stop hung 0.08 rad above it with the pin 8 cm apart, and its speed built up unchecked because it was
never "at" the limit. My first reading of the trace blamed the limit row; printing the errors pass by pass
showed them alternating. The pass now measures again after every row. The same sweep in 2D (108 arms) failed
56 times before and none after; that fix is on `main`.

Things measured, not assumed:
- A sphere on a ball joint swings with the compound-pendulum period `2 pi sqrt((I + m L^2)/(m g L))` to 0.5%.
- Swung in a circle it keeps its angular momentum about the vertical to 1% (at 960 Hz).
- Spin about an **off-centre** pin is slowly lost: 6% of the speed in 5 s at 240 Hz, a quarter of that at
  960 Hz, and the 2D engine loses the same. Each step moves the body along the tangent, the position pass
  pulls it back to the circle, and projecting the velocity onto the new tangent costs about `(w dt)^2` of the
  energy. Spin about a pin through the centre of mass loses nothing. My first tests demanded exact
  conservation and were wrong to.
- A hinged arm follows the 2D engine's revolute joint to 2 mm over 3 s.
- A torque-limited motor spins a wheel up at `torque / I`; it lifts an arm at 1.1 times the torque the weight
  needs and not at 0.9.
- An 8-link chain carrying a ball 26 times a link's mass, dropped from horizontal: worst pin gap 13.5 mm
  (10 sweeps), and its energy never rises. Heavier is worse (128x: 24 mm; 511x: 12 cm); more sweeps help.
- The undamped spring is not lossless: implicit Euler leaves about half the speed after 2 s at 240 Hz. Its
  period is right to 2%, and with damping ratio 0.2 successive peaks shrink by the textbook factor 0.277.

Tests (28, 144 in total): every row against a finite difference of the error it claims to be the rate of; the
measurements above; momentum of two free jointed bodies; the stored impulse of a resting pendulum equal to
`m g dt`; a pulled-apart, twisted hinge closed by the position pass with exactly zero velocity created;
limits, including the arm resting on its stop; islands, sleeping, waking by motor; a jointed arm lying on the
floor. Of 34 deliberate breakages 31 were caught at first. Two showed gaps (doubling the damping; dropping
the 8 degree cap on one correction) and now have tests; the third flips the sign of one hinge direction,
which changes nothing.

In the demo, scene `6` is a swinging chain, a plank bridge and a spring; scene `7` is a door with limits, a
two-part flail on crossed hinges, and a motor driving a paddle through crates. Joints are drawn in orange,
hinge axes in yellow. Checked by headless renders only.

Not done: a slider (prismatic) joint, a fixed (weld) joint, cone and twist limits for the ball joint, a mouse
joint. Hinge limits must lie inside (-pi, pi).

## Step 9: GJK and EPA (`gjk.hpp`, `src/gjk.cpp`)

Distance and penetration between any two convex shapes, from one function, `closest(a, b)`. It is not used
by the world yet: this step is the two algorithms on their own, checked against known answers.

The box and sphere code knows its shapes' geometry. These algorithms know one thing about a shape, its
**support function**: which of your points is farthest in direction `d`? Both work on the Minkowski
difference `A - B` (every point of A minus every point of B), which is never built, only sampled: its
farthest point along `d` is A's farthest along `d` minus B's farthest along `-d`. A and B overlap exactly
when `A - B` contains the origin.

- **GJK** finds the point of `A - B` nearest the origin. It keeps a simplex of up to four samples, finds
  the simplex's nearest point to the origin (a closest-point-on-segment/triangle/tetrahedron problem), and
  samples again in the direction of the origin from there. It stops when the new sample is no nearer than
  the simplex already gets (apart), or when the simplex holds the origin (overlapping).
- **EPA** takes over when they overlap. It grows GJK's tetrahedron outward, always pushing on the face
  nearest the origin, until that face cannot move: it is then on the surface of `A - B`, and its normal and
  distance are the direction and depth of the shortest way apart.
- A shape is a **core plus a rounding radius**: a sphere is a point, a capsule a segment. The algorithms run
  on the cores. If the cores are apart, the answer for the rounded shapes is the cores' distance less the two
  radii, exact, with no EPA even when the rounded shapes overlap.

Shapes so far: sphere, capsule, box, the convex hull of a point set, and any of those rounded. The result is
a signed distance (negative = depth), a normal from a to b, and the closest (or deepest) point of each.

**How it was checked.** Known answers first (spheres, box faces, a corner over a face, two crossed edges,
capsules, a small box inside a big one, identical boxes in one place). Then against the code that knows the
shapes: 3000 box pairs against a plain 15-axis SAT (worst difference under 1 mm) and 2000 sphere-box pairs
against `collide()`. Then an **audit that uses only support functions**, so it applies to every pairing: the
two points lie on their shapes; they are `distance` apart along the normal; flattened onto the normal the
shapes are apart by exactly `distance`; and for an overlapping pair no nearby direction gets out sooner.

A stress program runs that audit on two million random pairs: five kinds of shape, at scales 0.01, 1 and 50,
some 300 units from the origin, one pair in thirteen slid into exact contact first. All pass but one, where
the deepest points are about 5 mm out (the depth and normal are right); the code detects that case and
returns `converged = false`. Timing, release build: 0.2 microseconds for a separated pair, 1.4 when EPA runs.
GJK averages 4.7 iterations (worst 22), EPA 4.8 (worst 22).

**What the stress runs found that the unit tests had not.** The first version passed 18 of 19 tests at once.
Then:
1. *Deepest points wrong on flat sides.* A flat side of `A - B` (two box faces together) is covered by
   several triangles in one plane, and EPA's "nearest face" picks any of them; only one contains the point
   wanted. The deepest points now come from a second GJK run on the pair pulled just clear along the normal.
2. *Overlap reported for shapes a hair apart.* Rounding let a nearly flat tetrahedron pass the "origin is
   inside" test. A tetrahedron flatter than a ten-thousandth of its size now encloses nothing.
3. *A normal a third of a degree out* for capsules whose axes nearly cross: the closest points were a tiny
   difference of large numbers. For a triangle simplex the normal is now the triangle's own.
4. *A pair in exact contact reported 14 cm apart.* See the next paragraph.

**A mistake of method, corrected.** Following my rule (if breaking some code changes no test and no stress
result, delete it), I deleted five pieces after a 500000-pair stress showed them idle. Two of them were not
idle. The stress used random pairs, and random pairs are never exactly touching, which is exactly when those
two pieces act. With touching pairs added, about 1% of them failed. One piece went back: GJK treats anything
closer than a ten-thousandth of the shapes' size as touching and hands over to EPA, because single precision
cannot give a direction for a separation that small (so distances under 0.1 mm on a 1 m object read as zero).
The other, a guard for zero-area triangles, went back and was measured again with touching pairs, still
changed nothing in two million, and is deleted. The rule stands, but the stress has to contain the cases the
code is for.

Tests (28, 172 in total): the known answers; the comparisons with SAT and `collide()`; the audit on 6000
random pairs, 4000 pairs slid into exact contact and 6000 a hair either side of it; order of the shapes,
position in the world and scale; the closest-point routines against a slow independent version in all seven
regions of a triangle; a ball centred exactly on a box's corner, edge or face; a rod ending exactly on a box
corner; the capsule and deepest-point regressions; iteration counts. Of about 40 deliberate breakages of the
final code, all are caught except one line: the "still overlapping after being pulled clear" half of the
`converged` check has never been seen to happen, so nothing tests it.

Not done: wiring this into `collide()` and the world, which needs contact MANIFOLDS (GJK/EPA give one point;
a box resting on a face needs four) and the new shapes as bodies with mass and inertia.

## Planned next

Contact manifolds from GJK/EPA, then capsules and convex hulls as bodies in the world.
