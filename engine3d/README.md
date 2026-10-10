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

*Correction (step 11):* points 1 and 3 were both wrong in the end. The farthest-point rule has a hole of its
own, and the impulse I removed as doing nothing is needed; neither showed until capsules gave the stress
long thin bodies. See step 11 for the rule as it now stands.

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

## Step 10: contact manifolds from GJK and EPA (`collide_convex` in `collision.hpp`)

GJK and EPA answer with a normal, a depth and ONE point. A box resting on one point rocks, so the solver
needs the whole patch. `collide_convex(a, b, manifold)` keeps the `collide()` contract (normal from a to b,
up to four points, each halfway through the overlap, each with an id) for any two `Convex` shapes, and
`World::narrowphase = NarrowPhase::Convex` sends every pair, and the continuous-collision sweep, through it.
The default is still the per-shape code.

How the patch is built: each shape offers the **feature** it presents to the other along the normal. A box
offers the face that faces that way most squarely, a capsule its axis, a sphere its centre.
- If a box face is what is touching (the normal is within 0.8 degrees of the face's own), that face is the
  reference: the other shape's feature (a face, a segment or a point) is clipped to the face's outline, and
  what has sunk below the face is the patch. This is the box code's face clipping with the incident "face"
  allowed to be a segment or a point, and with both shapes' rounding radii. A capsule lying on a box gets a
  point under each end; tilted, only the lower end is below the face, so it gets one.
- Two capsules with (nearly) parallel axes get the two ends of the stretch where they run side by side,
  each end measured across for itself.
- Anything else (a corner, crossed edges, anything round) meets at a point: EPA's.

Two things that had to be got right:
- **The depth goes with the normal.** A face manifold uses the face's own normal in place of EPA's, and the
  two can differ by a fraction of a degree. Keeping EPA's depth with the face's normal understated the
  overlap, by up to 2.5 cm in 14600 of 300000 nearly-aligned pairs. The depth is now recomputed along the
  normal actually used, from the support functions. My first tests did not notice: I had the face tolerance
  at 2.5 degrees, tightened it for another reason, and the symptom shrank below the test's tolerance while
  the cause stayed. A targeted stress found it and the audit now checks it directly.
- **How loosely "a face is touching" is judged.** At 2.5 degrees a face that was not the real contact could
  take over and report a deeper overlap than the true one. When a face really is touching, EPA's normal
  is that face's normal to rounding, so the tolerance can be tight.

Measured against the code that knows the shapes:
- 20000 random box pairs: the two always agree on hit or miss. Where they choose the same normal (77% of
  hits), the points and depths match in all but 3 of 10089. Where they do not, it is the box code's
  preference for face axes against EPA's true minimum.
- Sphere-sphere and sphere-box: the same point, depth and normal.
- A tower of 8 crates ends in the same place to four decimals on either path and sleeps. A closed room of 120
  mixed boxes and balls settles alike; the convex path costs 0.112 ms a step against 0.084 (release).

Tests (17, 189 in total): the patch for each kind of contact (box on box, overhanging an edge, capsule
flat, tilted, upright and overhanging, capsules side by side, crossed and end to end, a rounded box); the
comparisons above; an audit of 6000 shallow contacts between spheres, boxes, capsules and rounded boxes (each
point inside both shapes, depths consistent, ids distinct, moving out by the depth separates the pair); ids
that survive a small movement and change when a box rests on a different face; and in the world: the tower,
a settling pile, friction on a slope, a fast box stopped at a thin wall, and the switch itself. Of 37
deliberate breakages 30 were caught at first. Five showed gaps, now tested (among them the depth above);
two changed nothing even on 300000 targeted pairs and the code they guarded is gone.

In the demo `N` switches the narrow phase. Not tried by hand.

Not done, and the reason hulls come next: a hull here is a cloud of points with no record of its faces, so
it offers only its farthest vertex and gets a one-point contact. A hull resting on a face would rock.
Capsules have full manifolds but are not bodies yet (no mass, inertia, bounding box or drawing).

## Step 11: capsules as bodies

A capsule is every point within `radius` of a segment along the body's own y axis: `Shape::capsule`,
`Body::solid_capsule`, `Body::fixed_capsule`. It needed a mass and inertia, a bounding box, `contains`, and a
place in `collide()`; there is no hand-written capsule collision, every pair with a capsule in it goes
through `collide_convex` (step 10).

**Inertia.** Split the mass by volume into the cylinder's share `mc` and the two end caps' `ms` (together one
sphere). About the axis it is a cylinder plus a sphere, `mc r^2/2 + ms 2r^2/5`. Across the axis the cylinder
gives `mc (L^2/12 + r^2/4)` and the caps, moved out to the ends by the parallel axis theorem (a hemisphere's
centre of mass is 3r/8 from its flat face), `ms (2r^2/5 + h^2 + 3hr/4)` with `h = L/2`. Checked against a
brute-force integral to 0.4% for three proportions, and against the two limits: no cylinder is a sphere, no
radius is a thin rod.

Measured in the world:
- Dropped flat, it rests on two contact points and sleeps; dropped at an angle it ends lying down.
- It rolls down a slope at `g sin(theta) / (1 + I/(m r^2))` to 2%, without slipping.
- Hung from a hinge by its tip, it swings with the compound-pendulum period to 0.5%.
- Three layers of logs stacked crosswise, each contact a single point between two round surfaces, stay put
  and sleep. A pile of 30 settles in a closed room.
- In collisions with spheres, boxes and other capsules, linear momentum is kept to a part in ten thousand
  and energy never rises. Angular momentum drifts by up to 5%. I measured where: mostly in free flight after
  an off-centre hit sets a long thin body tumbling at up to 19 rad/s, which is the gyroscopic step, not the
  contact.

**A bug in the sweep, older than capsules, that only capsules showed.** The rule for a body already touching
a wall (step 7) followed the body's farthest point along the direction it entered by. Of 20000 harsh capsule
shots, 22 went through the wall. A thin capsule that arrives end-on and turns side-on as it goes in reaches
no deeper with its farthest point, since the end swings back as fast as the body advances, while all of it
passes through. Watching the centre as well left one shot in 60000: its touching end was swinging away from
the wall at 340 rad/s, so the solver saw a separating contact and did nothing, and the sweep, which held the
pose but not the velocity, let it creep in a quarter radius per step. That is the case for the impact impulse
I had tried and deleted twice because no box stress needed it.

With the impulse in, I broke each part of the rule in turn and ran every stress I have: 200000 harsh capsule
shots, 100000 thin boxes, the original 3000 box shots (which also demand that every shot is heading back),
and 50000 bodies that start against a thin wall with up to 80 rad/s of spin and no speed toward it. What
they say:
- Watching the farthest point is not needed once the centre is watched; nothing changes without it.
- Watching the centre is needed (3 capsules in 100000 get through without it).
- The impulse is needed, and so is its bounce (one box shot stops dead without it).
- Giving the impulse an arm and a torque changes nothing.
- The cap of a quarter of the inscribed radius is needed: at a full radius, 58 of the 50000 spinners get
  their centre into the wall.
So the rule is now: a body already in contact may move its **centre** no more than a quarter of its
inscribed radius further in per step, and if it is held there, the part of its velocity that was carrying
it in is reversed as in any bounce. It is shorter than what it replaced, a ratio instead of a search, and
the function that found a shape's farthest point is gone. After it: all four stresses clean, no shot
through, none stuck, none gaining energy.

Tests (18, 207 in total): the inertia integral and limits; mass, bounding box against the support function,
`contains` against the distance to the axis, `collide()` against `collide_convex`; the measurements above;
120 fast spinning capsules at a 4 cm wall with the sweep on and off; 500 harsh shots; the shots that got
through earlier rules, kept by number; the spinners that fail with a full-radius cap; and a held body
bouncing by the bouncier of the two materials. 21 deliberate breakages of the capsule code: 20 caught, and
the one that was not led to everything in the paragraph above. Of the final sweep rule, 6 breakages, all
caught after two more tests.

In the demo capsules fall in scene 3 and a log rolls down the slope in scene 4 (headless renders only).

Not checked: the 2D sweep. It uses a different rule for a body in contact (a cap on the overlap depth) and
passed 300000 harsh shots, but none of those shapes was as thin for its length as these capsules.

## Step 12: a convex hull builder (`hull.hpp`, `src/hull.cpp`)

GJK and EPA are content with a convex shape as a cloud of points. A contact patch is not: to rest on a face
a shape has to know it has one. `build_hull(points, count, hull)` turns the cloud into the corners that
matter and the flat polygons between them, each with an outward normal and its corners in order. Points
inside, on a face or along an edge are dropped. A cloud that encloses no volume is refused.

The method is EPA's polytope growing, run to completion (Quickhull without the bookkeeping that makes it
fast on big inputs):
1. Start with a tetrahedron of four of the points, spread as widely as possible.
2. Take the point farthest outside the current surface. If none is outside by more than the tolerance, stop.
3. Remove every triangle that point can see, spreading out from one across shared edges. Close the hole
   with a fan of triangles from its rim to the point. Go to 2.
4. Merge triangles that lie in one plane into polygons, so a cube is six squares and not twelve triangles.
   Faces are grown from a seed triangle, and a neighbour joins if its corners lie in the SEED's plane.
5. A point left on only two faces is a point along an edge, not a corner; drop it.

The tolerance, a hundred-thousandth of the cloud's size, decides two things only: how far outside a point
must be to be worth adding, and how level two triangles must be to be one face.

**What went wrong on the way, all found by stress and none by the first tests.** The first version built
every solid I knew the answer for and then failed on one very flat cloud in forty.
- *Visibility with a tolerance is wrong whichever way it leans.* "Can this point see this triangle" has to
  be answered for triangles the point is almost level with. Count level as seen, and triangles the point is
  slightly BELOW are removed too: the surface shrinks away from corners it already had. Count level as not
  seen, and the triangle that closes the hole can cut under a corner by far more than the tolerance (the
  error is levered by how close the point is to the rim). Use zero, in single precision, and rounding decides
  for each triangle of one flat face separately: some go, some stay, and new triangles are laid on top of
  old ones.
- *The fix is to make the question exact.* The builder does its own arithmetic in double precision (input
  and output stay in the engine's `Real`). Single precision cannot even place the plane of a thin triangle
  to within the tolerance: one a hundredth as wide as it is long has a normal good to a part in 100000.
  In double, "in front of" means in front of, the removed patch is always a disc with one rim, and the
  tolerance goes back to deciding only what is worth adding and what counts as flat.
- *Merging neighbours pairwise lets a face creep round a curve*, each step within tolerance. Hence the seed.
- *A plane computed from positions far from the origin* (Newell's sum of cross products) loses everything
  to cancellation; the face now simply takes its seed triangle's plane.

**Then, deleting.** After the move to double precision, several pieces written for the single-precision
version did nothing in 200000 stress clouds, and went:
- special handling for a point exactly in line with a rim edge. It cannot arise: such a point is at least
  twice as far outside every triangle as the edge's far corner ever was, so farthest-first brings it in
  before that corner exists;
- a fallback for clouds whose points all share one x (such a cloud is flat, and refused anyway);
- remembering which points were already corners; growing faces from their largest triangle; requiring that
  a joined triangle face the same way as its seed (I could not construct a case, and a hand-built thin wedge
  showed why: a corner is only added if it is more than the tolerance outside, so a solid cannot have two
  opposite faces within the tolerance of each other); a fallback margin and a rim check in step 3.
What is left as a net: a triangle with no area, a merged face that is not convex, fewer than four faces, or
`V - E + F != 2` all refuse the build. Only the convexity refusal has ever been seen to fire since
(12 flat clouds in 150000, kept as a test); the others overlap each other and breaking any one of them
alone changes nothing, so they are untested as individuals.

**Checked how.** An audit that knows nothing of how the hull was built: closed (every edge has its twin),
every face flat, convex and wound counter-clockwise, no input point in front of any face, every vertex an
input point on at least three faces, `V - E + F = 2`, no two neighbouring faces squarely coplanar. A stress
of a million clouds of eight kinds (random, on a sphere, random lattice points, prisms, boxes with points
on their faces, pyramids, slabs a thousand times wider than thick, clouds with repeats) at scales 0.01, 1
and 100, half rotated, a fifth far from the origin: every hull built passes; 1074 are refused, all flat
lattice subsets or slabs. And one independent opinion: the faces enclose exactly what GJK says the raw
points enclose.

Tests (12, 219 in total): the cube (and a 1 x 2 x 3 box: six rectangles at the right distances); a cube
buried in interior, on-face, on-edge and repeated points, in 30 orders; a 4 x 3 x 5 lattice, rotated 40
ways, that must come out as 8 corners and 6 faces; tetrahedron, octahedron, icosahedron; prisms of 3 to 12
sides keeping their ends as single faces; points on a sphere; 400 random clouds whose hull must reach
exactly as far as the cloud in every direction; the comparison with GJK; independence of point order; 600
flat clouds and the 12 dented ones; clouds with no volume. Cost: 14 microseconds for 8 points, 0.25 ms for
100, 24 ms for 1000 (it is quadratic and more; fine for shapes, not for point-cloud scans).

Not done: nothing uses the faces yet. Hulls still collide through their point cloud and get a one-point
contact.

## Planned next

Face contacts for hulls in `collide_convex`, then hulls as bodies (mass and inertia of a polyhedron).
