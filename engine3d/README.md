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

## Step 13: face contacts for hulls

A `Convex` can now carry a built hull (`Convex::hull(built, pos, q)`), and `collide_convex` uses it: the
feature a hull presents along the contact normal is the face that faces that way most squarely, a convex
polygon of up to 16 corners. The clipping of step 10 was written for a box's rectangle; it now takes any
convex polygon as the reference (one clipping plane through each side, facing outward) and any polygon,
segment or point as the thing clipped. Boxes go through the same code.

- A tetrahedron rests on three points. A hexagonal prism on its end gets four of its six corners, spread
  round it; lying on a side, that rectangle's four.
- A cube overhanging a slanted side of a hexagon is cut along that side, and the patch is the same whichever
  of the two is asked first.
- A capsule lying on a hull's face gets a point under each end; a ball, one; a rounded hull is its radius
  taller.
- A cloud of points that was never built into a hull still works, with one point.

Checked against boxes, which are the shape with a known answer: a hull built from a box's eight corners
gives the box's normal, depth and patch on 4000 shallow contacts, alone and against a real box. And the box
code and the convex path still agree on 20000 random box pairs: 10087 of the 10088 contacts with the same
normal have the same points. An audit of 6000 shallow contacts over ten hulls (prisms, a tetrahedron, random
hulls) against hulls, boxes, capsules and balls, half of them pressed face to face, holds every manifold to:
depth equal to the overlap along its normal, every point no deeper than that, ids distinct, and moving out
by the depth separates the pair. A stress of 200000 such contacts over 30 hulls is clean.

Three things that changed outside the new code, each because a hull showed it up:
- **The clip margin overshoots at sharp corners.** The 2 mm margin (step 3) lets a corner lying exactly on a
  side stay a corner. On a rectangle a point can end up 2.8 mm outside the face. On a thin triangular face
  of a random hull the two margins meet 14 cm beyond the tip. Points the margin let through are now brought
  back to the nearest point of the face's outline. The box code got the same correction, so that the two
  paths stay identical; its contact points near a face's edge moved by up to 2 mm.
- **Which four points are kept was decided by rounding.** A hexagon resting flat has six equally deep
  corners. An exact "pick the deepest, then the farthest..." chose a different four after a nudge of a
  millimetre, and the solver would have lost its stored impulses every time. A later candidate now wins only
  if it is clearly better (0.1 mm of depth, a thousandth in distance or area). The choice is greedy, not a
  search: on stretched and lopsided faces it spans at least 65% of the best area any four corners could.
- **"A contact point is inside both shapes" is only true of boxes.** A point halfway through the overlap,
  under a corner of a face whose neighbouring sides slope inward (a tetrahedron's, most hulls'), lies just
  outside that sloping side, by less than the overlap. The hull audit allows exactly that much.

Tests (9, 228 in total): the above. Of 24 deliberate breakages 20 were caught at first. Two showed gaps, now
tested (the spread of the four points; stability under a tilt too slight to see). Two changed nothing and
their code is gone: a term in the contact id for which shape is the reference, and a redundant check.

Not done: hulls are not bodies, so no simulation uses any of this yet. A hull face with more than 16
corners falls back to a one-point contact.

## Step 14: hulls as bodies

`Body::solid_hull(points, count, density, origin, q)` builds the hull of the points, works out what a solid of
that shape weighs and how it resists turning, and returns a body (or nothing, if the points enclose no
volume). `Body::fixed_hull` is the immovable one. `Shape::Type::Hull` holds the hull by shared pointer, so a
hundred copies of one shape share one hull.

**Where the body is.** A body's `pos` is its centre of mass; everything in the engine assumes it. So the
points may be given in any convenient frame, the caller says where that frame's origin goes, and the factory
moves the hull's corners to be measured from the centroid and sets `pos` to where the centroid lands. The
corners end up in the world exactly where they were asked to be.

**Mass, centre and inertia** (`mass_properties` in `hull.hpp`). Cut the solid into tetrahedra, one from the
origin to each triangle of each face. Their volumes are signed, so wherever the origin is, the parts outside
the solid cancel. For a tetrahedron with corners 0, a, b, c:
- volume `V = a . (b x c) / 6`, centre `(a + b + c) / 4`;
- the integral of `x x^T` over it is `V / 20 * (a a^T + b b^T + c c^T + (a + b + c)(a + b + c)^T)`.
Summing gives the second moment `S`; the inertia tensor is `trace(S) E - S`, moved to the centroid by taking
away `V c c^T`. In double precision, because the pieces are large and mostly cancel. The tensor is generally
not diagonal (a lopsided solid spun about one axis is thrown about the others), and the body keeps it as it
is: nothing in the engine needed it diagonal.

The rest is what every shape needs: a bounding box (the corners, turned), the radii the sweep uses (to the
nearest face, to the farthest corner), `contains` (behind every face), and collision through `collide_convex`
with the face contacts of step 13.

Measured:
- A box-shaped hull has the box's volume and inertia, at the origin and 75 units from it.
- Prisms, an icosahedron, a tetrahedron and six random lumps match a count of 343000 little cells to 2% in
  volume, centre and every entry of the tensor.
- A tower of eight box-shaped hulls ends within 2 cm of a tower of boxes and sleeps.
- A dropped regular tetrahedron rests on a face on three points with its centre `a sqrt(6) / 12` up; a
  twenty-sided die thrown spinning ends on a face, its centre `g^2 / sqrt(3)` half-edges up.
- A hexagonal prism stands on its end, and lying on a side does not roll.
- 36 lumps, prisms and dice settle in a closed room. 150 fast spinning hulls at a 4 cm wall: none through or
  stuck with the sweep, most through without.
- A free lopsided hull's centre keeps to a straight line within a millimetre over four seconds.

Tests (15, 243 in total). 26 deliberate breakages: 25 caught at first; the one that was not (density left
out of the inertia) showed that the test only checked the mass at a density other than 1, and now it checks
both. This step went in without a bug found, which after the last four is worth saying, and worth not
trusting too far: the stress programs of steps 9 to 13 already exercise everything underneath it.

In the demo, lumps and prisms fall in scene 3, drawn by their edges (one headless render).

Not done: hull faces with more than 16 corners get one contact point; nothing reduces a hull with very many
corners (support is a search through all of them).

## Step 15: where the time goes

Nothing so far was written to be fast, and nobody knew what it cost. `bench3d/phys3d_bench` drops bodies
into a closed room, lets them settle with sleeping off (so every body is simulated every step), and reports
the time per step and per phase; `StepStats` now carries the phase times. Four scenes: crates, mixed shapes
(spheres, boxes, capsules, hulls), hulls only, and chains of ball-jointed links.

First numbers, 1000 bodies, release build, ms per step at 120 Hz:

| scene | step | broad | narrow | solver |
|---|---|---|---|---|
| crates, per-shape code | 5.1 | 3.4 | 0.55 | 1.1 |
| crates, GJK/EPA | 5.3 | 1.6 | 2.6 | 1.1 |
| mixed | 6.0 | 3.3 | 1.7 | 1.0 |
| hulls | 9.8 | 3.3 | 5.3 | 1.1 |
| chains (1000 links) | 1.9 | 0.45 | 0 | 0.56, and 0.86 in the position pass |

**The broad phase was the biggest cost, and worse than doing nothing clever.** At 1000 crates the tree
visited 694 nodes per body: more than brute force's 500 box tests. Three causes, three changes.
- *A lived-in tree is a bad tree.* Bodies fall, and each is taken out and re-inserted many times on the way.
  Every re-insertion is placed as well as it can be in the tree as it stands, and nothing revisits the
  choices made around it. The same boxes inserted fresh needed 200 visits. The tree now has `rebuild()`:
  keep the leaves, split them in halves at the median along the direction they are most spread, and again
  for each half. The broad phase calls it once a quarter of the leaves have been re-inserted (waiting for
  all of them cost two to four times the work; rebuilding on every single re-insertion gained under a
  tenth more). A pile at rest re-inserts nothing and never pays.
- *One question per body repeats the same walk.* Asking "who overlaps this box?" for each body starts at the
  root every time. `for_each_overlapping_pair` asks the tree about itself once: the pairs inside a node are
  the pairs inside each child plus those between the two children, and two subtrees whose boxes do not
  overlap are dismissed with one comparison however many leaves they hold. 694 visits per body became about
  60 pair tests.
- *The padding's extra pairs went to the narrow phase.* They are now checked against the real boxes first, so
  the tree returns exactly the brute-force list.

**The convex narrow phase did work it threw away.** Profiled on a settled pile of hulls: 0.94 of the 1.06
microseconds per candidate pair was in `closest()`, and about half of that for touching pairs was the second
GJK run that finds the deepest points (step 9). A face contact never looks at those points. `closest` can now
leave them out, and `collide_convex` asks for them only when it ends up with a one-point contact.

After, same scenes and machine:

| scene | step | broad | narrow | solver |
|---|---|---|---|---|
| crates, per-shape code | 1.9 | 0.33 | 0.42 | 1.2 |
| crates, GJK/EPA | 3.3 | 0.40 | 1.7 | 1.2 |
| mixed | 2.8 | 0.45 | 1.3 | 1.1 |
| hulls | 6.2 | 0.62 | 4.3 | 1.2 |
| chains | 1.9 | 0.32 | 0 | 0.58, and 0.93 in the position pass |

The machine was not idle for either set, and two runs of one scene differed by up to a factor of two in the
broad phase before the changes; the improvement there is well outside that. At 2000 bodies: crates 4.4 ms,
mixed 6.7, hulls 12.5.

What is left, in order of size:
- The GJK/EPA narrow phase, 0.9 microseconds per candidate pair, is most of a hull scene. It starts from
  nothing every step; starting from last step's answer for the pair is the usual remedy.
- The solver, 33 ns per contact point per sweep, is most of a crate scene.
- The joint position pass costs more than the joint solver. It rebuilds all of a joint's rows for each row
  it corrects, which is what made it correct in step 8.

Tests (3, 246 in total): the tree asked about itself returns exactly the brute-force pairs, on sparse, dense
and piled-up boxes, before and after a rebuild; a rebuild keeps every leaf's id and box, halves the work on
a tree that 600 boxes have settled through, and gives the same tree twice; the broad phase rebuilds only
after a quarter of the bodies have moved. 21 deliberate breakages: 19 caught. One showed a line that could
never matter and is gone. The other is the rule that breaks ties by leaf id when splitting at the median: it
is there so that the tree is the same on every machine, which a test on one machine cannot show.

## Step 16: more joints (weld, slider, a cone for the ball joint)

Three more, on the same rows and the same position pass as step 8.

- **Weld**: two bodies held as one rigid piece. The ball joint's three point rows, and three more that
  forbid any relative turning. The joint remembers b's orientation as seen from a when it was made
  (`conj(q_a) q_b`). Whatever rotation takes that to what it is now is the error: a turn of some angle about
  some axis, read off the quaternion, turned into the world frame, and shared out over three rows, one per
  world axis.
- **Slider**: b slides along an axis fixed in a, and nothing else. The weld's three orientation rows, two
  rows across the axis, and along it a motor (a speed, with a force limit) and two end stops. A row along a
  direction `e` that rides with body a measures `C = e . (pb - pa)`, and because `e` turns with a,
  `dC/dt = e . (v_b + w_b x rb - v_a - w_a x ra) + (w_a x e) . d`: a's angular part has `ra + d` where a pin
  joint's has `ra`. (The 2D prismatic joint has the same term.) That uses all eight row slots.
- **Cone** on a ball joint (`limit_cone`): each body carries a copy of an axis, and the angle between the
  copies may not exceed `cone_angle`. That angle grows when b turns relative to a about `a1 x b1`, so it is
  one one-sided row about that direction, built like a hinge's limit. Turning about the axis itself is not
  touched.

A motor told to move wakes a sleeping slider, as it does a hinge.

Measured:
- Two welded bodies tumbling free: pin within 1 mm, relative orientation within a tenth of a degree, one
  angular velocity between them, momentum kept.
- A beam welded to the world at one end: the weld's stored impulses are the beam's weight (to 0.2%) and its
  weight times the distance to its centre (to 0.5%).
- A block on a rail tilted 0.6 rad, struck with velocity and spin mostly off the rail: it moves only along
  it (2 mm, a fifth of a degree), at `g cos(0.6)` to 0.5%.
- End stops hold, the lower one carrying exactly the block's weight. A motor at its force limit gives
  `v = F t / m`, and lifts a load at 1.1 times its weight and not at 0.9.
- The slider follows the 2D prismatic joint to a millimetre on a tilted rail with end stops.
- An arm swung hard on a ball joint with a 0.5 rad cone goes out to 0.5 and no further, and the stop adds
  no energy; swung gently it ends in the same place to 0.01 mm with the cone or without. Leaning on the cone
  and spinning about its own length, it keeps its spin.

Tests (12, 258 in total): the above, and every new row against a finite difference of the error it claims
to be the rate of. The first run failed only where my tests were wrong (a step count; an energy tolerance
tighter than the integrator's own wobble), which is exactly when a first-run pass means least, so: 30
deliberate breakages, 29 caught. The one that was not is the line that reads a quaternion "the short way
round": `q` and `-q` are the same orientation, and without it a body whose quaternion happens to be written
the other way would look a full turn out of place. No test wrote one that way; one does now.

Not done: a twist limit for the ball joint (the cone leaves twist free), a mouse joint, a fixed-distance
"rope" with slack. No demo scene uses the new joints.

## Step 17: starting GJK from last step's answer

Two bodies resting on each other are in almost the same place every step, and GJK was finding the same
simplex from scratch every time. `GjkCache` keeps, for one pair, the directions in which GJK sampled the
corners of its final simplex; the next call samples A - B in those directions again. For a pair that was
touching those four points usually still enclose the origin, which is all GJK has to establish before EPA
takes over; for a pair that was apart they are already next to the closest point. The world keeps one cache
per candidate pair. Its pair list and last step's are both sorted, so finding last step's entry is a walk
along it, not a lookup.

Measured before building it, on settled piles of 1000 crates and 1000 hulls:
- starting only from last step's NORMAL saved 3 to 9% per pair, too little to be worth carrying around;
- starting from last step's SIMPLEX saved a third (0.36 to 0.24 microseconds per pair for crates, 0.61 to
  0.41 for hulls), with every answer the same.

In the world, narrow phase per step at 1000 bodies: hulls 4.3 to 3.2 ms, mixed 1.26 to 1.06, crates on the
GJK path 1.74 to 1.35. At 2000 hulls, 8.7 to 6.6. A cache is only a place to start: one from a different pair
altogether still gives the right answer, and a million random pairs each started from a cache warmed on a
slightly different pose pass the step 9 audit without a single failure.

**A bug this found, older than the cache.** A test that compares the world's warm-started contacts with
fresh calls disagreed on one pair: a capsule resting on the 60 m floor. The fresh answer was the wrong one,
its normal tilted by 0.8 degrees and its depth inconsistent with that normal by 0.87 m. GJK's closest point on
a triangle (step 9) decides which region of the triangle the origin is in by multiplying dot products, whose
size goes as the fourth power of the triangle's. For samples of a 60 m floor against something 20 cm above
it that is 10^7, all single precision can resolve, and it picked the wrong region. That routine now works in
double precision, and 500 boxes and capsules resting on the 60 m floor must get a level normal.

**A limit it also showed, not fixed here.** GJK's tolerances are relative to the size of the numbers it works
with, and against a big floor those are the floor's corners. A stress of small shapes on floors of three sizes:

| floor half-width | distances wrong by over 1 mm | worst |
|---|---|---|
| 5 m | 0 of 200000 | under 1 mm |
| 30 m | 1167 | 4.5 mm |
| 300 m | 26044 | 11 cm |

Normals stay level; it is distances near contact that suffer (within a ten-thousandth of 42 m, about 4 mm,
everything reads as touching). The next step is to do GJK and EPA's own arithmetic in double precision, as
the hull builder does.

`StepStats` now also counts the pairs that ran GJK in a step and their iterations: on a settled pile each
takes one (started afresh, the same pairs took 3.3 each).

Tests (5, 263 in total): a warm cache gives the same answers in under 70% of the iterations; a cache from
anywhere else does no harm; the world's warm-started contacts equal fresh calls, and settled pairs average
under two iterations; the level floor, both for the 500 resting shapes and for the eight floor-stress cases
that went wrong. Deliberate breakages: the world not reading the cache, not writing it, or forgetting it each
step, and the triangle routine back in single precision, are all caught. The world reading the wrong pair's
entry is not: it costs time, never correctness, which the cache-from-anywhere test is there to show.

## Step 18: GJK and EPA in double precision

GJK and EPA now do their own arithmetic in double precision, as the hull builder does; the engine's Real
is still what goes in and what comes out. The double vector they share is in `src/vec3d.hpp`.

What single precision had cost was not the arithmetic itself but the tolerances it forced. GJK's answers are
differences of support points, and next to a big shape those are big numbers, placed to a part in ten
million. So "touching" had to mean "within a ten-thousandth of the shapes' size", which next to a 300 m
floor is 4 cm. In double the same points are good to 10^-16 of their size, and three tolerances tighten:

| tolerance | single | double |
|---|---|---|
| GJK: closer than this counts as touching | 10^-4 of the size | 10^-10 |
| GJK: a tetrahedron flatter than this encloses nothing | 10^-4 | 10^-9 |
| EPA: a face this close to the surface is on it | 2 * 10^-5 | 10^-9 |

The floor stress of step 17, now:

| floor half-width | distances wrong by over 1 mm (of 400000) | worst |
|---|---|---|
| 5 m | 0 | 1.2 * 10^-7 m |
| 30 m | 0 (was 1167) | 1.2 * 10^-7 m (was 4.5 mm) |
| 300 m | 0 (was 26044) | 1.2 * 10^-7 m (was 11 cm) |

The 10^-7 m left over is the shapes' own: a float placing the bottom of a crate. The million-pair audit of
step 9 used to pass at a tolerance of 2 * 10^-3 of the shapes' size; it now passes at 10^-6 without a single
failure (step 17's GJK fails 1489 pairs there).

**What did not need to change.** Several other tolerances were tightened too, and putting each back made no
difference anyone could measure, so they were put back: GJK's convergence test (one part in a million of the
distance squared; the cores are polytopes, and GJK nearly always ends first by sampling a corner it already
has), EPA's "two samples are the same point", and how far apart find_deepest_points pulls the shapes.
Choosing a support corner in double precision made no difference either: a wrong choice is between two corners
within 10^-7 of a tie.

**What did.** Turning that corner and moving it into place in double. A crate resting on the floor, tilted by
10^-7 radians by the solver, has its lowest corner 9 * 10^-8 m into the floor. Rounded to single precision and
then measured in double it came out 1.5 * 10^-8 m above, so the contact began a step later, and a tower that
had stood exactly as it does on the box code (a test since step 10) ended 2 mm to one side of it. The answer
had been within a rounding step all along; double precision moved it to the other side of zero. So
`Convex::core_corner` chooses the corner in the shape's own frame, and GJK turns and places it in double. That
also makes the answer the same wherever in the world the pair is: added in single precision, the shape's
position rounded each support point to the float spacing there, a millimetre at 10 km.

Cost: the narrow phase on the GJK paths is 5 to 16% slower (crates on GJK/EPA: 250 0.26 to 0.30 ms, 1000
1.32 to 1.44; 1000 lumps 3.10 to 3.25), about a third of what step 17 saved. Whole steps are about 3% slower.

One test changed: the pile of boxes and balls had required every body to be below 0.2 m/s at the end, "give
or take a ball still rolling". With this step a ball was: resting on a crate tilted by 0.1 degrees (which the
pile leaves it at), it rolled at 0.005 m/s^2, as 5/7 g sin(tilt) less the friction of the wall beside it
predicts, until it rolled off the crate's edge. Which ball does that changes with any change to the engine.
The test now requires the crates to be at rest and the balls under 2 m/s.

Tests (3, 266 in total): small shapes near 30 m and 300 m floors get their distance to 10^-5 m, converged,
on the right side of zero and with level normals; a crate 10^-2 to 10^-5 m above a 300 m floor is found apart
by GJK, with its distance and both points to a micrometre; 200 pairs give the same answer to a micrometre
8 km from the origin as at it. Deliberate breakages: each of the three tolerances back to single precision,
triangle weights rounded to float, the pulled-apart pair not pulled apart, the corner turned in single
precision, the position added in single precision, the whole support back in single precision: all caught.
Which of two tied hull points a support returns is not, and need not be: either is right.

## Step 19: solver speed

After the narrow phase work of steps 15 to 18, the solver was the biggest phase in the crate scenes: 1.1 ms
of a 1.9 ms step at 1000 crates. Timed on its own on a snapshot of a settled pile (the same bodies and
contacts every run), 90% of that was the ten sweeps, about 30 ns per contact point per sweep, and arithmetic
turned out to be the smaller part of it.

Four changes, each measured on 1000 crates (solver alone):

| change | solver |
|---|---|
| step 18 | 1.18 ms |
| rows in Jacobian form: r x dir and I^-1 (r x dir) worked out once per step, so a sweep never multiplies by an inertia matrix or takes a cross product (as the joint rows already did) | 1.04 ms |
| each contact's two velocities read into locals at its start and written back at its end, so its rows work in registers and not through memory | 0.86 ms |
| ... and written back only to bodies that can move | 0.75 ms |
| contacts solved in colour order (below) | 0.56 ms |

The last two are the same discovery. Gauss-Seidel is a chain: each contact starts from the velocities the
last one left, and the processor cannot begin one until it knows them. Every contact with the floor used to
write the floor's velocity (adding zero to it) and the next floor contact had to wait to read it: a chain
through one memory location from one end of the pile to the other. And contacts sorted by body pair nearly
always share a body with the one before ((5, 7), (5, 9), (5, 12)). So `solving_order` colours them greedily:
each takes the lowest colour no earlier contact on either of its moving bodies has, and they are solved colour
by colour. Within a colour no two share a moving body, and the processor overlaps them. The joints get the
same: a joint's rows (which share its bodies) are solved in registers, and joints are coloured too, which more
than halves the joint sweeps of the chains scene (a chain's joints all share a link with the next).

The first three changes give the same answers to rounding (1.6 * 10^-7 m/s apart, on velocities of 0.1). The
colour order changes the order, and with ten sweeps that are not run to convergence the order shows: one step's
velocities differ by up to 0.07 m/s. What matters is whether piles settle as well, and they do:

| 1000 bodies, after 8 s | mean speed (old / new) | mean depth | deepest |
|---|---|---|---|
| crates | 8 * 10^-9 / 3 * 10^-7 m/s | 2.57 / 2.56 mm | 14.7 / 14.2 mm |
| lumps | 1.2 * 10^-5 / 2.3 * 10^-8 | 2.12 / 2.10 | 25.3 / 17.9 |
| 4000 crates | 1.5 * 10^-8 / 1.6 * 10^-8 | 2.58 / 2.57 | 15.6 / 15.5 |

**The joint position pass** (the "move" column of the chains scene) prepared each joint once more than it has
rows, and rebuilt its stand-in for the world each time: 0.89 to 0.68 ms at 1000 joints, the same corrections.

**The benchmark's chains scene was broken.** A chain is 7 m long and leans up to 55 degrees, and the room ended
2 m beyond the outermost chains: the corner chains started through the walls, were pushed out the far side
while their joints pulled them back, and their last links flew off at up to 276 m/s and hung metres from
their chains. Every chains number in the tables above (steps 15 to 18) timed that. The room now reaches 7 m
beyond them, and no joint opens more than 2 mm.

In the world, step 18 against step 19, both on the fixed chains scene:

| scene | solver | step |
|---|---|---|
| crates 1000 | 1.11 to 0.59 ms | 1.86 to 1.45 ms |
| lumps 1000 | 1.16 to 0.61 | 4.93 to 4.34 |
| chains 1000 | 0.58 to 0.25 (and position pass 0.95 to 0.74) | 1.90 to 1.38 |
| crates 2000 | 2.29 to 1.21 | 4.26 to 3.66 |

**Something this showed, not fixed here.** The broad phase reads slower in step 19 (2000 crates 0.94 to 1.41 ms)
for no reason of its own: its tree is only laid out afresh after a quarter of its leaves have been re-inserted,
so a pile that comes to rest keeps whatever tree it had, and a different history leaves a different tree. Laid
out afresh, the settled trees of either step take the same 123000 to 137000 node visits per step; as left, they
took 130000 to 296000 depending on luck.

Tests (1, 267 in total): three pendulums hung from the world move bit for bit as each does alone (the order
joints are solved in must not matter to bodies that share nothing, and joints to the world, all with first body
-1, must not be mixed up). One test changed: the broad-phase equivalence test counted awake contacts and needed
more than 10; with this step its pyramid is asleep by then, so it counts stored ones. Deliberate breakages: 17,
15 caught, one of them (joint rows grouped by their first body alone) only once the pendulum test was written
for it. Not caught: writing back to static bodies, and everything in one colour; both cost time and never
change a result.

## Step 20: when to lay the tree out afresh

Step 15's rule was: rebuild the broad phase's tree once a quarter of its leaves have been re-inserted. It never
looks at the tree. A pile that comes to rest a little short of a quarter keeps that tree for as long as it
rests, and settled piles were measured searching trees from 1.0 to 2.2 times as dear as fresh ones, by luck of
their histories alone (step 19).

Now the broad phase also watches its own search. Each search costs some number of node visits; divided by the
pairs it finds plus the leaves it has (what even a perfect tree pays for), that is a measure of the tree. The
first search after a rebuild sets the mark. When a search costs a tenth more than the mark, the tree is laid out
afresh before the next one. A scene that changes (a pile forming) moves the measure too; that costs a rebuild and
sets a new mark, nothing more. Why a tenth: a rebuild costs about a third of one search (0.15 ms against
0.45 ms for 2000 leaves), so one bought at a tenth's wear has paid for itself after three quiet steps.
Thresholds tried on the benchmark piles, 250 to 2000 bodies:

| rule | settled trees, times as dear as fresh | broad phase while 2000 crates / lumps settle (8 s) |
|---|---|---|
| a quarter re-inserted (step 15) | 1.00 to 2.16 | 1603 / 2069 ms |
| ... or a third dearer | 1.01 to 1.27 | 1078 / 1352 |
| ... or a fifth dearer | 1.00 to 1.17 | 1116 / 1310 |
| ... or a tenth dearer | 0.94 to 1.05 (chains, never still, 1.14) | 989 / 1347 |

(Under 1 means the tree as it stood was better than a fresh median split, which happens.) Settled, per step:

| scene | broad phase, step 19 | step 20 |
|---|---|---|
| crates 1000 | 0.36 ms | 0.28 ms |
| crates 2000 | 1.48 | 0.97 |
| lumps 2000 | 2.33 | 1.32 |
| mixed 2000 | 1.53 | 1.08 |

Two shapes of wear that one might expect turned out not to wear the tree at all: bodies creeping just past their
padding in a dense pile, and bodies re-inserted anywhere in a sparse scene. What wears it is how a pile forms.
`StepStats::broadphase_rebuilt` now says when a step rebuilt.

Tests (1, 268 in total): a tree whose searches have become dear is rebuilt once, and then left alone. (Its
search is made dearer without touching the tree: 512 boxes overlapping their neighbours by 2 cm all shrink by
2 cm a side, inside their padding, so the same work finds no pairs.) Deliberate breakages: 6; caught: the wear
never calling for a rebuild, the mark set only at the first rebuild, cost per leaf alone, and the re-insertion
rule dropped. Not caught: cost per pair alone, without the leaves. That only differs in a scene with no pairs,
where it is infinite; the leaves are there to keep the measure finite, and no scene was found where it shows.
A guard against a zero mark went the other way: its breakage survived because the mark can only be zero for a
tree with no leaves, so it was deleted.

## Step 21: the joint position pass

After step 19 the chains scene spent more time in the joint position pass (0.74 ms for 1000 joints) than in
the whole velocity solver. The pass corrects a joint one row at a time, measuring each afresh from where the
last correction left the bodies (step 7: corrected together from one measurement, independent rows
overshoot each other and never settle). Measuring meant calling `prepare` for the whole joint, every row of
it, once per row: twelve thousand times a step for a thousand ball joints. Profiled, about 30 ns a time, of
which turning the inertia into world space was only 4: the rest was measuring itself.

Two small savings first: `prepare` can be asked for just one row (`only_row`), and a correction no longer
turns a body that cannot move (it was renormalising the world stand-in's orientation). 0.67 to 0.62 ms.

Then the real one. A ball, hinge or weld holds its two anchor points together with three rows, along the
world axes, and those three are now corrected as one: measure the gap C = pb - pa once, and find the push P
that closes it from the 3 x 3 matrix that says how the gap answers a push,

    K = (1/ma + 1/mb) E - [ra] IA^-1 [ra] - [rb] IB^-1 [rb]        ([r] the matrix of r x),    P = -K^-1 C,

applied to both bodies equally and oppositely (so no centre of mass moves), at most 20 cm at a time as
before. This is not "correcting all together from one measurement" in step 7's sense: K couples the three
rows exactly, where the rows one at a time see each other only through what they have already moved.
Measured, two free bodies pulled 5 mm apart:

| joint | one pass, row by row | one pass, as one | four passes, row by row | four, as one |
|---|---|---|---|---|
| ball | 2.9% of the gap left | 0.31% | 2.4 * 10^-7 m | 8 * 10^-9 m |
| weld | 2.9% | 0.31% | 2.9 * 10^-5 m | 6 * 10^-8 m |
| hinge | 23% | 23% | 5.2 * 10^-5 m | 5.2 * 10^-5 m |

(A hinge's two turning rows come after its anchor and open it again as they turn the bodies back, whichever
way the anchor was closed.) In the chains scene, over seven seconds, the joints are held about twice as
tightly: mean gap 0.20 to 0.11 mm, 99th percentile 0.77 to 0.50 mm, worst 3.1 to 2.1 mm. And faster:

| chains | position pass and movement | step |
|---|---|---|
| 250 | 0.19 to 0.08 ms | 0.31 to 0.21 ms |
| 1000 | 0.74 to 0.34 | 1.37 to 0.98 |
| 2000 | 1.48 to 0.66 | 3.01 to 2.16 |

Tests (3, 271 in total): one correction closes a ball's or a weld's anchor gap to under 1% and moves no centre
of mass; an anchor a metre open closes 20 cm per correction; a sleeping joint is left alone even when moved
open. Deliberate breakages: 12, 9 caught (all three of the new tests were written for ones that at first were
not). Not caught: a weld's anchor corrected before its turning rows rather than after (an order), a hinge's
anchor corrected row by row as before (no difference measured, above), and corrections moving bodies that
cannot move (the same result).

## Step 22: a hinge's turning rows as a block (tried, and not kept)

Step 21 left a hinge's two turning rows one at a time, after its anchor, and they open the anchor again (23% of
a 5 mm gap left after one pass). The obvious next move was to solve more of a joint at once. Built and measured:
every row of a joint that holds both ways (all but limits, the cone and motors) corrected together, K = J M^-1 J^T
solved by Cholesky, with the one-sided rows after it one at a time, a fallback to rows one at a time when K is
singular, and the ball joint's closed-form anchor kept as a fast path. On two bodies it does what it promises:

| one pass, 5 mm gap | step 21 | all two-way rows as one |
|---|---|---|
| hinge | 23% left | 0.25% |
| weld | 0.31% | 0.18% |

and on a single hinge chain swinging alone it is the same as step 21. But in the benchmark's chains scene with
hinges (or welds) in place of ball joints, where chains strike each other, it holds the joints worse:

| 1000 links, seconds 1 to 8 | step 20 | step 21 | all as one | anchor, then the two turning rows as a 2 x 2 block |
|---|---|---|---|---|
| hinge chains: mean gap | 0.28 mm | 0.23 | 0.78 | 0.62 |
| hinge chains: joint-steps over 1 cm | 0 | 0 | 0.9% | 0.36% |
| weld chains: mean gap | 3.9 mm | 3.7 | 7.8 | |
| weld chains: worst | 5.5 cm | 4.2 | 20 | |

Re-closing the anchor after the block did not help (0.77 mm), so it is not the anchor: it is correcting a joint's
turning rows together. Every solve was well conditioned (no fallbacks in 3.8 million) and fewer than 0.2% of them
had an error past its limit. The likeliest reading is that the coupled solve turns a thin link through the
rotation that is cheapest for it, about its own long axis, and the next joint's axis turns with it; one row at a
time the turns are smaller and disturb the neighbours less. Whatever the reason, it is measured, so step 21's pass
stays. Step 21 itself was re-checked against step 20 on the same hinge and weld chains (above): as good or
better, at half the cost.

Two tests written for it are kept, because they hold for step 21's pass as well and cover what nothing did: an
anchor a metre open closes 20 cm per correction for hinges and welds too, not only balls; and a ball joint bent
0.6 rad past a 0.3 rad cone comes back 0.14 rad per correction and settles inside it (leaving the cone's error out
of the position pass was caught before only by a test of the row's derivative). 272 tests in total.

Also seen: welded chains are held poorly whatever the pass (a mean gap of 4 mm, the worst 4 cm, when ten links
welded end to end and leaning collide with their neighbours). That is the case for a different method (solving a
whole chain, not a joint, at once) and is left open.

## Step 23: joints solved exactly, a tree at a time

Step 22 ended with welded chains held poorly when they collide. Measured first, the cause is not the collisions:
one chain of ten links welded end to end and to the world, alone, sagged and opened by 2.3 mm on average and 14 mm
at worst. Gauss-Seidel passes a correction along a chain one joint per sweep, and ten sweeps cannot make a 7 m
welded rod rigid (160 sweeps: 0.2 mm). More position passes barely helped.

For joints that form a tree (chains, ragdolls: anything without a loop) the equations can be solved exactly in
time proportional to their number (Baraff, "Linear-Time Dynamics using Lagrange Multipliers", 1996): take the bodies'
velocity changes and the rows' impulses together,

    [ M  J^T ] [ dv ]   [ 0 ]
    [ J   0  ] [ mu ] = [ c ]        c = -(J v + bias),  the impulse is -mu,

a matrix with a block per body and per joint whose graph is the joints' own tree, and eliminate it leaves first:
each elimination changes only the parent's diagonal block, so nothing fills in. `JointTrees` (src/joint_trees.hpp)
finds the trees among the joints each step, factors them (in double precision; Cholesky for each block, which is
positive definite for a body and negative definite for a joint), and solves them in one pass down and one up.
Only the rows that hold both ways rigidly are taken: anchors, hinge and weld axes, slider rails. Limits, motors, the
cone and springs, and joints in a loop (a structure fixed to the world at two places is a loop through the world)
or on a body that cannot turn about some axis, are left to Gauss-Seidel as before.

Two things were tried on the way and are worth knowing:

- A first estimate, sweeping the joints 30 times within each sweep, made colliding welded chains blow up (links at
  96 m/s). The chains scene with welds in place of ball joints welds each chain to the world as well: rigid 7 m
  cantilevers, started at random leans, some through each other, with no physical way apart. Soft joints had hidden
  that by bending. The exact solve does not blow up even there, and on a fair version (welded rods hung from a ball
  joint, so they can swing into each other) neither method does.
- Solving the trees exactly in every sweep held the joints no better than solving them exactly in the last sweep
  only, with every row swept one at a time before it (the mean gap was the same; the worst moved by a millimetre or
  so either way), and cost 0.2 ms more at a thousand joints. So the last sweep does it.

What it costs was mostly not arithmetic. Clearing a kilobyte of node per joint and per body every step cost more
than factoring; the per-node loops of unknown length (3, 5 or 6 rows) were slow until each size had its own
unrolled copy; and the world now keeps its `JointSolver` between steps so that the trees' storage is kept too.

| 1000 links, seconds 1 to 8 | before (step 22) | trees |
|---|---|---|
| one welded cantilever: mean / worst gap | 2.3 / 14 mm | 0 / 0 (it does not move) |
| welded cantilevers colliding: mean / worst | 3.5 mm / 4.2 cm | 0.05 / 2.5 mm |
| welded rods on a ball joint, colliding: mean / worst | 0.88 mm / 3.3 cm | 0.10 / 5.1 mm |
| hinge chains colliding: mean / worst | 0.22 / 5.7 mm | 0.11 / 2.0 mm |
| ball chains colliding: mean / worst | 0.11 / 2.1 mm | 0.11 / 1.6 mm |

and in the benchmark (ball chains, which needed it least): solver 0.25 to 0.55 ms at 1000 links, step 0.96 to 1.26;
the crate scenes, without joints, unchanged.

Tests (5, 277 in total): a chain of every rigid joint, held or free, and a ragdoll-shaped tree are each solved
exactly in one sweep; a free tree keeps its momentum, linear and angular; a welded cantilever holds its shape; and a
ring, a bridge held at both ends, and a chain with a body that cannot turn about one axis are still solved, by
sweeping, and that body is never turned about that axis. Deliberate breakages: 14; caught: no trees, a body that
cannot turn taken into a tree (once that test turned the body about an axis a joint could actually push on), the
parent's block not reduced, either pass dropping its link to the parent, the impulse's sign, the mass matrix without
its inertia's off-diagonal, soft rows taken into trees, the exact solve never asked for, and tree rows swept and
skipped at once (which hung, also a catch). Not caught: the rows' bias left out of the right-hand side (every row
the trees take has none; kept, being the right formula), and the trees not interleaved by height, or solved in every
sweep (both cost time only). A check that a tree has at most one joint held to the world was deleted: a second one is
left with an empty diagonal block, which the factoring already refuses.

## Step 24: four contacts at a time

Contacts of one colour share no moving body (step 19), so they can be solved side by side. `WideConstraint` holds four
contacts, one per lane of a 4-wide vector (`float __attribute__((vector_size(16)))`, which GCC and Clang turn into NEON or
SSE), and does exactly the arithmetic `ContactConstraint` does, on four of everything: the bodies' velocities are
gathered into lanes, friction and then the normals are solved, and they are scattered back to the bodies that can move.
A lane without a contact, or a point slot a contact does not have, is given no mass and moves nothing. The last colour,
shared by whatever found the first 63 taken (a body with more than 63 contacts), may repeat a body within four, and is
still solved one contact at a time.

Four things it took, measured on the solver alone:

- Plain loops over four lanes were unrolled by the compiler but not vectorised (its cost model said not worth it), and the
  sweeps were no faster. Written with the vector type, they were 40% faster.
- Building each contact as a scalar constraint, keeping them all, and then copying them into lanes, cost more than the
  sweeps saved. Each is now prepared on the stack and dropped straight into its lane.
- A colour's order is free, so its contacts are sorted by number of points: four are solved for as many points as the
  most any of them has, and a sphere's single point beside a crate's four wasted three slots.
- `emplace_back` clears a new element, and a WideConstraint is three kilobytes; it now has an empty constructor, and
  packing clears only what is read.

| solver, ms per step | step 23 | step 24 |
|---|---|---|
| crates 250 | 0.138 | 0.087 |
| crates 1000 | 0.603 | 0.455 |
| crates 2000 | 1.262 | 0.924 |
| mixed 1000 | 0.566 | 0.481 |
| lumps 1000 | 0.617 | 0.504 |

Whole steps: crates 1000 1.40 to 1.24 ms, crates 2000 3.21 to 2.83; the chains scene (joints) unchanged. On a snapshot,
the velocities after a solve are those of step 23 to 2 * 10^-7.

Tests (2, 279 in total): a slab carrying 80 cubes (81 contacts on one body) holds them all level, and a heavy ball with 70
small ones closing on it in weightlessness keeps the total momentum, which it would not if the contacts sharing it beyond
the 63rd colour were solved four at a time (only the last lane's change to its velocity would count). Deliberate
breakages: 10; caught: the shared colour never solved, or solved in fours; friction not limited; normals that pull;
missing points not cleared (caught once the automatic clearing went); the warm start, or the kept impulses, dropped; b
pushed with a's inverse mass; b's spin not written back. Not caught: contacts not grouped by point count (time only).

## Step 25: a compact array of velocities (measured, not built)

After step 24 the contact solver's largest remaining cost is moving velocities between the bodies and the lanes: each
four contacts read and write eight bodies' velocity and spin, out of 232-byte bodies. Engines often solve on a compact
array of just those (32 bytes a body), copied in before the solve and out after. Doing that properly here means the
joint solver and the joint trees moving to the array too, so first an experiment on contacts alone (exact for scenes
without joints): it made the whole contact solve 5 to 9% faster (1000 crates: 0.40 to 0.38 ms), about 2% of a step.
What remains is the shuffling into lanes itself, which an array does not remove. Not worth two homes for every
velocity; not built.

## Step 26: joint trees in the position pass

Step 23 solved joint trees exactly for velocities; the position pass still corrected them a joint at a time. Step 22
had found that correcting a joint's rows together, joint by joint, made colliding chains worse. A whole tree at once is
different: each joint's correction is made knowing what it does to all the others. Now each pass first measures every
rigid two-way row of every tree (one `prepare` of all the joints), solves the same tree system as step 23 with the
errors on the right (each limited as a row's own correction would be), and moves and turns the bodies by the result.
Then, as before, one row at a time: the rows of joints in no tree, and the one-sided rows (limits, the cone) of joints
in one.

A tree correction is exact to first order; what it leaves is about the square of what it found, so a second pass finds
little unless the first found a lot. Measured in the chains scene: one tree pass always left gaps of 6 * 10^-7 m on
average but a few of 1 to 3 mm where a joint had been far off; four passes for every tree cost twice as much, and so did
a second pass for all of them whenever any needed one. So a tree gets another pass only if its own largest error was
over a millimetre (or a milliradian), and only those trees are measured again.

Seven seconds of chains striking each other, 1000 links:

| | step 24, mean / worst gap | step 26 |
|---|---|---|
| ball chains | 0.11 mm / 1.6 mm | 0.6 / 4.8 micrometres |
| hinge chains | 0.11 / 2.0 mm | 0.6 / 11 micrometres |
| welded cantilevers | 0.05 / 2.5 mm | 0.03 / 1.4 micrometres |
| welded rods hung from a ball joint | 0.10 / 4.7 mm | 0.5 / 3.9 micrometres |

It costs the chains benchmark 0.2 ms a step at 1000 links (movement 0.33 to 0.55 ms; step 1.26 to 1.51), measured on a
busy machine; scenes without joints are unchanged. (A ring of joints, a loop, is still corrected row by row, and closes
as slowly as before: a 2.4 cm gap is 4 mm after twenty passes.)

Tests (3, 282 in total): one pass corrects a twelve-link chain of every rigid joint, knocked 4 mm and 4 mrad out, to
10^-5 (row by row it left a tenth of the error); a chain knocked 5 cm out needs more than one pass and gets them; a ring
is still corrected, a little more each pass. Two first versions of these were wrong, not the code: a slider's anchors
may part along its rail, so the tests measure the rows' errors, not the anchors' gap. Deliberate breakages: 8; caught:
no tree correction, no second pass, tree joints' limits and cones left alone, errors not limited, bodies turned the
wrong way or not at all. Not caught: every tree given every pass, and tree joints' two-way rows corrected again row by
row on the way to their limits (both time only). An index one past the end, found by the chains measurement crashing,
not by a test (the tests passed over it), is fixed; the sanitizer build, which would have caught it, was run after.

## Step 27: loops of joints

Steps 23 and 26 solve joints exactly when they form a tree. A loop (a ring of joints, or anything fixed to the world in
two places, a rope bridge, which is a loop through the world) was still corrected a row at a time, and slowly: a ring
opened 2.4 cm took twenty position passes to get to 4 mm, and a 40-plank rope bridge swung with its joints open by 1.6
mm on average and 9 mm at worst.

Each loop is now cut at one joint and the cut put back exactly through the tree. Which joint: the joints are taken in
turn, each joining two sets of bodies (the world one more body); one that would join a set to itself closes a loop and
is cut. Joints between two moving bodies are taken before joints to the world, so that a structure fixed at both ends
is one tree fixed at one end and cut at the other (taken in any order, both ends could be fixed first and the cut fall
in the middle, between two trees, where neither can carry it). Putting the cut back: an impulse on the cut rows moves
the tree's bodies by A lambda, A found once per step by one tree solve per cut row with that row's unit impulse on the
right (and nothing asked of the tree's own rows, so they still hold); the cut rows' speeds then change by W lambda,
W = J_cut A. After each tree solve, lambda = W^-1 (what the cut rows still need), and A lambda is applied: for the
velocities in the last sweep, for the positions in every pass (less, there, what the tree's own correction already did
to them, to first order). A tree with more than 24 cut rows is left to sweeping.

W can be singular, and the first test found one: a straight bridge cannot be pulled along its length, so its cut row
along the bridge is implied by the others. The Cholesky factoring of W now drops such a row (it gets no impulse) and
solves the rest exactly. A joint given twice is the same case, every row of the second dropped.

| | step 26, mean / worst gap | step 27 |
|---|---|---|
| rope bridge, 12 planks | 0.24 / 2.0 mm | 1.1 / 5.7 micrometres |
| rope bridge, 40 planks | 1.6 / 9.0 mm | 1.2 / 5.7 micrometres |
| ring, one position pass from 2.4 cm open | (twenty passes: 4 mm) | 0.05 mm; two passes: 10^-7 m |

For a few hundredths of a millisecond. The chains benchmark, which has no loops, is unchanged within the noise of a busy
machine.

**A bookkeeping fault this found, in step 23.** A new test checks that each row's kept impulse is what it really did: from
a cold start, one sweep, then every body's velocity must be its start plus M^-1 J^T lambda over all the rows' impulses.
For a chain it was off by 4 * 10^-4: the trees' answer has a part for the bodies and a part for the rows, and step 23
took the bodies' velocities from the first and the impulses from the second, which agree only to single precision
(a few parts in ten thousand over light links). The velocity solves now give each row's impulse to its bodies through
the row itself, as Gauss-Seidel does, so velocities and kept impulses agree by construction, and momentum is kept exactly;
and the solves are done in double precision, without which the rows were then off instead (a chain solved "exactly" in one
sweep was 10^-4 out). About 0.03 ms more per exact solve at a thousand joints.

Tests (4 new, 2 rewritten, 286 in total): a ring and a bridge held at both ends are solved exactly in one sweep (the body
that cannot turn is still swept); a ring opened 2.4 cm is closed by one position pass to 0.2 mm and by two to a micrometre;
a rope bridge started as a V with slack sags and holds together within 0.1 mm; a resting rope bridge's two ends carry its
weight between them; the impulses kept are the impulses applied, for a tree, a ring and a bridge; and a joint given twice
is no trouble. (The first rope bridge was built straight and taut, which no chain can hold up: under its own weight its
tension is unbounded. That was the test.) Deliberate breakages: 12; caught: cuts never put back, fixed joints not taken
last, every row dropped, the cut's impulse pushed the wrong way, the tree's own correction not counted in the position
pass, the cut rows' impulses or the tree rows' changes not applied, an impulse kept without its velocity or the other way
round. Not caught: no row ever dropped (only an exactly singular row was met, which a pivot of exactly zero drops anyway;
the threshold is not pinned by a test), and cut joints not marked as done (time only).

## Step 28: demo scenes for the joints

Steps 21 to 27 changed how joints are solved, and had only been seen as numbers. Three scenes in `demo3d` show them:

- `8` a rope bridge: sixteen planks hinged end to end and hung with slack from posts at both ends (a loop through the
  world, step 27), started as a V; it swings into its sag and five crates are dropped on it.
- `9` ragdolls: three, pushed off the top of a flight of six steps. Each is one tree (step 23): neck, shoulders and
  hips are ball joints held inside cones, elbows and knees hinges with limits that bend one way.
- `0` welded structures: five rods of six to ten welded links, each hung from a ball joint and started leaning the
  same way, so that, swinging at different rates, they drift out of step and knock together; and a beam of ten welded
  links sticking out of a post with a heavy ball welded to its tip, which stays straight (row by row it sagged).

The text now also shows the largest gap between a joint's two anchors (balls, hinges, welds): in all three it stays
under a micrometre or so. Checked headless (`--capture`) at the start and after 4 and 10 seconds; the rods were first
started through each other, and one ragdoll fell off the back of the steps, both changed. Two things seen are not
faults: the bridge is still swinging, slowly, after ten seconds, and a ragdoll's head rocks back and forth on the
ground on its neck, because nothing in the engine damps a swing or a roll (there is no joint damping or rolling
resistance), so neither falls asleep for a long while.

## Step 29: joint damping

Nothing in the engine slowed a joint's swing: a ragdoll's head rocked on its neck and a bridge swayed for as long as one
cared to watch, and so neither fell asleep. `Joint::damping` now adds viscous damping to the motion a joint allows: a
torque of -damping times the relative angular velocity for a ball (about all three axes) or a hinge (about its axis), a
force of -damping times the sliding speed for a slider. It is taken implicitly, as a soft row with target speed zero and
softness gamma = 1 / (damping dt) (the form the distance spring already uses): such a row settles at lambda =
-(J v) / gamma = -damping dt (J v), the damping impulse over the step at the speed the step ends with, and so any amount
is stable (a great deal simply holds the joint still). Being soft, the rows are left to Gauss-Seidel by the trees, and
have nothing for the position pass to correct. A hinge and a slider already had eight rows' worth of impulse slots, so
there are now nine.

In the demo, the ragdolls' joints are damped (0.02 N m s at neck, shoulders and hips, 0.01 at elbows and knees: the parts
weigh grams at density 1) and they come to rest and fall asleep about six seconds after being pushed (twenty without:
damped three times less, ten); the rope bridge's hinges are damped (0.012), and it still sways for a long while: its
slowest swing carries much more than the hinges take out each time, and its peak speed comes down only slowly through
the sleep threshold. The welded rods are left undamped, swinging being the point of them.

Tests (6, 292 in total): a wheel on a damped hinge, a ball on a damped ball joint (spun about a tilted axis) and a carriage
on a damped slider each slow exactly as implicit damping says, w0 / (1 + c dt / I)^n; damping between two free bodies, by
a ball, a hinge or a slider anchored away from both centres, keeps their momentum and angular momentum; a motor against
damping settles where its torque runs out (w = torque / c); a damped pendulum loses its swing and stays on its pin.
Deliberate breakages: 7; caught: no damping rows, gamma upside down, rows rigid not soft, a ball damped about one axis,
a slider's row without its turning part, damping pushing both bodies the same way (the last two only once the momentum
test had joints between free bodies, anchored off centre; with one end on the world they changed nothing). Not caught: a
hinge's damping row in the motor's impulse slot (they share what is kept for the warm start; the motor test settles the
same either way).

**Seen while setting it up, and next.** The worst gap between a joint's anchors over a whole run is not the steady value
the demo shows: in the first seconds, while crates land on the bridge and ragdolls hit the steps, joints open by up to
2 cm on the bridge and 10 cm on a ragdoll, and close again. Two causes, measured: in the last sweep the contacts are solved
after the joints and have the last word (solving the joint trees after them instead takes the bridge's worst to 0.1 mm);
and continuous collision pulls a fast body back along its path without regard to its joints, after the joint position
pass (with it off, a ragdoll's worst was 2.6 cm, not 8).

## Planned next

The two orderings above.
