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

## Planned next

A 3D broad phase (the world tests every pair, so it is limited to a few hundred bodies), GJK/EPA for general
convex shapes, and joints, sleeping and continuous collision ported from the 2D designs.
