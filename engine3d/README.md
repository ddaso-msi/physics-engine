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

## Planned next

Step 3: shapes (sphere, box) and 3D contact generation: sphere cases, then box-box by SAT over 15 axes with face
clipping and the edge-edge case. Then the sequential-impulse solver in 3D, a wireframe SDL demo, and GJK/EPA for
general convex shapes.
