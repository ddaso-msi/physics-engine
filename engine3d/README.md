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

## Planned next

Step 2: rigid body state in 3D (`Body` with `Vec3` velocities, `Quat` orientation, world-frame inverse inertia,
the gyroscopic term for torque-free spin). Step 3: shapes (sphere, box) and 3D contact generation: sphere cases,
then box-box by SAT over 15 axes with face clipping and the edge-edge case. Then the sequential-impulse solver in
3D, a wireframe SDL demo, and GJK/EPA for general convex shapes.
