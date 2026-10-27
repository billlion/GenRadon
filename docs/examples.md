# Toy acquisition schemes

All examples use the cube `[-1,1]^D` and a nonnegative density made from two Gaussian
blobs sampled at voxel centers, then represented as piecewise constant values.
Units are arbitrary. Output data are bin integrals with natural coarea weight.

The 2D configurations use equally spaced angles on a circle of radius 3. In 3D,
Fibonacci directions give approximately uniform coverage of a sphere of radius 3.
There are no duplicate angular endpoints or repeated pole samples. These are
small, reproducible acquisition schemes for exercising a projector; they are not
claims of sampling optimality, complete data or unique reconstruction. The 3D
configurations deliberately sample a small subset of each full surface family.

For each exterior direction `n`, choose a unit tangent `t`, set
`a=3*n`, `b=3*(cos(pi/3)*n + sin(pi/3)*t)`, and `axis=-n`.
The focal separation is 3. In 2D the tangent is the positive 90-degree rotation;
in 3D it is the azimuthal tangent recorded through the resulting `b` coordinates.

| Case | Toy interpretation | Level bins |
|---|---|---|
| Cone / V-line | Source density projected onto emission-Compton cones with vertex `a` and inward axis | 48 uniform bins on [0.7,1.01), including a guard margin |
| Spheroid / ellipse | Reflectivity integrated over equal total path-length surfaces from `a` to `b` | 48 uniform bins on [2,10) |
| Hyperboloid / hyperbola | Source density integrated over equal path-length-difference surfaces for receivers `a,b` | 48 uniform bins on [-3.1,3.1) |

These bin ranges cover the reconstruction cube. The final cone bin extends beyond
the physical cosine limit 1 by 0.01. This guard margin contains the affine
surrogate's small overshoot near the axis: the selected tolerance and bin widths
give a maximum accepted level error below 0.01. The exact generating function
never exceeds 1, and there is no physical surface at levels above 1. Padding the
outer measurement bin permits a conservation check without renormalising the
operator; the overshoot shrinks with affine refinement. A nonuniform angular or time response can later be
represented by another bin grid or detector-response integration; the current
cone bins are specifically uniform in cosine, not angle or energy.

Wave speed is omitted, so path length rather than time is the level coordinate.
The spheroid example has no coherent phase, spreading factor or antenna pattern.
The cone example has no scattering cross-section, attenuation or detector
efficiency. Physical acquisition weights would generally modify the operator.

Run without `--output` for summaries only. With `--output DIRECTORY`, each case
has its own folder:

* `object.csv`: linear cell index, spatial coordinates, input density and `A^T A f`.
* `data.csv`: configuration index, bin index, edges and the bin-integrated data.
* `geometry.csv`: the source/receiver positions and cone axis for each configuration.
* `metadata.txt`: dimensions, grid shape, normalisation, integration settings and
  work counters.

The `transpose_of_data` field is an unnormalised backprojection, not a reconstructed
image. The program checks that summing each configuration's data reproduces the
voxelised object's mass. Midpoint and affine runs write the same layouts so their
forward data can be compared directly in Python or MATLAB before bindings exist.

Future examples should include properly specified Cormack curve families in 2D,
and spindle/apple/lemon surfaces in 3D motivated by Compton scattering work.
Those geometries may require branch masks and domains avoiding singularities;
they are not included as built-ins yet.
