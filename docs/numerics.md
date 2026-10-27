# Numerical conventions

## Discrete quantities

`Grid<D>` holds piecewise constant **density values**, not voxel masses. Cell
volumes are included in the forward coefficients. Physical coordinates are given
by the lower and upper grid corners; no unit-voxel assumption is made.

Coordinate zero is fastest: in 3D, `(i,j,k)` maps to `i + nx*(j + ny*k)`.
Data are configuration-major: `(view,bin)` maps to `view * number_of_bins + bin`.
Bindings must preserve these conventions explicitly, including NumPy/MATLAB
array strides. Bins can be nonuniform; each is half open `[lower,upper)`, including
the last bin. Values outside the supplied bin range contribute zero. No
renormalisation compensates for truncated acquisition ranges.

For integrated data, `A_ki` is the voxel-band intersection volume. For averaged
data, it is that volume divided by the bin width. The spatial density remains
constant within a voxel even when its integration geometry is subdivided.

## Adjointness

`backproject(y)` is `A^T y` in ordinary Euclidean array inner products. Both
directions invoke exactly the same coefficient generator and geometry-dependent
subdivision. They never adapt to image or data values.

`weighted_adjoint(y, data_weights)` returns `A^T (data_weights * y) / voxel_volume`.
For **bin-average** data and acquisition quadrature weights `q_view`, choose
`data_weights[view,bin] = q_view * bin_width[bin]`. This approximates the adjoint
with respect to spatial volume and `dtheta ds` measures. For **bin-integrated**
data interpreted as a piecewise constant data density, the corresponding data
inner-product weights are `q_view / bin_width[bin]`.

The choice of data inner product is explicit; a statistical noise precision matrix
or a least-squares weighting need not equal these geometric quadrature weights.

## Midpoint quadrature

`points_per_axis = m` places `m^D` equally weighted composite midpoint samples in
each voxel. The total sample weight is its volume. With full bin coverage, the
sum over bins conserves object mass for each configuration, up to roundoff.
Increasing `m` improves resolution of cut regions but does not provide the
high-order accuracy of Gaussian quadrature on a smooth integrand: the band
indicator is discontinuous. A thin band can miss every sample.

## Adaptive affine integration

For a subcell centered at `c`, approximate

\[
\phi(c+u)=\phi(c)+g\cdot u+R(u).
\]

If `M` bounds the Hessian operator norm throughout the subcell and `h` is its
vector of half-widths, then `|R| <= E = M |h|^2 / 2`. Accept the approximation when
`E <= affine_tolerance * smallest_bin_width`. Otherwise bisect the longest spatial
edge and recurse. These subcells currently start at each voxel; sharing one affine
approximation across several voxels is a later optimisation.

The tolerance controls **error in the level coordinate**, not a guaranteed
relative error in each matrix coefficient or reconstructed image. A coefficient
near a grazing intersection can be tiny. Tightening the tolerance should be
assessed through a convergence study appropriate to the acquisition.

An affine surrogate can extend slightly beyond the true level range. For example,
near a cone axis its value can exceed cosine 1. At a tightly truncated bin range
this can lose a small amount of mass even if the exact range is covered. The toy
acquisitions pad their outer cone bin by more than the accepted level error.
Neither the generic projector nor backprojector silently clamps or renormalises
the surrogate. Full-coverage conservation tests must cover the approximate range
as well as the exact one.

The implementation intersects a normalised square or cube with the two affine
half-spaces. It computes polygon areas or convex polyhedron volumes directly;
it does not subtract two almost equal cumulative volumes. Clipping uses extended
precision where `long double` provides it. Geometry remains approximate for
nonlinear `phi`; ordinary floating point roundoff still applies.

Built-in geometries provide conservative analytic Hessian norm bounds. These are
mathematical bounds evaluated in floating point, not interval-certified bounds
including every rounding error. Before excluding a subcell, the code includes
the Taylor remainder in its level range. Same-sign corner values alone are not
used to exclude nonlinear intersections.

For a custom geometry without a regional bound, users must explicitly set
`allow_estimated_bounds = true`. The Frobenius norm of the center Hessian is then
used as a heuristic; it can miss variation elsewhere in the subcell. Statistics
report the number of regions treated this way. This mode offers no certified
remainder or intersection-detection guarantee.

Reaching `max_depth` before the tolerance is met raises an error. There is no silent
switch to a lower-accuracy approximation. Exceptions from invalid or nonfinite
geometry evaluations are propagated. The adaptive recursion is bounded in depth,
but an excessively tight tolerance may still require considerable work.

## Singularities and nonregular levels

The cone generating function is undefined at the vertex; focal distance
derivatives are undefined at a focus. The examples keep these points outside the
closed reconstruction box. A custom geometry must be valid on its integration
domain. Regions containing built-in singularities have an infinite Hessian bound
and eventually trigger the depth limit if they cannot be excluded from the domain.

For distance sums, the shortest level degenerates to the focal segment; distance
differences have limiting degenerate levels. The band-integral model still has a
volume interpretation where `phi` is defined, but the regular surface/coarea
interpretation is asserted only where its spatial gradient is nonzero. Cone
levels `s=+-1` are degenerate as well. The core does not test Bolker conditions,
injectivity, or completeness of acquisition.

## Current computational scope

The core is serial, header-only C++17, and matrix-free. It recomputes coefficients
in forward and transpose applications instead of caching a potentially large
system matrix. The 3D clipping implementation prioritises clarity and validation;
its small polygon allocations and repeated per-voxel geometry evaluations are
future optimisation targets. No speed advantage over a specialised projector is
claimed from these initial tests.

Physical weights, complex amplitudes and intersections of two level functions are
not implemented in this version. A later weighted implementation must integrate
the weight over each cut region and conjugate complex coefficients in its adjoint.

## Verification

The standalone C++ test checks:

* 2D and 3D affine volumes against an independent inclusion-exclusion distribution
  formula, with positive and negative gradient components;
* narrow slabs, grazing simplex volumes, positivity and partition of voxel volume;
* gradients and Hessians against finite differences, and sampled directional
  checks of the analytic Hessian bounds;
* Euclidean and weighted adjoint identities for all built-in geometries in both
  dimensions and with both integration and normalisation choices;
* mass conservation on full acquisition ranges and nonuniform bins;
* convergence to analytic annulus and spherical-shell volumes for a custom
  squared-radius generating function;
* indexing, invalid inputs, value-only geometries, explicit heuristic-bound
  selection, and failure on an exhausted subdivision depth.

The optional generated-code test compares SymPy-generated spheroid jets with the
hand-written analytic implementation. Toy acquisition tests exercise the six
examples and check conservation. These verify the projector, not an inverse
problem solution or a complete physical measurement model.
