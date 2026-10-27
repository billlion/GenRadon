# GenRadon

A small C++17 library for forward and adjoint calculations of binned generalized
Radon transforms defined by generating functions in two and three dimensions.

This first implementation concentrates on **cones, spheroids and hyperboloids of
revolution**, with ordinary planes as a validation case. User-defined C++ types
can supply other generating functions. The core is header-only, uses the C++
standard library, and has no Python or MATLAB runtime dependency.

## Mathematical convention

For a configuration `theta` and a scalar level coordinate `s`, the continuous model is

\[
T f(\theta,s)=\int_\Omega f(x)\,\delta(s-\phi(x,\theta))\,dx.
\]

On regular level sets this is the surface integral of
`f / |gradient_x phi|`. For a piecewise constant voxel density `f_i`, bin-integrated
data have coefficients

\[
A_{(\theta,k),i}
=\left|V_i\cap\{s_k^-\leq\phi(x,\theta)<s_k^+\}\right|.
\]

The default is **bin integrals with the natural coarea weight**. Bin averages are
also available. `backproject` is exactly the Euclidean transpose of the chosen
discrete projector, up to roundoff; `weighted_adjoint` explicitly implements
`M_X^{-1} A^T M_Y`. Backprojection alone is not an inverse.

See [the numerical conventions](docs/numerics.md) for indexing, inner products,
error tolerances, singularities and the limitations of this initial version.

## Build and test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
```

The default build needs only CMake 3.16 or later and a C++17 compiler. Without CMake:

```sh
c++ -std=c++17 -O2 -Iinclude tests/test_core.cpp -o genradon_tests
./genradon_tests
c++ -std=c++17 -O2 -Iinclude examples/toy_transforms.cpp -o genradon_toy
./genradon_toy --smoke
```

To use an installed CMake package:

```sh
cmake --install build --prefix /your/install/prefix
```

Then use `find_package(GenRadon CONFIG REQUIRED)` and link the interface target
`GenRadon::genradon`. Alternatively, add `include/` to your include path.

## Small example

```cpp
#include <genradon/genradon.hpp>

int main() {
    using namespace genradon;
    Grid<3> grid({12, 12, 12}, {-1, -1, -1}, {1, 1, 1});
    auto bins = Bins::uniform(2.0, 10.0, 48);
    std::vector<Spheroid<3>> views{
        Spheroid<3>({-3, 0, 0}, {0, 3, 0}),
        Spheroid<3>({0, -3, 0}, {3, 0, 0})
    };
    Options options;
    options.integration = Integration::adaptive_affine;
    options.affine_tolerance = 0.05;
    Projector<3, Spheroid<3>> projector(grid, bins, views, options);
    std::vector<double> image(grid.size(), 1.0);
    auto data = projector.project(image);
    auto backprojection = projector.backproject(data);
}
```

Two integration methods share the coefficient-generation path:

* **Composite midpoint quadrature:** one or more points per coordinate per voxel;
  only function evaluations are required.
* **Adaptive affine integration:** use gradients and Hessian bounds to subdivide
  voxels, then clip each subcell by two planes. The slab volume is calculated
  directly, including narrow bins and oblique cuts.

## Toy acquisitions

```sh
./build/genradon_toy --output output-midpoint
./build/genradon_toy --affine --output output-affine
```

For a multi-configuration build, the executable may be in `build/Release/`.
Each command runs six cases: cone/V-line, spheroid/ellipse, and
hyperboloid/hyperbola transforms in 2D and 3D. The default examples use a `24 x 24`
grid with 64 configurations in 2D, and a `12 x 12 x 12` grid with 48 configurations
in 3D, with 48 level bins. `--smoke` selects smaller cases for quick checks.

Configuration positions lie on an exterior circle or Fibonacci sphere. Outputs
contain the phantom, forward data, transpose backprojection, acquisition
parameters and metadata. See [the acquisition descriptions](docs/examples.md).
They are geometrical toy versions of emission Compton, travel-time and
time-difference problems, not complete detector or coherent radar models.

## Custom geometries and generated derivatives

A geometry is a C++ type providing `double value(const Vec<D>&) const`.
Adaptive affine integration additionally uses `Jet<D> jet(...) const`, containing
the value, spatial gradient and Hessian, and preferably
`double hessian_bound(const Box<D>&) const`.

Each geometry object contains one configuration's parameters. The same compiled
type therefore handles different source positions without recompilation.
See [the geometry interface](docs/geometry.md).

An optional SymPy example generates the spheroid derivative code:

```sh
python -m pip install sympy
python tools/generate_spheroid.py --output generated/spheroid3.hpp
```

To generate, compile and compare that code against the analytic implementation:

```sh
cmake -S . -B build-codegen -DGENRADON_TEST_CODEGEN=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-codegen --parallel
ctest --test-dir build-codegen --output-on-failure
```

This is a demonstration of build-time generation, not yet a general symbolic
expression frontend. The default build does not need SymPy or Python.

## Scope and next steps

The current core supports **one scalar level function**, real voxel densities,
Cartesian grids, nonuniform level bins shared across configurations, and serial
matrix-free application. It does not yet implement physical amplitude weights,
complex phasors, two-level-function intersections, Python or MATLAB bindings,
GPU execution, or an inversion algorithm.

Next development targets are Python/NumPy and MATLAB/MEX bindings; reuse of affine
approximations across larger voxel blocks; weighted transforms; and further toy
families including Cormack curves and Compton spindle/apple/lemon surfaces.
Their branch and parameter conventions will be specified before implementation.

## Background

* V. P. Palamodov, *A uniform reconstruction formula in integral geometry*,
  Inverse Problems 28 (2012), 065014.
  [Paper](https://arxiv.org/abs/1111.6514).
* R. I. Saye, *High-order quadrature methods for implicitly defined surfaces and
  volumes in hyperrectangles*, SIAM J. Sci. Comput. 37 (2015), A993–A1019.
  [DOI](https://doi.org/10.1137/140966290).
* J. W. Webber and E. T. Quinto, *Microlocal properties of seven-dimensional lemon
  and apple Radon transforms with applications in Compton scattering tomography*,
  Inverse Problems 38 (2022), 064001.
  [Paper](https://arxiv.org/abs/2201.03793).
* J. W. Webber, S. Holman and E. T. Quinto, *Surface of revolution Radon transforms
  with centers on generalized surfaces in R^n*.
  [Paper](https://arxiv.org/abs/2312.15635).

These references motivate the framework and future examples. The numerical
projector does not require, or establish, the hypotheses of an inversion theorem.
