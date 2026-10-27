# Geometry interface

The geometry is a C++ value type, copied into the projector's list of acquisition
configurations. Its evaluators are `const` and should be deterministic. There is
no required base class or virtual dispatch per point.

```cpp
struct MyGeometry {
    double value(const genradon::Vec<2>& x) const;
    genradon::Jet<2> jet(const genradon::Vec<2>& x) const; // optional for midpoint
    double hessian_bound(const genradon::Box<2>& region) const; // preferred
};
```

`jet` returns derivatives with respect to **spatial coordinates**, holding the
acquisition parameters fixed. Its Hessian is a full symmetric matrix. The bound
method returns an upper bound on the spatial Hessian operator norm throughout the
entire region. A pointwise Hessian norm is not a regional bound.

A minimal custom geometry:

```cpp
struct CircleSquaredRadius {
    genradon::Vec<2> center;
    double value(const genradon::Vec<2>& x) const {
        auto r = genradon::subtract(x, center);
        return genradon::dot(r, r);
    }
    genradon::Jet<2> jet(const genradon::Vec<2>& x) const {
        auto r = genradon::subtract(x, center);
        return {value(x), {2*r[0], 2*r[1]}, {{{2, 0}, {0, 2}}}};
    }
    double hessian_bound(const genradon::Box<2>&) const { return 2; }
};
```

Its level coordinate is **squared radius**. Changing this to radius changes the
natural weight and the interpretation of uniform bins. A generating function is
part of the measurement specification, not just a representation of its zero set.

## Built-ins

| Type | Generating function | Parameters |
|---|---|---|
| `Plane<D>` | `n.x - offset`, with normal normalised to unit length | normal, offset |
| `Spheroid<D>` | `|x-a| + |x-b|` | two foci |
| `Hyperboloid<D>` | `|x-a| - |x-b|` | two distinct foci |
| `Cone<D>` | `axis.(x-vertex) / |x-vertex|`, unit axis | vertex, axis |

`D` is 2 or 3. `Spheroid<2>` gives ellipses. `Hyperboloid` describes the signed
distance-difference family: a nonzero admissible level selects one sheet; zero
gives the perpendicular bisector. This is **not** a one-sheet hyperboloid family.
`Cone<2>` gives V-lines. Cone levels are cosines of opening angles; both signs are
available, with no squaring that would silently join opposite cones.

The focal Hessians use `(I - nn^T)/r`, with signs chosen for sums or differences.
The norm bounds use the minimum distance from each focus to the integration box.
The cone uses the analytic gradient and Hessian, with a conservative bound
`sqrt(2) / minimum_distance_to_vertex^2`. To see this, write `c=axis.n`,
`t=axis-c*n` and `P=I-n*n^T`; the Hessian times squared radius is
`-(t*n^T+n*t^T)-c*P`, whose norm is at most `|t|+|c| <= sqrt(2)`.
The bound is infinite when the box contains
the singular point.

## Symbolic generation

`tools/generate_spheroid.py` illustrates the build-time route:

1. Define `phi(x; parameters)` symbolically with real coordinates and parameters.
2. Differentiate only in `x`.
3. Eliminate common subexpressions jointly across value, gradient and the six
   independent Hessian entries.
4. Emit C++ into an ordinary geometry class.

The generated spheroid class delegates its cheap value evaluation and regional
bound to the analytic implementation; its fused jet is generated. The bound is
not inferred automatically from the symbolic Hessian. Focal coordinates remain
runtime arguments. An optional CMake target compiles and validates this example.

The production built-ins are hand-written and remain usable without SymPy.
Future generated geometries should retain explicit domain assumptions, singularity
handling and any necessary branch restrictions. A full user-expression compiler
and Python/MATLAB construction interfaces are subsequent work.
