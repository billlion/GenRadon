#pragma once

#include "types.hpp"

namespace genradon {

namespace detail {
template <std::size_t D> Jet<D> distance_jet(const Vec<D>& x, const Vec<D>& focus) {
    const auto r = subtract(x, focus);
    const double length = norm(r);
    if (!(length > 0) || !std::isfinite(length))
        throw std::domain_error("Distance derivatives are undefined at the focus");
    Jet<D> j;
    j.value = length;
    for (std::size_t d = 0; d < D; ++d) j.gradient[d] = r[d] / length;
    for (std::size_t d = 0; d < D; ++d)
        for (std::size_t e = 0; e < D; ++e)
            j.hessian[d][e] = ((d == e ? 1.0 : 0.0) - j.gradient[d] * j.gradient[e]) / length;
    return j;
}
}

template <std::size_t D> class Plane {
public:
    explicit Plane(Vec<D> normal, double offset = 0) : n_(unit(normal)), offset_(offset) {
        if (!std::isfinite(offset)) throw std::invalid_argument("Nonfinite plane offset");
    }
    double value(const Vec<D>& x) const { return dot(n_, x) - offset_; }
    Jet<D> jet(const Vec<D>& x) const { return {value(x), n_, {}}; }
    double hessian_bound(const Box<D>&) const { return 0; }
private:
    Vec<D> n_;
    double offset_;
};

// Sign +1: ellipses / prolate spheroids. Sign -1: signed distance-difference
// hyperbolae / hyperboloids (one sheet per nonzero signed level).
template <std::size_t D, int Sign> class FocalGeometry {
    static_assert(Sign == 1 || Sign == -1, "Invalid focal sign");
public:
    FocalGeometry(Vec<D> a, Vec<D> b) : a_(a), b_(b) {
        require_finite(a);
        require_finite(b);
        if constexpr (Sign == -1) {
            if (norm(subtract(a, b)) == 0)
                throw std::invalid_argument("Distance-difference foci must be distinct");
        }
    }
    double value(const Vec<D>& x) const {
        return norm(subtract(x, a_)) + Sign * norm(subtract(x, b_));
    }
    Jet<D> jet(const Vec<D>& x) const {
        auto a = detail::distance_jet(x, a_);
        const auto b = detail::distance_jet(x, b_);
        a.value += Sign * b.value;
        for (std::size_t d = 0; d < D; ++d) {
            a.gradient[d] += Sign * b.gradient[d];
            for (std::size_t e = 0; e < D; ++e) a.hessian[d][e] += Sign * b.hessian[d][e];
        }
        return a;
    }
    double hessian_bound(const Box<D>& box) const {
        const double da = distance_to_box(a_, box), db = distance_to_box(b_, box);
        if (da == 0 || db == 0) return std::numeric_limits<double>::infinity();
        return 1 / da + 1 / db;
    }
private:
    Vec<D> a_, b_;
};

template <std::size_t D> using Spheroid = FocalGeometry<D, 1>;
template <std::size_t D> using Hyperboloid = FocalGeometry<D, -1>;

// Level s = cos(opening angle). In 2D the level is a pair of rays (a V-line).
template <std::size_t D> class Cone {
public:
    Cone(Vec<D> vertex, Vec<D> axis) : vertex_(vertex), axis_(unit(axis)) {
        require_finite(vertex);
    }
    double value(const Vec<D>& x) const {
        auto r = subtract(x, vertex_);
        const double length = norm(r);
        if (!(length > 0) || !std::isfinite(length))
            throw std::domain_error("Cone generating function is undefined at its vertex");
        for (double& v : r) v /= length;
        return std::clamp(dot(axis_, r), -1.0, 1.0);
    }
    Jet<D> jet(const Vec<D>& x) const {
        auto n = subtract(x, vertex_);
        const double length = norm(n);
        if (!(length > 0) || !std::isfinite(length))
            throw std::domain_error("Cone derivatives are undefined at its vertex");
        for (double& v : n) v /= length;
        Jet<D> j;
        j.value = std::clamp(dot(axis_, n), -1.0, 1.0);
        for (std::size_t d = 0; d < D; ++d) {
            j.gradient[d] = (axis_[d] - j.value * n[d]) / length;
            for (std::size_t e = 0; e < D; ++e)
                j.hessian[d][e] = (-axis_[d] * n[e] - n[d] * axis_[e]
                    - (d == e ? j.value : 0.0) + 3 * j.value * n[d] * n[e]) / (length * length);
        }
        return j;
    }
    double hessian_bound(const Box<D>& box) const {
        const double d = distance_to_box(vertex_, box);
        // Let c=axis.n, t=axis-c*n, P=I-n*n^T. Then r^2 H is
        // -(t*n^T+n*t^T)-c*P, whose norm is <= |t|+|c| <= sqrt(2).
        return d > 0 ? std::sqrt(2.0) / (d * d) : std::numeric_limits<double>::infinity();
    }
private:
    Vec<D> vertex_, axis_;
};

} // namespace genradon
