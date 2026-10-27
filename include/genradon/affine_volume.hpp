#pragma once

#include "types.hpp"

namespace genradon {
namespace detail {

using Point = std::array<long double, 3>;
using Polygon = std::vector<Point>;

inline Point minus(const Point& a, const Point& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}
inline Point cross(const Point& a, const Point& b) {
    return {a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]};
}
inline long double inner(const Point& a, const Point& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
inline bool same_point(const Point& a, const Point& b) {
    constexpr long double eps = 64 * std::numeric_limits<long double>::epsilon();
    return std::abs(a[0]-b[0]) <= eps && std::abs(a[1]-b[1]) <= eps && std::abs(a[2]-b[2]) <= eps;
}
inline void append_unique(Polygon& p, const Point& x) {
    for (const auto& y : p) if (same_point(x, y)) return;
    p.push_back(x);
}

// Keep n.x <= t. Coordinates are scaled to [-1,1]^D before clipping.
inline Polygon clip_polygon(const Polygon& p, const Point& n, long double t, Polygon* cap = nullptr) {
    Polygon out;
    if (p.empty()) return out;
    auto a = p.back();
    long double da = inner(n, a) - t;
    for (const auto& b : p) {
        const long double db = inner(n, b) - t;
        if ((da <= 0) != (db <= 0)) {
            const long double fraction = da / (da - db);
            Point v{};
            for (std::size_t d = 0; d < 3; ++d) v[d] = a[d] + fraction * (b[d] - a[d]);
            append_unique(out, v);
            if (cap) append_unique(*cap, v);
        }
        if (db <= 0) append_unique(out, b);
        a = b;
        da = db;
    }
    return out;
}

inline void clip_polyhedron(std::vector<Polygon>& faces, const Point& n, long double t) {
    Polygon cap;
    std::vector<Polygon> kept;
    for (const auto& face : faces) {
        auto clipped = clip_polygon(face, n, t, &cap);
        if (clipped.size() >= 3) kept.push_back(std::move(clipped));
    }
    if (cap.size() >= 3) {
        Point center{};
        for (const auto& p : cap)
            for (std::size_t d = 0; d < 3; ++d) center[d] += p[d] / cap.size();
        std::size_t axis = 0;
        for (std::size_t d = 1; d < 3; ++d) if (std::abs(n[d]) < std::abs(n[axis])) axis = d;
        Point basis{};
        basis[axis] = 1;
        auto u = cross(n, basis);
        const auto length = std::sqrt(inner(u, u));
        for (auto& value : u) value /= length;
        const auto v = cross(n, u);
        std::sort(cap.begin(), cap.end(), [&](const Point& a, const Point& b) {
            const auto da = minus(a, center), db = minus(b, center);
            return std::atan2(inner(da, v), inner(da, u)) < std::atan2(inner(db, v), inner(db, u));
        });
        kept.push_back(std::move(cap));
    }
    faces = std::move(kept);
}

inline long double polyhedron_volume(const std::vector<Polygon>& faces) {
    Point center{};
    std::size_t count = 0;
    for (const auto& face : faces) for (const auto& p : face) {
        for (std::size_t d = 0; d < 3; ++d) center[d] += p[d];
        ++count;
    }
    if (count == 0) return 0;
    for (auto& value : center) value /= count;
    long double volume = 0;
    // The mean of all face vertices lies inside this convex polyhedron.
    // Absolute tetrahedron volumes avoid relying on face orientation.
    for (const auto& f : faces)
        for (std::size_t j = 1; j + 1 < f.size(); ++j)
            volume += std::abs(inner(minus(f[0], center),
                cross(minus(f[j], center), minus(f[j+1], center)))) / 6;
    return volume;
}
}

// Volume of lo <= center_value + gradient.(x-box.center) < hi in box.
// Directly clip by both planes instead of subtracting nearly equal cumulative
// volumes. This is exact for affine geometry, up to floating point roundoff.
template <std::size_t D>
double affine_band_volume(const Box<D>& box, double center_value, const Vec<D>& gradient,
                          double lo, double hi) {
    require_finite(box.center);
    require_finite(box.half_width);
    require_finite(gradient);
    if (!std::isfinite(center_value) || !std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo))
        throw std::invalid_argument("Invalid affine band");
    for (double h : box.half_width)
        if (!(h > 0)) throw std::invalid_argument("Box half widths must be positive");
    detail::Point n{};
    long double span = 0;
    for (std::size_t d = 0; d < D; ++d) {
        n[d] = static_cast<long double>(gradient[d]) * box.half_width[d];
        span += std::abs(n[d]);
    }
    if (span == 0) return lo <= center_value && center_value < hi ? box.volume() : 0;
    long double lower = static_cast<long double>(lo) - center_value;
    long double upper = static_cast<long double>(hi) - center_value;
    if (upper <= -span || lower >= span) return 0;
    if (lower <= -span && upper >= span) return box.volume();
    // Scale the plane equation to improve conditioning of cap construction.
    for (auto& v : n) v /= span;
    lower /= span;
    upper /= span;
    detail::Point negative{-n[0], -n[1], -n[2]};
    long double fraction = 0;
    if constexpr (D == 2) {
        detail::Polygon p{{-1,-1,0}, {1,-1,0}, {1,1,0}, {-1,1,0}};
        if (upper < 1) p = detail::clip_polygon(p, n, upper);
        if (lower > -1) p = detail::clip_polygon(p, negative, -lower);
        long double area = 0;
        for (std::size_t j = 1; j + 1 < p.size(); ++j)
            area += std::abs(detail::cross(detail::minus(p[j], p[0]), detail::minus(p[j+1], p[0]))[2]) / 2;
        fraction = area / 4;
    } else {
        std::vector<detail::Polygon> faces{
            {{-1,-1,-1},{-1,1,-1},{1,1,-1},{1,-1,-1}},
            {{-1,-1,1},{1,-1,1},{1,1,1},{-1,1,1}},
            {{-1,-1,-1},{1,-1,-1},{1,-1,1},{-1,-1,1}},
            {{-1,1,-1},{-1,1,1},{1,1,1},{1,1,-1}},
            {{-1,-1,-1},{-1,-1,1},{-1,1,1},{-1,1,-1}},
            {{1,-1,-1},{1,1,-1},{1,1,1},{1,-1,1}}};
        if (upper < 1) detail::clip_polyhedron(faces, n, upper);
        if (lower > -1) detail::clip_polyhedron(faces, negative, -lower);
        fraction = detail::polyhedron_volume(faces) / 8;
    }
    return box.volume() * static_cast<double>(std::clamp(fraction, 0.0L, 1.0L));
}

} // namespace genradon
