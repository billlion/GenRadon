#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace genradon {

template <std::size_t D> using Vec = std::array<double, D>;
template <std::size_t D> using Matrix = std::array<Vec<D>, D>;

template <std::size_t D> struct Jet {
    double value{};
    Vec<D> gradient{};
    Matrix<D> hessian{};
};

template <std::size_t D> double dot(const Vec<D>& a, const Vec<D>& b) {
    double result = 0;
    for (std::size_t d = 0; d < D; ++d) result += a[d] * b[d];
    return result;
}

template <std::size_t D> double norm(const Vec<D>& a) {
    double result = 0;
    for (double v : a) result = std::hypot(result, v);
    return result;
}

template <std::size_t D> Vec<D> subtract(const Vec<D>& a, const Vec<D>& b) {
    Vec<D> result{};
    for (std::size_t d = 0; d < D; ++d) result[d] = a[d] - b[d];
    return result;
}

template <std::size_t D> void require_finite(const Vec<D>& x) {
    for (double v : x)
        if (!std::isfinite(v)) throw std::invalid_argument("Coordinates must be finite");
}

template <std::size_t D> Vec<D> unit(Vec<D> x) {
    require_finite(x);
    const double length = norm(x);
    if (!(length > 0) || !std::isfinite(length))
        throw std::invalid_argument("Direction must have finite nonzero length");
    for (double& v : x) v /= length;
    return x;
}

template <std::size_t D> struct Box {
    static_assert(D == 2 || D == 3, "GenRadon supports 2D and 3D");
    Vec<D> center{};
    Vec<D> half_width{};

    double volume() const {
        double v = 1;
        for (double h : half_width) v *= 2 * h;
        return v;
    }
};

template <std::size_t D> double distance_to_box(const Vec<D>& p, const Box<D>& b) {
    Vec<D> delta{};
    for (std::size_t d = 0; d < D; ++d)
        delta[d] = std::max(0.0, std::abs(p[d] - b.center[d]) - b.half_width[d]);
    return norm(delta);
}

// Cartesian cells; x (coordinate zero) is the fastest varying index.
template <std::size_t D> class Grid {
public:
    Grid(std::array<std::size_t, D> shape, Vec<D> lower, Vec<D> upper)
        : shape_(shape), lower_(lower), upper_(upper) {
        require_finite(lower);
        require_finite(upper);
        size_ = 1;
        for (std::size_t d = 0; d < D; ++d) {
            if (shape[d] == 0 || !(upper[d] > lower[d]))
                throw std::invalid_argument("Grid needs positive shape and extent");
            if (shape[d] > std::numeric_limits<std::size_t>::max() / size_)
                throw std::overflow_error("Grid size overflow");
            size_ *= shape[d];
            spacing_[d] = (upper[d] - lower[d]) / static_cast<double>(shape[d]);
            if (!(spacing_[d] > 0) || !std::isfinite(spacing_[d]))
                throw std::invalid_argument("Grid spacing is not representable");
        }
        if (!(cell_volume() > 0) || !std::isfinite(cell_volume()))
            throw std::invalid_argument("Cell volume is not representable");
    }

    std::size_t size() const { return size_; }
    const auto& shape() const { return shape_; }
    const auto& spacing() const { return spacing_; }
    double cell_volume() const {
        double v = 1;
        for (double h : spacing_) v *= h;
        return v;
    }
    Box<D> cell(std::size_t index) const {
        if (index >= size_) throw std::out_of_range("Cell index");
        Box<D> b;
        for (std::size_t d = 0; d < D; ++d) {
            const auto k = index % shape_[d];
            index /= shape_[d];
            b.center[d] = lower_[d] + (static_cast<double>(k) + 0.5) * spacing_[d];
            b.half_width[d] = spacing_[d] / 2;
        }
        return b;
    }

private:
    std::array<std::size_t, D> shape_;
    Vec<D> lower_, upper_, spacing_{};
    std::size_t size_{};
};

class Bins {
public:
    explicit Bins(std::vector<double> edges) : edges_(std::move(edges)) {
        if (edges_.size() < 2) throw std::invalid_argument("Need at least two bin edges");
        min_width_ = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < edges_.size(); ++i) {
            if (!std::isfinite(edges_[i])) throw std::invalid_argument("Bin edges must be finite");
            if (i) {
                const double w = edges_[i] - edges_[i - 1];
                if (!(w > 0) || !std::isfinite(w))
                    throw std::invalid_argument("Bin edges must strictly increase with finite widths");
                min_width_ = std::min(min_width_, w);
            }
        }
    }
    static Bins uniform(double lo, double hi, std::size_t count) {
        if (count == 0 || count == std::numeric_limits<std::size_t>::max() ||
            !std::isfinite(lo) || !std::isfinite(hi) || !(hi > lo) || !std::isfinite(hi - lo))
            throw std::invalid_argument("Invalid uniform bins");
        std::vector<double> edges(count + 1);
        for (std::size_t i = 0; i <= count; ++i)
            edges[i] = lo + (hi - lo) * (static_cast<double>(i) / static_cast<double>(count));
        edges.back() = hi;
        return Bins(std::move(edges));
    }
    std::size_t size() const { return edges_.size() - 1; }
    const auto& edges() const { return edges_; }
    double width(std::size_t k) const { return edges_.at(k + 1) - edges_.at(k); }
    double min_width() const { return min_width_; }
    // Every bin is half open; the final upper edge is excluded too.
    std::size_t locate(double s) const {
        if (!std::isfinite(s)) throw std::domain_error("Generating function returned a nonfinite value");
        if (s < edges_.front() || s >= edges_.back()) return size();
        return static_cast<std::size_t>(std::upper_bound(edges_.begin(), edges_.end(), s)
                                       - edges_.begin() - 1);
    }

private:
    std::vector<double> edges_;
    double min_width_{};
};

} // namespace genradon
