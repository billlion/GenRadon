#pragma once

#include "affine_volume.hpp"
#include <type_traits>
#include <utility>

namespace genradon {

enum class Integration { midpoint, adaptive_affine };
enum class Normalization { bin_integral, bin_average };

struct Options {
    Integration integration = Integration::midpoint;
    Normalization normalization = Normalization::bin_integral;
    std::size_t points_per_axis = 2;
    double affine_tolerance = 0.05; // allowed level error / smallest bin width
    unsigned max_depth = 15;       // binary subdivisions of a voxel
    bool allow_estimated_bounds = false;
};

struct Statistics {
    std::size_t point_evaluations = 0;
    std::size_t jet_evaluations = 0;
    std::size_t affine_regions = 0;
    std::size_t subdivisions = 0;
    std::size_t estimated_bound_regions = 0;
};

namespace detail {
template <class G, std::size_t D, class = void> struct has_jet : std::false_type {};
template <class G, std::size_t D> struct has_jet<G, D,
    std::void_t<decltype(std::declval<const G&>().jet(std::declval<const Vec<D>&>()))>> : std::true_type {};
template <class G, std::size_t D, class = void> struct has_bound : std::false_type {};
template <class G, std::size_t D> struct has_bound<G, D,
    std::void_t<decltype(std::declval<const G&>().hessian_bound(std::declval<const Box<D>&>()))>> : std::true_type {};
}

// Geometry is a C++ value type: value(x), optionally jet(x), hessian_bound(box).
// Each element of views is one acquisition configuration. No interpreter callback
// occurs in the integration loop, and no system matrix is stored.
template <std::size_t D, class Geometry> class Projector {
public:
    Projector(Grid<D> grid, Bins bins, std::vector<Geometry> views, Options options = {})
        : grid_(std::move(grid)), bins_(std::move(bins)), views_(std::move(views)), options_(options) {
        if (views_.empty()) throw std::invalid_argument("Need at least one acquisition configuration");
        if (views_.size() > std::numeric_limits<std::size_t>::max() / bins_.size())
            throw std::overflow_error("Data size overflow");
        if (options_.points_per_axis == 0 || options_.points_per_axis > 1024)
            throw std::invalid_argument("points_per_axis must be in [1,1024]");
        if (!std::isfinite(options_.affine_tolerance) || !(options_.affine_tolerance > 0) ||
            options_.max_depth > 30)
            throw std::invalid_argument("Invalid affine tolerance or depth (maximum 30)");
        if (options_.integration == Integration::adaptive_affine) {
            if constexpr (!detail::has_jet<Geometry, D>::value)
                throw std::invalid_argument("Adaptive affine integration requires jet(x)");
            if constexpr (!detail::has_bound<Geometry, D>::value)
                if (!options_.allow_estimated_bounds)
                    throw std::invalid_argument("Supply hessian_bound(box) or explicitly allow estimated bounds");
        }
    }

    const Grid<D>& grid() const { return grid_; }
    const Bins& bins() const { return bins_; }
    std::size_t data_size() const { return views_.size() * bins_.size(); }

    std::vector<double> project(const std::vector<double>& object, Statistics* statistics = nullptr) const {
        check_vector(object, grid_.size());
        std::vector<double> result(data_size(), 0);
        coefficients([&](std::size_t row, std::size_t col, double a) {
            result[row] += a * object[col];
        }, statistics);
        return result;
    }

    // Euclidean matrix transpose, not an implicitly weighted continuous adjoint.
    std::vector<double> backproject(const std::vector<double>& data, Statistics* statistics = nullptr) const {
        check_vector(data, data_size());
        std::vector<double> result(grid_.size(), 0);
        coefficients([&](std::size_t row, std::size_t col, double a) {
            result[col] += a * data[row];
        }, statistics);
        return result;
    }

    // M_X^{-1} A^T M_Y, M_X = voxel volume * I. Callers explicitly choose
    // the data measure, e.g. view_weight * bin_width for bin-average data.
    std::vector<double> weighted_adjoint(const std::vector<double>& data,
                                         const std::vector<double>& data_weights) const {
        check_vector(data, data_size());
        check_vector(data_weights, data_size());
        auto weighted = data;
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (!(data_weights[i] > 0)) throw std::invalid_argument("Data inner-product weights must be positive");
            weighted[i] *= data_weights[i];
        }
        auto result = backproject(weighted);
        for (double& v : result) v /= grid_.cell_volume();
        return result;
    }

private:
    static void check_vector(const std::vector<double>& x, std::size_t size) {
        if (x.size() != size) throw std::invalid_argument("Array size mismatch");
        for (double v : x) if (!std::isfinite(v)) throw std::invalid_argument("Arrays must contain finite values");
    }

    template <class Emit> void coefficients(Emit emit, Statistics* statistics) const {
        Statistics stats;
        for (std::size_t view = 0; view < views_.size(); ++view) {
            for (std::size_t col = 0; col < grid_.size(); ++col) {
                auto add = [&](std::size_t bin, double volume) {
                    if (options_.normalization == Normalization::bin_average) volume /= bins_.width(bin);
                    emit(view * bins_.size() + bin, col, volume);
                };
                const auto box = grid_.cell(col);
                if (options_.integration == Integration::midpoint) midpoint(views_[view], box, add, stats);
                else {
                    if constexpr (detail::has_jet<Geometry, D>::value)
                        affine(views_[view], box, 0, add, stats);
                }
            }
        }
        if (statistics) *statistics = stats;
    }

    template <class Add> void midpoint(const Geometry& geometry, const Box<D>& box,
                                       Add& add, Statistics& stats) const {
        std::size_t count = 1;
        for (std::size_t d = 0; d < D; ++d) count *= options_.points_per_axis;
        const double weight = box.volume() / static_cast<double>(count);
        for (std::size_t q = 0; q < count; ++q) {
            auto index = q;
            Vec<D> x{};
            for (std::size_t d = 0; d < D; ++d) {
                const auto k = index % options_.points_per_axis;
                index /= options_.points_per_axis;
                x[d] = box.center[d] + box.half_width[d] *
                    (2 * (static_cast<double>(k) + 0.5) / static_cast<double>(options_.points_per_axis) - 1);
            }
            const auto bin = bins_.locate(geometry.value(x));
            ++stats.point_evaluations;
            if (bin != bins_.size()) add(bin, weight);
        }
    }

    template <class Add> void affine(const Geometry& geometry, const Box<D>& box,
                                     unsigned depth, Add& add, Statistics& stats) const {
        double bound = 0;
        if constexpr (detail::has_bound<Geometry, D>::value) {
            bound = geometry.hessian_bound(box);
            if (std::isnan(bound) || bound < 0) throw std::domain_error("Invalid Hessian bound");
            // Avoid differentiating at a singular center while it can be refined.
            if (!std::isfinite(bound)) {
                split(geometry, box, depth, add, stats);
                return;
            }
        }
        const auto j = geometry.jet(box.center);
        ++stats.jet_evaluations;
        if (!std::isfinite(j.value)) throw std::domain_error("Nonfinite generating function");
        for (std::size_t d = 0; d < D; ++d) {
            if (!std::isfinite(j.gradient[d])) throw std::domain_error("Nonfinite gradient");
            for (std::size_t e = 0; e < D; ++e) {
                if (!std::isfinite(j.hessian[d][e])) throw std::domain_error("Nonfinite Hessian");
                if constexpr (!detail::has_bound<Geometry, D>::value)
                    bound = std::hypot(bound, j.hessian[d][e]);
            }
        }
        if constexpr (!detail::has_bound<Geometry, D>::value) ++stats.estimated_bound_regions;
        const double error = 0.5 * bound * dot(box.half_width, box.half_width);
        double span = 0;
        for (std::size_t d = 0; d < D; ++d) span += std::abs(j.gradient[d]) * box.half_width[d];
        if (!std::isfinite(error) || !std::isfinite(span))
            throw std::domain_error("Affine error or span is not representable");
        // Bound the true range before excluding a region; do not use corner signs.
        if (j.value + span + error < bins_.edges().front() ||
            j.value - span - error >= bins_.edges().back()) return;
        if (error > options_.affine_tolerance * bins_.min_width()) {
            split(geometry, box, depth, add, stats);
            return;
        }
        ++stats.affine_regions;
        if (span == 0) {
            const auto k = bins_.locate(j.value);
            if (k != bins_.size()) add(k, box.volume());
            return;
        }
        const auto& edges = bins_.edges();
        const auto start = std::upper_bound(edges.begin(), edges.end(), j.value - span);
        std::size_t k = start == edges.begin() ? 0 : static_cast<std::size_t>(start - edges.begin() - 1);
        for (; k < bins_.size() && edges[k] < j.value + span; ++k) {
            const double volume = affine_band_volume(box, j.value, j.gradient, edges[k], edges[k+1]);
            if (volume > 0) add(k, volume);
        }
    }

    template <class Add> void split(const Geometry& geometry, const Box<D>& box,
                                    unsigned depth, Add& add, Statistics& stats) const {
        if (depth >= options_.max_depth)
            throw std::runtime_error("Affine tolerance not reached: increase max_depth, relax tolerance, or exclude singularities");
        const auto axis = static_cast<std::size_t>(std::max_element(box.half_width.begin(), box.half_width.end())
                                                   - box.half_width.begin());
        auto left = box, right = box;
        left.half_width[axis] /= 2;
        right.half_width[axis] /= 2;
        left.center[axis] -= left.half_width[axis];
        right.center[axis] += right.half_width[axis];
        if (left.center[axis] == box.center[axis] || right.center[axis] == box.center[axis])
            throw std::runtime_error("Subdivision is below coordinate precision");
        ++stats.subdivisions;
        affine(geometry, left, depth + 1, add, stats);
        affine(geometry, right, depth + 1, add, stats);
    }

    Grid<D> grid_;
    Bins bins_;
    std::vector<Geometry> views_;
    Options options_;
};

} // namespace genradon
