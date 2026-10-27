#include <genradon/genradon.hpp>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <string>

using namespace genradon;

namespace {
std::size_t checks = 0;
void check(bool ok, const std::string& message) {
    ++checks;
    if (!ok) throw std::runtime_error(message);
}
void near(double actual, double expected, double tolerance, const std::string& message) {
    check(std::isfinite(actual) && std::abs(actual - expected) <= tolerance,
          message + ": actual=" + std::to_string(actual) + ", expected=" + std::to_string(expected));
}
template <class F> void throws(F f, const std::string& message) {
    bool caught = false;
    try { f(); } catch (const std::exception&) { caught = true; }
    check(caught, message);
}
double inner(const std::vector<double>& a, const std::vector<double>& b) {
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

// Independent inclusion-exclusion formula for the CDF of sum of uniforms.
// Used only as a test oracle, with well-separated coefficients.
template <std::size_t D> long double cdf(long double s, const Vec<D>& widths) {
    long double answer = 0, denominator = 1;
    for (std::size_t d = 0; d < D; ++d) denominator *= (d + 1) * widths[d];
    for (unsigned mask = 0; mask < (1U << D); ++mask) {
        long double t = s;
        int sign = 1;
        for (std::size_t d = 0; d < D; ++d) if (mask & (1U << d)) {
            t -= widths[d];
            sign = -sign;
        }
        if (t > 0) answer += sign * std::pow(t, static_cast<int>(D));
    }
    return std::clamp(answer / denominator, 0.0L, 1.0L);
}

template <std::size_t D> void test_clipping() {
    Box<D> box;
    box.half_width.fill(0.5);
    box.center.fill(0.5);
    std::mt19937 rng(41 + D);
    std::uniform_real_distribution<double> random(0.2, 2.0);
    for (int trial = 0; trial < 80; ++trial) {
        Vec<D> gradient{}, widths{};
        double shift = 0, center = 0, span = 0;
        for (std::size_t d = 0; d < D; ++d) {
            widths[d] = random(rng);
            gradient[d] = trial & (1 << d) ? -widths[d] : widths[d];
            center += gradient[d] / 2;
            if (gradient[d] < 0) shift += gradient[d];
            span += widths[d];
        }
        double total = 0;
        for (int k = 0; k < 17; ++k) {
            const double a = shift + span * k / 17;
            const double b = shift + span * (k + 1) / 17;
            const double v = affine_band_volume(box, center, gradient, a, b);
            const double expected = static_cast<double>(cdf(b - shift, widths) - cdf(a - shift, widths));
            near(v, expected, 2e-12, "Affine clipping vs independent distribution formula");
            check(v >= 0, "Positive affine volume");
            total += v;
        }
        near(total, 1, 3e-13, "Affine volume partition");
    }
    Vec<D> axis{};
    axis[0] = 1;
    near(affine_band_volume(box, 0.5, axis, 0.5, 0.5 + 1e-10),
         (0.5 + 1e-10) - 0.5, 1e-16, "Narrow slab");
    Vec<D> diagonal;
    diagonal.fill(1);
    const double width = 1e-5;
    double expected = std::pow(width, D);
    for (std::size_t d = 1; d <= D; ++d) expected /= d;
    near(affine_band_volume(box, D / 2.0, diagonal, 0, width), expected,
         expected * 1e-8, "Grazing simplex sliver");
    Vec<D> zero{};
    near(affine_band_volume(box, 0, zero, 0, 1), 1, 0, "Constant level at lower edge");
    near(affine_band_volume(box, 1, zero, 0, 1), 0, 0, "Constant level at final edge");
}

template <std::size_t D, class G> void test_derivatives(const G& g) {
    Vec<D> x{};
    for (std::size_t d = 0; d < D; ++d) x[d] = 0.12 + 0.17 * d;
    const auto j = g.jet(x);
    near(j.value, g.value(x), 2e-14, "Jet value");
    const double h = 1e-5;
    for (std::size_t d = 0; d < D; ++d) {
        auto p = x, m = x;
        p[d] += h;
        m[d] -= h;
        near(j.gradient[d], (g.value(p) - g.value(m)) / (2*h), 2e-8, "Gradient finite difference");
        const auto jp = g.jet(p), jm = g.jet(m);
        for (std::size_t e = 0; e < D; ++e) {
            near(j.hessian[e][d], (jp.gradient[e] - jm.gradient[e]) / (2*h), 2e-8, "Hessian finite difference");
            near(j.hessian[e][d], j.hessian[d][e], 1e-14, "Hessian symmetry");
        }
    }
    Box<D> box{x, {}};
    box.half_width.fill(0.1);
    const double bound = g.hessian_bound(box);
    std::mt19937 rng(83);
    std::uniform_real_distribution<double> random(-1, 1);
    for (int k = 0; k < 40; ++k) {
        Vec<D> y{}, v{};
        for (std::size_t d = 0; d < D; ++d) {
            y[d] = x[d] + random(rng) * box.half_width[d];
            v[d] = random(rng);
        }
        v = unit(v);
        const auto jy = g.jet(y);
        Vec<D> hv{};
        for (std::size_t d = 0; d < D; ++d) hv[d] = dot(jy.hessian[d], v);
        check(norm(hv) <= bound + 1e-12, "Hessian bound on sample directions");
    }
}

template <std::size_t D, class G> void test_operator(const std::vector<G>& views) {
    std::array<std::size_t,D> shape;
    shape.fill(3);
    Vec<D> lower{}, upper{};
    lower.fill(-0.7);
    upper.fill(0.9);
    Grid<D> grid(shape, lower, upper);
    // Nonuniform, full coverage for these tests.
    Bins bins({-12,-4,-1,-0.2,0.15,0.5,0.8,1.01,2,4,6,8,12});
    std::mt19937 rng(108);
    std::uniform_real_distribution<double> random(-1, 1);
    std::vector<double> x(grid.size()), y(views.size()*bins.size()), w(y.size());
    for (auto& v : x) v = random(rng);
    for (std::size_t i = 0; i < y.size(); ++i) {
        y[i] = random(rng);
        w[i] = 0.7 * bins.width(i % bins.size());
    }
    for (auto mode : {Integration::midpoint, Integration::adaptive_affine}) {
        for (auto normalization : {Normalization::bin_integral, Normalization::bin_average}) {
            Options options;
            options.integration = mode;
            options.normalization = normalization;
            options.points_per_axis = 3;
            options.affine_tolerance = 0.15;
            Projector<D,G> a(grid, bins, views, options);
            const auto ax = a.project(x), aty = a.backproject(y);
            near(inner(ax,y), inner(x,aty), 3e-12, "Euclidean adjoint identity");
            auto wax = ax;
            for (std::size_t i = 0; i < ax.size(); ++i) wax[i] *= w[i];
            near(inner(wax,y), grid.cell_volume()*inner(x,a.weighted_adjoint(y,w)), 3e-12,
                 "Weighted adjoint identity");
            const double mass = grid.cell_volume()*std::accumulate(x.begin(),x.end(),0.0);
            for (std::size_t v = 0; v < views.size(); ++v) {
                double sum = 0;
                for (std::size_t k = 0; k < bins.size(); ++k)
                    sum += ax[v*bins.size()+k] * (normalization == Normalization::bin_average ? bins.width(k) : 1);
                near(sum,mass,3e-12,"Mass conservation across bins");
            }
            const auto zero = a.project(std::vector<double>(grid.size(),0));
            check(std::all_of(zero.begin(),zero.end(),[](double v){return v == 0;}),"Zero projects to zero");
        }
    }
}

template <std::size_t D> struct SquaredRadius {
    double value(const Vec<D>& x) const { return dot(x,x); }
    Jet<D> jet(const Vec<D>& x) const {
        Jet<D> j;
        j.value = value(x);
        for (std::size_t d = 0; d < D; ++d) { j.gradient[d] = 2*x[d]; j.hessian[d][d] = 2; }
        return j;
    }
    double hessian_bound(const Box<D>&) const { return 2; }
};

template <std::size_t D> void test_shell_convergence() {
    std::array<std::size_t,D> shape;
    shape.fill(4);
    Vec<D> lo{},hi{};
    lo.fill(-1); hi.fill(1);
    Grid<D> grid(shape,lo,hi);
    Bins bins({0.1,0.25,0.45,0.65});
    std::vector<double> errors;
    for (double tolerance : {0.15, 0.015}) {
        Options opt;
        opt.integration = Integration::adaptive_affine;
        opt.affine_tolerance = tolerance;
        opt.max_depth = 21;
        Projector<D,SquaredRadius<D>> p(grid,bins,{SquaredRadius<D>{}},opt);
        auto data = p.project(std::vector<double>(grid.size(),1));
        double error = 0;
        const double ball = D == 2 ? std::acos(-1.0) : 4*std::acos(-1.0)/3;
        for (std::size_t k=0;k<bins.size();++k) {
            const double exact = ball*(std::pow(bins.edges()[k+1],D/2.0)-std::pow(bins.edges()[k],D/2.0));
            error += std::abs(data[k]-exact);
        }
        errors.push_back(error);
    }
    check(errors[1] < errors[0] * 0.55, "Tighter affine tolerance improves analytic shell volumes");
    check(errors[1] < 0.02, "Shell volume absolute error");
}

void test_interfaces_and_errors() {
    throws([]{ Bins b({0,0,1}); },"Reject repeated edges");
    throws([]{ Bins b({0,std::numeric_limits<double>::infinity()}); },"Reject infinite edges");
    throws([]{ Grid<2> g({0,2},{0,0},{1,1}); },"Reject empty grid");
    throws([]{ Plane<2> p({0,0}); },"Reject zero normal");
    throws([]{ Cone<3> c({0,0,0},{1,0,0}); (void)c.value({0,0,0}); },"Reject cone apex");
    Bins b({0,1,2});
    check(b.locate(0)==0 && b.locate(1)==1 && b.locate(2)==b.size(),"Half-open bins");
    Grid<2> grid({2,3},{-1,-3},{1,3});
    near(grid.cell(1).center[0],0.5,0,"x fastest index");
    near(grid.cell(2).center[1],0,0,"next coordinate stride");
    struct ValueOnly { double value(const Vec<2>& x) const { return x[0]; } };
    Projector<2,ValueOnly> p(grid,Bins({-1,0,1}),{ValueOnly{}});
    const auto data = p.project(std::vector<double>(6,1));
    near(data[0],6,1e-14,"Value-only geometry");
    near(data[1],6,1e-14,"Value-only geometry second bin");
    throws([&]{ p.project({1}); },"Reject wrong input length");
    Options opt;
    opt.integration = Integration::adaptive_affine;
    throws([&]{ Projector<2,ValueOnly> a(grid,b,{ValueOnly{}},opt); },"Reject missing derivatives");
    struct Estimated {
        double value(const Vec<2>& x) const { return x[0]; }
        Jet<2> jet(const Vec<2>& x) const { return {x[0],{1,0},{}}; }
    };
    throws([&]{ Projector<2,Estimated> a(grid,b,{Estimated{}},opt); },"Reject implicit heuristic bounds");
    opt.allow_estimated_bounds = true;
    Projector<2,Estimated> estimated(grid,Bins({-1,0,1}),{Estimated{}},opt);
    Statistics stats;
    estimated.project(std::vector<double>(6,1),&stats);
    check(stats.estimated_bound_regions>0,"Report heuristic bound use");
    opt.max_depth=0;
    Projector<2,SquaredRadius<2>> limited(grid,b,{SquaredRadius<2>{}},opt);
    throws([&]{ limited.project(std::vector<double>(6,1)); },"Depth limit never silently ignores tolerance");
}
}

int main() {
    try {
        test_clipping<2>(); test_clipping<3>();
        test_derivatives<2>(Plane<2>({1,2}));
        test_derivatives<2>(Spheroid<2>({-2,0},{2,1}));
        test_derivatives<2>(Hyperboloid<2>({-2,0},{2,1}));
        test_derivatives<2>(Cone<2>({-2,0},{1,1}));
        test_derivatives<3>(Plane<3>({1,2,3}));
        test_derivatives<3>(Spheroid<3>({-2,0,0},{2,1,0}));
        test_derivatives<3>(Hyperboloid<3>({-2,0,0},{2,1,0}));
        test_derivatives<3>(Cone<3>({-2,0,0},{1,1,1}));
        test_operator<2>(std::vector<Plane<2>>{Plane<2>({1,0}),Plane<2>({1,1})});
        test_operator<2>(std::vector<Spheroid<2>>{Spheroid<2>({-2,0},{2,1})});
        test_operator<2>(std::vector<Hyperboloid<2>>{Hyperboloid<2>({-2,0},{2,1})});
        test_operator<2>(std::vector<Cone<2>>{Cone<2>({-2,0},{1,0})});
        test_operator<3>(std::vector<Plane<3>>{Plane<3>({1,0,0}),Plane<3>({1,1,1})});
        test_operator<3>(std::vector<Spheroid<3>>{Spheroid<3>({-2,0,0},{2,1,0})});
        test_operator<3>(std::vector<Hyperboloid<3>>{Hyperboloid<3>({-2,0,0},{2,1,0})});
        test_operator<3>(std::vector<Cone<3>>{Cone<3>({-2,0,0},{1,0,0})});
        test_shell_convergence<2>(); test_shell_convergence<3>();
        test_interfaces_and_errors();
        std::cout << "Passed " << checks << " checks\n";
    } catch (const std::exception& e) {
        std::cerr << "FAILED: " << e.what() << '\n';
        return 1;
    }
}
