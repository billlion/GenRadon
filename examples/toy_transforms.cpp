#include <genradon/genradon.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

using namespace genradon;
namespace fs = std::filesystem;

template <std::size_t D> struct Configuration { Vec<D> a, b, axis; };

template <std::size_t D> std::vector<Configuration<D>> acquisition(std::size_t count) {
    std::vector<Configuration<D>> result;
    const double pi = std::acos(-1.0);
    for (std::size_t k = 0; k < count; ++k) {
        Vec<D> n{}, tangent{};
        if constexpr (D == 2) {
            const double angle = 2*pi*static_cast<double>(k)/static_cast<double>(count);
            n = {std::cos(angle),std::sin(angle)};
            tangent = {-n[1],n[0]};
        } else {
            // Fibonacci sphere: equal weights, approximately uniform directions.
            const double z = 1-2*(static_cast<double>(k)+0.5)/static_cast<double>(count);
            const double angle = pi*(3-std::sqrt(5.0))*static_cast<double>(k);
            const double radial = std::sqrt(1-z*z);
            n = {radial*std::cos(angle),radial*std::sin(angle),z};
            tangent = {-std::sin(angle),std::cos(angle),0};
        }
        Configuration<D> c;
        for (std::size_t d = 0; d < D; ++d) {
            c.a[d] = 3*n[d];
            // Receiver 60 degrees away on the same exterior circle/sphere.
            c.b[d] = 3*(0.5*n[d]+std::sqrt(3.0)/2*tangent[d]);
            c.axis[d] = -n[d];
        }
        result.push_back(c);
    }
    return result;
}

std::ofstream output(const fs::path& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open output: " + path.string());
    f.exceptions(std::ios::badbit | std::ios::failbit);
    f << std::setprecision(17);
    return f;
}

template <std::size_t D, class G>
void run_case(const std::string& name, const Grid<D>& grid, const Bins& bins,
              const std::vector<G>& views, const std::vector<Configuration<D>>& configurations,
              const fs::path& root, Options options) {
    Projector<D,G> projector(grid,bins,views,options);
    std::vector<double> image(grid.size());
    for (std::size_t i = 0; i < grid.size(); ++i) {
        const auto x = grid.cell(i).center;
        double a = 0,b = 0;
        for (std::size_t d = 0; d < D; ++d) {
            const double ca = d==0 ? -0.3 : 0.1;
            const double cb = d==0 ? 0.35 : -0.25;
            a += std::pow((x[d]-ca)/0.25,2);
            b += std::pow((x[d]-cb)/0.18,2);
        }
        image[i] = std::exp(-a/2)+0.7*std::exp(-b/2);
    }
    const auto start = std::chrono::steady_clock::now();
    Statistics stats;
    const auto data = projector.project(image,&stats);
    const auto bp = projector.backproject(data);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    const double mass = grid.cell_volume()*std::accumulate(image.begin(),image.end(),0.0);
    double mass_error = 0;
    for (std::size_t v = 0; v < views.size(); ++v) {
        const double sum = std::accumulate(data.begin()+v*bins.size(),data.begin()+(v+1)*bins.size(),0.0);
        mass_error = std::max(mass_error,std::abs(sum-mass));
    }
    if (mass_error > 1e-10*std::max(1.0,mass))
        throw std::runtime_error(name+": acquisition lost mass (absolute error "+std::to_string(mass_error)+")");
    std::cout << name << ": " << grid.size() << " cells, " << views.size() << " views, "
              << bins.size() << " bins; max mass error " << mass_error << "; " << seconds << " s\n";
    if (root.empty()) return;
    const auto path = root/name;
    fs::create_directories(path);
    auto object_file = output(path/"object.csv");
    object_file << "cell";
    for (std::size_t d=0;d<D;++d) object_file << ",x" << d;
    object_file << ",density,transpose_of_data\n";
    for (std::size_t i=0;i<grid.size();++i) {
        object_file << i;
        for (double x:grid.cell(i).center) object_file << ',' << x;
        object_file << ',' << image[i] << ',' << bp[i] << '\n';
    }
    auto data_file = output(path/"data.csv");
    data_file << "view,bin,lower,upper,bin_integral\n";
    for (std::size_t v=0;v<views.size();++v)
        for (std::size_t k=0;k<bins.size();++k)
            data_file << v << ',' << k << ',' << bins.edges()[k] << ',' << bins.edges()[k+1]
                      << ',' << data[v*bins.size()+k] << '\n';
    auto geometry_file = output(path/"geometry.csv");
    geometry_file << "view";
    for (const auto* label:{"a","b","axis"})
        for (std::size_t d=0;d<D;++d) geometry_file << ',' << label << d;
    geometry_file << '\n';
    for (std::size_t i=0;i<configurations.size();++i) {
        geometry_file << i;
        for (const auto& x:{configurations[i].a,configurations[i].b,configurations[i].axis})
            for (double v:x) geometry_file << ',' << v;
        geometry_file << '\n';
    }
    auto metadata = output(path/"metadata.txt");
    metadata << "dimension=" << D << "\nshape=";
    for (auto n:grid.shape()) metadata << n << ' ';
    metadata << "\ndomain=[-1,1]^dimension\nindex_order=x0_fastest\nnormalization=bin_integral"
             << "\nweight=natural_coarea\nintegration="
             << (options.integration==Integration::midpoint?"midpoint":"adaptive_affine")
             << "\npoints_per_axis=" << options.points_per_axis
             << "\naffine_tolerance=" << options.affine_tolerance
             << "\nmax_depth=" << options.max_depth
             << "\npoint_evaluations=" << stats.point_evaluations
             << "\naffine_regions=" << stats.affine_regions
             << "\nsubdivisions=" << stats.subdivisions
             << "\nbackprojection=Euclidean_transpose_not_reconstruction\n";
}

template <std::size_t D> void run_dimension(const fs::path& path, Options options, bool smoke) {
    std::array<std::size_t,D> shape;
    shape.fill(smoke?4:(D==2?24:12));
    Vec<D> lower{},upper{};
    lower.fill(-1); upper.fill(1);
    Grid<D> grid(shape,lower,upper);
    const std::size_t views = smoke?6:(D==2?64:48), bins = smoke?12:48;
    const auto configurations = acquisition<D>(views);
    std::vector<Cone<D>> cones;
    std::vector<Spheroid<D>> spheroids;
    std::vector<Hyperboloid<D>> hyperboloids;
    for (const auto& c:configurations) {
        cones.emplace_back(c.a,c.axis);
        spheroids.emplace_back(c.a,c.b);
        hyperboloids.emplace_back(c.a,c.b);
    }
    const std::string suffix=std::to_string(D)+"d";
    // A guard margin beyond cosine=1 also contains the affine surrogate near
    // the cone axis. For these settings its maximum allowed error is < 0.01.
    run_case("cone"+suffix,grid,Bins::uniform(0.7,1.01,bins),cones,configurations,path,options);
    run_case("spheroid"+suffix,grid,Bins::uniform(2.0,10.0,bins),spheroids,configurations,path,options);
    run_case("hyperboloid"+suffix,grid,Bins::uniform(-3.1,3.1,bins),hyperboloids,configurations,path,options);
}

int main(int argc,char** argv) {
    try {
        fs::path path;
        Options options;
        options.points_per_axis=3;
        options.affine_tolerance=0.15;
        options.max_depth=24;
        bool smoke=false;
        for (int i=1;i<argc;++i) {
            const std::string arg=argv[i];
            if (arg=="--affine") options.integration=Integration::adaptive_affine;
            else if (arg=="--smoke") smoke=true;
            else if (arg=="--output" && i+1<argc) path=argv[++i];
            else throw std::invalid_argument("Usage: genradon_toy [--affine] [--smoke] [--output DIRECTORY]");
        }
        run_dimension<2>(path,options,smoke);
        run_dimension<3>(path,options,smoke);
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
