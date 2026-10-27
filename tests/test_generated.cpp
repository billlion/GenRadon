#include "spheroid3.hpp"
#include <iostream>
#include <random>

int main() {
    using namespace genradon;
    const GeneratedSpheroid3 generated({-2,1,0},{3,0,1});
    const Spheroid<3> reference({-2,1,0},{3,0,1});
    std::mt19937 random(19);
    std::uniform_real_distribution<double> uniform(-1,1);
    for (int k=0;k<100;++k) {
        Vec<3> x{uniform(random),uniform(random),uniform(random)};
        const auto a=generated.jet(x),b=reference.jet(x);
        if (std::abs(a.value-b.value)>1e-12) return 1;
        for (std::size_t i=0;i<3;++i) {
            if (std::abs(a.gradient[i]-b.gradient[i])>1e-12) return 1;
            for (std::size_t j=0;j<3;++j)
                if (std::abs(a.hessian[i][j]-b.hessian[i][j])>1e-12) return 1;
        }
    }
    std::cout << "Generated spheroid agrees with analytic derivatives at 100 points\n";
}
