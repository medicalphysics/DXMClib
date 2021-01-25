

#include "dxmc/interpolation.h"
#include "dxmc/dxmcrandom.h"
#include <iostream>
template<dxmc::Floating T>
bool testRITA()
{

    auto pdf = [](T x)->T { return 7 * x * std::exp(-4.0 * x) + 0.6 * std::exp(-12.5 * (x - 3.5) * (x - 3.5)); };

    const T min = 0;
    const T max = 4.5;


    dxmc::RITA<T, 100> rita(min, max, pdf);
    
    
    std::vector<std::size_t> hist(100, 0);
    const T step = (max - min) / hist.size();

    dxmc::RandomState state;
    for (std::size_t i = 0; i < 1e6; ++i) {
        const T x = rita(state, max/2);
        std::size_t ind = (x - min) / step;
        hist[ind]++;
    }

    for (std::size_t i = 0; i < hist.size(); ++i) {
        auto x = min + step * i;
        std::cout << x << ", ";
        std::cout << pdf(x) << ", ";
        std::cout << hist[i] << ", ";
        std::cout << "\n";
    }
    return true;
}



int main()
{


    testRITA<float>();

}
