

#include <algorithm>
#include <iostream>
#include <numeric>

#include "dxmc/tube.hpp"

using namespace dxmc;

constexpr double DOUBLEERRF = 1E-6;

template <typename T>
bool testHalfLayerCalculation()
{
    Tube<T> t;
    t.setAlFiltration(9.0);
    auto e = t.getEnergy();
    auto s = t.getSpecter(e);
    Material al(13);
    std::vector<T> att(e.size());
    std::transform(e.cbegin(), e.cend(), att.begin(), [&](auto e) -> T {
        return static_cast<T>(al.standardDensity()) * static_cast<T>(al.getTotalAttenuation(e));
    });

    const auto mmHVL = t.mmAlHalfValueLayer();

    auto I0 = std::reduce(s.cbegin(), s.cend(), T { 0 });
    auto I1 = std::transform_reduce(
        s.cbegin(), s.cend(), att.cbegin(), T { 0 }, std::plus<T>(),
        [=](auto i, auto a) -> T { return i * std::exp(-a * mmHVL * T { .1 }); });
    if ((std::abs(I1 / I0) - 0.5) < 0.01)
        return true;
    return false;
}

template <typename T>
void printSpecter()
{
    Tube<T> t;
    t.setAnodeAngleDeg(30);
    t.setVoltage(100);
    const auto s = t.getSpecter();
    std::cout << "Energy [keV], Intensity [A.U]\n";
    for (const auto [e, n] : s) {
        std::cout << e << ", " << n << '\n';
    }
}

int main(int argc, char* argv[])
{
    printSpecter<float>();
    bool success = testHalfLayerCalculation<float>();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
