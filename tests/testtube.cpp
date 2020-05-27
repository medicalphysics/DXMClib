

#include <algorithm>
#include <iostream>
#include <numeric>

#include "dxmc/tube.h"

constexpr double DOUBLEERRF = 1E-6;

bool testHalfLayerCalculation()
{
    Tube t;
    t.setAlFiltration(9.0);
    auto e = t.getEnergy();
    auto s = t.getSpecter(e);
    Material al(13);
    std::vector<double> att(e.size());
    std::transform(e.cbegin(), e.cend(), att.begin(), [&](double e) {
        return al.standardDensity() * al.getTotalAttenuation(e);
    });

    const double mmHVL = t.mmAlHalfValueLayer();

    auto I0 = std::reduce(s.cbegin(), s.cend(), 0.0);
    auto I1 = std::transform_reduce(
        s.cbegin(), s.cend(), att.cbegin(), 0.0, std::plus<double>(),
        [=](double i, double a) { return i * std::exp(-a * mmHVL * 0.1); });
    if ((std::abs(I1 / I0) - 0.5) < 0.01)
        return true;
    return false;
}

int main(int argc, char* argv[])
{
    bool success = testHalfLayerCalculation();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
