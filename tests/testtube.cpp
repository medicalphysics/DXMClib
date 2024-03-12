

#include <algorithm>
#include <iostream>
#include <numeric>

#include "dxmc/beams/tube/tube.hpp"

using namespace dxmc;

constexpr double DOUBLEERRF = 1E-6;

bool testHalfLayerCalculation()
{

    Tube t;
    t.setAlFiltration(2.0);
    auto e = t.getEnergy();
    auto s = t.getSpecter(e);
    const auto& al = AtomHandler::Atom(13);

    auto p = interpolate(al.photoel, e);
    auto i = interpolate(al.incoherent, e);
    auto c = interpolate(al.coherent, e);
    auto att = addVectors(p, i, c);

    const auto al_dens = al.standardDensity;

    std::transform(att.cbegin(), att.cend(), att.begin(), [al_dens](auto a) {
        return a * al_dens;
    });

    const auto mmHVL = t.mmAlHalfValueLayer();

    auto I0 = std::reduce(s.cbegin(), s.cend(), 0.0);
    auto I1 = std::transform_reduce(
        s.cbegin(), s.cend(), att.cbegin(), 0.0, std::plus<>(),
        [=](auto i, auto a) { return i * std::exp(-a * mmHVL * .1); });

    bool success = (std::abs(I1 / I0) - 0.5) < 0.01;
    if (success)
        std::cout << "SUCCESS ";
    else
        std::cout << "FAILURE ";

    return success;
}

void printSpecter()
{
    Tube t;
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
    std::cout << "Testing tube\n";
    bool success = testHalfLayerCalculation();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
