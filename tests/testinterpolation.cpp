

#include "dxmc/constants.h"
#include "dxmc/dxmcrandom.h"
#include "dxmc/interpolation.h"
#include <iostream>
#include <assert.h>
bool testSpline()
{

    std::vector<float> x(30);
    std::vector<float> y(x.size());
    for (int i = 0; i < x.size(); ++i) {
        x[i] = i * (360.0 / x.size()) * dxmc::DEG_TO_RAD<float>();
        y[i] = std::sin(x[i]);
    }
    dxmc::CubicSplineInterpolator S(x, y);
    float diff = 0;
    for (int i = 0; i < x.size() - 1; ++i) {
        std::cout << x[i] << ", ";
        std::cout << y[i] << ", ";
        std::cout << (x[i] + x[i + 1]) / 2 << ", ";
        std::cout << S((x[i] + x[i + 1]) / 2) << ", ";
        std::cout << std::sin((x[i] + x[i + 1]) / 2) << "\n";
        const auto d = S((x[i] + x[i + 1]) / 2) - std::sin((x[i] + x[i + 1]) / 2);
        diff += d * d;
    }

    const auto rms = std::sqrt(diff) / x.size();
    return rms < 0.001;
}

int main()
{
    bool splinetest = testSpline();
    assert(splinetest);
    return splinetest;
}
