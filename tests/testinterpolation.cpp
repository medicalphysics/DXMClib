

#include "dxmc/attenuationinterpolator.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material.hpp"

#include "xraylib.h"

#include <assert.h>
#include <iostream>

template <dxmc::Floating T>
bool testSpline()
{

    const T min = 0;
    const T max = 500 / 12.4;

    dxmc::Material mat(97);

    auto fun = [&](const T q) { return mat.getComptonNormalizedScatterFactor(q); };

    dxmc::CubicSplineInterpolator<T, 15> S(min, max, fun);

    std::vector<T> q(200);
    std::vector<T> sf(q.size());
    for (std::size_t i = 0; i < q.size(); ++i) {
        q[i] = min + (max - min) / (q.size() - 1) * i;
        sf[i] = fun(q[i]);
    }

    std::cout << "q, xraylib, spline\n";
    for (std::size_t i = 0; i < q.size(); ++i) {
        std::cout << q[i] << ", ";
        std::cout << sf[i] << ", ";
        std::cout << S(q[i]) << std::endl;
    }
    return true;
}

int main()
{
    bool splinetest = testSpline<float>();
    assert(splinetest);
    return 1;
}
