

#include "dxmc/attenuationinterpolator.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/interpolation.hpp"

#include <assert.h>
#include <iostream>

/* template <dxmc::Floating T>
bool testSpline()
{

    const T min = 0;
    const T max = 500 / 12.4;

    dxmc::Material mat(16);

    auto fun = [&](const T q) { return mat.getComptonNormalizedScatterFactor(q); };

    dxmc::CubicSplineInterpolator<T, 60> S(min, max, fun);

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
}*/

template <typename T>
void printMatrix(const dxmc::Matrix<T>& mat)
{
    const auto r = mat.nrows();
    const auto c = mat.ncolumns();

    for (std::size_t i = 0; i < r; ++i) {
        for (std::size_t j = 0; j < c; ++j) {
            std::cout << mat.getElement(i, j) << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
template <typename T>
bool testLSSplines()
{
    std::vector<T> x, y, t;
    for (int i = 0; i < 50; ++i) {
        x.push_back(i / T { 20 }-0.5);
        y.push_back(std::exp(-x[i] * x[i]));
    }
    t = { x[0], .1, .5, .9, x.back() };
    dxmc::CubicLSInterpolator<T> s(x, y, t);
    for (int i = 0; i < x.size(); ++i) {
        std::cout << x[i] << ", " << y[i] << ", " << s(x[i]) << std::endl;
    }
    return true;
}

int main()
{
    testLSSplines<double>();
    // bool splinetest = testSpline<float>();
    // assert(splinetest);
    return 1;
}
