

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
        x.push_back(i / T { 50 });
        y.push_back(1 - std::exp(-x[i]));
        if (i == 20) {
            x.push_back(i / T { 50 });
            y.push_back(1 - std::exp(-x[i]));
            y.back() += .3;
        }
        if (i > 20)
            y.back() += 0.3;
    }
    t = { x[0], .1, .5, .9, x.back(),   };
    std::sort(t.begin(), t.end());
    dxmc::CubicLSInterpolator<T> s(x, y, t);
    for (int i = 0; i < x.size(); ++i) {
        std::cout << x[i] << ", " << y[i] << ", " << s(x[i]) << std::endl;
    }
    return true;
}

bool testMatrix()
{
    std::vector<double> a { 1, 1, 1, 0, 2, 5, 2, 5, -1 };
    dxmc::Matrix<double> m(3, 3, a);
    std::vector<double> b { 6, -4, 27 };

    std::cout << "A:" << std::endl;
    auto res = m.solve(b);

    std::vector<double> test { 5, 3, -2 };
    auto y = m * res;
    auto y2 = m * test;

    return true;
}

int main()
{

    testLSSplines<double>();
    // bool splinetest = testSpline<float>();
    // assert(splinetest);
    return 1;
}
