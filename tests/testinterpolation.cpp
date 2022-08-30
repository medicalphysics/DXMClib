

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
template <dxmc::Floating T>
bool testMatrix()
{
    // std::vector<T> test { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
    // dxmc::Matrix m(4, 4, test);

    //std::vector<T> test { 1, 2, -1, -4, 2, 3, -1, -11, -2, 0, -3, 22 };
    std::vector<T> test { 5, -6, -7, 7, 3, -2, 5, -17, 2, 4, -3, 29 };
    dxmc::Matrix m(3, 4, test);
    printMatrix(m);
    printMatrix(m.transpose());
    // m.swapRows(0, 1);
    // printMatrix(m);

    printMatrix(m.reducedRowEchelonForm());

    return true;
}

int main()
{
    testMatrix<double>();
    // bool splinetest = testSpline<float>();
    // assert(splinetest);
    return 1;
}
