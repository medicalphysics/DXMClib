

#include "dxmc/constants.h"
#include "dxmc/dxmcrandom.h"
#include "dxmc/interpolation.h"
#include <iostream>

bool testSpline()
{

    std::vector<float> x(30);
    std::vector<float> y(x.size());
    for (int i = 0; i < x.size(); ++i) {
        x[i] = i * (360.0 / x.size()) * dxmc::DEG_TO_RAD<float>();
        y[i] = std::sin(x[i]);
    }
    dxmc::CubicSplineInterpolator S(x, y);
    for (int i = 0; i < x.size()-1; ++i) {
        std::cout << x[i] << ", ";
        std::cout << y[i] << ", ";
        std::cout << (x[i] + x[i+1])/2 << ", ";
        std::cout << S((x[i] + x[i + 1]) / 2) << "\n";
    }
    
return true;
}

int main()
{

    testSpline();
}
