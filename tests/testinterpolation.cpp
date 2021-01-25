

#include "dxmc/interpolation.h"

template<dxmc::Floating T>
bool testRITA()
{

    auto pdf = [](T x)->T { return 7 * x * std::exp(-4.0 * x) + 0.6 * std::exp(-12.5 * (x - 3.5) * (x - 3.5)); };

    T min = 0;
    T max = 4.5;


    dxmc::RITA rita(min, max, pdf);



    return true;
}



int main()
{


    testRITA<float>();

}
