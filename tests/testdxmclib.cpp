
#include "dxmc/transport.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world.hpp"
#include "dxmc/source.hpp"
#include "dxmc/floating.hpp"

#include <array>
#include <fstream>
#include <iostream>
#include <string>



void testInit()
{
    dxmc::PencilSource ps;
    dxmc::IsotropicSource is;
    dxmc::IsotropicCTSource isc;
    dxmc::DXSource dx;
    dxmc::CBCTSource cb;
    dxmc::CTAxialSource<T> ca;
    dxmc::CTSpiralSource<T> cs;
    dxmc::CTAxialDualSource<T> cad;
    dxmc::CTSpiralDualSource<T> csd;
    dxmc::CTTopogramSource<T> top;
    top.getExposure(0);
}


int main(int argc, char* argv[])
{
    testInit<float>();
    testInit<double>();

    return EXIT_SUCCESS;
}
