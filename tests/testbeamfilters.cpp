




#include "dxmc/beamfilters.h"
#include "dxmc/tube.h"

int main(int argc, char* argv[])
{
    constexpr double deg2rad = 3.14159265359 / 180.0;
    Tube t;


    t.setAnodeAngle(12 * deg2rad);


    HeelFilter f(t, 10.0 * deg2rad);

    auto w = f.sampleIntensityWeight(0*deg2rad, 50);


    return 0;
}