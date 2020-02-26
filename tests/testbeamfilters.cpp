


#include "dxmc/beamfilters.h"
#include "dxmc/tube.h"
#include "dxmc/dxmcrandom.h"

bool testUniformWeights(const HeelFilter& filter, double angle_span, double energy)
{
    std::uint64_t s[2];
    randomSeed(s);

    std::size_t N = 1e7;
    double acc = 0;
    for (std::size_t i = 0; i < N; ++i)
    {
        double angle = randomUniform(s, -angle_span * 0.5, angle_span * 0.5);
        acc += filter.sampleIntensityWeight(angle, energy);

    }

    auto mean = acc / N;

    return std::abs(mean - 1.0) < 0.001;
}


bool testWeightsSum(const HeelFilter& filter)
{
    auto as = filter.angleSize();
    auto es = filter.energySize();
    auto w = filter.weights();

    for (std::size_t i = 0; i < es; ++i)
    {
        auto ind = i * as;
        double sum = 0.0;
        for (std::size_t j = 0; j < as; ++j)
        {
            sum += w[ind + j];
        }
        sum = sum / as;
        if (std::abs(sum - 1.0) > 0.001)
            return false;
    }
    return true;
}


int main(int argc, char* argv[])
{
    constexpr double deg2rad = 3.14159265359 / 180.0;
    Tube t;


    t.setAnodeAngle(12 * deg2rad);





    HeelFilter f(t, 10.0 * deg2rad);

    auto w = f.sampleIntensityWeight(0*deg2rad, 90);
    bool weights = testWeightsSum(f);
    bool test = testUniformWeights(f, 10.0 * deg2rad, 90);

    return weights && test;
}