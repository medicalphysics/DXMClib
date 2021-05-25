

#include "dxmc/beamfilters.h"
#include "dxmc/dxmcrandom.h"
#include "dxmc/tube.h"

using namespace dxmc;
template <typename T>
bool testUniformWeights(const HeelFilter<T>& filter, T angle_span, T energy)
{
    
    RandomState s;

    std::size_t N = 1e7;
    T acc = 0;
    for (std::size_t i = 0; i < N; ++i) {
        const auto angle = s.randomUniform(-angle_span /2, angle_span /2);
        acc += filter.sampleIntensityWeight(angle, energy);
    }

    auto mean = acc / N;

    return std::abs(mean - 1.0) < 0.01;
}
template <typename T>
bool testWeightsSum(const HeelFilter<T>& filter)
{
    auto as = filter.angleSize();
    auto es = filter.energySize();
    auto w = filter.weights();

    for (std::size_t i = 0; i < es; ++i) {
        auto ind = i * as;
        double sum = 0.0;
        for (std::size_t j = 0; j < as; ++j) {
            sum += w[ind + j];
        }
        sum = sum / as;
        if (std::abs(sum - 1.0) >= 0.01)
            return false;
    }
    return true;
}

int main(int argc, char* argv[])
{
    constexpr float deg2rad = 3.14159265359 / 180.0;
    Tube<float> t;
    t.setVoltage(140.0);
    t.setAnodeAngle(12 * deg2rad);
    t.setAlFiltration(9.0);
    t.setSnFiltration(2.0);

    HeelFilter<float> f(t, 10.0 * deg2rad);

    auto w = f.sampleIntensityWeight(0 * deg2rad, 90);
    bool weights = testWeightsSum(f);
    bool test = testUniformWeights(f, 10 * deg2rad, float{90});

    return weights && test;
}