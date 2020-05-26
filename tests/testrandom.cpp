
#include "dxmc/dxmcrandom.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

constexpr double ERRF = 1e-4;
inline bool isEqual(double a, double b)
{
    return std::abs(a - b) < ERRF;
}

void testUniform()
{
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> v;
    v.fill(0);

    std::uint64_t seed[2];
    randomSeed(seed);

    for (std::size_t i = 0; i < 1e6; ++i) {
        const auto t = static_cast<std::size_t>(randomUniform<double>(seed, 100.0));
        ++v[t];
    }
    for (std::size_t i = 0; i < N; ++i) {
        std::cout << i << " " << v[i] << std::endl;
    }
}

void testRandomDistributionSingle()
{
    std::vector<double> v({ 1.0 });
    RandomDistribution rnd(v);

    for (std::size_t i = 0; i < 1e6; ++i) {
        auto ind = rnd.sampleIndex();
        assert(ind == 0);
    }
}

void testRandomDistribution()
{
    testRandomDistributionSingle();

    std::vector<double> weights = { 1.423E-04, 2.157E-04, 3.102E-04, 4.324E-04, 5.840E-04, 7.644E-04, 9.784E-04, 1.222E-03, 1.491E-03, 1.803E-03, 2.129E-03, 2.490E-03, 2.863E-03, 3.263E-03, 3.658E-03, 4.093E-03, 4.504E-03, 4.912E-03, 5.347E-03, 5.769E-03, 6.168E-03, 6.582E-03, 6.965E-03, 7.360E-03, 7.710E-03, 8.067E-03, 8.368E-03, 8.671E-03, 8.975E-03, 9.213E-03, 9.476E-03, 9.694E-03, 9.903E-03, 1.009E-02, 1.025E-02, 1.040E-02, 1.053E-02, 1.063E-02, 1.073E-02, 1.081E-02, 1.087E-02, 1.092E-02, 1.096E-02, 1.099E-02, 1.100E-02, 1.100E-02, 1.099E-02, 1.098E-02, 1.095E-02, 1.091E-02, 1.086E-02, 1.081E-02, 1.076E-02, 1.069E-02, 1.063E-02, 1.055E-02, 1.048E-02, 1.039E-02, 1.031E-02, 1.022E-02, 1.012E-02, 1.003E-02, 9.933E-03, 9.828E-03, 9.732E-03, 9.628E-03, 9.516E-03, 9.412E-03, 9.302E-03, 9.193E-03, 9.084E-03, 8.970E-03, 8.862E-03, 8.749E-03, 8.637E-03, 8.526E-03, 8.409E-03, 8.300E-03, 8.185E-03, 8.072E-03, 7.959E-03, 7.847E-03, 7.737E-03, 2.568E-02, 7.513E-03, 7.405E-03, 3.920E-02, 7.181E-03, 7.071E-03, 6.962E-03, 6.854E-03, 6.746E-03, 6.640E-03, 6.530E-03, 6.425E-03, 6.321E-03, 6.214E-03, 6.107E-03, 6.006E-03, 5.901E-03, 5.797E-03, 1.673E-02, 5.592E-03, 5.491E-03, 5.390E-03, 8.223E-03, 5.055E-03, 4.296E-03, 4.236E-03, 4.171E-03, 4.110E-03, 4.048E-03, 3.982E-03, 3.919E-03, 3.852E-03, 3.787E-03, 3.719E-03, 3.654E-03, 3.585E-03, 3.516E-03, 3.449E-03, 3.379E-03, 3.308E-03, 3.240E-03, 3.169E-03, 3.098E-03, 3.026E-03, 2.954E-03, 2.882E-03, 2.809E-03, 2.736E-03, 2.665E-03, 2.592E-03, 2.519E-03, 2.445E-03, 2.370E-03, 2.296E-03, 2.222E-03, 2.148E-03, 2.073E-03, 1.999E-03, 1.925E-03, 1.850E-03, 1.776E-03, 1.700E-03, 1.625E-03, 1.550E-03, 1.476E-03, 1.400E-03, 1.326E-03, 1.251E-03, 1.177E-03, 1.101E-03, 1.027E-03, 9.529E-04, 8.781E-04, 8.041E-04, 7.302E-04, 6.559E-04, 5.823E-04, 5.089E-04, 4.353E-04, 3.623E-04, 2.892E-04, 2.166E-04, 1.441E-04, 7.193E-05, 5.990E-06 };
    std::size_t N = 1e7;
    RandomDistribution rnd(weights);
    std::uint64_t seed[2];
    randomSeed(seed);
    std::vector<std::size_t> hist(weights.size(), 0);
    for (std::size_t i = 0; i < N; ++i) {
        const auto idx = rnd.sampleIndex();
        hist[idx] += 1;
    }

    std::cout << "probability, numberof, hist prob\n";
    for (std::size_t i = 0; i < weights.size(); ++i) {
        const double p = static_cast<double>(hist[i]) / N;
        std::cout << weights[i] << ", " << hist[i] << ", " << p << "\n";
        assert(isEqual(weights[i], p));
    }
}

void testWeights()
{

    std::uint64_t seed[2];
    randomSeed(seed);
    //std::uint64_t u1 = xoroshiro128plus(seed);
    double r1 = randomUniform<double>(seed);
    double r2 = randomUniform<double>(seed, 0.0, 3.14);

    std::cout << seed[0] << " " << seed[1] << " " << r1 << " " << r2 << std::endl;

    std::vector<double> weights({ 0, 1, 2, 3, 4, 5, 4, 3, 2, 1 });

    auto rd = std::make_unique<RandomDistribution>(weights);
    std::vector<std::size_t> res(weights.size());
    for (std::size_t i = 0; i < 1e7; ++i) {
        std::uint64_t ind = rd->sampleIndex();
        res[ind] += 1;
    }

    for (auto k : weights) {
        std::cout << k << " ";
    }
    std::cout << std::endl;

    double max_res = static_cast<double>(*(std::max_element(res.begin(), res.end())));
    double max_weights = *(std::max_element(weights.begin(), weights.end()));

    std::vector<double> res_norm;
    std::transform(res.begin(), res.end(), std::back_inserter(res_norm), [max_res, max_weights](auto d) { return d * max_weights / max_res; });

    for (auto k : res_norm) {
        std::cout << k << " ";
    }
    std::cout << std::endl;

    std::vector<bool> is_equal;
    std::transform(weights.begin(), weights.end(), res_norm.begin(), std::back_inserter(is_equal), [](double w, double r) { return std::abs(w - r) < 0.01; });

    auto success = std::all_of(is_equal.begin(), is_equal.end(), [](auto v) { return v == true; });
    assert(success);
}

int main(int argc, char* argv[])
{
    testUniform();
    testWeights();
    testRandomDistribution();
    return 0;
}
