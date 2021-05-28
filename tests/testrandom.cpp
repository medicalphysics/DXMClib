
#include "dxmc/dxmcrandom.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

using namespace dxmc;

constexpr double ERRF = 1e-4;
inline bool isEqual(double a, double b)
{
    return std::abs(a - b) < ERRF;
}

template <typename T>
void testUniform()
{
    constexpr std::size_t S = 1e8;
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> pcg;
    pcg.fill(0);

    RandomState state;

    T maxVal = 0;
    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto r = state.randomUniform<T>();
        maxVal = std::max(r, maxVal);
        const auto t = static_cast<std::size_t>(r * N);
        if (t == N)
            pcg[t - 1]++;
        else
            pcg[t]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Testing random uniform\nNominal count, PCG count, Difference\n";
    std::cout << "Max value is 1, obtained max value is " << maxVal << "\n";
    for (auto t : pcg) {
        const auto diff = t * N * 100.0 / S - 100.0;
        assert(std::abs(diff) < 5.0);
        std::cout << S / N << ", " << t << ", " << diff << "\n";
    }

    std::cout << "Counted random numbers: " << std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 }) << "\n";
    assert(std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 }) == S);

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "PCG time: " << pcg_time.count() << std::endl;
}

void testUniformRange()
{
    constexpr std::size_t S = 1e8;
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> pcg;
    pcg.fill(0);

    RandomState state;

    constexpr double r1 = 5.0;
    constexpr double r2 = 10.0;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto t = state.randomUniform(r1, r2);
        const auto idx = static_cast<std::size_t>(N * ((t - r1) / (r2 - r1)));
        pcg[idx]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Testing random uniform double range\nNominal count, PCG count, Difference\n";
    for (auto t : pcg) {
        const auto diff = t * N * 100.0 / S - 100.0;
        //assert(std::abs(diff) < 5.0);
        std::cout << S / N << ", " << t << ", " << diff << "\n";
    }

    std::cout << "Counted random numbers: " << std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 }) << "\n";
    assert(std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 }) == S);

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "PCG time: " << pcg_time.count() << std::endl;
}

void testrandomInteger(std::uint8_t max = 128)
{
    std::vector<std::size_t> vals(max, 0);
    RandomState state;

    std::size_t test = 4;
    state.randomInteger32BitCapped(test);

    for (std::size_t i = 0; i < max * 1e6; ++i) {
        const auto ind = state.randomInteger(max);
        vals[ind]++;
    }
    for (std::size_t i = 0; i < vals.size(); ++i) {
        std::cout << i << ", " << vals[i] << "\n";
    }
    auto mean = std::accumulate(vals.cbegin(), vals.cend(), 0.0) / max;
    auto var = std::transform_reduce(vals.cbegin(), vals.cend(), 0.0, std::plus<>(), [=](auto v) {
        const auto vd = static_cast<double>(v);
        const auto diff = vd - mean;
        return diff * diff;
    }) / (max - 1);
    std::cout << mean << ", " << var << ", " << std::sqrt(var) << "\n";
    assert(mean > var / 2);
}
void testUniformIndex()
{
    constexpr std::size_t S = 1e8;
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> pcg;
    pcg.fill(0);

    RandomState state;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto t = state.randomUniform(N);
        assert(t < N);
        pcg[t]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Testing random index\nNominal count, PCG count, Difference\n";
    for (auto t : pcg) {
        const auto diff = t * N * 100.0 / S - 100.0;
        assert(std::abs(diff) < 5.0);
        std::cout << S / N << ", " << t << ", " << diff << "\n";
    }

    std::cout << "Counted random numbers: " << std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 }) << "\n";
    assert(std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 }) == S);

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

    std::cout << "PCG time: " << pcg_time.count() << std::endl;
}

template <typename T>
void testRita()
{

    auto pdf = [](const T e) -> T { return e * e; };

    T min = 0;
    T max = 50;
    constexpr std::size_t steps = 150;
    std::vector<T> vals;
    for (std::size_t i = 0; i < steps; ++i) {
        vals.push_back(min + (max - min) * i / (steps - 1));
    }

    RITA<T, steps * 2> r(min, max, pdf);
    RandomState state;
    std::vector<T> hist(steps, 0);
    std::size_t nsamp = 1E7;
    for (std::size_t i = 0; i < nsamp; ++i) {
        const T v = r(state);
        const auto ind = (v - min) * (steps - 1) / (max - min);
        ++hist[static_cast<std::size_t>(ind)];
    }
    T rms = 0;
    std::cout << "RITA test, RITA, ANA" << std::endl;
    for (int i = 0; i < steps; ++i) {
        const T ana = i * i / (T { steps } * steps * steps) * 3;
        std::cout << i << ", " << hist[i] / nsamp << ", " << ana << std::endl;
        if (i > 0 && i < steps - 2) {
            const T diff = hist[i] / nsamp - ana;
            rms += diff * diff;
        }
    }
    rms /= steps - 2;
    rms = std::sqrt(rms);
    std::cout << "RMS: " << rms << std::endl;
    assert(rms < 0.001);
}

int main(int argc, char* argv[])
{
    //testrandomInteger();
    testUniform<double>();
    testUniform<float>();
    //testUniformRange();
    //testUniformIndex();
    //testRita<float>();
    return 0;
}
