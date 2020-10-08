
#include "dxmc/dxmcrandom.h"

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

void testUniform()
{
    constexpr std::size_t S = 1e8;
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> pcg;
    pcg.fill(0);

    RandomState state;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto t = static_cast<std::size_t>(state.randomUniform<double>() * N);
        pcg[t]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();

    std::cout << "Testing random uniform\nNominal count, PCG count, Difference\n";
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

int main(int argc, char* argv[])
{
    testUniform();
    testUniformRange();
    testUniformIndex();

    return 0;
}
