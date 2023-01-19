
#include "dxmc/dxmcrandom.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numbers>
#include <numeric>
#include <vector>

using namespace dxmc;

constexpr double ERRF = 1e-4;
inline bool isEqual(double a, double b)
{
    return std::abs(a - b) < ERRF;
}

template <typename T>
bool testUniform()
{
    constexpr std::size_t S = 1e8;
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> pcg;
    pcg.fill(0);

    RandomState state;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto r = state.randomUniform<T>();
        const auto t = static_cast<std::size_t>(std::nextafter(r * N, T { 0 }));
        pcg[t]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();

    // check total number
    const auto pcg_sum = std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 });
    bool success = pcg_sum == S;

    T max_deviation = 0;
    for (auto t : pcg) {
        const auto diff = t * N * 100.0 / S - 100.0;
        // deviation i percent
        const auto dev = (static_cast<T>(t) * N) / static_cast<T>(S) - T { 1 };
        max_deviation = std::max(std::abs(dev), max_deviation);
    }
    // max deviation is less than 0.5 %
    success = success && max_deviation * T { 100 } < T { 0.5 };

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::string msg = success ? "SUCCSESS" : "FAILURE";
    std::cout << msg << " Test random uniform, time: " << pcg_time.count() << "ms" << std::endl;
    return success;
}

template <typename T>
bool testUniformRange()
{
    constexpr std::size_t S = 1e8;
    constexpr std::size_t N = 100;
    std::array<std::size_t, N> pcg;
    pcg.fill(0);

    RandomState state;

    constexpr T r1 = 5;
    constexpr T r2 = 10;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto t = state.randomUniform(r1, r2);
        const auto idx = static_cast<std::size_t>(std::nextafter(N * ((t - r1) / (r2 - r1)), T { 0 }));
        pcg[idx]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();

    // check total number
    const auto pcg_sum = std::accumulate(pcg.cbegin(), pcg.cend(), std::size_t { 0 });
    bool success = pcg_sum == S;

    T max_deviation = 0;
    for (auto t : pcg) {
        const auto diff = t * N * 100.0 / S - 100.0;
        // deviation i percent
        const auto dev = (static_cast<T>(t) * N) / static_cast<T>(S) - T { 1 };
        max_deviation = std::max(std::abs(dev), max_deviation);
    }
    // max deviation is less than 0.5 %
    success = success && max_deviation * T { 100 } < T { 0.5 };

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::string msg = success ? "SUCCSESS" : "FAILURE";
    std::cout << msg << " Test random uniform range, time: " << pcg_time.count() << "ms" << std::endl;
    return success;
}

bool testrandomInteger(std::uint8_t max = 128)
{
    std::vector<std::size_t> vals(max, 0);
    RandomState state;
    const std::size_t S = 1E6 * max;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < S; ++i) {
        const auto ind = state.randomInteger(max);
        vals[ind]++;
    }
    const auto end = std::chrono::high_resolution_clock::now();
    auto sum = std::accumulate(vals.cbegin(), vals.cend(), std::size_t { 0 });

    bool success = sum == S;

    double max_diff = 0;
    for (auto n : vals) {
        const auto per_n = static_cast<double>(S) / max;
        const auto diff = static_cast<double>(n) / per_n - 1;
        max_diff = std::max(std::abs(diff), max_diff);
    }
    success = success && max_diff * 100 < 1.0;

    auto pcg_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::string msg = success ? "SUCCSESS" : "FAILURE";
    std::cout << msg << " Test random uniform index, time: " << pcg_time.count() << "ms" << std::endl;
    return success;
}

template <typename T>
bool testCPDFdist(bool print = false)
{
    auto f = [](T x) { return std::exp(-(x * x) / 2); };
    const T xmin = -3;
    const T xmax = 3;
    const std::size_t N = 50;
    dxmc::CPDFSampling samp(xmin, xmax, f);
    dxmc::RandomState state;
    std::vector<T> hist(N, T { 0 });
    const auto step = (xmax - xmin) / (N);
    constexpr std::size_t Nsamp = 1E7;

    const auto start = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < Nsamp; ++i) {
        const T x = samp(state);
        const auto idx = static_cast<std::size_t>(std::nextafter((x - xmin) / step, T { 0 }));
        hist[idx] += 1;
    }
    const auto end = std::chrono::high_resolution_clock::now();
    std::for_each(hist.begin(), hist.end(), [](auto& h) { h /= Nsamp; });

    auto func = hist;
    T func_sum = 0;
    for (std::size_t i = 0; i < hist.size(); ++i) {
        const auto x = xmin + step * i + step / 2;
        const auto p = f(x);
        func_sum += p;
        func[i] = p;
    }
    std::for_each(func.begin(), func.end(), [&](auto& p) { p /= func_sum; });

    std::vector<T> error(hist.size(), 0);

    if (print)
        std::cout << "x, sampled PDF, analytical pdf" << std::endl;

    T maxError = 0;
    for (std::size_t i = 0; i < hist.size(); ++i) {
        const auto x = xmin + step * i + step / 2;
        error[i] = std::abs(hist[i] - func[i]);
        if (print)
            std::cout << x << ", " << hist[i] << ", " << func[i] << std::endl;
        maxError = std::max(error[i], maxError);
    }
    auto success = maxError < T { 0.01 };

    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::string msg = success ? "SUCCSESS" : "FAILURE";
    std::cout << msg << " Test Gaussian analytical pdf, time: " << time.count() << "ms" << std::endl;
    return success;
}

int main(int argc, char* argv[])
{
    std::cout << "Testing random number generator:" << std::endl;

    bool success = true;
    success = success && testrandomInteger();
    success = success && testUniform<double>();
    success = success && testUniform<float>();
    success = success && testUniformRange<double>();
    success = success && testUniformRange<float>();
    success = success && testCPDFdist<double>();
    success = success && testCPDFdist<float>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
