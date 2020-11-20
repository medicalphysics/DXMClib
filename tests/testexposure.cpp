

#include "dxmc/exposure.h"
#include "dxmc/vectormath.h"

#include <iostream>

template <typename T>
bool testExposure()
{
    std::array<T, 6> cos = { 1, 0, 0, 0, 1, 0 };
    std::array<T, 3> pos = { 0, 0, 0 };
    const T angle = 60;
    std::array<T, 2> ang = { angle * dxmc::DEG_TO_RAD<T>(), 0 };
    dxmc::Exposure<T> exp(pos,
        cos,
        ang);
    std::array<T, 3> dir0;
    dxmc::vectormath::cross(cos.data(), dir0.data());

    std::vector<T> angs(90, 0);
    const auto ang_step = (angle) / angs.size();

    dxmc::RandomState state;

    for (std::uint64_t i = 0; i < 1e7; ++i) {
        auto p = exp.sampleParticle(state);
        const auto angle_s = std::atan2(p.dir[0], p.dir[2]);
        std::size_t idx = angs.size() / 2 + angle_s / (ang_step * dxmc::DEG_TO_RAD<T>());
        if (idx < angs.size())
            angs[idx]++;
    }
    T min = angs[0];
    T max = angs[0];
    for (std::size_t i = 0; i < angs.size(); ++i) {
        std::cout << -angle / 2 + i * ang_step << ", ";
        std::cout << angs[i] << std::endl;
        min = std::min(min, angs[i]);
        max = std::max(max, angs[i]);
    }
    auto test = (max - min) / min;
    if (test < 0.05)
        return true;
    return false;
}

int main()
{
    return testExposure<float>();
}