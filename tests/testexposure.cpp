

#include "dxmc/exposure.h"
#include "dxmc/vectormath.h"

#include <iostream>

template <typename T>
bool testExposure()
{
    std::array<T, 6> cos = { 1, 0, 0, 0, 1, 0 };
    std::array<T, 3> pos = { 0, 0, 0 };
    std::array<T, 3> dir0;
    dxmc::vectormath::cross(cos.data(), dir0.data());

    const T angleX = std::atan(T { 160 } / 600) * 2;
    const T angleY = std::atan(T { 40 } / 600) * 2;
    std::array<T, 2> ang = { angleX, angleY };

    dxmc::vectormath::rotate(&cos[0], &dir0[0], 165 * dxmc::DEG_TO_RAD<T>());
    dxmc::vectormath::rotate(&cos[3], &dir0[0], 165 * dxmc::DEG_TO_RAD<T>());

    std::array<T, 3> direction;
    dxmc::vectormath::cross(&cos[0], &direction[0]);

    dxmc::Exposure<T> exp(pos,
        cos,
        ang);

    std::vector<std::size_t> angsX(90, 0);
    std::vector<std::size_t> angsY(angsX.size(), 0);

    const auto ang_stepX = angleX / angsX.size();
    const auto ang_stepY = angleY / angsX.size();

    dxmc::RandomState state;

    for (std::uint64_t i = 0; i < 1e7; ++i) {
        auto p = exp.sampleParticle(state);
        //const auto angle1 = std::atan2(p.dir[0], p.dir[2]);
        const auto angle1 = dxmc::vectormath::angleBetweenOnPlane(p.dir.data(), direction.data(), &cos[3]);
        const std::size_t idx1 = (angleX / 2 + angle1) / ang_stepX;
        ++angsX[idx1];
        //const auto angle2 = std::atan2(p.dir[1], p.dir[2]);
        const auto angle2 = dxmc::vectormath::angleBetweenOnPlane(p.dir.data(), direction.data(), &cos[0]);
        const std::size_t idx2 = (angleY / 2 + angle2) / ang_stepY;
        ++angsY[idx2];
    }

    for (std::size_t i = 0; i < angsX.size(); ++i) {
        std::cout << (-angleX / 2 + ang_stepX * (i + .5)) * dxmc::RAD_TO_DEG<T>() << ", " << angsX[i] << ", ";
        std::cout << (-angleY / 2 + ang_stepY * (i + .5)) * dxmc::RAD_TO_DEG<T>() << ", " << angsY[i] << "\n";
    }
    const auto rmsX = std::sqrt(std::transform_reduce(angsX.cbegin(), angsX.cend(), T { 0 }, std::plus<>(), [](const auto s) -> T { return s * s; }) / angsX.size());
    const auto rmsY = std::sqrt(std::transform_reduce(angsY.cbegin(), angsY.cend(), T { 0 }, std::plus<>(), [](const auto s) -> T { return s * s; }) / angsY.size());

    const T lim = T { 1e7 } / angsX.size();
    if (std::abs((rmsX - lim) / lim) > 0.01)
        return false;
    if (std::abs((rmsY - lim) / lim) > 0.01)
        return false;

    return true;
}

int main()
{
    return testExposure<float>();
}