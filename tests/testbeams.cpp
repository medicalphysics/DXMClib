/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2022 Erlend Andersen
*/

#include "dxmc/beams/beamtype.hpp"
#include "dxmc/beams/dxbeam.hpp"
#include "dxmc/beams/isotropicbeam.hpp"
#include "dxmc/beams/isotropicmonoenergybeam.hpp"
#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/beams/ctspiralbeam.hpp"

#include <iostream>

template <typename T, dxmc::BeamType<T> B>
bool initiateBeam(B beam)
{
    auto e = beam.exposure(0);
    dxmc::RandomState state;
    auto p = e.sampleParticle(state);
    return true;
}

template <typename T>
bool testDXBeam()
{
    dxmc::DXBeam<T> beam;

    auto& tube = beam.tube();
    return initiateBeam<T>(beam);
}

template <typename T>
bool testpencilbeam()
{

    dxmc::PencilBeam<T> beam;

    auto e = beam.exposure(0);

    dxmc::RandomState state;

    auto p = e.sampleParticle(state);

    return initiateBeam<T, dxmc::PencilBeam<T>>(beam);
}

template <typename T>
bool testIsotropicMonoEnergyBeam()
{
    using Beam = dxmc::IsotropicMonoEnergyBeam<T>;

    Beam beam({ -60, 0, 0 }, { 0, 1, 0, 0, 0, 1 });

    const T collangle_y = std::atan(T { 16 } / T { 60 });
    const auto collangle_z = std::atan(T { 4 } / T { 60 });
    beam.setCollimationAngles({ -collangle_y, -collangle_z, collangle_y, collangle_z });

    const std::array<T, 3> co_x = { 0, 1, 0 };
    const std::array<T, 3> co_y = { 0, 0, 1 };
    const std::array<T, 3> pos = { -60, 0, 0 };

    dxmc::RandomState state;

    std::vector<std::vector<std::size_t>> ang_histy;

    for (std::size_t angInt = 0; angInt < 360; angInt = angInt + 45) {
        const T angle = static_cast<T>(angInt) * dxmc::DEG_TO_RAD<T>();
        auto x = dxmc::vectormath::rotate(co_x, { 0, 0, 1 }, angle);
        auto p_ang = dxmc::vectormath::rotate(pos, { 0, 0, 1 }, angle);
        beam.setPosition(p_ang);
        beam.setDirectionCosines(x, co_y);
        auto exposure = beam.exposure(0);
        const auto dir = dxmc::vectormath::cross(x, co_y);
        dxmc::RandomState state;

        std::vector<std::size_t> h(100, 0);

        for (std::size_t i = 0; i < 1E6; ++i) {
            auto p = exposure.sampleParticle(state);
            auto pdir = p.dir;
            dxmc::vectormath::normalize(pdir);
            auto ang = std::acos(dxmc::vectormath::dot(dir, pdir));
            int ind = static_cast<int>((100 * ang) / collangle_z);
            if (ind < h.size())
                h[ind]++;
        }
        ang_histy.push_back(h);
    }
    for (int j = 0; j < ang_histy[0].size(); ++j) {
        std::cout << (collangle_z * j) / 100 << ", ";
        for (int i = 0; i < ang_histy.size(); ++i) {
            std::cout << ang_histy[i][j] << ", ";
        }
        std::cout << std::endl;
    }

    return true;
}

template <typename T>
bool testCTSpiralSource()
{

    using Beam = dxmc::CTSpiralBeam<T>;
    Beam beam;
    beam.setStartStop({ 0, 0, 0 }, { 0, 0, 1 });
    auto e = beam.exposure(0);

    return true;
}

int main()
{
    std::cout << "Testing beams\n";

    bool success = true;
    success = success && testDXBeam<float>();
    success = success && testDXBeam<double>();
    success = success && testIsotropicMonoEnergyBeam<double>();
    success = success && testIsotropicMonoEnergyBeam<float>();
    success = success && testpencilbeam<float>();
    success = success && testpencilbeam<double>();
    success = success && testCTSpiralSource<float>();
    success = success && testCTSpiralSource<double>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}