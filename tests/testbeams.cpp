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
#include "dxmc/beams/ctspiralbeam.hpp"
#include "dxmc/beams/dxbeam.hpp"
#include "dxmc/beams/isotropicbeam.hpp"
#include "dxmc/beams/isotropicmonoenergybeam.hpp"
#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/vectormath.hpp"

#include <iostream>

template <typename T, dxmc::BeamType<T> B>
bool initiateBeam(B& beam)
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

    Beam beam;
    return initiateBeam<T, Beam>(beam);
}

template <typename T>
bool testCTSpiralBeam()
{
    using Beam = dxmc::CTSpiralBeam<T>;
    Beam beam;
    beam.setStartStop({ 0, 0, 0 }, { 0, 0, 10 });
    beam.setNumberOfParticlesPerExposure(1E6);
    beam.setStepAngleDeg(5);
    beam.setSourceDetectorDistance(115);
    auto f = beam.calibrationFactor();

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
    success = success && testCTSpiralBeam<float>();
    success = success && testCTSpiralBeam<double>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}