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
#include "dxmc/beams/pencilbeam.hpp"

#include <iostream>

template <dxmc::Floating T, dxmc::BeamType<T> B>
bool initiateBeam(B beam)
{
    auto e = beam.exposure(0);
    dxmc::RandomState state;
    auto p = e.sampleParticle(state);
    return true;
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

int main()
{
    std::cout << "Testing beams\n";

    bool success = true;
    success = success && testpencilbeam<float>();
    success = success && testpencilbeam<double>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}