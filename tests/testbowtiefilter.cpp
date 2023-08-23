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

Copyright 2023 Erlend Andersen
*/

#include "dxmc/beams/filters/bowtiefilter.hpp"
#include "dxmc/dxmcrandom.hpp"

#include <iostream>

template <typename T>
bool testBowtie()
{
    std::vector<T> angle, weight;
    angle.push_back(0.166511074);
    weight.push_back(3.53208);
    angle.push_back(0);
    weight.push_back(13.9167);
    angle.push_back(0.041992107);
    weight.push_back(12.5868);
    angle.push_back(0.083836642);
    weight.push_back(9.41943);
    angle.push_back(0.246954945);
    weight.push_back(1.96665);
    angle.push_back(0.324269441);
    weight.push_back(1.27605);
    angle.push_back(0.390607044);
    weight.push_back(0.947716);

    dxmc::BowtieFilter<T> filter(angle, weight);

    // test sampling
    T w = 0;
    constexpr std::size_t N = 1E6;
    dxmc::RandomState state;
    for (std::size_t i = 0; i < N; ++i) {
        const T ang = state.randomUniform(T { 0.390607044 });
        w += filter(ang);
    }
    w /= N;

    auto success = std::abs(w - 1) < T { 0.01 };

    return success;
}

int main()
{

    bool success = true;
    success = success && testBowtie<float>();
    success = success && testBowtie<double>();
    if (success) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}