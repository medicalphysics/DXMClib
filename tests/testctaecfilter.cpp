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

#include "dxmc/beams/filters/ctaecfilter.hpp"

template <typename T>
bool testAECCTFilter()
{

    dxmc::CTAECFilter<T> f;

    std::vector<T> w = { 5, 6, 7, 6, 3, 5, 6, 2, 6 };

    T start = -5;
    T stop = 4;

    f.setData({ 0, 0, start }, { 0, 0, stop }, w);

    constexpr std::size_t N = 100;

    std::vector<T> res;
    T K = 0;
    for (std::size_t i = 0; i < N; ++i) {
        T k = start + (stop - start) / N;
        res.push_back(f({ 0, 0, k }));
        K += res.back();
    }

    auto test = f.integrate();
    T kk = K / N;
    return kk > 0.95 && kk < 1.05;
}

int main()
{

    bool success = true;
    success = success && testAECCTFilter<float>();
    success = success && testAECCTFilter<double>();
    if (success) {
        return EXIT_SUCCESS;
    } else {
        return EXIT_FAILURE;
    }
}