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

#include "atomicshell.hpp"




void AtomicShell::setPhotoelectricData(const std::vector<double>& data, double minEnergy, double maxEnergy, double barnToAtt)
{
    const auto N = data.size() / 2;
    m_photoel.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * 0.001;
        if (e >= minEnergy && e <= maxEnergy) {
            const double a = data[ind + 1] * barnToAtt;
            m_photoel.push_back(std::make_pair(e, a));
        }
    }
    m_photoel.shrink_to_fit();
}