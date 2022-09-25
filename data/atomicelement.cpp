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

#include "atomicelement.hpp"

#include "dxmc/interpolation.hpp"

AtomicElement::AtomicElement(std::uint8_t Z)
    : m_Z(Z)
{
}

void AtomicElement::setCoherentData(const std::vector<double>& data)
{

    const auto N = data.size() / 2;
    m_coherent.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * 0.001;
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * barnToAtt();
            m_coherent.push_back(std::make_pair(e, a));
        }
    }
    m_coherent.shrink_to_fit();
}
void AtomicElement::setIncoherentData(const std::vector<double>& data)
{

    const auto N = data.size() / 2;
    m_incoherent.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * 0.001;
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * barnToAtt();
            m_incoherent.push_back(std::make_pair(e, a));
        }
    }
    m_incoherent.shrink_to_fit();
}

void AtomicElement::setShellBindingEnergy(std::uint8_t shell, double bindingEnergy)
{
}

void AtomicElement::setShellPhotoelectricData(std::uint8_t shell, const std::vector<double>& data)
{
}

void AtomicElement::setPhotoelectricData(const std::vector<double>& data)
{

    const auto N = data.size() / 2;
    m_photoel.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * 0.001;
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * barnToAtt();
            m_photoel.push_back(std::make_pair(e, a));
        }
    }
    m_photoel.shrink_to_fit();
}