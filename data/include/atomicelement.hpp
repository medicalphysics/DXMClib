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

#pragma once

#include "atomicshell.hpp"

#include <vector>

class AtomicElement {
public:
    AtomicElement() {};
    AtomicElement(std::uint8_t Z);

    void setZ(std::uint8_t Z) { m_Z = Z; }
    void setAtomicWeight(double AW) { m_atomicWeight = AW; }

    void setPhotoelectricData(const std::vector<double>& data);
    void setCoherentData(const std::vector<double>& data);
    void setIncoherentData(const std::vector<double>& data);

    constexpr double maxPhotonEnergy() { return 500.0; }
    constexpr double minPhotonEnergy() { return 1.0; }
    double barnToAtt() {
        constexpr double u = 1.6605402;
        return 1.0 / (u * m_atomicWeight);            
    }

private:
    std::uint8_t m_Z = 0;
    double m_atomicWeight = 0.0;
    std::vector<std::pair<double, double>> m_coherent;
    std::vector<std::pair<double, double>> m_incoherent;
    std::vector<std::pair<double, double>> m_photoel;


    std::vector<AtomicShell> m_shells;



};
