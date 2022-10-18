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

#include <utility>
#include <vector>

class AtomicShell {
    enum class Shell {};

public:
    AtomicShell() {};
    AtomicShell(std::uint8_t shell)
        : m_shell(shell)
    {
    }
    void setShell(std::uint8_t shell) { m_shell = shell; }
    void setBindingEnergy(double en) { m_bindingEnergy = en; }
    void setNumberOfElectrons(double N) { m_numberOfElectrons = N; }

    void setPhotoelectricData(const std::vector<std::pair<double, double>>& data) { m_photoel = data; }
    const auto& photoelectricData() const { return m_photoel; }

private:
    std::uint8_t m_shell = 0;
    double m_numberOfElectrons = 0;
    double m_bindingEnergy = 0;
    double m_HartreeFockOrbital_0 = 0;
    std::vector<std::pair<double, double>> m_photoel;
};