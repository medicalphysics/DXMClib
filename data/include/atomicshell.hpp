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

#include <algorithm>
#include <iterator>
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
    void setNumberOfPhotonsPerInitVacancy(double N) { m_numberOfPhotonsPerInitVacancy = N; }
    void setEnergyOfPhotonsPerInitVacancy(double keV) { m_energyOfPhotonsPerInitVacancy = keV; }

    std::uint64_t shell() const { return m_shell; }
    double numberOfElectrons() const { return m_numberOfElectrons; }
    double bindingEnergy() const { return m_bindingEnergy; }
    double hartreeFockOrbital_0() const { return m_HartreeFockOrbital_0; }
    double numberOfPhotonsPerInitVacancy() const { return m_numberOfPhotonsPerInitVacancy; }
    double energyOfPhotonsPerInitVacancy() const { return m_energyOfPhotonsPerInitVacancy; }

    void setPhotoelectricData(const std::vector<std::pair<double, double>>& data) { m_photoel = data; }
    const auto& photoelectricData() const { return m_photoel; }

    std::vector<char> toBinary() const;
    char* fromBinary(std::vector<char>& data, char* begin);

private:
    std::uint64_t m_shell = 0;
    double m_numberOfElectrons = 0;
    double m_bindingEnergy = 0;
    double m_HartreeFockOrbital_0 = 0;
    double m_numberOfPhotonsPerInitVacancy = 0;
    double m_energyOfPhotonsPerInitVacancy = 0;
    std::vector<std::pair<double, double>> m_photoel;
};