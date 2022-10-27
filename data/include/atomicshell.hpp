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

#include "serialize.hpp"

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

    std::uint8_t shell() const { return m_shell; }
    double numberOfElectrons() const { return m_numberOfElectrons; }
    double bindingEnergy() const { return m_bindingEnergy; }
    double hartreeFockOrbital_0() const { return m_HartreeFockOrbital_0; }
    double numberOfPhotonsPerInitVacancy() const { return m_numberOfPhotonsPerInitVacancy; }
    double energyOfPhotonsPerInitVacancy() const { return m_energyOfPhotonsPerInitVacancy; }

    void setPhotoelectricData(const std::vector<std::pair<double, double>>& data) { m_photoel = data; }
    const auto& photoelectricData() const { return m_photoel; }

    std::vector<char> toBinary() const
    {
        std::vector<char> buffer;
        serialize(&m_shell, buffer);
        serialize(&m_numberOfElectrons, buffer);
        serialize(&m_bindingEnergy, buffer);
        serialize(&m_HartreeFockOrbital_0, buffer);
        serialize(&m_numberOfPhotonsPerInitVacancy, buffer);
        serialize(&m_energyOfPhotonsPerInitVacancy, buffer);
        serialize(m_photoel, buffer);
        const std::uint64_t buffer_size = buffer.size();
        auto buffer_size_ptr = &buffer_size;
        // std::copy(buffer_size_ptr, buffer_size_ptr + sizeof(buffer_size), std::front_inserter(buffer));
        return buffer;
    }
    std::vector<char>::iterator fromBinary(std::vector<char>::iterator begin, std::vector<char>::iterator end);

private:
    std::uint8_t m_shell = 0;
    double m_numberOfElectrons = 0;
    double m_bindingEnergy = 0;
    double m_HartreeFockOrbital_0 = 0;
    double m_numberOfPhotonsPerInitVacancy = 0;
    double m_energyOfPhotonsPerInitVacancy = 0;
    std::vector<std::pair<double, double>> m_photoel;
};