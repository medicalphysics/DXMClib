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
#include "serialize.hpp"

bool AtomicShell::operator==(const AtomicShell& other) const
{
    bool eq = m_shell == other.m_shell;
    eq = eq && m_bindingEnergy == other.m_bindingEnergy;
    eq = eq && m_energyOfPhotonsPerInitVacancy == other.m_energyOfPhotonsPerInitVacancy;
    eq = eq && m_HartreeFockOrbital_0 == other.m_HartreeFockOrbital_0;
    eq = eq && m_numberOfElectrons == other.m_numberOfElectrons;
    eq = eq && m_numberOfPhotonsPerInitVacancy == other.m_numberOfPhotonsPerInitVacancy;
    eq = eq && m_photoel == other.m_photoel;

    return eq;
}

std::vector<char> AtomicShell::toBinary() const
{
    std::vector<char> buffer(sizeof(std::uint64_t));
    serialize(m_shell, buffer);
    serialize(m_numberOfElectrons, buffer);
    serialize(m_bindingEnergy, buffer);
    serialize(m_HartreeFockOrbital_0, buffer);
    serialize(m_numberOfPhotonsPerInitVacancy, buffer);
    serialize(m_energyOfPhotonsPerInitVacancy, buffer);
    serialize(m_photoel, buffer);
    std::uint64_t buffer_size = buffer.size() - sizeof(std::uint64_t);
    auto buffer_size_addr = reinterpret_cast<char*>(&buffer_size);
    // writing over first value to contain size of serialized data
    std::copy(buffer_size_addr, buffer_size_addr + sizeof(std::uint64_t), buffer.begin());
    return buffer;
}

char* AtomicShell::fromBinary(std::vector<char>& data, char* begin)
{
    std::uint64_t size;
    auto start = deserialize(size, begin);
    start = deserialize(m_shell, start);
    start = deserialize(m_numberOfElectrons, start);
    start = deserialize(m_bindingEnergy, start);
    start = deserialize(m_HartreeFockOrbital_0, start);
    start = deserialize(m_numberOfPhotonsPerInitVacancy, start);
    start = deserialize(m_energyOfPhotonsPerInitVacancy, start);
    start = deserialize(m_photoel, start);

    return start;
}
