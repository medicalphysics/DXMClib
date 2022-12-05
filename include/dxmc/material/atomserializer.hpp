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

#include "dxmc/material/atomicelement.hpp"
#include "dxmc/material/atomicshell.hpp"

#include <algorithm>
#include <concepts>
#include <iterator>
#include <map>
#include <vector>

namespace dxmc {
template <typename T>
concept Number = std::is_integral<T>::value || std::is_floating_point<T>::value;

class AtomSerializer {
public:
    static std::vector<char> serializeAtoms(const std::map<std::uint64_t, AtomicElement<double>>& elements)
    {
        std::vector<char> buffer;
        std::uint64_t n_elements = elements.size();
        serialize(n_elements, buffer);
        for (const auto& [Z, atom] : elements) {
            serializeAtomicElement(atom, buffer);
        }
        return buffer;
    }
    static std::map<std::uint64_t, AtomicElement<double>> deserializeAtoms(std::vector<char>& buffer)
    {
        std::map<std::uint64_t, AtomicElement<double>> elements;
        if (buffer.size() == 0)
            return elements;
        std::uint64_t number_elements;
        auto start = deserialize(number_elements, &(buffer[0]));

        for (std::uint64_t i = 0; i < number_elements; i++) {
            AtomicElement<double> atom;
            start = deserializeAtomicElement(atom, start);
            elements[atom.Z] = atom;
        }
        return elements;
    }

protected:
    template <Number T>
    static void serialize(T in, std::vector<char>& buffer)
    {
        auto dest = std::back_inserter(buffer);
        auto in_c = reinterpret_cast<char*>(&in);
        std::copy(in_c, in_c + sizeof(T), dest);
    }

    static void appendToBuffer(const std::vector<char>& data, std::vector<char>& buffer)
    {
        std::copy(data.cbegin(), data.cend(), std::back_inserter(buffer));
    }

    template <typename T>
        requires(!std::is_same<T, char>::value)
    static void serialize(const std::vector<T>& data, std::vector<char>& buffer)
    {
        std::uint64_t size = data.size() * sizeof(T);
        serialize(size, buffer);
        auto in_c = reinterpret_cast<const char*>(data.data());
        auto dest = std::back_inserter(buffer);
        std::copy(in_c, in_c + size, dest);
    }

    static void serializeAtomicShell(const AtomicShell<double>& shell, std::vector<char>& buffer)
    {
        std::vector<char> data;
        serialize(shell.shell, data);
        serialize(shell.numberOfElectrons, data);
        serialize(shell.bindingEnergy, data);
        serialize(shell.kineticEnergy, data);
        serialize(shell.HartreeFockOrbital_0, data);
        serialize(shell.numberOfPhotonsPerInitVacancy, data);
        serialize(shell.energyOfPhotonsPerInitVacancy, data);
        serialize(shell.photoel, data);

        std::uint64_t data_size = data.size();
        serialize(data_size, buffer); // adding data size to buffer
        appendToBuffer(data, buffer); // adding data
    }

    static void serializeAtomicElement(const AtomicElement<double>& atom, std::vector<char>& buffer)
    {
        std::vector<char> data;
        serialize(atom.Z, data);
        serialize(atom.atomicWeight, data);
        serialize(atom.standardDensity, data);
        serialize(atom.coherent, data);
        serialize(atom.incoherent, data);
        serialize(atom.photoel, data);
        serialize(atom.formFactor, data);
        serialize(atom.incoherentSF, data);

        // adding shells
        // adding number of shells
        std::uint64_t n_shells = atom.shells.size();
        serialize(n_shells, data);
        // adding each shell
        for (const auto& [id, shell] : atom.shells) {
            serializeAtomicShell(shell, data);
        }
        std::uint64_t data_size = data.size();
        serialize(data_size, buffer); // adding data size to buffer
        appendToBuffer(data, buffer); // adding data to buffer
    }

    template <Number T>
    static char* deserialize(T& value, char* begin)
    {
        auto val_ptr = reinterpret_cast<T*>(begin);
        value = *val_ptr;
        return begin + sizeof(T);
    }
    template <typename T>
    static char* deserialize(std::vector<T>& val, char* begin, std::size_t size)
    {
        auto n_elements = size / (sizeof(T));
        val.clear();
        val.reserve(n_elements);
        auto start = reinterpret_cast<T*>(begin);
        std::copy(start, start + n_elements, std::back_inserter(val));
        return begin + n_elements * sizeof(T);
    }

    template <typename T>
    static char* deserialize(std::vector<T>& val, char* begin)
    {
        std::uint64_t size;
        auto start = deserialize(size, begin);
        return deserialize(val, start, size);
    }
    static char* deserializeAtomicShell(AtomicShell<double>& shell, char* begin)
    {
        std::uint64_t size { 0 };
        auto start = deserialize(size, begin);
        start = deserialize(shell.shell, start);
        start = deserialize(shell.numberOfElectrons, start);
        start = deserialize(shell.bindingEnergy, start);
        start = deserialize(shell.kineticEnergy, start);
        start = deserialize(shell.HartreeFockOrbital_0, start);
        start = deserialize(shell.numberOfPhotonsPerInitVacancy, start);
        start = deserialize(shell.energyOfPhotonsPerInitVacancy, start);
        start = deserialize(shell.photoel, start);
        return start;
    }
    static char* deserializeAtomicElement(AtomicElement<double>& atom, char* begin)
    {
        std::uint64_t size;
        auto start = deserialize(size, begin);
        start = deserialize(atom.Z, start);
        start = deserialize(atom.atomicWeight, start);
        start = deserialize(atom.standardDensity, start);
        start = deserialize(atom.coherent, start);
        start = deserialize(atom.incoherent, start);
        start = deserialize(atom.photoel, start);
        start = deserialize(atom.formFactor, start);
        start = deserialize(atom.incoherentSF, start);

        std::uint64_t n_shells { 0 };
        atom.shells.clear();
        start = deserialize(n_shells, start);
        for (std::uint64_t i = 0; i < n_shells; ++i) {
            AtomicShell<double> shell;
            start = deserializeAtomicShell(shell, start);
            atom.shells[shell.shell] = shell;
        }
        return start;
    }
};

}