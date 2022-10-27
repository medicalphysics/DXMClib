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
#include "serialize.hpp"

#include <map>
#include <utility>
#include <vector>

class AtomicElement {
public:
    AtomicElement() {};
    AtomicElement(std::uint8_t Z);

    void setZ(std::uint8_t Z) { m_Z = Z; }
    void setAtomicWeight(double AW) { m_atomicWeight = AW; }

    std::uint8_t Z() const { return m_Z; }
    double atomicWeight() const { return m_atomicWeight; }

    void setPhotoelectricData(const std::vector<double>& data);
    void setCoherentData(const std::vector<double>& data);
    void setIncoherentData(const std::vector<double>& data);

    void setFormFactor(const std::vector<double>& data);
    void setImaginaryAnomalousSF(const std::vector<double>& data);
    void setRealAnomalousSF(const std::vector<double>& data);
    void setIncoherentSF(const std::vector<double>& data);

    const auto& photoelectricData() const { return m_photoel; }
    const auto& coherentData() const { return m_coherent; }
    const auto& incoherentData() const { return m_incoherent; }

    const auto& formFactor() const { return m_formFactor; }
    const auto& imaginaryAnomalousSF() const { return m_imagAnomSF; }
    const auto& realAnomalousSF() const { return m_realAnomSF; }
    const auto& incoherentSF() const { return m_incoherentSF; }

    void setShellBindingEnergy(const std::vector<double>& data);
    void setShellPhotoelectricData(const std::uint8_t shell, const std::vector<double>& data);
    void setShellNumberOfElectrons(const std::vector<double>& data);
    void setShellNumberOfPhotonsPerInitVacancy(const std::vector<double>& data);
    void setShellEnergyOfPhotonsPerInitVacancy(const std::vector<double>& data);

    const std::map<std::uint8_t, AtomicShell>& shells() const { return m_shells; }

    std::vector<char> toBinary() const
    {
        std::vector<char> buffer;

        serialize(&m_Z, buffer);
        serialize(&m_atomicWeight, buffer);
        serialize(m_coherent, buffer);
        serialize(m_incoherent, buffer);
        serialize(m_photoel, buffer);
        serialize(m_formFactor, buffer);
        serialize(m_incoherentSF, buffer);

        // adding shells
        // adding number of shells
        std::uint64_t n_shells = m_shells.size();
        serialize(&n_shells, buffer);
        // adding each shell
        for (const auto& [id, shell] : m_shells) {
            auto s_buffer = shell.toBinary();
            serialize(s_buffer, buffer);
        }

        const std::uint64_t buffer_size = buffer.size();
        auto buffer_size_ptr = &buffer_size;

        // std::copy(buffer_size_ptr, buffer_size_ptr + sizeof(buffer_size), std::front_inserter(buffer));
        return buffer;
    }

    std::vector<char>::iterator fromBinary(std::vector<char>::iterator begin, std::vector<char>::iterator end);

    static double momentumTransfer(double energy, double angle);

    static constexpr double maxPhotonEnergy()
    {
        return 500.0;
    }
    static constexpr double minPhotonEnergy() { return 1.0; }
    static constexpr double MeVTokeV() { return 1000; }
    double barnToAtt()
    {
        constexpr double u = 1.6605402;
        return 1.0 / (u * m_atomicWeight);
    }

protected:
private:
    std::uint8_t m_Z = 0;
    double m_atomicWeight = 0.0;
    std::vector<std::pair<double, double>> m_coherent;
    std::vector<std::pair<double, double>> m_incoherent;
    std::vector<std::pair<double, double>> m_photoel;
    std::vector<std::pair<double, double>> m_formFactor;
    std::vector<std::pair<double, double>> m_imagAnomSF;
    std::vector<std::pair<double, double>> m_realAnomSF;
    std::vector<std::pair<double, double>> m_incoherentSF;

    std::map<std::uint8_t, AtomicShell> m_shells;
};
