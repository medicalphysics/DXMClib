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
#include "serialize.hpp"

#include <cmath>

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
        const double e = data[ind] * MeVTokeV();
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
        const double e = data[ind] * MeVTokeV();
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * barnToAtt();
            m_incoherent.push_back(std::make_pair(e, a));
        }
    }
    m_incoherent.shrink_to_fit();
}

void AtomicElement::setFormFactor(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_formFactor.clear();
    m_formFactor.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E6; // for units of per Å
        const auto f = data[i + 1];
        m_formFactor.push_back(std::make_pair(x, f));
    }
}

void AtomicElement::setImaginaryAnomalousSF(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_imagAnomSF.clear();
    m_imagAnomSF.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E3;
        const auto f = data[i + 1];
        m_imagAnomSF.push_back(std::make_pair(x, f));
    }
}

void AtomicElement::setRealAnomalousSF(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_realAnomSF.clear();
    m_realAnomSF.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E3;
        const auto f = data[i + 1];
        m_realAnomSF.push_back(std::make_pair(x, f));
    }
}

void AtomicElement::setIncoherentSF(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_incoherentSF.clear();
    m_incoherentSF.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E6;
        const auto f = data[i + 1];
        m_incoherentSF.push_back(std::make_pair(x, f));
    }
}

void AtomicElement::setPhotoelectricData(const std::vector<double>& data)
{
    const auto c = barnToAtt();
    const auto N = data.size() / 2;
    m_photoel.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * MeVTokeV();
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * c;
            m_photoel.push_back(std::make_pair(e, a));
        }
    }
    m_photoel.shrink_to_fit();
}
void AtomicElement::setShellPhotoelectricData(std::uint8_t shell, const std::vector<double>& data)
{
    if (!m_shells.contains(shell)) {
        m_shells[shell] = AtomicShell(shell);
    }

    const auto N = data.size() / 2;
    std::vector<std::pair<double, double>> photoel;
    photoel.reserve(N);

    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * MeVTokeV();
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * barnToAtt();
            photoel.push_back(std::make_pair(e, a));
        }
    }
    photoel.shrink_to_fit();
    m_shells[shell].setPhotoelectricData(photoel);
}

void AtomicElement::setShellBindingEnergy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint8_t>(data[i]);
        const auto E = data[i + 1] * MeVTokeV();
        if (!m_shells.contains(shell)) {
            m_shells[shell] = AtomicShell(shell);
        }
        m_shells[shell].setBindingEnergy(E);
    }
}

void AtomicElement::setShellNumberOfElectrons(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint8_t>(data[i]);
        const auto N = data[i + 1];
        if (!m_shells.contains(shell)) {
            m_shells[shell] = AtomicShell(shell);
        }
        m_shells[shell].setNumberOfElectrons(N);
    }
}

void AtomicElement::setShellNumberOfPhotonsPerInitVacancy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint8_t>(data[i]);
        const auto N = data[i + 1];
        if (!m_shells.contains(shell)) {
            m_shells[shell] = AtomicShell(shell);
        }
        m_shells[shell].setNumberOfPhotonsPerInitVacancy(N);
    }
}

void AtomicElement::setShellEnergyOfPhotonsPerInitVacancy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint8_t>(data[i]);
        const auto E = data[i + 1] * MeVTokeV();
        if (!m_shells.contains(shell)) {
            m_shells[shell] = AtomicShell(shell);
        }
        m_shells[shell].setEnergyOfPhotonsPerInitVacancy(E);
    }
}

std::vector<char> AtomicElement::toBinary() const
{
    std::vector<char> buffer(sizeof(std::uint64_t));
    serialize(m_Z, buffer);
    serialize(m_atomicWeight, buffer);
    serialize(m_coherent, buffer);
    serialize(m_incoherent, buffer);
    serialize(m_photoel, buffer);
    serialize(m_formFactor, buffer);
    serialize(m_incoherentSF, buffer);

    // adding shells
    // adding number of shells
    std::uint64_t n_shells = m_shells.size();
    serialize(n_shells, buffer);
    // adding each shell
    for (const auto& [id, shell] : m_shells) {
        auto s_buffer = shell.toBinary();
        serialize(s_buffer, buffer);
    }

    std::uint64_t buffer_size = buffer.size() - sizeof(std::uint64_t);
    auto buffer_size_addr = reinterpret_cast<char*>(&buffer_size);
    // writing over first value to contain size of serialized data
    std::copy(buffer_size_addr, buffer_size_addr + sizeof(std::uint64_t), buffer.begin());

    return buffer;
}
char* AtomicElement::fromBinary(std::vector<char>& data, char* begin)
{
    std::uint64_t size;
    auto start = deserialize(size, begin);
    start = deserialize(m_Z, start);
    start = deserialize(m_atomicWeight, start);
    start = deserialize(m_coherent, start);
    start = deserialize(m_incoherent, start);
    start = deserialize(m_photoel, start);
    start = deserialize(m_formFactor, start);
    start = deserialize(m_incoherentSF, start);

    std::uint64_t n_shells;
    m_shells.clear();
    start = deserialize(n_shells, start);
    for (std::size_t i = 0; i < n_shells; ++i) {
        AtomicShell shell;
        start = shell.fromBinary(data, start);
        m_shells[shell.shell()] = shell;
    }

    return start;
}
double AtomicElement::momentumTransfer(double energy, double angle)
{
    constexpr double hc_si = 1.239841193E-6; // ev*m
    constexpr double m2A = 1E10; // meters to Ångstrøm
    constexpr double eV2keV = 1E-3; // eV to keV
    constexpr double hc = hc_si * m2A * eV2keV; // kev*Å
    constexpr double hc_inv = 1.0 / hc;
    return energy * std::sin(angle * 0.5) * hc_inv; // per Å
}