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

#include "atomicelementhandler.hpp"

#include <cmath>

AtomicElementHandler::AtomicElementHandler(std::uint64_t Z)
{
    m_atom.Z = Z;
}

void AtomicElementHandler::setCoherentData(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_atom.coherent.clear();
    m_atom.coherent.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * MeVTokeV();
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * barnToAtt();
            m_atom.coherent.push_back(std::make_pair(e, a));
        }
    }
    m_atom.coherent.shrink_to_fit();
}
void AtomicElementHandler::setIncoherentData(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_atom.incoherent.clear();
    m_atom.incoherent.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * MeVTokeV();
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * barnToAtt();
            m_atom.incoherent.push_back(std::make_pair(e, a));
        }
    }
    m_atom.incoherent.shrink_to_fit();
}

void AtomicElementHandler::setFormFactor(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_atom.formFactor.clear();
    m_atom.formFactor.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E6; // for units of per Å
        const auto f = data[i + 1];
        m_atom.formFactor.push_back(std::make_pair(x, f));
    }
}
/*
void AtomicElementHandler::setImaginaryAnomalousSF(const std::vector<double>& data)
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

void AtomicElementHandler::setRealAnomalousSF(const std::vector<double>& data)
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
*/
void AtomicElementHandler::setIncoherentSF(const std::vector<double>& data)
{
    const auto N = data.size() / 2;
    m_atom.incoherentSF.clear();
    m_atom.incoherentSF.reserve(N);
    for (std::size_t i = 0; i < N * 2; i = i + 2) {
        const auto x = data[i] * 1E6; // for units of per Å
        const auto f = data[i + 1];
        m_atom.incoherentSF.push_back(std::make_pair(x, f));
    }
}

void AtomicElementHandler::setPhotoelectricData(const std::vector<double>& data)
{
    const auto c = barnToAtt();
    const auto N = data.size() / 2;
    m_atom.photoel.clear();
    m_atom.photoel.reserve(N);
    for (std::size_t i = 0; i < N; ++i) {
        const auto ind = i * 2;
        const double e = data[ind] * MeVTokeV();
        if (e >= minPhotonEnergy() && e <= maxPhotonEnergy()) {
            const double a = data[ind + 1] * c;
            m_atom.photoel.push_back(std::make_pair(e, a));
        }
    }
    m_atom.photoel.shrink_to_fit();
}
void AtomicElementHandler::setShellPhotoelectricData(std::uint64_t shell, const std::vector<double>& data)
{
    if (!m_atom.shells.contains(shell)) {
        m_atom.shells[shell] = dxmc::AtomicShell<double>(shell);
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
    m_atom.shells[shell].photoel = photoel;
}

void AtomicElementHandler::setShellBindingEnergy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto E = data[i + 1] * MeVTokeV();
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = dxmc::AtomicShell<double>(shell);
        }
        m_atom.shells[shell].bindingEnergy = E;
    }
}

void AtomicElementHandler::setShellNumberOfElectrons(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto N = data[i + 1];
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = dxmc::AtomicShell<double>(shell);
        }
        m_atom.shells[shell].numberOfElectrons = N;
    }
}

void AtomicElementHandler::setShellNumberOfPhotonsPerInitVacancy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto N = data[i + 1];
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = dxmc::AtomicShell<double>(shell);
        }
        m_atom.shells[shell].numberOfPhotonsPerInitVacancy = N;
    }
}

void AtomicElementHandler::setShellEnergyOfPhotonsPerInitVacancy(const std::vector<double>& data)
{
    for (std::size_t i = 0; i < data.size(); i = i + 2) {
        const auto shell = static_cast<std::uint64_t>(data[i]);
        const auto E = data[i + 1] * MeVTokeV();
        if (!m_atom.shells.contains(shell)) {
            m_atom.shells[shell] = dxmc::AtomicShell<double>(shell);
        }
        m_atom.shells[shell].energyOfPhotonsPerInitVacancy = E;
    }
}

bool AtomicElementHandler::operator==(const AtomicElementHandler& other) const
{
    bool eq = m_atom.Z == other.m_atom.Z;
    eq = eq && m_atom.atomicWeight == other.m_atom.atomicWeight;
    eq = eq && m_atom.photoel == other.m_atom.photoel;
    eq = eq && m_atom.coherent == other.m_atom.coherent;
    eq = eq && m_atom.incoherent == other.m_atom.incoherent;
    eq = eq && m_atom.incoherentSF == other.m_atom.incoherentSF;
    eq = eq && m_atom.formFactor == other.m_atom.formFactor;

    for (const auto& [key, shell] : m_atom.shells) {
        if (!other.m_atom.shells.contains(key))
            return false;
        const auto& lshell = m_atom.shells.at(key);
        const auto& rshell = other.m_atom.shells.at(key);
        eq = eq && lshell.shell == rshell.shell;
        eq = eq && lshell.bindingEnergy == rshell.bindingEnergy;
        eq = eq && lshell.energyOfPhotonsPerInitVacancy == rshell.energyOfPhotonsPerInitVacancy;
        eq = eq && lshell.HartreeFockOrbital_0 == rshell.HartreeFockOrbital_0;
        eq = eq && lshell.numberOfElectrons == rshell.numberOfElectrons;
        eq = eq && lshell.numberOfPhotonsPerInitVacancy == rshell.numberOfPhotonsPerInitVacancy;
        eq = eq && lshell.photoel == rshell.photoel;
    }
    return eq;
}
double AtomicElementHandler::momentumTransfer(double energy, double angle)
{
    constexpr double hc_si = 1.239841193E-6; // ev*m
    constexpr double m2A = 1E10; // meters to Ångstrøm
    constexpr double eV2keV = 1E-3; // eV to keV
    constexpr double hc = hc_si * m2A * eV2keV; // kev*Å
    constexpr double hc_inv = 1.0 / hc;
    return energy * std::sin(angle * 0.5) * hc_inv; // per Å
}