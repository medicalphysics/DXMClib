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

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/atomhandler.hpp"
#include "dxmc/material/nistmaterials.hpp"

#include <algorithm>
#include <execution>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace dxmc {

template <Floating T>
class MassEnergyTransfer {

public:
    MassEnergyTransfer(const std::map<std::size_t, T>& massFractions)
    {
        constructInterpolationTable(massFractions);
    }
    MassEnergyTransfer(const std::size_t Z)
    {
        std::map<std::size_t, T> massFractions;
        massFractions[Z] = T { 1 };
        constructInterpolationTable(massFractions);
    }
    MassEnergyTransfer(const std::string& nistName)
    {
        const auto& massFractions = NISTMaterials<T>::Composition(nistName);
        if (massFractions.empty())
            return;
        constructInterpolationTable(massFractions);
    }

    T operator()(const T energy)
    {
        return interpolate<T, false, false>(m_data, energy);
    }
    std::vector<T> operator()(const std::vector<T>& energy)
    {
        return interpolate(m_data, energy);
    }

protected:
    void constructInterpolationTable(const std::map<std::size_t, T>& massFractions)
    {
        std::vector<std::size_t> Zs;
        T sumMass = 0;
        for (const auto& [Z, massFrac] : massFractions) {
            Zs.push_back(Z);
            sumMass += massFrac;
        }
        auto energy = energyGrid(Zs);
        std::vector<T> data(energy.size(), T { 0 });

        for (const auto& [Z, massFrac] : massFractions) {
            appendMET(data, energy, Z, massFrac / sumMass);
        }

        m_data.resize(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), data.cbegin(), m_data.begin(), [](const T e, const T d) { return std::make_pair(e, d); });

        return;
    }

    static void addVector(std::vector<T>& to, const std::vector<T>& from, const T fraction = T { 1 })
    {
        std::transform(std::execution::par_unseq, from.cbegin(), from.cend(), to.cbegin(), to.begin(), [fraction](const auto from, const auto to) { return to + fraction * from; });
    }

    static void appendMET(std::vector<T>& data, const std::vector<T>& energy, std::size_t Z, T fraction = T { 1 })
    {
        const auto& atom = AtomHandler<T>::Atom(Z);

        // photoelectric ignoring fluro
        auto fpe = interpolate(atom.photoel, energy);
        addVector(data, fpe, fraction);

        // compton ignoring fluro
        auto finco_scatterEn = interpolate(atom.incoherentMeanScatterEnergy, energy);
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), finco_scatterEn.cbegin(), finco_scatterEn.begin(), [](const T e, const T scatter_e) {
            return T { 1 } - scatter_e / e;
        });
        auto finco = interpolate(atom.incoherent, energy);
        std::transform(std::execution::par_unseq, finco.cbegin(), finco.cend(), finco_scatterEn.cbegin(), finco.begin(), std::multiplies<T>());
        addVector(data, finco, fraction);
    }

    static std::vector<T> energyGrid(const std::vector<std::size_t>& atomicNumbers)
    {
        std::vector<std::vector<std::pair<T, T>>> data;
        for (const auto Z : atomicNumbers) {
            const auto& element = AtomHandler<T>::Atom(Z);
            data.push_back(element.incoherent);
            data.push_back(element.photoel);
            data.push_back(element.incoherentMeanScatterEnergy);
            for (const auto& [shIdx, shell] : element.shells) {
                data.push_back(shell.photoel);
            }
        }
        return energyGrid(data);
    }

    static std::vector<T> energyGrid(const std::size_t atomicNumber)
    {
        const auto& element = AtomHandler<T>::Atom(atomicNumber);

        std::vector<std::vector<std::pair<T, T>>> data;

        data.push_back(element.incoherent);
        data.push_back(element.photoel);
        data.push_back(element.incoherentMeanScatterEnergy);
        for (const auto& [shIdx, shell] : element.shells) {
            data.push_back(shell.photoel);
        }
        return energyGrid(data);
    }

    static std::vector<T> energyGrid(const std::vector<std::vector<std::pair<T, T>>>& data)
    {
        std::size_t N = 0;
        for (const auto& vec : data) {
            N += vec.size();
        }
        std::vector<T> energy(N);
        auto begin = energy.begin();
        for (const auto& vec : data) {
            begin = std::transform(std::execution::par_unseq, vec.cbegin(), vec.cend(), begin, [](const auto& pair) {
                return pair.first;
            });
        }
        std::sort(energy.begin(), energy.end());
        // removing duplicates

        auto erase_from = std::unique(std::execution::par_unseq, energy.begin(), energy.end(), [](const auto& lh, const auto& rh) {
            return std::abs(lh - rh) <= (std::numeric_limits<T>::epsilon() * 5);
        });
        if (std::distance(erase_from, energy.end()) != 0) {
            energy.erase(erase_from, energy.end());
        }
        return energy;
    }

private:
    std::vector<std::pair<T, T>> m_data;
};

}