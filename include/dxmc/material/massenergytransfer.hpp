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
         //return interpolate<T, false, false>(m_data, energy);
        const auto elog = std::log(energy);
        const auto r = CubicLSInterpolator<T>::evaluateSpline(elog, m_data2);
        return std::exp(r);
    }
    std::vector<T> operator()(const std::vector<T>& energy)
    {
         //return interpolate(m_data, energy);
        std::vector<T> r(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), r.begin(), [&](const auto e) {
            const auto elog = std::log(e);
            const auto r = CubicLSInterpolator<T>::evaluateSpline(elog, m_data2);
            return std::exp(r);
        });
        return r;
    }

    const std::vector<std::array<T, 3>>& getDataTable() const
    {
        return m_data2;
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

        std::vector<std::pair<T, T>> d_log(m_data.size());

        std::transform(std::execution::par_unseq, m_data.cbegin(), m_data.cend(), d_log.begin(), [](const auto& p) { const auto p_fix = p.second; return std::make_pair(std::log(p.first), std::log(p_fix)); });

        d_log.erase(std::remove_if(std::begin(d_log), std::end(d_log), [](const auto& d) { return std::isnan(d.second); }), std::end(d_log));

        auto inter = CubicLSInterpolator(d_log, 30, true);
        m_data2 = inter.getDataTable();
        return;
    }

    static void addVector(std::vector<T>& to, const std::vector<T>& from, const T fraction = T { 1 })
    {
        std::transform(std::execution::par_unseq, from.cbegin(), from.cend(), to.cbegin(), to.begin(), [fraction](const auto from, const auto to) { return to + fraction * from; });
    }

    static void appendMET(std::vector<T>& data, const std::vector<T>& energy, std::size_t Z, T fraction = T { 1 })
    {
        const auto& atom = AtomHandler<T>::Atom(Z);

        // photoelectric with fluro
        {
            const auto pe = interpolate(atom.photoel, energy);
            std::vector<T> fluro(pe.size(), T { 0 });
            for (const auto& [shIdx, shell] : atom.shells) {
                auto pe_shell = interpolate(shell.photoel, energy);
                std::transform(std::execution::par_unseq, pe.cbegin(), pe.cend(), pe_shell.cbegin(), pe_shell.begin(), [](const T p, const T psh) {
                    return psh / p;
                });
                std::transform(std::execution::par_unseq, pe_shell.cbegin(), pe_shell.cend(), energy.cbegin(), pe_shell.begin(), [&shell](const T p, const T e) {
                    if (e >= shell.bindingEnergy)
                        return p * shell.energyOfPhotonsPerInitVacancy;
                    else
                        return T { 0 };
                });
                addVector(fluro, pe_shell);
            }
            std::transform(
                std::execution::par_unseq, fluro.cbegin(), fluro.cend(), energy.cbegin(), fluro.begin(), [](const T f, const T e) {
                    return 1 - f / e;
                });
            std::transform(std::execution::par_unseq, pe.cbegin(), pe.cend(), fluro.cbegin(), fluro.begin(), std::multiplies<T>());
            addVector(data, fluro, fraction);
        }

        // compton with fluro
        {
            const auto Z = atom.Z;
            std::vector<T> inco_fluro(energy.size(), T { 0 });
            auto inco_scatterEn = interpolate(atom.incoherentMeanScatterEnergy, energy);
            for (const auto& [shIdx, shell] : atom.shells) {
                std::vector<T> inco_fluro_shell(energy.size());
                std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), inco_fluro_shell.begin(), [&shell, Z](const T e) {
                    if (e >= shell.bindingEnergy)
                        return shell.energyOfPhotonsPerInitVacancy / Z;
                    else
                        return T { 0 };
                });
                addVector(inco_fluro, inco_fluro_shell);
            }
            std::vector<T> inco_f(energy.size());
            // <E_scatter> + X_fluro
            std::transform(std::execution::par_unseq, inco_fluro.cbegin(), inco_fluro.cend(), inco_scatterEn.cbegin(), inco_f.begin(), std::plus<>());

            // 1-(<E_scatter> + X_fluro)/E
            std::transform(std::execution::par_unseq, inco_f.cbegin(), inco_f.cend(), energy.cbegin(), inco_f.begin(), [](const T f, const T e) { return 1 - f / e; });

            auto inco_u = interpolate(atom.incoherent, energy);
            // u * (1-(<E_scatter> + X_fluro)/E)
            std::transform(std::execution::par_unseq, inco_f.cbegin(), inco_f.cend(), inco_u.cbegin(), inco_f.begin(), std::multiplies<T>());

            addVector(data, inco_f, fraction);
        }
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
    std::vector<std::array<T, 3>> m_data2;
    std::vector<std::pair<T, T>> m_data;
};

}