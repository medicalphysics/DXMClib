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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include "dxmc/betheHeitlerCrossSection.h"
#include "dxmc/constants.h"
#include "dxmc/floating.h"
#include "dxmc/material.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <utility>
#include <vector>

namespace dxmc {

template <Floating T = double>
class Tube {
public:
    Tube(T tubeVoltage = 120.0, T anodeAngleDeg = 12.0, T energyResolution = 1.0)
        : m_voltage(tubeVoltage)
        , m_energyResolution(energyResolution)
    {
        setAnodeAngleDeg(anodeAngleDeg);
        m_hasCachedHVL = false;
    }
    Tube(const Tube<T>& other)
    {
        m_anodeAngle = other.anodeAngle();
        m_voltage = other.voltage();
        m_energyResolution = other.energyResolution();
        m_filtrationMaterials = other.filtrationMaterials();
        m_hasCachedHVL = false;
    }
    static constexpr T maxVoltage() { return T { 150 }; }
    static constexpr T minVoltage() { return T { 50 }; }

    T voltage() const { return m_voltage; }
    void setVoltage(T voltage)
    {
        m_voltage = std::min(std::max(voltage, minVoltage()), maxVoltage());
        m_hasCachedHVL = false;
    }

    T anodeAngle() const { return m_anodeAngle; }
    T anodeAngleDeg() const { return m_anodeAngle * RAD_TO_DEG<T>(); }
    void setAnodeAngle(T angle)
    {
        auto a = std::abs(angle);
        if (a > PI_VAL<T>() * T { 0.5 })
            a = PI_VAL<T>() * T { 0.5 };
        m_anodeAngle = a;
        m_hasCachedHVL = false;
    }
    void setAnodeAngleDeg(T angle)
    {
        setAnodeAngle(angle * DEG_TO_RAD<T>());
    }

    void addFiltrationMaterial(const Material& filtrationMaterial, T mm)
    {
        m_filtrationMaterials.push_back(std::make_pair(filtrationMaterial, std::abs(mm)));
        m_hasCachedHVL = false;
    }
    std::vector<std::pair<Material, T>>& filtrationMaterials() { return m_filtrationMaterials; }
    const std::vector<std::pair<Material, T>>& filtrationMaterials() const { return m_filtrationMaterials; }

    void setAlFiltration(T mm)
    {
        bool hasAlFiltration = false;
        bool alFiltrationSet = false;
        for (auto& [material, thickness] : m_filtrationMaterials)
            if (material.name().compare("Al") == 0)
                if (!alFiltrationSet) {
                    thickness = std::abs(mm);
                    hasAlFiltration = true;
                    alFiltrationSet = true;
                }
        if (!hasAlFiltration) {
            Material al(13);
            addFiltrationMaterial(al, std::abs(mm));
        }
    }
    void setCuFiltration(T mm)
    {
        bool hasCuFiltration = false;
        bool cuFiltrationSet = false;
        for (auto& [material, thickness] : m_filtrationMaterials)
            if (material.name().compare("Cu") == 0)
                if (!cuFiltrationSet) {
                    thickness = std::abs(mm);
                    hasCuFiltration = true;
                    cuFiltrationSet = true;
                }
        if (!hasCuFiltration) {
            Material cu(29);
            addFiltrationMaterial(cu, std::abs(mm));
        }
    }
    void setSnFiltration(T mm)
    {
        bool hasSnFiltration = false;
        bool snFiltrationSet = false;
        for (auto& [material, thickness] : m_filtrationMaterials)
            if (material.name().compare("Sn") == 0)
                if (!snFiltrationSet) {
                    thickness = std::abs(mm);
                    hasSnFiltration = true;
                    snFiltrationSet = true;
                }
        if (!hasSnFiltration) {
            Material sn(50);
            addFiltrationMaterial(sn, std::abs(mm));
        }
    }
    T AlFiltration() const
    {
        for (auto& [material, thickness] : m_filtrationMaterials) {
            if (material.name().compare("Al") == 0)
                return thickness;
        }
        return T { 0 };
    }
    T CuFiltration() const
    {
        for (auto& [material, thickness] : m_filtrationMaterials) {
            if (material.name().compare("Cu") == 0)
                return thickness;
        }
        return T { 0 };
    }
    T SnFiltration() const
    {
        for (auto& [material, thickness] : m_filtrationMaterials) {
            if (material.name().compare("Sn") == 0)
                return thickness;
        }
        return T { 0 };
    }

    void clearFiltrationMaterials() { m_filtrationMaterials.clear(); }

    void setEnergyResolution(T energyResolution) { m_energyResolution = energyResolution; }

    T energyResolution() const { return m_energyResolution; }

    std::vector<T> getEnergy() const
    {
        std::vector<T> energies;
        T hv = m_energyResolution;
        const auto n_elem = static_cast<std::size_t>(std::ceil(m_voltage / m_energyResolution));
        energies.reserve(n_elem);
        while (hv <= m_voltage) {
            energies.push_back(hv);
            hv = hv + m_energyResolution;
        }
        return energies;
    }

    std::vector<std::pair<T, T>> getSpecter(bool normalize = true) const
    {
        auto energies = getEnergy();

        auto specter = this->getSpecter(energies, normalize);

        std::vector<std::pair<T, T>> map;
        map.reserve(specter.size());
        for (std::size_t i = 0; i < specter.size(); ++i)
            map.push_back(std::make_pair(energies[i], specter[i]));
        return map;
    }

    std::vector<T> getSpecter(const std::vector<T>& energies, const T anodeAngle, bool normalize = true) const
    {
        std::vector<T> specter;
        specter.resize(energies.size());
        std::transform(std::execution::par_unseq, energies.begin(), energies.end(), specter.begin(), [&](auto hv) -> T {
            const auto vd = this->voltage();
            const auto bh_t = BetheHeitlerCrossSection::betheHeitlerSpectra(vd, hv, anodeAngle);
            return bh_t;
        });

        //adding characteristic radiation
        addCharacteristicEnergy(energies, specter);
        filterSpecter(energies, specter);
        if (normalize) {
            normalizeSpecter(specter);
        }
        return specter;
    }
    std::vector<T> getSpecter(const std::vector<T>& energies, bool normalize = true) const
    {
        return getSpecter(energies, m_anodeAngle, normalize);
    }

    T mmAlHalfValueLayer()
    {
        if (!m_hasCachedHVL) {
            m_cachedHVL = calculatemmAlHalfValueLayer();
            m_hasCachedHVL = true;
        }
        return m_cachedHVL;
    }
    [[nodiscard]] T mmAlHalfValueLayer() const
    {
        if (!m_hasCachedHVL) {
            return calculatemmAlHalfValueLayer();
        }
        return m_cachedHVL;
    }

protected:
    void addCharacteristicEnergy(const std::vector<T>& energy, std::vector<T>& specter) const
    {
        auto energyBegin = energy.begin();
        auto energyEnd = energy.end();

        const auto voltage_d = voltage();
        const auto anode_d = anodeAngle();

        const auto kEdge = BetheHeitlerCrossSection::characteristicTungstenKedge(voltage_d, anode_d);
        for (const auto& [e, n] : kEdge) {
            //find closest energy
            auto eIdx = std::lower_bound(energyBegin, energyEnd, e);
            if (eIdx != energyEnd) {

                if (std::abs(e - *eIdx) <= T { 2.0 }) { // we only add characteristic radiation if specter energy is closer than 2 keV from the K edge
                    auto nIdx = specter.begin();
                    std::advance(nIdx, std::distance(energyBegin, eIdx));
                    *nIdx = *nIdx + n; // adding characteristic k edge intensity to specter
                }
            }
        }
    }
    void filterSpecter(const std::vector<T>& energies, std::vector<T>& specter) const
    {
        for (auto const& [material, mm] : m_filtrationMaterials) {
            const T cm = mm * T { 0.1 }; //for mm -> cm
            std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), energies.cbegin(), specter.begin(),
                [&, material = material](const auto n, const auto e) -> T { return n * std::exp(-material.getTotalAttenuation(e) * material.standardDensity() * cm); });
        }
    }
    void normalizeSpecter(std::vector<T>& specter) const
    {
        const auto sum = std::reduce(std::execution::par_unseq, specter.begin(), specter.end());
        std::for_each(std::execution::par_unseq, specter.begin(), specter.end(), [=](auto& n) { n = n / sum; });
    }

    T calculatemmAlHalfValueLayer() const
    {
        auto energy = getEnergy();
        auto specter = getSpecter(energy);

        Material al(13);
        std::vector<T> att(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), att.begin(), [&](auto e) -> T {
            const T dens = static_cast<T>(al.standardDensity());
            const T att = static_cast<T>(al.getTotalAttenuation(e));
            return att * dens;
        });

        T x { 0.5 };
        T step { 1 };
        T g;
        do {
            g = std::transform_reduce(std::execution::par_unseq, specter.cbegin(), specter.cend(), att.cbegin(), T { 0 }, std::plus<T>(), [=](auto s, auto a) -> T { return s * std::exp(-a * x); });
            step = g - T { 0.5 };
            x = x + step;

        } while (std::abs(step) > T { 0.01 });

        return x * T { 10.0 }; // cm -> mm
    }

private:
    T m_voltage, m_energyResolution, m_anodeAngle;
    T m_cachedHVL = 0;
    bool m_hasCachedHVL = false;
    std::vector<std::pair<Material, T>> m_filtrationMaterials;
};
}