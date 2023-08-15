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

#include "dxmc/beams/tube/betheHeitlerCrossSection.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/atomhandler.hpp"

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

    static constexpr T maxVoltage() { return T { 150 }; }
    static constexpr T minVoltage() { return T { 50 }; }

    T voltage() const { return m_voltage; }
    void setVoltage(T voltage)
    {
        m_voltage = std::clamp(voltage, minVoltage(), maxVoltage());
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

    void addFiltrationMaterial(const std::size_t Z, const T mm)
    {
        m_filtrationMaterials[Z] = std::abs(mm);
        m_hasCachedHVL = false;
    }

    const std::map<std::size_t, T>& filtrationMaterials() const { return m_filtrationMaterials; }

    void setAlFiltration(T mm)
    {
        addFiltrationMaterial(13, mm);
    }
    void setCuFiltration(T mm)
    {
        addFiltrationMaterial(29, mm);
    }
    void setSnFiltration(T mm)
    {
        addFiltrationMaterial(50, mm);
    }

    T filtration(std::size_t Z) const
    {
        if (m_filtrationMaterials.contains(Z))
            return m_filtrationMaterials.at(Z);
        return T { 0 };
    }

    T AlFiltration() const
    {
        return filtration(13);
    }
    T CuFiltration() const
    {
        return filtration(29);
    }
    T SnFiltration() const
    {
        return filtration(50);
    }

    void clearFiltrationMaterials()
    {
        m_filtrationMaterials.clear();
        m_hasCachedHVL = false;
    }

    void setEnergyResolution(T energyResolution)
    {
        m_energyResolution = std::clamp(energyResolution, T { 0.1 }, T { 10 });
        m_hasCachedHVL = false;
    }

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
        std::vector<T> specter(energies.size());
        const auto kVp = voltage();
        std::transform(std::execution::par_unseq, energies.begin(), energies.end(), specter.begin(), [kVp, anodeAngle](auto hv) -> T {
            return BetheHeitlerCrossSection::betheHeitlerSpectra(kVp, hv, anodeAngle);
        });

        // adding characteristic radiation
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
            // find closest energy
            auto eIdx = std::lower_bound(energyBegin, energyEnd, e);
            if (eIdx != energyEnd) {

                if (std::abs(e - *eIdx) <= T { 2.0 }) {
                    // we only add characteristic radiation if specter energy is closer than 2 keV from the K edge
                    auto nIdx = specter.begin();
                    std::advance(nIdx, std::distance(energyBegin, eIdx));
                    *nIdx = *nIdx + n; // adding characteristic k edge intensity to specter
                }
            }
        }
    }
    void filterSpecter(const std::vector<T>& energies, std::vector<T>& specter) const
    {
        std::vector<T> totAtt(energies.size(), T { 0 });
        for (auto const& [Z, mm] : m_filtrationMaterials) {
            const auto& atom = AtomHandler<T>::Atom(Z);
            const T cm = mm * T { 0.1 }; // for mm -> cm
            const auto dist = atom.standardDensity * cm;
            auto p = interpolate(atom.photoel, energies);
            auto in = interpolate(atom.incoherent, energies);
            auto co = interpolate(atom.coherent, energies);
            for (std::size_t i = 0; i < energies.size(); ++i) {
                totAtt[i] += (p[i] + in[i] + co[i]) * dist;
            }
        }
        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), totAtt.cbegin(), specter.begin(),
            [&](const auto s, const auto att) { return s * std::exp(-att); });
    }

    void normalizeSpecter(std::vector<T>& specter) const
    {
        const auto sum = std::reduce(std::execution::par_unseq, specter.cbegin(), specter.cend());
        std::for_each(std::execution::par_unseq, specter.begin(), specter.end(), [sum](auto& n) { n = n / sum; });
    }

    T calculatemmAlHalfValueLayer() const
    {
        const auto energy = getEnergy();
        const auto specter = getSpecter(energy);

        const auto& Al = AtomHandler<T>::Atom(13);

        const auto photo = interpolate(Al.photoel, energy);
        const auto incoherent = interpolate(Al.incoherent, energy);
        const auto coherent = interpolate(Al.coherent, energy);

        auto att = addVectors(photo, incoherent, coherent);

        const auto dens = Al.standardDensity;
        std::for_each(std::execution::par_unseq, att.begin(), att.end(), [dens](auto& a) { a *= dens; });

        // Boosted gradient decent for finding half value layer
        T x = T { 0.5 };
        T step = 0;
        int iter = 0;
        T g;

        do {
            x = x + step;
            g = std::transform_reduce(std::execution::par_unseq, specter.cbegin(), specter.cend(), att.cbegin(), T { 0 }, std::plus<T>(), [x](auto s, auto a) -> T { return s * std::exp(-a * x); });
            step = (g - T { 0.5 }) * std::max(5 - iter, 1);
        } while ((std::abs(g - T { 0.5 }) > T { 0.005 } && iter++ < 10));

        return x * T { 10.0 }; // cm -> mm
    }

private:
    T m_voltage = 120;
    T m_energyResolution = 1;
    T m_anodeAngle = 0.21f; // about 12 degrees
    T m_cachedHVL = 0;
    std::map<std::size_t, T> m_filtrationMaterials;
    bool m_hasCachedHVL = false;
};
}