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

#include "dxmc/floating.h"
#include "dxmc/material.h"

#include <array>
#include <cmath>
#include <vector>

namespace dxmc {
template <Floating T>
class AttenuationLutInterpolator {
private:
    constexpr T shell_offset()
    {
        return T { 0.0001 };
    }
    std::vector<T> m_x;
    std::vector<T> m_coefficients;

    std::size_t m_resolution = 40;
    std::size_t m_linearIndex = 0;
    T m_linearStep = 0;
    T m_linearEnergy = 0;

protected:
public:
    AttenuationLutInterpolator(const std::vector<Material>& materials, const T minEnergy, const T maxEnergy, std::size_t resolution = 40)
    {
        m_resolution = std::max(static_cast<std::size_t>(maxEnergy - minEnergy), std::size_t { 50 });
        const T logMinEnergy = std::log10(minEnergy);
        const T logMaxEnergy = std::log10(maxEnergy);
        //https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation
        // get Shell edge energies
        m_x.clear();
        std::vector<T> binding_energies;
        for (const auto& mat : materials) {
            const auto minEnergyD = static_cast<double>(minEnergy);
            const auto bindingEnergies = mat.getBindingEnergies(minEnergyD);
            for (const auto e : bindingEnergies) {
                const T bindingE = static_cast<T>(e);
                binding_energies.push_back(std::log10(bindingE - shell_offset()));
                binding_energies.push_back(std::log10(bindingE));
            }
        }
        std::sort(binding_energies.begin(), binding_energies.end());
        //remove duplicates
        binding_energies.erase(std::unique(binding_energies.begin(), binding_energies.end()), binding_energies.end());

        m_x.reserve(m_resolution + 1 + binding_energies.size());

        for (std::size_t i = 0; i <= m_resolution; ++i) {
            const T e = logMinEnergy + (i * (logMaxEnergy - logMinEnergy)) / (m_resolution - 1);
            m_x.push_back(e);
        }
        // adding binding energies
        for (const auto e : binding_energies) {
            m_x.push_back(e);
        }

        std::sort(m_x.begin(), m_x.end());
        //erasing duplicates
        m_x.erase(std::unique(m_x.begin(), m_x.end()), m_x.end());

        m_coefficients.resize(materials.size() * (m_x.size() - 1) * 2 * 3);
        std::fill(m_coefficients.begin(), m_coefficients.end(), T { 0 });

        for (std::size_t matIdx = 0; matIdx < materials.size(); ++matIdx) {
            for (std::size_t typeIdx = 0; typeIdx < 3; ++typeIdx) {
                std::vector<T> y(m_x.size());
                if (typeIdx == 0)
                    std::transform(m_x.cbegin(), m_x.cend(), y.begin(), [&](const auto le) -> T { return std::log10(materials[matIdx].getPhotoelectricAttenuation(std::pow(T { 10 }, le))); });
                else if (typeIdx == 1) {
                    std::transform(m_x.cbegin(), m_x.cend(), y.begin(), [&](const auto le) -> T { return std::log10(materials[matIdx].getComptonAttenuation(std::pow(T { 10 }, le))); });
                } else {
                    std::transform(m_x.cbegin(), m_x.cend(), y.begin(), [&](const auto le) -> T { return std::log10(materials[matIdx].getRayleightAttenuation(std::pow(T { 10 }, le))); });
                }

                for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
                    const auto offset = matIdx * (m_x.size() - 1) * 6 + typeIdx * 2 + i * 6;
                    const T w = (y[i + 1] - y[i]) / (m_x[i + 1] - m_x[i]);
                    const T a = w;
                    const T b = y[i] - m_x[i] * w;
                    m_coefficients[offset] = b;
                    m_coefficients[offset + 1] = a;
                }
            }
        }
        m_x.erase(m_x.begin());
        if (binding_energies.size() > 0) {
            const auto pos = std::upper_bound(m_x.cbegin(), m_x.cend(), binding_energies.back());
            if (pos == m_x.end()) {
                m_linearIndex = m_x.size();
                m_linearStep = m_x[m_x.size() - 1] - m_x[m_x.size() - 2];
                m_linearEnergy = m_x.back();
            } else {
                m_linearIndex = std::distance(m_x.cbegin(), pos);
                m_linearStep = m_x[m_linearIndex + 1] - m_x[m_linearIndex];
                m_linearEnergy = m_x[m_linearIndex];
            }
        } else {
            m_linearIndex = 0;
            m_linearStep = m_x[1] - m_x[0];
            m_linearEnergy = m_x.front();
        }
        m_resolution = m_x.size();
    }

    std::array<T, 3> operator()(const std::size_t materialIdx, const T energy) const
    {
        std::array<T, 3> res;
        const T logEnergy = std::log10(energy);
        if (logEnergy > m_linearEnergy)
            [[likely]]
            {
                const auto index = static_cast<std::size_t>((logEnergy - m_linearEnergy) / m_linearStep) + m_linearIndex;
                const auto offset = materialIdx * m_resolution * 6 + index * 6;
                for (std::size_t i = 0; i < 3; ++i) {
                    res[i] = std::pow(T { 10 }, m_coefficients[offset + 2 * i] + m_coefficients[offset + 2 * i + 1] * logEnergy);
                }
            }
        else
            [[unlikly]]
            {
                const auto pos = std::upper_bound(m_x.cbegin(), m_x.cend(), logEnergy);
                const auto index = pos != m_x.cend() ? std::distance(m_x.cbegin(), pos) : m_resolution - 1;
                const auto offset = materialIdx * m_resolution * 6 + index * 6;
                for (std::size_t i = 0; i < 3; ++i) {
                    res[i] = std::pow(T { 10 }, m_coefficients[offset + 2 * i] + m_coefficients[offset + 2 * i + 1] * logEnergy);
                }
            }

        return res;
    }
};
}
