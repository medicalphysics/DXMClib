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

#include "dxmc/floating.hpp"
#include "dxmc/material.hpp"
#include "dxmc/world.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace dxmc {
template <Floating T>
class AttenuationLutInterpolator {
private:
    constexpr T shell_offset()
    {
        return T { 0.001 };
    }
    std::vector<T> m_x;
    std::vector<T> m_coefficients;
    std::vector<T> m_maxCoefficients;

    std::size_t m_resolution = 40;
    std::size_t m_linearIndex = 0;
    T m_linearStep = 0;
    T m_linearEnergy = 0;
    T m_removedEnergy = 0;

protected:
    template <typename DensIter, typename MatIter>
    void generateMaxAttenuationLut(const std::vector<Material>& materials, const DensIter densBegin, const DensIter densEnd, const MatIter matBegin, const T firstEnergy)
    {
        std::vector<T> maxDens(materials.size());
        for (std::size_t matIdx = 0; matIdx < maxDens.size(); ++matIdx) {
            maxDens[matIdx] = std::transform_reduce(
                std::execution::par_unseq,
                densBegin, densEnd, matBegin, T { 0 }, [](const auto l, const auto r) { return std::max(r, l); },
                [=](const auto d, const auto mIdx) -> T {
                    return mIdx == matIdx ? d : T { 0 };
                });
        }

        auto x = m_x;
        x.insert(x.cbegin(), firstEnergy);
        std::vector<T> y(x.size());

        for (std::size_t xIdx = 0; xIdx < x.size(); ++xIdx) {
            T maxVal = 0;
            for (std::size_t matIdx = 0; matIdx < maxDens.size(); ++matIdx) {
                const auto vals = this->operator()(matIdx, std::pow(T { 10 }, x[xIdx]));
                const auto val = maxDens[matIdx] * std::reduce(vals.cbegin(), vals.cend());
                maxVal = std::max(maxVal, val);
            }
            y[xIdx] = std::log10(1 / maxVal);
        }

        m_maxCoefficients.resize(m_x.size() * 2);
        for (std::size_t i = 0; i < x.size() - 1; ++i) {
            const auto offset = i * 2;
            const T w = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
            const T a = w;
            const T b = y[i] - x[i] * w;
            m_maxCoefficients[offset] = b;
            m_maxCoefficients[offset + 1] = a;
        }
    }

public:
    AttenuationLutInterpolator() { }
    AttenuationLutInterpolator(const World<T>& world, T maxEnergy, T minEnergy)
    {
        minEnergy = std::max(T { 0.1 }, minEnergy);
        maxEnergy = std::max(maxEnergy, T { 50 });
        const auto& materials = world.materialMap();
        auto densBegin = world.densityArray()->cbegin();
        auto densEnd = world.densityArray()->cend();
        auto matBegin = world.materialIndexArray()->cbegin();
        generate(materials, densBegin, densEnd, matBegin, maxEnergy, minEnergy);
    }
    AttenuationLutInterpolator(const std::vector<Material>& materials, T maxEnergy, T minEnergy)
    {
        minEnergy = std::max(T { 0.1 }, minEnergy);
        maxEnergy = std::min(maxEnergy, T { 50 });
        std::vector<T> dens(materials.size());
        std::transform(materials.cbegin(), materials.cend(), dens.begin(), [](const auto& m) { return static_cast<T>(m.standardDensity()); });
        std::vector<std::uint8_t> materialIndex(materials.size());
        std::iota(materialIndex.begin(), materialIndex.end(), 0);
        auto densBegin = dens.cbegin();
        auto densEnd = dens.cend();
        auto matBegin = materialIndex.cbegin();
        generate(materials, densBegin, densEnd, matBegin, maxEnergy, minEnergy);
    }
    template <typename DensIter, typename MatIter>
    AttenuationLutInterpolator(const std::vector<Material>& materials, const DensIter densBegin, const DensIter densEnd, const MatIter matBegin, T maxEnergy, T minEnergy)
    {
        minEnergy = std::max(T { 0.1 }, minEnergy);
        maxEnergy = std::min(maxEnergy, T { 50 });
        generate(materials, densBegin, densEnd, matBegin, maxEnergy, minEnergy);
    }

    template <typename DensIter, typename MatIter>
    void generate(const std::vector<Material>& materials, const DensIter densBegin, const DensIter densEnd, const MatIter matBegin, const T maxEnergy, const T minEnergy)

    {
        m_resolution = std::max(static_cast<std::size_t>(maxEnergy - minEnergy), std::size_t { 10 });
        const T logMinEnergy = std::log10(minEnergy);
        const T logMaxEnergy = std::log10(maxEnergy);

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
            const T e = logMinEnergy + (i * (logMaxEnergy - logMinEnergy)) / (m_resolution);
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

        const T missingEnergy = m_x.front();
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

        generateMaxAttenuationLut(materials, densBegin, densEnd, matBegin, missingEnergy);
    }

    std::array<T, 3> operator()(const std::size_t materialIdx, const T energy) const
    {
        std::array<T, 3> res;
        const T logEnergy = std::log10(energy);
        if (logEnergy > m_linearEnergy)
            [[likely]] {
            const auto index = std::min(static_cast<std::size_t>((logEnergy - m_linearEnergy) / m_linearStep) + m_linearIndex, m_resolution - 1);
            const auto offset = materialIdx * m_resolution * 6 + index * 6;
            for (std::size_t i = 0; i < 3; ++i) {
                res[i] = std::pow(T { 10 }, m_coefficients[offset + 2 * i] + m_coefficients[offset + 2 * i + 1] * logEnergy);
            }
        } else
            [[unlikely]] {
            const auto pos = std::upper_bound(m_x.cbegin(), m_x.cend(), logEnergy);
            const auto index = pos != m_x.cend() ? std::distance(m_x.cbegin(), pos) : m_resolution - 1;
            const auto offset = materialIdx * m_resolution * 6 + index * 6;
            for (std::size_t i = 0; i < 3; ++i) {
                res[i] = std::pow(T { 10 }, m_coefficients[offset + 2 * i] + m_coefficients[offset + 2 * i + 1] * logEnergy);
            }
        }

        return res;
    }

    T maxAttenuationInverse(const T energy) const
    {
        const T logEnergy = std::log10(energy);
        if (logEnergy > m_linearEnergy)
            [[likely]] {
            const auto index = static_cast<std::size_t>((logEnergy - m_linearEnergy) / m_linearStep) + m_linearIndex;
            const auto offset = index + index;
            const auto res = std::pow(T { 10 }, m_maxCoefficients[offset] + m_maxCoefficients[offset + 1] * logEnergy);
            return res;
        } else
            [[unlikely]] {
            const auto pos = std::upper_bound(m_x.cbegin(), m_x.cend(), logEnergy);
            const auto index = pos != m_x.cend() ? std::distance(m_x.cbegin(), pos) : m_resolution - 1;
            const auto offset = index + index;
            const auto res = std::pow(T { 10 }, m_maxCoefficients[offset] + m_maxCoefficients[offset + 1] * logEnergy);
            return res;
        }
    }
};
}
