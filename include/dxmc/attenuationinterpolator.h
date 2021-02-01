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
#include <vector>

namespace dxmc {
template <Floating T>
class AttenuationLutInterpolator {
private:
    std::vector<double> m_x;
    std::vector<double> m_coefficients;
    double m_shell_offset = 0.0001;
    std::size_t m_resolution = 40;
    T maxBindingEnergy = 0;

protected:
    static std::vector<double> gaussSplineElimination(const std::vector<double>& h, const std::vector<double>& D)
    {
        const std::size_t m = h.size() - 1; // rows
        const std::size_t n = m + 1; // columns

        auto idx = [=](int i, int j) -> int { return j + i * m; }; //(row, column)
        //construct matrix
        std::vector<double> A(n * m, 0);
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                const auto index = idx(i, j);
                if (i == j)
                    A[index] = 2 * (h[i] + h[i + 1]);
                else if (j - i == 1)
                    A[index] = h[i];
                else if (j - i == -1)
                    A[index] = h[i];
                else if (j == n)
                    A[index] = D[j + 1];
            }
        }

        for (int i = 0; i < m - 1; i++) {
            //Partial Pivoting
            for (int k = i + 1; k < m; k++) {
                //If diagonal element(absolute vallue) is smaller than any of the terms below it
                if (std::abs(A[idx(i, i)]) < std::abs(A[idx(k, i)])) {
                    //Swap the rows
                    for (int j = 0; j < n; j++) {
                        const auto temp = A[idx(i, j)];
                        A[idx(i, j)] = A[idx(k, j)];
                        A[idx(k, j)] = temp;
                    }
                }
            }
            //Begin Gauss Elimination
            for (int k = i + 1; k < m; k++) {
                const auto term = A[idx(k, i)] / A[idx(i, i)];
                for (int j = 0; j < n; j++) {
                    A[idx(k, j)] = A[idx(k, j)] - term * A[idx(i, j)];
                }
            }
        }
        std::vector<double> sigma(h.size() + 1, 0);

        //Begin Back-substitution
        for (int i = m - 1; i >= 0; i--) {
            sigma[i + 1] = A[idx(i, n - 1)];
            for (int j = i + 1; j < n - 1; j++) {
                sigma[i + 1] = sigma[i + 1] - A[idx(i, j)] * sigma[j + 1];
            }
            sigma[i + 1] = sigma[i + 1] / A[idx(i, i)];
        }
        return sigma;
    }

public:
    AttenuationLutInterpolator(const std::vector<Material>& materials, const T minEnergy, const T maxEnergy, std::size_t resolution = 40)
    {
        // get Shell edge energies
        m_x.clear();
        for (const auto& mat : materials) {
            const auto minEnergyD = static_cast<double>(minEnergy);
            const auto bindingEnergies = mat.getBindingEnergies(minEnergyD);
            for (const auto e : bindingEnergies) {
                m_x.push_back(std::log10(e - m_shell_offset));
                m_x.push_back(std::log10(e));
            }
        }

        const double logMinEnergy = std::log10(minEnergy);
        const double logMaxEnergy = std::log10(maxEnergy);
        for (std::size_t i = 0; i <= m_resolution; ++i) {
            const T e = logMinEnergy + (i * (logMaxEnergy - logMinEnergy)) / (m_resolution - 1);
            m_x.push_back(e);
        }
        std::sort(m_x.begin(), m_x.end());

        m_coefficients.resize(4 * materials.size() * (m_x.size() - 1) * 3);

        for (std::size_t matIdx = 0; matIdx < materials.size(); ++matIdx) {
            const auto offset = 4 * 3 * (m_x.size() - 1) * matIdx;
            for (std::size_t typeIdx = 0; typeIdx < 3; ++typeIdx) {
                std::vector<double> y(m_x.size());
                if (typeIdx == 0)
                    std::transform(m_x.cbegin(), m_x.cend(), y.begin(), [&](const double le) { return materials[matIdx].getPhotoelectricAttenuation(std::pow(10, le)); });
                else if (typeIdx == 1) {
                    std::transform(m_x.cbegin(), m_x.cend(), y.begin(), [&](const double le) { return materials[matIdx].getComptonAttenuation(std::pow(10, le)); });
                } else {
                    std::transform(m_x.cbegin(), m_x.cend(), y.begin(), [&](const double le) { return materials[matIdx].getRayleightAttenuation(std::pow(10, le)); });
                }
                std::vector<double> h(m_x.size() - 1);
                std::vector<double> rho(m_x.size() - 1);
                std::vector<double> delta(m_x.size() - 1);
                for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
                    h[i] = m_x[i + 1] - m_x[i];
                    delta[i] = (y[i + 1] - y[i]) / h[i];
                }

                const double sigma0 = 1;
                const double sigmaN = 1;
                std::vector<double> D(m_x.size() - 1);
                for (std::size_t i = 1; i < D.size(); ++i) {
                    D[i] = 6 * (delta[i] - delta[i - 1]);
                }
                D[0] = 0;
                D[1] += -h[1] * sigma0;
                D[m_x.size() - 2] += -h[m_x.size() - 2] * sigmaN;

                auto sigma1N = gaussSplineElimination(h, D);
                sigma1N[0] = sigma0;
                sigma1N[m_x.size() - 1] = sigmaN;

                for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
                    m_coefficients[offset + i * 4 * typeIdx + 0] = (sigma1N[i] * m_x[i + 1] * m_x[i + 1] * m_x[i + 1] - sigma1N[i + 1] * m_x[i] * m_x[i] * m_x[i] + 6 * (y[i] * m_x[i + 1] - y[i + 1] * m_x[i])) / (h[i] * 6);
                    m_coefficients[offset + i * 4 * typeIdx + 0] += h[i] * (sigma1N[i + 1] * m_x[i] - sigma1N[i] * m_x[i + 1]) / 6;

                    m_coefficients[offset + i * 4 * typeIdx + 1] = (sigma1N[i + 1] * m_x[i] * m_x[i] - sigma1N[i] * m_x[i + 1] * m_x[i + 1] + 2 * (y[i + 1] - y[i])) / (h[i] * 2) + h[i] * (sigma1N[i] - sigma1N[i + 1]) / 6;
                    m_coefficients[offset + i * 4 * typeIdx + 2] = (sigma1N[i] * m_x[i + 1] - sigma1N[i + 1] * m_x[i]) / (2 * h[i]);
                    m_coefficients[offset + i * 4 * typeIdx + 3] = (sigma1N[i + 1] - sigma1N[i]) / (6 * h[i]);
                }
            }
        }
    }

    std::array<T, 3> getAttenuation(const std::size_t materialIdx, const T energy) const
    {
        const auto logEnergy = std::log10(energy);
        auto upper_iter = std::upper_bound(m_x.cbegin(), m_x.cend(), logEnergy);
        if (upper_iter == m_x.cend() || upper_iter == m_x.cbegin())
            return 0;
        const std::size_t i = std::distance(m_x.cbegin(), upper_iter) - 1;
        const auto offset = 4 * 3 * materialIdx * m_x.size() + 4 * 3 * i;
        std::array<T, 3> res;
        const auto logEnergy2 = logEnergy * logEnergy;
        for (std::size_t t = 0; t < 3; ++t) {
            res[t] = m_coefficients[offset] + m_coefficients[offset + 1] * logEnergy + m_coefficients[offset + 2] * logEnergy2 + m_coefficients[offset + 3] * logEnergy * logEnergy2;
        }
        return res;
    }
};
}
