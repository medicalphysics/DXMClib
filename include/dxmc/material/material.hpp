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

#include "dxmc/floating.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/atomhandler.hpp"
#include "dxmc/material/atomicelement.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <execution>

#include <iostream>

namespace dxmc {
template <Floating T>
class Material2Shell {
public:
    Material2Shell(const AtomicShell<T>& shell) { }

private:
    std::size_t m_shellIdx = 0;
    T m_nElectrons = 0;
    T m_bindingEnergy = 0;
    T m_hartreeFockProfile_0 = 0;
    T m_shellNumberOfPhotonsPerVacancy = 0;
    T m_shellEnergyOfPhotonsPerVacancy = 0;
    std::array<std::array<T, 3>, 30> m_photoProb;
};
template <Floating T>
struct AttenuationValues {
    T photoelectric;
    T incoherent;
    T coherent;
    T sum() { return photoelectric + incoherent + coherent; }
};
template <Floating T, int N = 5>
class Material2 {

public:
    static std::optional<Material2<T>> byZ(std::size_t Z)
    {
        auto a = AtomHandler<T>::Atom(Z);
        if (a.Z == Z) {
            std::map<std::size_t, T> w;
            w[Z] = T { 1 };
            return byWeight(w);
        }
        return std::nullopt;
    }
    static std::optional<Material2<T>> byWeight(const std::map<std::size_t, T>& weights)
    {
        auto m = constructMaterial(weights);
        return m;
    }

    static std::optional<Material2<T>> byChemicalFormula(const std::string& str)
    {
        auto numberDensCompound = parseCompoundStr(str);
        if (numberDensCompound.size() == 0)
            return std::nullopt;
        std::map<std::size_t, T> weight;
        for (const auto [Z, numDens] : numberDensCompound) {
            if (Z == 0 || Z > 99)
                return std::nullopt;
            const auto& atom = AtomHandler<T>::Atom(Z);
            weight[Z] = numDens * atom.atomicWeight;
        }
        return Material2<T>::byWeight(weight);
    }

    inline T formFactor(const T energy, const T angle) const
    {
        const auto mt = momentumTransfer(energy, angle);
        return CubicLSInterpolator<T>::evaluateSpline(mt, m_attenuationTableOffset[3], m_attenuationTableOffset[4]);
    }

    inline T scatterFactor(const T energy, const T angle) const
    {
        const auto mt = momentumTransfer(energy, angle);
        return CubicLSInterpolator<T>::evaluateSpline(mt, m_attenuationTableOffset[4], m_attenuationTableOffset[5]);
    }

    // photo, coherent, incoherent
    inline AttenuationValues<T> attenuationValues(const T energy) const
    {
        const auto logEnergy = std::log(energy);
        AttenuationValues<T> att {
            std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, m_attenuationTableOffset[0], m_attenuationTableOffset[1])),
            std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, m_attenuationTableOffset[1], m_attenuationTableOffset[2])),
            std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, m_attenuationTableOffset[2], m_attenuationTableOffset[3]))
        };
        return att;
    }
    inline T attenuationPhotoeletric(const T energy) const
    {
        const auto logEnergy = std::log(energy);
        return std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, m_attenuationTableOffset[0], m_attenuationTableOffset[1]));
    }
    static T momentumTransfer(T energy, T angle)
    {
        constexpr double hc_si = 1.239841193E-6; // ev*m
        constexpr double m2A = 1E10; // meters to �ngstr�m
        constexpr double eV2keV = 1E-3; // eV to keV
        constexpr double hc = hc_si * m2A * eV2keV; // kev*�
        constexpr double hc_inv = 1.0 / hc;
        return energy * std::sin(angle * 0.5) * hc_inv; // per �
    }

protected:
    Material2()
    {
    }

    /**
     * This function parses a chemical formula string and returns a map of elements
     * Z and number density (not normalized). It's kinda messy but supports parenthesis
     * in the expression.
     */
    static std::map<std::uint64_t, T> parseCompoundStr(const std::string& str)
    {
        const std::array<std::string, 100> S { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm" };
        //        const std::array<std::string, 11> N { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "." };

        std::map<std::uint64_t, T> parsed;

        auto begin = str.begin();
        std::vector<std::pair<std::string, std::string>> vals;
        while (begin != str.end()) {
            const char l = *begin;
            if (std::isupper(l)) {
                auto v = std::pair<std::string, std::string>();
                v.first = std::string(&l, 1);
                vals.push_back(v);
            } else if (std::islower(l)) {
                if (vals.size() > 0)
                    vals.back().first.append(&l, 1);
            } else if ((l == '.') || std::isdigit(l)) {
                if (vals.size() > 0)
                    vals.back().second.append(&l, 1);
            } else if (l == '(') {
                auto startP = begin + 1;
                auto endP = startP;
                while (endP != str.end()) {
                    const char l1 = *endP;
                    if (l1 == ')') {
                        std::string strP(startP, endP);
                        auto parseP = Material2<T>::parseCompoundStr(strP);
                        auto endPdig = endP + 1;
                        T w = 1;
                        std::string str_w;
                        while (endPdig != str.end()) {
                            const char l2 = *endPdig;
                            if (std::isdigit(l2) || l2 == '.') {
                                str_w.append(&l2, 1);
                                begin = endPdig;
                                endPdig++;
                            } else {
                                endPdig = str.end();
                            }
                        }
                        try {
                            w = static_cast<T>(std::stod(str_w));
                        } catch (const std::invalid_argument& e) {
                        }
                        for (const auto& [el, frac] : parseP) {
                            if (!parsed.contains(el)) {
                                parsed[el] = 0;
                            }
                            parsed[el] += (frac * w);
                        }
                        endP = str.end();
                    } else {
                        endP++;
                    }
                }
            }
            begin++;
        }
        for (const auto& [el, num] : vals) {
            auto pos = std::find(S.begin(), S.end(), el);
            if (pos != S.end()) {
                std::uint64_t Z = std::distance(S.begin(), pos) + 1;
                double w = 1.0;
                try {
                    w = std::stod(num);
                } catch (const std::invalid_argument& e) {
                }
                if (!parsed.contains(Z)) {
                    parsed[Z] = 0;
                }
                parsed[Z] += static_cast<T>(w);
            }
        }
        return parsed;
    }

    enum class LUTType {
        photoelectric,
        coherent,
        incoherent,
        formfactor,
        scatterfactor
    };

    static CubicLSInterpolator<T> constructSplineInterpolator(const std::vector<std::vector<std::pair<T, T>>>& data, const std::vector<T>& weights, bool loglog = false)
    {
        const auto Nvec = std::reduce(data.cbegin(), data.cend(), std::size_t { 0 }, [](const auto lh, const auto& rh) -> std::size_t { return rh.size() + lh; });

        std::vector<T> w(data.size());
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (data.size() != weights.size())
                w[i] = T { 1 };
            else
                w[i] = weights[i];
        }
        const auto sum_weights = std::reduce(w.cbegin(), w.cend(), T { 0 });
        std::transform(w.cbegin(), w.cend(), w.begin(), [sum_weights](const auto& ww) { return ww / sum_weights; });

        std::vector<std::pair<T, T>> arr(Nvec);
        std::vector<std::size_t> idx(Nvec);
        std::size_t start = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            std::copy(std::execution::par_unseq, data[i].cbegin(), data[i].cend(), arr.begin() + start);
            std::fill(idx.begin() + start, idx.begin() + start + data[i].size(), i);
            start += data[i].size();
        }
        // applying weights for first array in data
        std::transform(std::execution::par_unseq, arr.cbegin(), arr.cend(), idx.cbegin(), arr.begin(), [&](const auto& pair, const auto i) {
            return std::make_pair(pair.first, pair.second * w[i]);
        });

        for (std::size_t i = 0; i < data.size(); ++i) {
            std::transform(std::execution::par, arr.begin(), arr.end(), idx.cbegin(), arr.begin(), [&](const auto& pair, const auto index) {
                if (index == i)
                    return pair;
                else
                    return std::make_pair(pair.first, pair.second + w[i] * interpolate(data[i], pair.first));
            });
        }

        if (loglog) {
            std::for_each(std::execution::par_unseq, arr.begin(), arr.end(), [](auto& v) {
                v.first = std::log(v.first);
                v.second = std::log(v.second);
            });
        }
        // sorting for generating interpolation lut
        std::sort(arr.begin(), arr.end(), [](const auto& lh, const auto& rh) {
            if (lh.first < rh.first)
                return true;
            else if (lh.first == rh.first)
                return lh.second < rh.second;
            return false;
        });

        // erasing duplicate items
        auto erase_from = std::unique(arr.begin(), arr.end(), [](const auto& lh, const auto& rh) {
            constexpr auto e = std::numeric_limits<T>::epsilon() * 10;
            if (lh.first == rh.first)
                return std::abs(lh.second - rh.second) <= e;
            else
                return std::abs(lh.first - rh.first) <= e;
        });
        if (std::distance(erase_from, arr.end()) != 0)
            arr.erase(erase_from, arr.end());

        auto interpolator = CubicLSInterpolator(arr);
        return interpolator;
    }

    // constructs a least squares spline interpolator from attenuation data from a compound
    // This function is a bit convoluted since we want to preserve discontiuities in the interpolation
    // But is basicly an weighted average of attenuation coefficients for each element.
    static CubicLSInterpolator<T> constructSplineInterpolator(const std::map<std::size_t, T>& normalizedWeight, LUTType type)
    {
        auto getAtomArr = [&](const AtomicElement<T>& atom, LUTType type = LUTType::photoelectric) -> const std::vector<std::pair<T, T>>& {
            if (type == LUTType::photoelectric)
                return atom.photoel;
            else if (type == LUTType::incoherent)
                return atom.incoherent;
            else if (type == LUTType::coherent)
                return atom.coherent;
            else if (type == LUTType::formfactor)
                return atom.formFactor;
            else if (type == LUTType::scatterfactor)
                return atom.incoherentSF;
            return atom.photoel;
        };

        std::vector<std::vector<std::pair<T, T>>> data;
        std::vector<T> weights;
        for (const auto& [Z, w] : normalizedWeight) {
            const auto& a = AtomHandler<T>::Atom(Z);
            data.push_back(getAtomArr(a, type));
            weights.push_back(w);
        }
        const bool loglog = type == LUTType::photoelectric || type == LUTType::incoherent || type == LUTType::coherent;
        auto lut = constructSplineInterpolator(data, weights, loglog);
        return lut;
    }

    static std::optional<Material2<T>> constructMaterial(const std::map<std::size_t, T>& compositionByWeight)
    {
        for (const auto& [Z, w] : compositionByWeight) {
            const auto& a = AtomHandler<T>::Atom(Z);
            if (a.Z != Z)
                return std::nullopt;
        }

        auto weight = compositionByWeight;
        const auto totalWeight = std::reduce(weight.cbegin(), weight.cend(), T { 0 }, [](const auto& acc, const auto& right) -> T { return acc + right.second; });
        std::for_each(weight.begin(), weight.end(), [=](auto& w) { w.second /= totalWeight; });

        std::array<CubicLSInterpolator<T>, 5> attenuation = {
            constructSplineInterpolator(weight, LUTType::photoelectric),
            constructSplineInterpolator(weight, LUTType::incoherent),
            constructSplineInterpolator(weight, LUTType::coherent),
            constructSplineInterpolator(weight, LUTType::formfactor),
            constructSplineInterpolator(weight, LUTType::scatterfactor)
        };
        Material2<T> m;
        m.m_attenuationTable.clear();
        std::array<std::size_t, attenuation.size()> offset;
        for (std::size_t i = 0; i < attenuation.size(); ++i) {
            const auto& table = attenuation[i].getDataTable();
            auto begin = m.m_attenuationTable.insert(m.m_attenuationTable.end(), table.cbegin(), table.cend());
            offset[i] = std::distance(m.m_attenuationTable.begin(), begin);
        }
        for (std::size_t i = 0; i < attenuation.size(); ++i) {
            m.m_attenuationTableOffset[i] = m.m_attenuationTable.begin() + offset[i];
        }
        m.m_attenuationTableOffset[attenuation.size()] = m.m_attenuationTable.end();

        createMaterialAtomicShells(m, weight);
        return m;
    }
    static void createMaterialAtomicShells(Material2<T>& material, const std::map<std::size_t, T>& normalizedWeight)
    {
        struct Shell {
            std::uint64_t Z = 0;
            std::uint64_t S = 0;
            T weight = 0;
            T bindingEnergy = 0;
        };

        std::vector<Shell> shells;
        for (const auto& [Z, w] : normalizedWeight) {
            const auto& atom = AtomHandler<T>::Atom(Z);
            for (const auto& [S, shell] : atom.shells) {
                shells.emplace_back(Z, S, w, shell.bindingEnergy);
            }
        }
        std::sort(shells.begin(), shells.end(), [](const auto& lh, const auto& rh) {
            return lh.bindingEnergy > rh.bindingEnergy;
        });

        std::vector<T> energies;
        for (int i = 0; i < std::min(N, static_cast<int>(shells.size())); ++i) {
            const auto& shell = AtomHandler<T>::Atom(shells[i].Z).shells.at(shells[i].S);
            std::transform(shell.photoel.cbegin(), shell.photoel.cend(), std::back_inserter(energies), [](const auto& v) { return std::log(v.first); });
        }
        std::sort(std::execution::unseq, energies.begin(), energies.end());
        auto erase_from = std::unique(energies.begin(), energies.end(), [](const auto& lh, const auto& rh) {
            constexpr auto e = std::numeric_limits<T>::epsilon() * 10;
            return std::abs(lh - rh) <= e;
        });
        if (std::distance(erase_from, energies.end()) != 0)
            energies.erase(erase_from, energies.end());
        std::array<std::vector<T>, N> photo;
        for (int i = 0; i < std::min(N, static_cast<int>(shells.size())); ++i) {
            photo[i].resize(energies.size());
            const auto& shell = AtomHandler<T>::Atom(shells[i].Z).shells.at(shells[i].S);
            const auto w = normalizedWeight.at(shells[i].Z);
            std::transform(energies.cbegin(), energies.cend(), photo[i].begin(), [&](const auto& v) {
                const auto e = std::exp(v);
                return w * interpolate<T, false, true>(shell.photoel, e) / material.attenuationPhotoeletric(e);
            });
        }

        auto test = photo[0];
        std::fill(test.begin(), test.end(), T { 0 });
        for (int i = 0; i < test.size(); ++i)
            for (const auto& v : photo)
                test[i] += v[i];

        

        auto slett = true;
    }

private:
    std::vector<std::array<T, 3>> m_attenuationTable;
    std::array<typename std::vector<std::array<T, 3>>::iterator, 6> m_attenuationTableOffset;
    // std::array<Material2Shell<T>, N> m_shells;
};
}