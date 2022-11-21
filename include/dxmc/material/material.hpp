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

namespace dxmc {

template <Floating T>
class Material2 {

public:
    Material2(std::uint64_t Z)
    {
        auto a = AtomHandler<T>(Z);
    }
    Material2(const std::map<std::size_t, T>& compositionByWeight)
    {
        auto weight = compositionByWeight;
        const auto totalWeight = std::reduce(weight.cbegin(), weight.cend(), T { 0 }, [](const auto& acc, const auto& right) -> T { return acc + right.second; });
        std::for_each(weight.begin(), weight.end(), [=](auto& w) { w.second /= totalWeight; });
    }
    static std::optional<Material2<T>> byZ(std::size_t Z)
    {
        if (Z > 0 && Z <= 100) {
            std::map<std::size_t, T> w;
            w[Z] = T { 1 };
            return byWeight(w);
        }

        return std::nullopt;
    }
    static std::optional<Material2<T>> byWeight(const std::map<std::size_t, T>& weights)
    {
        auto m = constructMaterial(weights);
        return std::nullopt;
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

protected:
    /**
     * This function parses a chemical formula string and returns a map of elements
     * Z and number density (not normalized). It's kinda messy but supports parenthesis
     * in the expression.
     */
    static std::map<std::uint64_t, T> parseCompoundStr(const std::string& str)
    {
        const std::array<std::string, 100> S { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm" };
        const std::array<std::string, 11> N { "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "." };

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

private:
    static bool constructMaterial(const std::map<std::size_t, T>& compositionByWeight)
    {
        auto weight = compositionByWeight;
        const auto totalWeight = std::reduce(weight.cbegin(), weight.cend(), T { 0 }, [](const auto& acc, const auto& right) -> T { return acc + right.second; });
        std::for_each(weight.begin(), weight.end(), [=](auto& w) { w.second /= totalWeight; });

        // photoel
        std::vector<std::pair<T, T>> photoel;
        std::vector<std::size_t> Zs;
        for (const auto& [Z, w] : weight) {
            const auto& a = AtomHandler<T>::Atom(Z);
            auto begin = photoel.insert(photoel.cend(), a.photoel.cbegin(), a.photoel.cend());
            std::for_each(begin, photoel.end(), [=](auto& p) { p.second *= w; });
            Zs.insert(Zs.cend(), a.photoel.size(), Z);
        }
        for (const auto& [Z, w] : weight) {
            const auto& a = AtomHandler<T>::Atom(Z);
            for (std::size_t i = 0; i < photoel.size(); ++i) {
                if (Zs[i] != Z) {
                    photoel[i].second += w * interpolate(a.photoel, photoel[i].first);
                }
            }
        }
        std::sort(photoel.begin(), photoel.end(), [](const auto& lh, const auto& rh) {
            if (lh.first < rh.first)
                return true;
            else if (lh.first == rh.first)
                return lh.second < rh.second;
            return false;
        });
        std::for_each(photoel.begin(), photoel.end(), [](auto& v) { v.first = std::log(v.first); v.second=std::log(v.second); });
        auto erase_from = std::unique(photoel.begin(), photoel.end(), [](const auto& lh, const auto& rh) {
            constexpr auto e = std::numeric_limits<T>::epsilon() * 10;
            if (lh.first == rh.first)
                return std::abs(lh.second - rh.second) <= e;
            else
                return std::abs(lh.first - rh.first) <= e;
        });
        photoel.erase(erase_from, photoel.end());

        // remove duplicates
        auto test = CubicLSInterpolator(photoel);
        return true;
    }
};
}