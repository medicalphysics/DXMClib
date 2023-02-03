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

#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/material/atomhandler.hpp"
#include "dxmc/material/atomicelement.hpp"
#include "dxmc/material/nistmaterials.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <execution>

namespace dxmc {
template <Floating T>
struct Material2Shell {
    T numberOfElectronsFraction = 0;
    T numberOfElectrons = 0;
    T bindingEnergy = 0;
    T HartreeFockOrbital_0 = 0;
    T numberOfPhotonsPerInitVacancy = 0;
    T energyOfPhotonsPerInitVacancy = 0;
};
template <Floating T>
struct AttenuationValues {
    T photoelectric;
    T incoherent;
    T coherent;
    T sum() const { return photoelectric + incoherent + coherent; }
};
template <Floating T, std::size_t N = 5>
class Material2 {
public:
    static std::optional<Material2<T, N>> byZ(std::size_t Z)
    {
        auto a = AtomHandler<T>::Atom(Z);
        if (a.Z == Z) {
            std::map<std::size_t, T> w;
            w[Z] = T { 1 };
            return byWeight(w);
        }
        return std::nullopt;
    }
    static std::optional<Material2<T, N>> byWeight(const std::map<std::size_t, T>& weights)
    {
        auto m = constructMaterial(weights);
        return m;
    }

    static std::optional<Material2<T, N>> byChemicalFormula(const std::string& str)
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
        return Material2<T, N>::byWeight(weight);
    }

    static std::optional<Material2<T, N>> byNistName(const std::string& name)
    {
        const auto& w = NISTMaterials<T>::Composition(name);
        if (!w.empty())
            return Material2<T, N>::byWeight(w);

        return std::nullopt;
    }

    static std::vector<std::string> listNistCompoundNames()
    {
        return NISTMaterials<T>::listNames();
    }

    inline T effectiveZ() const { return m_effectiveZ; }

    inline T scatterFactor(const T momentumTransfer) const
    {
        const auto logmomt = std::log(momentumTransfer);
        return std::exp(CubicLSInterpolator<T>::evaluateSpline(logmomt, m_attenuationTableOffset[4], m_attenuationTableOffset[5]));
    }

    inline T formFactor(const T momentumTransfer) const
    {
        const auto logmomt = std::log(momentumTransfer);
        return std::exp(CubicLSInterpolator<T>::evaluateSpline(logmomt, m_attenuationTableOffset[3], m_attenuationTableOffset[4]));
    }

    inline T sampleSquaredMomentumTransferFromFormFactorSquared(T qsquared_max, RandomState& state) const
    {
        return m_formFactorInvSamp(qsquared_max, state);
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
    // photo, coherent, incoherent
    std::vector<AttenuationValues<T>> attenuationValues(const std::vector<T>& energy) const
    {
        std::vector<AttenuationValues<T>> att(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin, energy.cend(), att.begin(), [&](const auto e) {
            const auto logEnergy = std::log(energy);
            AttenuationValues<T> a {
                std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, m_attenuationTableOffset[0], m_attenuationTableOffset[1])),
                std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, m_attenuationTableOffset[1], m_attenuationTableOffset[2])),
                std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, m_attenuationTableOffset[2], m_attenuationTableOffset[3]))
            };
            return a; });
        return att;
    }
    inline std::vector<AttenuationValues<T>> totalAttenuationValue(const std::vector<T>& energy) const
    {
        auto att = attenuationValues(energy);
        std::vector<T> sum(att.size());
        std::transform(std::execution::par_unseq, att.cbegin(), att.cend(), sum.begin(), [](const auto& a) { return a.sum(); });
        return sum;
    }

    inline T attenuationPhotoeletric(const T energy) const
    {
        const auto logEnergy = std::log(energy);
        return std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, m_attenuationTableOffset[0], m_attenuationTableOffset[1]));
    }
    inline T attenuationPhotoelectricShell(std::uint8_t shell, const T energy) const
    {
        const std::uint8_t n_photo_shells = std::min(numberOfShells(), static_cast<std::uint8_t>(N));
        const auto start = shell < n_photo_shells ? m_attenuationTableOffset[5 + shell] : m_attenuationTableOffset[0];
        const auto stop = shell < n_photo_shells ? m_attenuationTableOffset[5 + shell + 1] : m_attenuationTableOffset[1];
        const auto logEnergy = std::log(energy);
        // prevents extrapolation to lower energies
        return logEnergy < (*start)[0] ? T { 0 } : std::exp(CubicLSInterpolator<T>::evaluateSpline(logEnergy, start, stop));
    }
    inline std::uint8_t numberOfShells() const
    {
        return m_numberOfShells;
    }
    inline const Material2Shell<T>& shell(std::size_t shell) const
    {
        return m_shells[shell];
    }
    inline const auto& shells() const { return m_shells; }

    static T momentumTransfer(T energy, T angle)
    {
        return momentumTransferMax(energy) * std::sin(angle * T { 0.5 }); // per Å
    }

    static T momentumTransferCosAngle(T energy, T cosAngle)
    {
        return momentumTransferMax(energy) * std::sqrt((1 - cosAngle) * T { 0.5 }); // per Å
    }

    static constexpr T momentumTransferMax(T energy)
    {
        constexpr T hc_si = 1.239841193E-6; // ev*m
        constexpr T m2A = 1E10; // meters to Ångstrøm
        constexpr T eV2keV = 1E-3; // eV to keV
        constexpr T hc = hc_si * m2A * eV2keV; // kev*Å
        constexpr T hc_inv = T { 1 } / hc;
        return energy * hc_inv; // per Å
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
                        auto parseP = Material2<T, N>::parseCompoundStr(strP);
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
        scatterfactor,
    };

    static CubicLSInterpolator<T> constructSplineInterpolator(const std::vector<std::vector<std::pair<T, T>>>& data, const std::vector<T>& weights, bool loglog = false, std::size_t nknots = 15)
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
            // removing zero and negative items
            auto last = std::remove_if(arr.begin(), arr.end(), [](const auto pair) -> bool {
                return pair.first <= T { 0 };
            });
            arr.erase(last, arr.end());

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
            constexpr auto e = std::numeric_limits<T>::epsilon() * 5;
            if (lh.first == rh.first)
                return std::abs(lh.second - rh.second) <= e;
            else
                return std::abs(lh.first - rh.first) <= e;
        });
        if (std::distance(erase_from, arr.end()) != 0)
            arr.erase(erase_from, arr.end());

        auto interpolator = CubicLSInterpolator(arr, nknots, true);
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
            if (type == LUTType::coherent) {
                // here we are finding the highest dip, energy wise, in coherent data and introduces an discontinuity
                auto arr = getAtomArr(a, type);
                auto max_binding_energy = MAX_ENERGY<T>();
                for (const auto& [sidx, shell] : a.shells) {
                    const auto binding_energy = shell.bindingEnergy;
                    if (max_binding_energy > MIN_ENERGY<T>() && binding_energy < max_binding_energy) {
                        auto pos = std::upper_bound(arr.begin(), arr.end(), shell.bindingEnergy, [](const auto v, const auto& el) { return v < el.first; });
                        if (pos != arr.end() && pos != arr.cbegin()) {
                            auto prev = pos - 1;
                            arr.insert(pos, std::make_pair(prev->first, pos->second));
                        }
                        max_binding_energy = binding_energy - 5;
                    }
                }
                data.push_back(arr);
            } else {
                data.push_back(getAtomArr(a, type));
            }

            weights.push_back(w);
        };

        std::size_t nknots = 10;
        if (type == LUTType::photoelectric) {
            nknots = 5;
        } else if (type == LUTType::coherent) {
            nknots = 7;
        } else if (type == LUTType::incoherent) {
            nknots = 10;
        } else if (type == LUTType::scatterfactor) {
            nknots = 20;
        } else if (type == LUTType::formfactor) {
            nknots = 20;
        }
        constexpr auto loglog = true;
        auto lut = constructSplineInterpolator(data, weights, loglog, nknots);
        return lut;
    }

    static std::optional<Material2<T, N>> constructMaterial(const std::map<std::size_t, T>& compositionByWeight)
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
            constructSplineInterpolator(weight, LUTType::scatterfactor),
        };
        Material2<T, N> m;
        m.m_effectiveZ = std::reduce(weight.cbegin(), weight.cend(), T { 0 }, [](const auto z, const auto& pair) { return z + pair.first * pair.second; });

        m.m_attenuationTable.clear();
        std::array<std::size_t, attenuation.size() + N> offset;
        for (std::size_t i = 0; i < attenuation.size(); ++i) {
            const auto& table = attenuation[i].getDataTable();
            auto begin = m.m_attenuationTable.insert(m.m_attenuationTable.end(), table.cbegin(), table.cend());
            offset[i] = std::distance(m.m_attenuationTable.begin(), begin);
        }
        createMaterialAtomicShells(m, weight, offset);

        for (std::size_t i = 0; i < offset.size(); ++i) {
            if (i < attenuation.size() + std::min(std::uint8_t { N }, m.numberOfShells()))
                m.m_attenuationTableOffset[i] = m.m_attenuationTable.begin() + offset[i];
            else
                m.m_attenuationTableOffset[i] = m.m_attenuationTable.end();
        }
        m.m_attenuationTableOffset[offset.size()] = m.m_attenuationTable.end();

        // creating lookuptable for inverse sampling of formfactor
        generateFormFactorInverseSampling(m);

        return m;
    }
    static void generateFormFactorInverseSampling(Material2<T, N>& material)
    {
        const T qmax = Material2<T, N>::momentumTransferMax(MAX_ENERGY<T>());
        constexpr T qmin = T { 0.001 };
        auto func = [&material](T qsquared) -> T {
            const auto q = std::sqrt(qsquared);
            const auto f = material.formFactor(q);
            return f * f;
        };

        material.m_formFactorInvSamp = CPDFSampling<T, 20>(qmin * qmin, qmax * qmax, func);
    }

    static void createMaterialAtomicShells(Material2<T, N>& material, const std::map<std::size_t, T>& normalizedWeight, std::array<std::size_t, 5 + N>& offset)
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

        const auto sum_weight = std::reduce(shells.cbegin(), shells.cend(), T { 0 }, [](T r, const auto& s) { return s.weight + r; });

        const auto Nshells = std::min(shells.size(), N);
        material.m_numberOfShells = static_cast<std::uint8_t>(Nshells);
        for (std::size_t i = 0; i < Nshells; ++i) {
            const auto& shell = AtomHandler<T>::Atom(shells[i].Z).shells.at(shells[i].S);
            std::vector<std::pair<T, T>> photolog(shell.photoel.size());
            std::transform(std::execution::par_unseq, shell.photoel.cbegin(), shell.photoel.cend(), photolog.begin(),
                [=](const auto& p) {
                    return std::make_pair(std::log(p.first), std::log(p.second * shells[i].weight));
                });
            CubicLSInterpolator<T> inter(photolog, 30, false);

            auto begin = inter.getDataTable().begin();
            auto end = inter.getDataTable().end();
            auto table_beg = material.m_attenuationTable.insert(material.m_attenuationTable.end(), begin, end);
            offset[i + 5] = std::distance(material.m_attenuationTable.begin(), table_beg);

            auto& materialshell = material.m_shells[i];

            materialshell.numberOfElectronsFraction = shells[i].weight * shell.numberOfElectrons / sum_weight;
            materialshell.numberOfElectrons = shell.numberOfElectrons;
            materialshell.bindingEnergy = shell.bindingEnergy;
            materialshell.HartreeFockOrbital_0 = shell.HartreeFockOrbital_0;
            materialshell.numberOfPhotonsPerInitVacancy = shell.numberOfPhotonsPerInitVacancy;
            materialshell.energyOfPhotonsPerInitVacancy = shell.energyOfPhotonsPerInitVacancy;
        }
        // Filling remainder shell
        if (shells.size() > Nshells) {
            material.m_numberOfShells++;
            auto& materialshell = material.m_shells[Nshells];
            const T mean_fac = T { 1 } / (shells.size() - Nshells);
            for (std::size_t i = Nshells; i < shells.size(); ++i) {
                const auto& shell = AtomHandler<T>::Atom(shells[i].Z).shells.at(shells[i].S);
                const auto w = shells[i].weight;

                materialshell.bindingEnergy += shell.bindingEnergy * mean_fac;
                materialshell.numberOfElectronsFraction += shell.numberOfElectrons * w / sum_weight;
                materialshell.numberOfElectrons += shell.numberOfElectrons;
                materialshell.HartreeFockOrbital_0 += shell.HartreeFockOrbital_0 * mean_fac;
                materialshell.numberOfPhotonsPerInitVacancy += shell.numberOfPhotonsPerInitVacancy * mean_fac;
                materialshell.energyOfPhotonsPerInitVacancy += shell.energyOfPhotonsPerInitVacancy * mean_fac;
            }
        }

        // normalize number og electrons fraction
        const auto sumElFraction = std::reduce(material.m_shells.cbegin(), material.m_shells.cend(), T { 0 }, [](T r, const auto& s) { return r + s.numberOfElectronsFraction; });
        std::for_each(material.m_shells.begin(), material.m_shells.end(), [sumElFraction](auto& s) { s.numberOfElectronsFraction /= sumElFraction; });
    }

private:
    T m_effectiveZ = 0;
    std::vector<std::array<T, 3>> m_attenuationTable;
    std::array<typename std::vector<std::array<T, 3>>::iterator, 5 + N + 1> m_attenuationTableOffset;
    CPDFSampling<T, 20> m_formFactorInvSamp;
    std::array<Material2Shell<T>, N + 1> m_shells;
    std::uint8_t m_numberOfShells = 0;
};
}