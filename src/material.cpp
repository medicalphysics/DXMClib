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

#include "dxmc/material.h"
#include "dxmc/interpolation.h"
#include "dxmc/vectormath.h"
#include "xraylib.h"

#include <algorithm>
#include <execution>
#include <numeric>

namespace dxmc {

double Material::getPhotoelectricAttenuation(double energy) const
{
    return CS_Photo_CP(m_name.c_str(), energy, nullptr);
}

double Material::getRayleightAttenuation(double energy) const
{
    return CS_Rayl_CP(m_name.c_str(), energy, nullptr);
}

double Material::getComptonAttenuation(double energy) const
{
    return CS_Compt_CP(m_name.c_str(), energy, nullptr);
}

double Material::getTotalAttenuation(double energy) const
{
    return CS_Total_CP(m_name.c_str(), energy, nullptr);
}

double Material::getTotalAttenuation(int atomicNumber, double energy)
{
    return CS_Total(atomicNumber, energy, nullptr);
}

double Material::getMassEnergyAbsorbtion(double energy) const
{
    return CS_Energy_CP(m_name.c_str(), energy, nullptr);
}

double Material::getAtomicWeight(int Z)
{
    return AtomicWeight(Z, nullptr);
}

std::string Material::getAtomicNumberToSymbol(int Z)
{
    return std::string(AtomicNumberToSymbol(Z, nullptr));
}

Material::Material(int atomicNumber)
{
    setByAtomicNumber(atomicNumber);
}

const std::string& Material::prettyName(void) const
{
    if (m_prettyName.empty())
        return m_name;
    return m_prettyName;
}

Material::Material(const std::string& xraylibMaterialNameOrCompound, const std::string& prettyName, const double density)
{
    setByMaterialName(xraylibMaterialNameOrCompound);
    if (!m_valid) {
        setByCompoundName(xraylibMaterialNameOrCompound);
    }
    if (prettyName.size() > 0) {
        m_prettyName = prettyName;
    } else {
        m_prettyName = m_name;
    }
    if (density > 0) {
        setStandardDensity(density);
    }
}

void Material::setStandardDensity(double density)
{
    if (density > 0.0) {
        m_density = density;
        m_hasDensity = true;
    }
}

std::vector<std::string> Material::getNISTCompoundNames(void)
{
    int n_strings;
    auto charArray = GetCompoundDataNISTList(&n_strings, nullptr);
    auto materialNames = std::vector<std::string>();
    for (int i = 0; i < n_strings; ++i) {
        if (charArray[i]) {
            std::string s(charArray[i]);
            materialNames.push_back(s);
        }
    }
    return materialNames;
}

int Material::getAtomicNumberFromSymbol(const std::string& symbol)
{
    return SymbolToAtomicNumber(symbol.c_str(), nullptr);
}
std::string Material::getSymbolFromAtomicNumber(int Z)
{
    char* chars = AtomicNumberToSymbol(Z, nullptr);
    std::string st(chars);
    xrlFree(chars);
    return st;
}

std::vector<double> Material::getRayleightFormFactorSquared(const std::vector<double>& momentumTransfer) const
{
    std::vector<double> fraction;
    std::vector<int> elements;

    struct compoundData* m = CompoundParser(m_name.c_str(), nullptr);
    if (m) {
        fraction.resize(m->nElements);
        elements.resize(m->nElements);
        for (int i = 0; i < m->nElements; i++) {
            elements[i] = m->Elements[i];
            fraction[i] = m->nAtoms[i] / (m->nAtomsAll);
        }
        FreeCompoundData(m);
        m = nullptr;
    }
    struct compoundDataNIST* n = GetCompoundDataNISTByName(m_name.c_str(), nullptr);
    if (n) {
        fraction.resize(n->nElements);
        elements.resize(n->nElements);
        for (int i = 0; i < n->nElements; i++) {
            elements[i] = n->Elements[i];
            fraction[i] = n->massFractions[i] / AtomicWeight(elements[i], nullptr);
        }
        const double weight = std::accumulate(fraction.cbegin(), fraction.cend(), 0.0);
        std::transform(fraction.cbegin(), fraction.cend(), fraction.begin(), [=](double w) { return w / weight; });
        FreeCompoundDataNIST(n);
        n = nullptr;
    }

    std::vector<double> formFactor(momentumTransfer.size(), 0.0);
    for (std::size_t i = 0; i < formFactor.size(); ++i) {
        for (std::size_t j = 0; j < fraction.size(); ++j) {
            formFactor[i] += fraction[j] * FF_Rayl(elements[j], momentumTransfer[i], nullptr);
        }
        formFactor[i] = formFactor[i] * formFactor[i];
    }
    return formFactor;
}

template <typename T>
    requires std::is_same<T, compoundData>::value || std::is_same<T, compoundDataNIST>::value std::array<ElectronShellConfiguration<double>, 12> electronConfiguration(T* compound)
{
    std::vector<ElectronShellConfiguration<double>> configs;
    std::vector<int> elements(compound->Elements, compound->Elements + compound->nElements);

    std::vector<double> massFractions(compound->massFractions, compound->massFractions + compound->nElements);
    std::vector<double> numberFractions(elements.size());
    std::transform(elements.cbegin(), elements.cend(), massFractions.cbegin(), numberFractions.begin(), [](const auto Z, const auto m) { return m / AtomicWeight(Z, nullptr); });
    auto const numberNormalization = std::accumulate(numberFractions.cbegin(), numberFractions.cend(), 0.0);
    std::transform(numberFractions.cbegin(), numberFractions.cend(), numberFractions.begin(), [=](const auto f) { return f / numberNormalization; });

    for (int i = 0; i < elements.size(); ++i) {
        const int Z = elements[i];
        const double numberFraction = numberFractions[i];
        // we group obitals together, i.e the K shell is one orbital, L1,L2,L3 is grouped into a <L> shell, all the way up to including N shells
        xrl_error* errorEdge = nullptr;
        int shell = 0;

        while (!errorEdge) {
            const double bindingEnergy = EdgeEnergy(Z, shell, &errorEdge); // binding energy
            if (!errorEdge) {
                const double totalElectrons = numberFraction * ElectronConfig(Z, shell, nullptr); // number of electrons in each shell
                const double CProfile = ComptonProfile_Partial(Z, shell, 0.0, &errorEdge); // Hartree Fock orbital for electron momentum =0

                std::array<double, 3> fluroProbabilities = { 0, 0, 0 };
                std::array<double, 3> fluroEnergy = { 0, 0, 0 };
                std::array<int, 3> lines { 0, 0, 0 };

                int line_start = 0;
                int line_stop = 0;
                double yield = 0;
                if (shell == 0) {
                    line_start = -1;
                    line_stop = -29;
                    yield = FluorYield(Z, K_SHELL, nullptr);
                } else if (shell == 1) {
                    line_start = -30;
                    line_stop = -58;
                    yield = FluorYield(Z, L1_SHELL, nullptr) + FluorYield(Z, L2_SHELL, nullptr) * CosKronTransProb(Z, FL12_TRANS, nullptr);
                    yield += (CosKronTransProb(Z, FL13_TRANS, nullptr) + CosKronTransProb(Z, FL12_TRANS, nullptr) * CosKronTransProb(Z, FL23_TRANS, nullptr)) * FluorYield(Z, L3_SHELL, nullptr);
                } else if (shell == 2) {
                    line_start = -59;
                    line_stop = -85;
                    yield = FluorYield(Z, L2_SHELL, nullptr) + FluorYield(Z, L3_SHELL, nullptr) * CosKronTransProb(Z, FL23_TRANS, nullptr);
                } else if (shell == 3) {
                    line_start = -86;
                    line_stop = -113;
                    yield = FluorYield(Z, L3_SHELL, nullptr);
                }
                for (int line = line_start; line >= line_stop; --line) {
                    const auto amin = vectormath::argmin3<int, double>(fluroProbabilities.data());
                    xrl_error* errorLine = nullptr;
                    const double prob = RadRate(Z, line, &errorLine);
                    if (prob > fluroProbabilities[amin] && !errorLine) {
                        fluroProbabilities[amin] = prob;
                        fluroEnergy[amin] = LineEnergy(Z, line, nullptr);
                    }
                }
                const auto fluroLineprobs = std::reduce(fluroProbabilities.cbegin(), fluroProbabilities.cend(), 0.0);
                if (fluroLineprobs > 0) {
                    std::transform(
                        fluroProbabilities.cbegin(), fluroProbabilities.cend(), fluroProbabilities.begin(), [=](auto p) { return p / fluroLineprobs; });
                }

                ElectronShellConfiguration<double> new_config(
                    bindingEnergy,
                    totalElectrons,
                    CProfile,
                    0.0,
                    yield,
                    fluroProbabilities,
                    fluroEnergy,
                    Z,
                    shell);
                configs.push_back(new_config);
            }
            ++shell;
        }
    }

    std::sort(configs.begin(), configs.end(), [](const auto& lh, const auto& rh) { return lh.bindingEnergy > rh.bindingEnergy; });

    std::array<ElectronShellConfiguration<double>, 12> configs_a;
    for (std::size_t i = 0; i < std::min(configs.size(), configs_a.size()); ++i) {
        configs_a[i] = configs[i];
    }

    const auto electrons_sum = std::transform_reduce(configs_a.cbegin(), configs_a.cend(), 0.0, std::plus<>(), [](const auto c) { return c.numberElectrons; });
    for (auto& c : configs_a) {
        c.numberElectrons /= electrons_sum;
    }

    //calculating shell probabilities;
    std::vector<double> energy(400);
    std::iota(energy.begin(), energy.end(), 0.5);

    for (std::size_t i = 0; i < configs_a.size(); ++i) {
        if (configs_a[i].Z > 0) {
            const auto c_prob = configs_a[i].numberElectrons * std::transform_reduce(std::execution::par_unseq, energy.cbegin(), energy.cend(), 0.0, std::plus<>(), [=](const auto e) { return CSb_Photo_Partial(configs_a[i].Z, configs_a[i].shell, e, nullptr); });
            double tot_prob = c_prob;

            for (std::size_t j = i + 1; j < configs_a.size(); ++j) {
                tot_prob += configs_a[j].numberElectrons * std::transform_reduce(std::execution::par_unseq, energy.cbegin(), energy.cend(), 0.0, std::plus<>(), [=](const auto e) { return CSb_Photo_Partial(configs_a[j].Z, configs_a[j].shell, e, nullptr); });
            }
            if (tot_prob > 0) {
                configs_a[i].photoIonizationProbability = c_prob / tot_prob;
            }
        }
    }
    return configs_a;
}

std::array<ElectronShellConfiguration<double>, 12> Material::getElectronConfiguration() const
{
    std::array<ElectronShellConfiguration<double>, 12> config;

    struct compoundData* m = CompoundParser(m_name.c_str(), nullptr);
    if (m) {
        config = electronConfiguration(m);
        FreeCompoundData(m);
        m = nullptr;
        return config;
    }
    struct compoundDataNIST* n = GetCompoundDataNISTByName(m_name.c_str(), nullptr);
    if (n) {
        config = electronConfiguration(n);
        FreeCompoundDataNIST(n);
        n = nullptr;
    }
    return config;
}
std::vector<double> Material::getComptonNormalizedScatterFactor(const std::vector<double>& momentumTransfer) const
{
    std::vector<double> fraction;
    std::vector<int> elements;

    struct compoundData* m = CompoundParser(m_name.c_str(), nullptr);
    if (m) {
        fraction.resize(m->nElements);
        elements.resize(m->nElements);
        for (int i = 0; i < m->nElements; i++) {
            elements[i] = m->Elements[i];
            fraction[i] = m->nAtoms[i] / (m->nAtomsAll);
        }
        FreeCompoundData(m);
        m = nullptr;
    }
    struct compoundDataNIST* n = GetCompoundDataNISTByName(m_name.c_str(), nullptr);
    if (n) {
        fraction.resize(n->nElements);
        elements.resize(n->nElements);
        for (int i = 0; i < n->nElements; i++) {
            elements[i] = n->Elements[i];
            fraction[i] = n->massFractions[i] / AtomicWeight(elements[i], nullptr);
        }
        const double weight = std::accumulate(fraction.cbegin(), fraction.cend(), 0.0);
        std::transform(fraction.cbegin(), fraction.cend(), fraction.begin(), [=](double w) { return w / weight; });
        FreeCompoundDataNIST(n);
        n = nullptr;
    }

    std::vector<double> formFactor(momentumTransfer.size(), 0.0);
    for (std::size_t i = 0; i < formFactor.size(); ++i) {
        for (std::size_t j = 0; j < fraction.size(); ++j) {
            formFactor[i] += fraction[j] * SF_Compt(elements[j], momentumTransfer[i], nullptr) / elements[j];
        }
    }
    return formFactor;
}

template <typename T>
    requires std::is_same<T, compoundData>::value || std::is_same<T, compoundDataNIST>::value std::vector<double> calculateBindingEnergies(T* compound, const double minEnergy)
{
    std::vector<int> elements(compound->Elements, compound->Elements + compound->nElements);
    std::vector<double> edges;

    for (const int Z : elements) {
        int shell = 0;
        xrl_error* errorEdge = nullptr;
        while (!errorEdge) {
            const double edge = EdgeEnergy(Z, shell, &errorEdge);
            if (!errorEdge && edge > minEnergy) {
                edges.push_back(edge);
            }
            ++shell;
        }
    }
    std::sort(edges.begin(), edges.end(), std::greater<double>());
    return edges;
}

std::vector<double> Material::getBindingEnergies(const double minValue) const
{
    std::vector<double> edges;
    struct compoundData* m = CompoundParser(m_name.c_str(), nullptr);
    if (m) {
        edges = calculateBindingEnergies(m, minValue);
        FreeCompoundData(m);
        m = nullptr;
    }
    struct compoundDataNIST* n = GetCompoundDataNISTByName(m_name.c_str(), nullptr);
    if (n) {
        edges = calculateBindingEnergies(n, minValue);
        FreeCompoundDataNIST(n);
        n = nullptr;
    }
    return edges;
}

void Material::setByAtomicNumber(int atomicNumber)
{
    char* raw_name = AtomicNumberToSymbol(atomicNumber, nullptr);
    if (raw_name) {
        m_name = raw_name;
        xrlFree(raw_name);
        raw_name = nullptr;
        m_density = ElementDensity(atomicNumber, nullptr);        
        m_hasDensity = true;
        m_valid = true;
    }
}

void Material::setByCompoundName(const std::string& name)
{
    compoundData* m = CompoundParser(name.c_str(), nullptr);
    if (m) {
        m_name = name;
        m_valid = true;
        m_hasDensity = false;        
        FreeCompoundData(m);
        m = nullptr;
    }
}

void Material::setByMaterialName(const std::string& name)
{
    struct compoundDataNIST* m = GetCompoundDataNISTByName(name.c_str(), nullptr);
    if (m) {
        m_name = m->name;
        m_valid = true;
        m_hasDensity = true;
        m_density = m->density;

        
        
        FreeCompoundDataNIST(m);
        m = nullptr;
    }
}
}