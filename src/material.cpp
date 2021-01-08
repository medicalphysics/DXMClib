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

double calculateBindingEnergy(int Z)
{
    std::vector<double> probs, energy;
    xrl_error* errorEdge = nullptr;
    int shell = 0;
    while (!errorEdge) {
        const double e = EdgeEnergy(Z, shell, &errorEdge);
        if (!errorEdge) {
            const double p = ElectronConfig(Z, shell, nullptr); // number of electrons in each shell
            probs.push_back(p);
            energy.push_back(e);
        }
        ++shell;
    }
    const double sum_probs = std::reduce(probs.cbegin(), probs.cend(), 0.0);
    const double mean_energy = std::transform_reduce(probs.cbegin(), probs.cend(), energy.cbegin(), 0.0, std::plus<>(), [=](auto p, auto e) { return e * p / sum_probs; });
    return mean_energy;
}

template <typename T>
    requires std::is_same<T, compoundData>::value || std::is_same<T, compoundDataNIST>::value double calculateMeanBindingEnergy(T* compound)
{
    std::vector<int> elements(compound->Elements, compound->Elements + compound->nElements);
    std::vector<double> massFractions(compound->massFractions, compound->massFractions + compound->nElements);
    std::vector<double> numberFractions(elements.size());
    std::transform(elements.cbegin(), elements.cend(), massFractions.cbegin(), numberFractions.begin(), [](const auto Z, const auto m) { return m / AtomicWeight(Z, nullptr); });
    auto const numberNormalization = std::accumulate(numberFractions.cbegin(), numberFractions.cend(), 0.0);
    std::transform(numberFractions.cbegin(), numberFractions.cend(), numberFractions.begin(), [=](const auto f) { return f / numberNormalization; });
    
    std::vector<double> bindingEnergy(elements.size());
    std::transform(elements.cbegin(), elements.cend(), bindingEnergy.begin(), [](const auto Z) { return calculateBindingEnergy(Z); });

    const auto meanBindingEnergy = std::transform_reduce(numberFractions.cbegin(), numberFractions.cend(), bindingEnergy.cbegin(), 0.0, std::plus<>(), std::multiplies<>());
    
    return meanBindingEnergy;
}

void Material::setByAtomicNumber(int atomicNumber)
{
    char* raw_name = AtomicNumberToSymbol(atomicNumber, nullptr);
    if (raw_name) {
        m_name = raw_name;
        xrlFree(raw_name);
        raw_name = nullptr;
        m_density = ElementDensity(atomicNumber, nullptr);
        m_meanBindingEnergy = calculateBindingEnergy(atomicNumber);
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
        //const double total_mass = std::reduce(m->massFractions, m->massFractions + m->nElements);
        //m_meanBindingEnergy = std::transform_reduce(m->massFractions, m->massFractions + m->nElements, m->Elements, 0.0, std::plus<>(), [=](double n, int z) -> double { return n * calculateMeanBindingEnergy(z) / total_mass; });
        m_meanBindingEnergy = calculateMeanBindingEnergy(m);
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

        //const double total_mass = std::reduce(m->massFractions, m->massFractions + m->nElements);
        //m_meanBindingEnergy = std::transform_reduce(m->massFractions, m->massFractions + m->nElements, m->Elements, 0.0, std::plus<>(), [=](double n, int z) -> double { return n * calculateMeanBindingEnergy(z) / total_mass; });
        m_meanBindingEnergy = calculateMeanBindingEnergy(m);
        FreeCompoundDataNIST(m);
        m = nullptr;
    }
}
}