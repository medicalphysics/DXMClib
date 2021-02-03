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

#include <algorithm>
#include <array>
#include <string>
#include <vector>

namespace dxmc {
template <Floating T>
struct ElectronShellConfiguration {
    T bindingEnergy = 0;
    T numberElectrons = 0;
    T hartreeFockOrbital_0 = 0;
    T photoIonizationProbability = 1;
    T fluorescenceYield = 0;
    std::array<T, 3> fluorLineProbabilities = { 1, 1, 1 };
    std::array<T, 3> fluorLineEnergies = { 0, 0, 0 };
    int Z = 0;
    int shell = 0;

    template <Floating U>
    ElectronShellConfiguration<U> cast() const
    {
        ElectronShellConfiguration<U> copy;
        copy.bindingEnergy = static_cast<U>(this->bindingEnergy);
        copy.numberElectrons = static_cast<U>(this->numberElectrons);
        copy.hartreeFockOrbital_0 = static_cast<U>(this->hartreeFockOrbital_0);
        copy.fluorescenceYield = static_cast<U>(this->fluorescenceYield);
        copy.photoIonizationProbability = static_cast<T>(this->photoIonizationProbability);
        for (std::size_t i = 0; i < 3; ++i) {
            copy.fluorLineProbabilities[i] = fluorLineProbabilities[i];
            copy.fluorLineEnergies[i] = fluorLineEnergies[i];
        }
        copy.Z = Z;
        copy.shell = shell;
        return copy;
    }
};

class Material {
public:
    Material(const std::string& xraylibMaterialNameOrCompound = "", const std::string& prettyName = "", const double density = -1.0);
    Material(int atomicNumber);
    bool isValid(void) const { return m_valid && m_hasDensity; }
    const std::string& name(void) const { return m_name; }
    const std::string& prettyName(void) const;

    bool hasStandardDensity(void) const { return m_hasDensity; }
    double standardDensity(void) const { return m_density; }
    void setStandardDensity(double density); // g/cm3

    double getRayleightFormFactorSquared(const double momentumTransfer) const;
    double getComptonNormalizedScatterFactor(const double momentumTransfer) const;

    template <Floating T>
    T getRayleightFormFactorSquared(const T momentumTransfer) const
    {
        return static_cast<T>(getRayleightFormFactorSquared(static_cast<double>(momentumTransfer)));
    }
    template <Floating T>
    T getComptonNormalizedScatterFactor(const T momentumTransfer) const
    {
        return static_cast<T>(getComptonNormalizedScatterFactor(static_cast<double>(momentumTransfer)));
    }

    double getPhotoelectricAttenuation(double energy) const;
    double getRayleightAttenuation(double energy) const;
    double getComptonAttenuation(double energy) const;
    double getTotalAttenuation(double energy) const;
    double getMassEnergyAbsorbtion(double energy) const;

    std::vector<double> getBindingEnergies(const double minValue = 1) const;

    template <Floating T>
    std::vector<T> getBindingEnergies(const T minValue = 1) const
    {
        const double dminValue = minValue;
        const auto vals = getBindingEnergies(dminValue);
        std::vector<T> res(vals.size());
        std::transform(vals.cbegin(), vals.cend(), res.begin(), [](const double e) -> T { return static_cast<T>(e); });
        return res;
    }

    std::array<ElectronShellConfiguration<double>, 12> getElectronConfiguration() const;

    template <Floating T>
    std::array<ElectronShellConfiguration<T>, 12> getElectronConfiguration() const
    {
        const auto conf_d = getElectronConfiguration();
        std::array<ElectronShellConfiguration<T>, 12> conf;
        for (std::size_t i = 0; i < conf.size(); ++i) {
            conf[i] = conf_d[i].cast<T>();
        }
        return conf;
    }

    static double getAtomicWeight(int Z);
    static std::string getAtomicNumberToSymbol(int Z);

    static std::vector<std::string> getNISTCompoundNames(void);
    static int getAtomicNumberFromSymbol(const std::string& symbol);
    static std::string getSymbolFromAtomicNumber(int Z);

    static double getTotalAttenuation(int atomicNumber, double energy);

protected:
    void setByCompoundName(const std::string& name);
    void setByMaterialName(const std::string& name);
    void setByAtomicNumber(int index);

private:
    std::string m_name;
    std::string m_prettyName;
    std::vector<int> m_elements;
    std::vector<double> m_elementNumberFraction;
    double m_density = -1.0;
    bool m_valid = false;
    bool m_hasDensity = false;
};
}