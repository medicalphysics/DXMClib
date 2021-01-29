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
    std::array<T, 3> fluorLineProbabilities;
    std::array<T, 3> fluorLineEnergies;

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

    std::vector<double> getRayleightFormFactorSquared(const std::vector<double>& momentumTransfer) const;
    std::vector<double> getComptonNormalizedScatterFactor(const std::vector<double>& momentumTransfer) const;

    template <Floating T>
    std::vector<T> getRayleightFormFactorSquared(const std::vector<T>& momentumTransfer) const
    {
        std::vector<double> in(momentumTransfer.cbegin(), momentumTransfer.cend());
        auto vec = getRayleightFormFactorSquared(in);
        std::vector<T> out(in.size());
        std::transform(vec.cbegin(), vec.cend(), out.begin(), [](double e) -> T { return static_cast<T>(e); });
        return out;
    }
    template <Floating T>
    std::vector<T> getComptonNormalizedScatterFactor(const std::vector<T>& momentumTransfer) const
    {
        std::vector<double> in(momentumTransfer.cbegin(), momentumTransfer.cend());
        auto vec = getComptonNormalizedScatterFactor(in);
        std::vector<T> out(in.size());
        std::transform(vec.cbegin(), vec.cend(), out.begin(), [](double e) -> T { return static_cast<T>(e); });
        return out;
    }

    double getPhotoelectricAttenuation(double energy) const;
    double getRayleightAttenuation(double energy) const;
    double getComptonAttenuation(double energy) const;
    double getTotalAttenuation(double energy) const;
    double getMassEnergyAbsorbtion(double energy) const;

    double getMeanBindingEnergy() const
    {
        return m_meanBindingEnergy;
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
    double m_density = -1.0;
    double m_meanBindingEnergy = 0; // binding energy in keV
    bool m_valid = false;
    bool m_hasDensity = false;
};
}