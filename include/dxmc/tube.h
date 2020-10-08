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

#include "dxmc/material.h"
#include <utility>
#include <vector>

namespace dxmc {

constexpr double TUBEMAXVOLTAGE = 150.0;
constexpr double TUBEMINVOLTAGE = 50.0;

class Tube {
public:
    Tube(double tubeVoltage = 120.0, double anodeAngleDeg = 12.0, double energyResolution = 1.0);
    Tube(const Tube& other);
    double voltage() const { return m_voltage; }
    void setVoltage(double voltage);
    static constexpr double maxVoltage() { return TUBEMAXVOLTAGE; }
    static constexpr double minVoltage() { return TUBEMINVOLTAGE; }

    double anodeAngle() const { return m_anodeAngle; }
    double anodeAngleDeg() const;
    void setAnodeAngle(double angle);
    void setAnodeAngleDeg(double angle);

    void addFiltrationMaterial(const Material& filtrationMaterial, double mm);
    std::vector<std::pair<Material, double>>& filtrationMaterials() { return m_filtrationMaterials; }
    const std::vector<std::pair<Material, double>>& filtrationMaterials() const { return m_filtrationMaterials; }

    void setAlFiltration(double mm);
    void setCuFiltration(double mm);
    void setSnFiltration(double mm);
    double AlFiltration() const;
    double CuFiltration() const;
    double SnFiltration() const;

    void clearFiltrationMaterials() { m_filtrationMaterials.clear(); }

    void setEnergyResolution(double energyResolution) { m_energyResolution = energyResolution; }
    double energyResolution() const { return m_energyResolution; }
    std::vector<double> getEnergy() const;

    std::vector<std::pair<double, double>> getSpecter(bool normalize = true) const;
    std::vector<double> getSpecter(const std::vector<double>& energies, bool normalize = true) const;
    std::vector<double> getSpecter(const std::vector<double>& energies, double anodeAngle, bool normalize = true) const;

    double mmAlHalfValueLayer();
    double mmAlHalfValueLayer() const;

protected:
    void addCharacteristicEnergy(const std::vector<double>& energy, std::vector<double>& specter) const;
    void filterSpecter(const std::vector<double>& energies, std::vector<double>& specter) const;
    void normalizeSpecter(std::vector<double>& specter) const;
    double calculatemmAlHalfValueLayer() const;

private:
    double m_voltage, m_energyResolution, m_anodeAngle;
    double m_cachedHVL = 0.0;
    bool m_hasCachedHVL = false;
    std::vector<std::pair<Material, double>> m_filtrationMaterials;
};
}