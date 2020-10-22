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

#pragma once // include guard

#include "dxmc/beamfilters.h"
#include "dxmc/constants.h"
#include "dxmc/dxmcrandom.h"
#include "dxmc/exposure.h"
#include "dxmc/floating.h"
#include "dxmc/progressbar.h"
#include "dxmc/transport.h"
#include "dxmc/tube.h"
#include "dxmc/world.h"

#include <array>
#include <memory>
#include <vector>

namespace dxmc {

template <Floating T>
class Source {
public:
    enum class Type { None,
        CTSpiral,
        CTAxial,
        DX,
        CTDual,
        Pencil,
        Isotropic };

private:
    void normalizeDirectionCosines(void)
    {
        vectormath::normalize(&m_directionCosines[0]);
        vectormath::normalize(&m_directionCosines[3]);
    }

protected:
    std::array<T, 3> m_position = { 0, 0, 0 };
    std::array<T, 6> m_directionCosines = { 1, 0, 0, 0, 1, 0 };
    std::uint64_t m_historiesPerExposure = 1;
    Type m_type = Type::None;

public:
    Source()
        : m_type(Type::None) {};
    virtual Exposure<T> getExposure(std::uint64_t i) const = 0;
    virtual T maxPhotonEnergyProduced() const { return Tube<T>::maxVoltage(); }

    void setPosition(const std::array<T, 3>& position) { m_position = position; }
    void setPosition(T x, T y, T z)
    {
        m_position[0] = x;
        m_position[1] = y;
        m_position[2] = z;
    };
    std::array<T, 3>& position(void) { return m_position; }
    const std::array<T, 3>& position(void) const { return m_position; }

    void setDirectionCosines(const std::array<T, 6>& cosines)
    {
        m_directionCosines = cosines;
        normalizeDirectionCosines();
    }
    const std::array<T, 6>& directionCosines(void) const { return m_directionCosines; }
    std::array<T, 6>& directionCosines(void) { return m_directionCosines; }

    void setHistoriesPerExposure(std::uint64_t histories) { m_historiesPerExposure = histories; }
    std::uint64_t historiesPerExposure(void) const { return m_historiesPerExposure; }

    virtual std::uint64_t totalExposures(void) const = 0;

    Source::Type type() const { return m_type; }

    virtual T getCalibrationValue(ProgressBar<T>* progress = nullptr) const = 0;

    virtual bool isValid(void) const = 0;
    virtual bool validate(void) = 0;
    virtual void updateFromWorld(const World<T>& world) { }
};

template <Floating T>
class PencilSource final : public Source<T> {
public:
    PencilSource()
    {
        m_type = Type::Pencil;
    }
    virtual ~PencilSource() = default;
    Exposure<T> getExposure(std::uint64_t i) const override
    {
        constexpr std::array<T, 2> collimationAngles { 0, 0 };
        Exposure<T> exposure(m_position,
            m_directionCosines,
            collimationAngles,
            m_historiesPerExposure);

        exposure.setMonoenergeticPhotonEnergy(m_photonEnergy);
        return exposure;
    }

    void setPhotonEnergy(T energy)
    {
        if (energy > ELECTRON_REST_MASS<T>())
            m_photonEnergy = ELECTRON_REST_MASS<T>();
        else if (energy < 1)
            m_photonEnergy = 1;
        else
            m_photonEnergy = energy;
    }
    T photonEnergy() const { return m_photonEnergy; }
    T maxPhotonEnergyProduced() const override { return m_photonEnergy; }
    std::uint64_t totalExposures(void) const override { return m_totalExposures; }
    void setTotalExposures(std::uint64_t exposures)
    {
        if (exposures > 0)
            m_totalExposures = exposures;
    }
    void setAirDose(T Gycm2)
    {
        if (Gycm2 > 0.0)
            m_airDose = Gycm2;
    }
    T airDose(void) const { return m_airDose; } // Gycm2

    T getCalibrationValue(ProgressBar<T>* = nullptr) const override
    {
        Material airMaterial("Air, Dry (near sea level)");
        const T nHistories = totalExposures() * historiesPerExposure();
        const T mea = static_cast<T>(airMaterial.getMassEnergyAbsorbtion(m_photonEnergy));
        const T calcOutput = nHistories * m_photonEnergy * mea * KEV_TO_MJ<T>();
        const T factor = m_airDose / calcOutput;
        return factor;
    }
    bool isValid(void) const override { return true; }
    bool validate(void) override { return true; }

protected:
    T m_photonEnergy = 100.0;
    T m_airDose = 1.0;
    std::uint64_t m_totalExposures = 10;
};

template <Floating T>
class IsotropicSource final : public Source<T> {
public:
    IsotropicSource()
        : Source<T>()
    {
        m_type = Type::Isotropic;
        // initializing a specterdistribution
        const std::vector<T> energies = { 60.0 };
        const std::vector<T> weights = { 1.0 };
        m_maxPhotonEnergy = 60.0;
        m_specterDistribution = SpecterDistribution<T>(weights, energies);
        //m_specterDistribution = std::make_unique<SpecterDistribution>(weights, energies);
    }
    virtual ~IsotropicSource() = default;
    Exposure<T> getExposure(std::uint64_t exposureNumber) const override
    {
        constexpr T weight { 1 };
        Exposure<T> exposure(m_position,
            m_directionCosines,
            m_collimationAngles,
            m_historiesPerExposure,
            weight,
            &m_specterDistribution);

        return exposure;
    }
    T maxPhotonEnergyProduced() const override
    {
        return m_maxPhotonEnergy;
    };
    void setTotalExposures(std::uint64_t nExposures)
    {
        m_totalExposures = nExposures;
    };
    std::uint64_t totalExposures() const override
    {
        return m_totalExposures;
    };

    T getCalibrationValue(ProgressBar<T>* progress = nullptr) const override
    {
        return T { 1.0 };
    };

    bool isValid() const override { return true; }
    bool validate() override { return true; }

    void setSpecter(const std::vector<T>& weights, const std::vector<T>& energies)
    {
        const auto maxEnergyElement = std::max_element(energies.cbegin(), energies.cend());
        m_maxPhotonEnergy = *maxEnergyElement;
        m_specterDistribution = SpecterDistribution<T>(weights, energies);
    }
    void setCollimationAngles(T xRad, T yRad)
    {
        if (xRad < -PI_VAL<T>())
            m_collimationAngles[0] = -PI_VAL<T>();
        else if (xRad > PI_VAL<T>())
            m_collimationAngles[0] = PI_VAL<T>();
        else
            m_collimationAngles[0] = xRad;

        constexpr T PI_H = PI_VAL<T>() * T { 0.5 };

        if (yRad < -PI_H)
            m_collimationAngles[1] = -PI_H;
        else if (yRad > PI_VAL<T>())
            m_collimationAngles[1] = PI_H;
        else
            m_collimationAngles[1] = yRad;
    }
    std::pair<T, T> collimationAngles() const
    {
        return std::pair<T, T>(m_collimationAngles[0], m_collimationAngles[1]);
    }

protected:
private:
    std::uint64_t m_totalExposures = 1;
    std::array<T, 2> m_collimationAngles = { 0, 0 };
    SpecterDistribution<T> m_specterDistribution;
    T m_maxPhotonEnergy = 1.0;
};

template <Floating T>
class DXSource final : public Source<T> {
public:
    DXSource()
        : Source<T>()
    {
        m_type = Type::DX;
        m_sdd = 1000.0;
        m_fieldSize[0] = 100.0;
        m_fieldSize[1] = 100.0;
        updateFieldSize(m_fieldSize);
        m_tube.setAlFiltration(2.0);
        setDirectionCosines(zeroDirectionCosines());
    }
    virtual ~DXSource() = default;
    Exposure<T> getExposure(std::uint64_t i) const override
    {
        constexpr T weight { 1 };
        Exposure<T> exposure(tubePosition(),
            m_directionCosines,
            m_collimationAngles,
            m_historiesPerExposure,
            weight,
            m_specterDistribution.get(),
            m_heelFilter.get());
        return exposure;
    }

    Tube<T>& tube(void)
    {
        m_specterValid = false;
        return m_tube;
    }

    const Tube<T>& tube(void) const { return m_tube; }

    T maxPhotonEnergyProduced() const override { return m_tube.voltage(); }

    std::uint64_t totalExposures(void) const override
    {
        return m_totalExposures;
    }
    void setTotalExposures(std::uint64_t exposures)
    {
        m_totalExposures = std::max(exposures, 1);
    }

    void setCollimationAngles(const std::array<T, 2>& radians)
    {
        std::array<double, 2> ang = { std::abs(angles[0]), std::abs(angles[1]) };
        updateCollimationAngles(ang);
    }
    const std::array<T, 2>& collimationAngles(void) const
    {
        return m_collimationAngles;
    }

    void setCollimationAnglesDeg(const std::array<T, 2>& degrees)
    {
        std::array<T, 2> ang = { std::abs(angles[0]), std::abs(angles[1]) };
        updateCollimationAngles(ang);
    }
    const std::array<T, 2> collimationAnglesDeg(void) const
    {
        std::array<T, 2> a { m_collimationAngles[0] * RAD_TO_DEG<T>(),
            m_collimationAngles[1] * RAD_TO_DEG<T>() };
        return a;
    }

    void setFieldSize(const std::array<T, 2>& mm)
    {
        std::array<T, 2> amm = { std::abs(mm[0]), std::abs(mm[1]) };
        updateFieldSize(amm);
    }

    const std::array<T, 2>& fieldSize(void) const { return m_fieldSize; }

    void setSourceDetectorDistance(T mm)
    {
        m_sdd = std::abs(mm);
        updateFieldSize(m_fieldSize);
    }
    T sourceDetectorDistance() const { return m_sdd; }

    void setSourceAngles(T primaryAngle, T secondaryAngle)
    {
        constexpr T ANGLE_ERRF = 1E-6;
        if (secondaryAngle > PI_VAL<T>() * 0.5 - ANGLE_ERRF)
            secondaryAngle = PI_VAL<T>() * 0.5 - ANGLE_ERRF;
        if (secondaryAngle < -PI_VAL<T>() * 0.5 + ANGLE_ERRF)
            secondaryAngle = -PI_VAL<T>() * 0.5 + ANGLE_ERRF;
        while (primaryAngle > PI_VAL<T>())
            primaryAngle -= PI_VAL<T>();
        while (primaryAngle < -PI_VAL<T>())
            primaryAngle += PI_VAL<T>();

        std::array<T, 6> cos = zeroDirectionCosines();
        std::array<T, 3> z = { .0, .0, 1.0 };
        vectormath::rotate(cos.data(), z.data(), primaryAngle);
        vectormath::rotate(&cos[3], z.data(), primaryAngle);
        std::array<T, 3> x = { 1.0, .0, .0 };
        vectormath::rotate(cos.data(), x.data(), -secondaryAngle);
        vectormath::rotate(&cos[3], x.data(), -secondaryAngle);
        //handling tube rotation
        std::array<T, 3> beam_dir;
        vectormath::cross(cos.data(), beam_dir.data());
        vectormath::rotate(cos.data(), beam_dir.data(), m_tubeRotationAngle);
        vectormath::rotate(&cos[3], beam_dir.data(), m_tubeRotationAngle);

        setDirectionCosines(cos);
    }
    void setSourceAngles(const std::array<T, 2>& angles)
    {
        setSourceAngles(angles[0], angles[1]);
    }
    std::array<T, 2> sourceAngles() const
    {
        constexpr T ANGLE_ERRF = 1E-6;
        auto cos = directionCosines();
        std::array<T, 3> beam_direction;
        vectormath::cross(cos.data(), beam_direction.data());
        vectormath::rotate(cos.data(), beam_direction.data(), -m_tubeRotationAngle);
        vectormath::rotate(&cos[3], beam_direction.data(), -m_tubeRotationAngle);
        vectormath::cross(cos.data(), beam_direction.data());

        //handling floating point
        T primAng, secAng;
        if (beam_direction[0] < -1.0 + ANGLE_ERRF) {
            primAng = PI_VAL<T>() / 2.0;
            secAng = 0.0;
        } else if (beam_direction[0] > 1.0 - ANGLE_ERRF) {
            primAng = -PI_VAL<T>() / 2.0;
            secAng = 0.0;
        } else {
            primAng = std::asin(-beam_direction[0]);
            secAng = -std::atan(beam_direction[2] / beam_direction[1]);
        }

        std::array<T, 2> angles = { primAng, secAng };

        T y[3] = { 0, 1, 0 };
        const auto dot_dir = vectormath::dot(y, beam_direction.data());
        if (dot_dir < 0.0) {
            if (angles[0] < 0) {
                angles[0] = -PI_VAL<T>() - angles[0];
            } else {
                angles[0] = PI_VAL<T>() - angles[0];
            }
        }

        return angles;
    }
    void setSourceAnglesDeg(T primaryAngle, T secondaryAngle)
    {
        setSourceAngles(primaryAngle * DEG_TO_RAD<T>(), secondaryAngle * DEG_TO_RAD<T>());
    }

    void setSourceAnglesDeg(const std::array<T, 2>& angles)
    {
        setSourceAnglesDeg(angles[0], angles[1]);
    }

    std::array<T, 2> sourceAnglesDeg() const
    {
        auto angles = sourceAngles();
        angles[0] *= RAD_TO_DEG<T>();
        angles[1] *= RAD_TO_DEG<T>();
        return angles;
    }

    void setTubeRotation(T angle)
    {
        const T diff_angle = angle - m_tubeRotationAngle;
        std::array<T, 3> beam_direction;
        auto cos = directionCosines();
        vectormath::cross(cos.data(), beam_direction.data());
        vectormath::rotate(cos.data(), beam_direction.data(), diff_angle);
        vectormath::rotate(&cos[3], beam_direction.data(), diff_angle);
        setDirectionCosines(cos);
        m_tubeRotationAngle = angle;
    }
    T tubeRotation() const { return m_tubeRotationAngle; }
    void setTubeRotationDeg(T angle) { setTubeRotation(angle * DEG_TO_RAD<T>()); }
    T tubeRotationDeg() const { return tubeRotation() * RAD_TO_DEG<T>(); }

    void setDap(T Gycm2)
    {
        if (Gycm2 > 0.0)
            m_dap = Gycm2;
    }
    T dap(void) const { return m_dap; } // Gycm2

    T getCalibrationValue(ProgressBar<T>* = nullptr) const override
    {
        auto specter = tube().getSpecter();
        std::vector<T> massAbsorb(specter.size(), T { 0 });
        Material airMaterial("Air, Dry (near sea level)");
        std::transform(specter.begin(), specter.end(), massAbsorb.begin(), [&](const auto& el) -> T { return airMaterial.getMassEnergyAbsorbtion(el.first); });

        const T nHist = totalExposures() * historiesPerExposure();

        T calcOutput = 0.0; // Air KERMA [keV/g]
        for (std::size_t i = 0; i < specter.size(); ++i) {
            auto& [keV, weight] = specter[i];
            calcOutput += keV * weight * nHist * massAbsorb[i];
        }
        calcOutput *= KEV_TO_MJ<T>() * 1000.0; // Air KERMA [mJ / kg] = [mGy]

        // (m_dap * 1000.0): converting from Gycm2 to mGycm2
        const T output = (m_dap * T { 1000.0 }) / (m_fieldSize[0] * m_fieldSize[1] * T { 0.01 }); //mm->cm
        const T factor = output / calcOutput; // mGy/mGy
        return factor;
    }

    const std::array<T, 3> tubePosition(void) const
    {
        std::array<T, 3> beamDirection;
        vectormath::cross(m_directionCosines.data(), beamDirection.data());
        std::array<T, 3> pos;
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] = m_position[i] - beamDirection[i] * m_sdd;
        }
        return pos;
    }

    bool isValid(void) const override { return m_specterValid; }
    bool validate(void) override
    {
        updateSpecterDistribution();
        return m_specterValid;
    }

    void setModelHeelEffect(bool on) { m_modelHeelEffect = on; };
    bool modelHeelEffect() const { return m_modelHeelEffect; };

protected:
    void updateFieldSize(const std::array<T, 2>& fieldSize)
    {
        for (std::size_t i = 0; i < 2; ++i) {
            m_fieldSize[i] = fieldSize[i];
            m_collimationAngles[i] = std::atan(m_fieldSize[i] * T { 0.5 } / m_sdd) * T { 2 };
        }
        m_specterValid = false;
    }

    void updateCollimationAngles(const std::array<T, 2>& fieldSize);
    void updateSpecterDistribution()
    {
        if (!m_specterValid) {
            auto energies = m_tube.getEnergy();
            auto n_obs = m_tube.getSpecter(energies);
            m_specterDistribution = std::make_unique<SpecterDistribution<T>>(n_obs, energies);
            if (m_modelHeelEffect)
                m_heelFilter = std::make_unique<HeelFilter<T>>(m_tube, m_collimationAngles[1]);
            else
                m_heelFilter = nullptr;
            m_specterValid = true;
        }
    }
    std::array<T, 6> zeroDirectionCosines() const
    {
        // changing this will break setSourceAngles
        std::array<T, 6> cos = { -1.0, .0, .0, .0, .0, 1.0 };
        return cos;
    }

private:
    T m_sdd = 1000.0;
    T m_dap = 1.0; // Gycm2
    std::array<T, 2> m_fieldSize;
    std::array<T, 2> m_collimationAngles;
    std::uint64_t m_totalExposures = 1000;
    Tube<T> m_tube;
    T m_tubeRotationAngle = 0.0;
    bool m_specterValid = false;
    std::shared_ptr<SpecterDistribution<T>> m_specterDistribution = nullptr;
    std::shared_ptr<HeelFilter<T>> m_heelFilter = nullptr;
    bool m_modelHeelEffect = true;
};

template <Floating T>
class CTSource : public Source<T> {
public:
    CTSource()
        : Source<T>()
    {
        m_type = Type::None;
        m_sdd = 1190.0;
        m_collimation = 38.4;
        m_fov = 500.0;
        m_startAngle = 0.0;
        m_exposureAngleStep = DEG_TO_RAD<T>();
        m_scanLenght = 100.0;
        auto& t = tube();
        t.setAlFiltration(7.0);
        const std::array<T, 6> ct_cosines { -1, 0, 0, 0, 0, 1 };
        setDirectionCosines(ct_cosines);
    }
    virtual ~CTSource() = default;
    virtual Exposure<T> getExposure(std::uint64_t i) const override = 0;

    Tube<T>& tube(void)
    {
        m_specterValid = false;
        return m_tube;
    };
    const Tube<T>& tube(void) const { return m_tube; }

    virtual T maxPhotonEnergyProduced() const override { return m_tube.voltage(); }

    void setBowTieFilter(std::shared_ptr<BowTieFilter<T>> filter) { m_bowTieFilter = filter; }
    std::shared_ptr<BowTieFilter<T>> bowTieFilter(void) { return m_bowTieFilter; }
    const std::shared_ptr<BowTieFilter<T>> bowTieFilter() const { return m_bowTieFilter; }

    void setSourceDetectorDistance(T sdd)
    {
        m_sdd = std::abs(sdd);
        m_specterValid = false;
    }
    T sourceDetectorDistance(void) const { return m_sdd; }

    void setCollimation(T mmCollimation)
    {
        m_collimation = std::abs(collimation);
        m_specterValid = false;
    }
    T collimation(void) const { return m_collimation; }

    void setFieldOfView(T fov)
    {
        m_fov = std::abs(fov);
    }
    T fieldOfView(void) const
    {
        return m_fov;
    }

    void setGantryTiltAngle(T angle)
    {
        if (angle < -PI_VAL<T>())
            m_gantryTiltAngle = -PI_VAL<T>();
        else if (angle > PI_VAL<T>())
            m_gantryTiltAngle = PI_VAL<T>();
        else
            m_gantryTiltAngle = angle;
    }
    T gantryTiltAngle() const
    {
        return m_gantryTiltAngle;
    }

    void setGantryTiltAngleDeg(T angle)
    {
        setGantryTiltAngle(angle * DEG_TO_RAD<T>());
    }

    T gantryTiltAngleDeg() const
    {
        return m_gantryTiltAngle * RAD_TO_DEG<T>();
    }

    void setStartAngle(T angle)
    {
        m_startAngle = angle;
    }

    T startAngle(void) const
    {
        return m_startAngle;
    }
    T startAngleDeg(void) const
    {
        return RAD_TO_DEG<T>() * m_startAngle;
    }
    void setStartAngleDeg(T angle)
    {
        m_startAngle = DEG_TO_RAD<T>() * angle;
    }

    void setExposureAngleStep(T angleStep)
    {
        const auto absAngle = std::abs(angleStep);
        if (absAngle < PI_VAL<T>())
            m_exposureAngleStep = absAngle;
        else if (absAngle > T { 0.1 } * DEG_TO_RAD<T>())
            m_exposureAngleStep = absAngle;
    }

    T exposureAngleStep(void) const
    {
        return m_exposureAngleStep;
    }

    void setExposureAngleStepDeg(T angleStep)
    {
        setExposureAngleStep(angleStep * DEG_TO_RAD<T>());
    }

    T exposureAngleStepDeg(void) const
    {
        return m_exposureAngleStep * RAD_TO_DEG<T>();
    }

    virtual void setScanLenght(T scanLenght)
    {
        m_scanLenght = std::abs(scanLenght);
    }
    T scanLenght(void) const
    {
        return m_scanLenght;
    }

    void setCtdiVol(T ctdivol)
    {
        if (ctdivol > 0.0)
            m_ctdivol = ctdivol;
    }
    T ctdiVol(void) const { return m_ctdivol; }

    bool useXCareFilter() const { return m_useXCareFilter; }
    void setUseXCareFilter(bool use) { m_useXCareFilter = use; }
    XCareFilter<T>& xcareFilter() { return m_xcareFilter; }
    const XCareFilter<T>& xcareFilter() const { return m_xcareFilter; }

    virtual std::uint64_t totalExposures(void) const override = 0;

    void setCtdiPhantomDiameter(std::uint64_t mm)
    {
        if (mm > 160)
            m_ctdiPhantomDiameter = mm;
        else
            m_ctdiPhantomDiameter = 160;
    }
    std::uint64_t ctdiPhantomDiameter(void) const { return m_ctdiPhantomDiameter; }

    virtual T getCalibrationValue(ProgressBar<T>* = nullptr) const = 0;

    virtual std::uint64_t exposuresPerRotatition() const
    {
        constexpr T pi_2 = T { 2 } * PI_VAL<T>();
        return static_cast<std::size_t>(pi_2 / m_exposureAngleStep);
    }

    bool isValid(void) const override { return m_specterValid; };
    virtual bool validate(void) override
    {
        updateSpecterDistribution();
        return m_specterValid;
    }

    void setAecFilter(std::shared_ptr<AECFilter<T>> filter) { m_aecFilter = filter; }
    std::shared_ptr<AECFilter<T>> aecFilter(void) { return m_aecFilter; }

    virtual void updateFromWorld(const World<T>& world) override
    {
        if (m_aecFilter)
            m_aecFilter->updateFromWorld(world);
    }

    void setModelHeelEffect(bool on) { m_modelHeelEffect = on; };
    bool modelHeelEffect() const { return m_modelHeelEffect; };

protected:
    template <typename U>
    requires std::is_base_of<CTSource<T>, U>::value T ctCalibration(
        const U& source, ProgressBar<T>* progressBar = nullptr) const
    {

        U sourceCopy(source);

        CTDIPhantom<T> world(sourceCopy.ctdiPhantomDiameter());

        sourceCopy.updateFromWorld(world);

        T meanWeight = 0;
        for (std::size_t i = 0; i < sourceCopy.totalExposures(); ++i) {
            Exposure<T> dummy;
            sourceCopy.getExposure(dummy, i);
            meanWeight += dummy.beamIntensityWeight();
        }
        meanWeight /= static_cast<double>(sourceCopy.totalExposures());

        const std::array<T, 6> cosines({ -1, 0, 0, 0, 0, 1 });
        sourceCopy.setDirectionCosines(cosines);

        sourceCopy.setUseXCareFilter(false); // we need to disable organ aec for ctdi statistics, this should be ok
        std::size_t statCounter = world.ctdiMinHistories() / 2 / (sourceCopy.exposuresPerRotatition() * sourceCopy.historiesPerExposure());
        if (statCounter < 1)
            statCounter = 1;

        const auto histories = sourceCopy.historiesPerExposure();
        sourceCopy.setHistoriesPerExposure(histories * statCounter); // ensuring enough histories for ctdi measurement
        sourceCopy.validate();

        Transport<T> transport;
        auto result = transport(world, &sourceCopy, progressBar);

        typedef CTDIPhantom<T>::HolePosition holePosition;
        std::array<CTDIPhantom<T>::HolePosition, 5> position = { holePosition::Center, holePosition::West, holePosition::East, holePosition::South, holePosition::North };

        std::array<T, 5> measureDose;
        measureDose.fill(0.0);
        for (std::size_t i = 0; i < 5; ++i) {
            auto holeIndices = world.holeIndices(position[i]);
            for (auto idx : holeIndices)
                measureDose[i] += result.dose[idx];
            measureDose[i] /= static_cast<T>(holeIndices.size());
        }

        std::array<std::uint32_t, 5> measureTally;
        measureTally.fill(0);
        for (std::size_t i = 0; i < 5; ++i) {
            auto holeIndices = world.holeIndices(position[i]);
            for (auto idx : holeIndices)
                measureTally[i] += result.nEvents[idx];
        }

        const T ctdiPher = (measureDose[1] + measureDose[2] + measureDose[3] + measureDose[4]) / T { 4.0 };
        const T ctdiCent = measureDose[0];
        const T ctdiw = (ctdiCent + T { 2 } * ctdiPher) / T { 3 } / static_cast<T>(statCounter);
        const T factor = sourceCopy.ctdiVol() / ctdiw / meanWeight;
        return factor;
    }
    virtual void updateSpecterDistribution()
    {
        if (!m_specterValid) {
            auto energies = m_tube.getEnergy();
            auto n_obs = m_tube.getSpecter(energies);
            m_specterDistribution = std::make_shared<SpecterDistribution<T>>(n_obs, energies);
            const T heel_span_angle = std::atan(m_collimation * T { 0.5 } / m_sdd) * T { 2.0 };
            if (m_modelHeelEffect)
                m_heelFilter = std::make_shared<HeelFilter<T>>(m_tube, heel_span_angle);
            else
                m_heelFilter = nullptr;
            m_specterValid = true;
        }
    }

    T m_sdd;
    T m_collimation;
    T m_fov;
    T m_startAngle;
    T m_exposureAngleStep;
    T m_scanLenght;
    T m_ctdivol = 1.0;
    T m_gantryTiltAngle = 0.0;
    std::shared_ptr<AECFilter<T>> m_aecFilter = nullptr;
    std::uint64_t m_ctdiPhantomDiameter = 320;
    std::shared_ptr<BowTieFilter<T>> m_bowTieFilter = nullptr;
    XCareFilter<T> m_xcareFilter;
    bool m_useXCareFilter = false;
    bool m_specterValid = false;
    Tube<T> m_tube;
    std::shared_ptr<SpecterDistribution<T>> m_specterDistribution = nullptr;
    std::shared_ptr<HeelFilter<T>> m_heelFilter = nullptr;
    bool m_modelHeelEffect = true;
};

template <Floating T>
class CTSpiralSource final : public CTSource<T> {
public:
    CTSpiralSource()
        : CTSource<T>()
    {
        m_type = Type::CTSpiral;
        m_pitch = 1.0;
    }
    Exposure<T> getExposure(std::uint64_t exposureIndex) const override
    {
        std::array<T, 3> pos = { 0, -m_sdd / 2.0, 0 };

        const auto angle = m_startAngle + m_exposureAngleStep * exposureIndex;

        auto directionCosines = m_directionCosines;
        T* rotationAxis = &directionCosines[3];
        T* otherAxis = &directionCosines[0];
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(rotationAxis, tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(otherAxis, tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis, angle);
        pos[2] += (exposureIndex * m_exposureAngleStep) * m_collimation * m_pitch / PI_2 + tiltCorrection[2];

        vectormath::rotate(otherAxis, rotationAxis, angle);
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] += m_position[i];
        }
        /*
        std::array<T, 3> rotationAxis, otherAxis;
        for (std::size_t i = 0; i < 3; ++i) {
            rotationAxis[i] = m_directionCosines[i + 3];
            otherAxis[i] = m_directionCosines[i];
        }
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(rotationAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(otherAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis.data(), angle);
        
        //adding transverse step
        pos[2] += (exposureIndex * m_exposureAngleStep) * m_collimation * m_pitch / PI_2 + tiltCorrection[2];

        vectormath::rotate(otherAxis.data(), rotationAxis.data(), angle);
        for (std::size_t i = 0; i < 3; ++i)
            pos[i] += m_position[i];
        */

        const std::array<2, T> collimationAngles = {
            std::atan(m_fov / m_sdd) * T { 2.0 },
            std::atan(m_collimation / m_sdd) * T { 2.0 }
        };

        T weight { 1.0 };
        if (m_aecFilter)
            weight *= m_aecFilter->sampleIntensityWeight(pos);
        if (m_useXCareFilter)
            weight *= m_xcareFilter.sampleIntensityWeight(angle);

        Exposure<T> exposure(pos,
            directionCosines,
            collimationAngles,
            m_historiesPerExposure,
            weight,
            m_specterDistribution.get(),
            m_heelFilter.get(),
            m_bowTieFilter.get());

        return exposure;
        /*
        exposure.setPosition(pos);
        exposure.setDirectionCosines(otherAxis, rotationAxis);
        exposure.setCollimationAngles(std::atan(m_fov / m_sdd) * 2.0, std::atan(m_collimation / m_sdd) * 2.0);
        exposure.setBeamFilter(m_bowTieFilter.get());
        exposure.setSpecterDistribution(m_specterDistribution.get());
        exposure.setHeelFilter(m_heelFilter.get());
        exposure.setNumberOfHistories(m_historiesPerExposure);
        T weight { 1.0 };
        if (m_aecFilter)
            weight *= m_aecFilter->sampleIntensityWeight(pos);
        if (m_useXCareFilter)
            weight *= m_xcareFilter.sampleIntensityWeight(angle);
        exposure.setBeamIntensityWeight(weight);
        return true;
        */
    }

    void setPitch(T pitch)
    {
        m_pitch = std::max(0.01, pitch);
    }

    T pitch(void) const { return m_pitch; }
    void setScanLenght(T scanLenght) override
    {
        m_scanLenght = std::max(std::abs(scanLenght), m_collimation * m_pitch * 0.5);
    }
    std::uint64_t totalExposures(void) const override
    {
        return static_cast<std::uint64_t>(m_scanLenght * PI_2 / (m_collimation * m_pitch * m_exposureAngleStep));
    }

    T getCalibrationValue(ProgressBar<T>* = nullptr) const override
    {
        const auto& me = *this;
        return ctCalibration(me, progressBar) * m_pitch;
    }

private:
    T m_pitch;
};

template <Floating T>
class CTAxialSource final : public CTSource<T> {
public:
    CTAxialSource()
        : CTSource<T>()
    {
        m_type = Type::CTAxial;
        m_step = m_collimation;
        m_scanLenght = m_step;
    }
    bool getExposure(Exposure<T>& exposure, std::uint64_t exposureIndex) const override
    {
        //calculating position
        std::array<T, 3> pos = { 0, -m_sdd / T { 2.0 }, 0 };

        constexpr T PI_2 = 2 * PI_VAL<T>();

        const std::uint64_t anglesPerRotation = static_cast<std::uint64_t>(PI_2 / m_exposureAngleStep);
        const std::uint64_t rotationNumber = exposureIndex / anglesPerRotation;

        const auto angle = m_startAngle + m_exposureAngleStep * (exposureIndex - (rotationNumber * anglesPerRotation));

        auto directionCosines = m_directionCosines;
        T* rotationAxis = &directionCosines[3];
        T* otherAxis = &directionCosines[0];
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(rotationAxis, tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(otherAxis, tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis, angle);
        pos[2] += m_step * rotationNumber + tiltCorrection[2];

        vectormath::rotate(otherAxis, rotationAxis, angle);
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] += m_position[i];
        }
        const std::array<2, T> collimationAngles = {
            std::atan(m_fov / m_sdd) * T { 2.0 },
            std::atan(m_collimation / m_sdd) * T { 2.0 }
        };

        T weight { 1.0 };
        if (m_aecFilter)
            weight *= m_aecFilter->sampleIntensityWeight(pos);
        if (m_useXCareFilter)
            weight *= m_xcareFilter.sampleIntensityWeight(angle);

        Exposure<T> exposure(pos,
            directionCosines,
            collimationAngles,
            m_historiesPerExposure,
            weight,
            m_specterDistribution.get(),
            m_heelFilter.get(), 
            m_bowTieFilter.get());

        return exposure;
        /*
        std::array<T, 3> rotationAxis, otherAxis;
        for (std::size_t i = 0; i < 3; ++i) {
            rotationAxis[i] = m_directionCosines[i + 3];
            otherAxis[i] = m_directionCosines[i];
        }
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(rotationAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(otherAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis.data(), angle);
        pos[2] += m_step * rotationNumber + tiltCorrection[2];
        vectormath::rotate(otherAxis.data(), rotationAxis.data(), angle);
        for (std::size_t i = 0; i < 3; ++i)
            pos[i] += m_position[i];

        exposure.setPosition(pos);
        exposure.setDirectionCosines(otherAxis, rotationAxis);
        exposure.setCollimationAngles(std::atan(m_fov / m_sdd) * 2.0, std::atan(m_collimation / m_sdd) * 2.0);
        exposure.setBeamFilter(m_bowTieFilter.get());
        exposure.setHeelFilter(m_heelFilter.get());
        exposure.setSpecterDistribution(m_specterDistribution.get());
        exposure.setNumberOfHistories(m_historiesPerExposure);
        T weight { 1.0 };
        if (m_aecFilter)
            weight *= m_aecFilter->sampleIntensityWeight(pos);
        if (m_useXCareFilter)
            weight *= m_xcareFilter.sampleIntensityWeight(angle);
        exposure.setBeamIntensityWeight(weight);
        return true;
        */
    }

    void setStep(T step)
    {
        const auto absStep = std::abs(step);
        const auto n_steps = m_scanLenght / m_step;
        m_step = absStep > 0.01 ? absStep : 0.01;
        setScanLenght(m_step * n_steps);
    }

    T step(void) const { return m_step; }
    void setScanLenght(T scanLenght) override
    {
        m_scanLenght = std::max(m_step * std::ceil(std::abs(scanLenght) / m_step), m_step);
    }

    std::uint64_t totalExposures(void) const override
    {
        const std::uint64_t anglesPerRotation = static_cast<std::uint64_t>(PI_VAL<T>() * T { 2 } / m_exposureAngleStep);
        const std::uint64_t rotationNumbers = static_cast<std::uint64_t>(std::round(m_scanLenght / m_step));
        return anglesPerRotation * rotationNumbers;
    }

    T getCalibrationValue(ProgressBar<T>* progressBar) const override
    {
        const auto& me = *this;
        return ctCalibration(me, progressBar);
    }

private:
    T m_step;
};

template <Floating T>
class CTDualSource final : public CTSource<T> {
public:
    //Source overrides
    CTDualSource()
        : CTSource<T>()
    {
        m_type = Type::CTDual;
        m_sddB = m_sdd;
        m_fovB = m_fov;
        m_startAngleB = m_startAngle + PI_VAL<T>() * T { 0.5 };

        m_pitch = 1.0;
        m_tube.setAlFiltration(7.0);
        m_tubeB.setAlFiltration(7.0);
    }

    bool getExposure(Exposure<T>& exposure, std::uint64_t exposureIndexTotal) const override
    {
        T sdd, startAngle, fov;
        uint64_t exposureIndex = exposureIndexTotal / 2;
        BeamFilter* bowTie = nullptr;
        SpecterDistribution* specterDistribution = nullptr;
        HeelFilter* heelFilter = nullptr;
        T weight { 1 };
        if (exposureIndexTotal % 2 == 0) {
            sdd = m_sdd;
            startAngle = m_startAngle;
            fov = m_fov;
            bowTie = m_bowTieFilter.get();
            specterDistribution = m_specterDistribution.get();
            heelFilter = m_heelFilter.get();
            weight = m_tubeAweight;
        } else {
            sdd = m_sddB;
            startAngle = m_startAngleB;
            fov = m_fovB;
            bowTie = m_bowTieFilterB.get();
            specterDistribution = m_specterDistributionB.get();
            heelFilter = m_heelFilterB.get();
            weight = m_tubeBweight;
        }

        std::array<T, 3> pos = { 0, -m_sdd / 2.0, 0 };
        const auto angle = startAngle + m_exposureAngleStep * exposureIndex;

        auto directionCosines = m_directionCosines;
        T* rotationAxis = &directionCosines[3];
        T* otherAxis = &directionCosines[0];
        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(rotationAxis, tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(otherAxis, tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis, angle);
        pos[2] += (exposureIndex * m_exposureAngleStep) * m_collimation * m_pitch / PI_2 + tiltCorrection[2];

        vectormath::rotate(otherAxis, rotationAxis, angle);
        for (std::size_t i = 0; i < 3; ++i) {
            pos[i] += m_position[i];
        }

        const std::array<2, T> collimationAngles = {
            std::atan(fov / sdd) * T { 2.0 },
            std::atan(m_collimation / sdd) * T { 2.0 }
        };

        if (m_aecFilter)
            weight *= m_aecFilter->sampleIntensityWeight(pos);
        if (m_useXCareFilter)
            weight *= m_xcareFilter.sampleIntensityWeight(angle);

        Exposure<T> exposure(pos,
            directionCosines,
            collimationAngles,
            m_historiesPerExposure,
            weight,
            specterDistribution,
            heelFilter,
            bowTie);

        return exposure;


        /*std::array<T, 3> pos = { 0, -sdd / 2.0, 0 };

        const T angle = startAngle + m_exposureAngleStep * exposureIndex;

        std::array<T, 3> rotationAxis, otherAxis;
        for (std::size_t i = 0; i < 3; ++i) {
            rotationAxis[i] = m_directionCosines[i + 3];
            otherAxis[i] = m_directionCosines[i];
        }

        std::array<T, 3> tiltAxis = { 1, 0, 0 };
        auto tiltCorrection = pos;
        vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(rotationAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(otherAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
        vectormath::rotate(pos.data(), rotationAxis.data(), angle);
        pos[2] += (exposureIndex * m_exposureAngleStep) * m_collimation * m_pitch / PI_2 + tiltCorrection[2];

        vectormath::rotate(otherAxis.data(), rotationAxis.data(), angle);
        for (std::size_t i = 0; i < 3; ++i)
            pos[i] += m_position[i];

        exposure.setPosition(pos);
        exposure.setDirectionCosines(otherAxis, rotationAxis);
        exposure.setCollimationAngles(std::atan(fov / sdd) * 2.0, std::atan(m_collimation / sdd) * 2.0);
        exposure.setBeamFilter(bowTie);
        exposure.setSpecterDistribution(specterDistribution);
        exposure.setHeelFilter(heelFilter);
        exposure.setNumberOfHistories(m_historiesPerExposure);
        if (m_aecFilter)
            weight *= m_aecFilter->sampleIntensityWeight(pos);
        if (m_useXCareFilter)
            weight *= m_xcareFilter.sampleIntensityWeight(angle);
        exposure.setBeamIntensityWeight(weight);
        return true;
        */
    }

    T tubeAmas() const { return m_tubeAmas; }
    T tubeBmas() const { return m_tubeBmas; }
    void setTubeAmas(T mas)
    {
        m_specterValid = false;
        m_tubeAmas = std::max(T { 0.0 }, mas);
    }
    void setTubeBmas(T mas)
    {
        m_specterValid = false;
        m_tubeBmas = std::max(T { 0.0 }, mas);
    }

    Tube<T>& tubeB(void)
    {
        m_specterValid = false;
        return m_tubeB;
    };
    const Tube<T>& tubeB(void) const { return m_tubeB; }

    T maxPhotonEnergyProduced() const override { return std::max(m_tube.voltage(), m_tubeB.voltage()); }

    std::uint64_t totalExposures(void) const override
    {
        auto singleSourceExposure = static_cast<std::uint64_t>(m_scanLenght * T { 2 } * PI_VAL<T>() / (m_collimation * m_pitch * m_exposureAngleStep));
        return singleSourceExposure * 2;
    }

    std::uint64_t exposuresPerRotatition() const override
    {
        return 2 * static_cast<std::size_t>(T { 2 } * PI_VAL<T>() / m_exposureAngleStep);
    }

    T getCalibrationValue(ProgressBar<T>* = nullptr) const override
    {
        const auto& me = *this;
        return ctCalibration(me, progressBar) * m_pitch;
    }

    void setBowTieFilterB(std::shared_ptr<BowTieFilter<T>> filter) { m_bowTieFilterB = filter; }
    std::shared_ptr<BowTieFilter<T>> bowTieFilterB(void) { return m_bowTieFilterB; }
    const std::shared_ptr<BowTieFilter<T>> bowTieFilterB() const { return m_bowTieFilterB; }

    void setSourceDetectorDistanceB(T sdd)
    {
        m_specterValid = false;
        m_sddB = std::abs(sdd);
    }
    T sourceDetectorDistanceB(void) const { return m_sddB; }

    void setFieldOfViewB(T fov) { m_fovB = std::abs(fov); }
    T fieldOfViewB(void) const { return m_fovB; }

    T pitch(void) const { return m_pitch; }
    void setPitch(T pitch)
    {
        m_pitch = std::max(T { 0.01 }, pitch);
    }

    void setScanLenght(T scanLenght) override
    {
        m_scanLenght = std::max(std::abs(scanLenght), m_collimation * m_pitch * T { 0.5 });
    }

    void setStartAngleB(T angle) { m_startAngleB = angle; }
    T startAngleB(void) const { return m_startAngleB; }
    void setStartAngleDegB(T angle)
    {
        m_startAngleB = DEG_TO_RAD<T>() * angle;
    }

    T startAngleDegB(void) const
    {
        return RAD_TO_DEG<T>() * m_startAngleB;
    }

    bool validate(void) override
    {
        updateSpecterDistribution();
        return m_specterValid;
    }

protected:
    void updateSpecterDistribution() override
    {
        if (!m_specterValid) {
            const auto energyA = m_tube.getEnergy();
            const auto energyB = m_tubeB.getEnergy();
            auto specterA = m_tube.getSpecter(energyA, false);
            auto specterB = m_tubeB.getSpecter(energyB, false);

            const auto sumA = std::accumulate(specterA.cbegin(), specterA.cend(), T { 0.0 });
            const auto sumB = std::accumulate(specterB.cbegin(), specterB.cend(), T { 0.0 });
            const auto weightA = m_tubeAmas * sumA;
            const auto weightB = m_tubeBmas * sumB;

            std::transform(std::execution::par_unseq, specterA.cbegin(), specterA.cend(), specterA.begin(), [=](auto i) { return i / sumA; });
            std::transform(std::execution::par_unseq, specterB.cbegin(), specterB.cend(), specterB.begin(), [=](auto i) { return i / sumB; });

            m_tubeAweight = weightA * T { 2.0 } / (weightA + weightB);
            m_tubeBweight = weightB * T { 2.0 } / (weightA + weightB);

            m_specterDistribution = std::make_shared<SpecterDistribution>(specterA, energyA);
            m_specterDistributionB = std::make_shared<SpecterDistribution>(specterB, energyB);

            const auto heel_span_angle = std::atan(m_collimation * T { 0.5 } / m_sdd) * T { 2.0 };
            m_heelFilter = std::make_shared<HeelFilter>(m_tube, heel_span_angle);
            m_heelFilterB = std::make_shared<HeelFilter>(m_tubeB, heel_span_angle);

            m_specterValid = true;
        }
    }

private:
    Tube<T> m_tubeB;
    std::shared_ptr<SpecterDistribution<T>> m_specterDistributionB = nullptr;
    T m_sddB;
    T m_fovB;
    T m_startAngleB;
    T m_pitch = 1.0;
    T m_tubeAmas = 100.0;
    T m_tubeBmas = 100.0;
    T m_tubeBweight = -1.0;
    T m_tubeAweight = -1.0;
    std::shared_ptr<BowTieFilter<T>> m_bowTieFilterB = nullptr;
    std::shared_ptr<HeelFilter<T>> m_heelFilterB = nullptr;
};
}
