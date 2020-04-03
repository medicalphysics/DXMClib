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


#include "dxmc/dxmcrandom.h"
#include "dxmc/beamfilters.h"
#include "dxmc/exposure.h"
#include "dxmc/progressbar.h"
#include "dxmc/tube.h"
#include "dxmc/world.h"

#include <vector>
#include <memory>
#include <array>


class Source
{
public:
	enum class Type { None, CTSpiral, CTAxial, DX, CTDual, Pencil, Isotropic };
    Source();
    virtual bool getExposure(Exposure& exposure, std::uint64_t i) const = 0;
	virtual double maxPhotonEnergyProduced() const {return  Tube::maxVoltage(); }

	void setPosition(const std::array<double, 3>& position) { m_position = position; }
	void setPosition(double x, double y, double z) {
		m_position[0] = x; m_position[1] = y; m_position[2] = z;
	};
	std::array<double, 3>& position(void) { return m_position; }
	const std::array<double, 3>& position(void) const { return m_position; }
	
	void setDirectionCosines(const std::array<double, 6>& cosines);
	const std::array<double, 6>& directionCosines(void) const { return m_directionCosines; }
	std::array<double, 6>& directionCosines(void) { return m_directionCosines; }

	void setHistoriesPerExposure(std::uint64_t histories);
	std::uint64_t historiesPerExposure(void) const;

	virtual std::uint64_t totalExposures(void) const = 0;

	Source::Type type() const { return m_type; }

	virtual double getCalibrationValue(ProgressBar* progress=nullptr) const = 0;
	
	virtual bool isValid(void) const = 0;
	virtual bool validate(void) = 0;
	virtual void updateFromWorld(const World& world) { }

protected:
	std::array<double, 3> m_position = { 0,0,0 };
	std::array<double, 6> m_directionCosines = { 1,0,0,0,1,0 };
	std::uint64_t m_historiesPerExposure = 1;
	Type m_type = Type::None;
private:
	void normalizeDirectionCosines(void);
};


class PencilSource final : public Source
{
public:
	PencilSource();
	virtual ~PencilSource() = default;
	bool getExposure(Exposure& exposure, std::uint64_t i) const override;

	void setPhotonEnergy(double energy);
	double photonEnergy() const { return m_photonEnergy; }

	double maxPhotonEnergyProduced() const override { return m_photonEnergy; }

	std::uint64_t totalExposures(void) const override { return m_totalExposures; }
	void setTotalExposures(std::uint64_t exposures) { if (exposures > 0) m_totalExposures = exposures; }

	void setAirDose(double Gycm2) { if (Gycm2 > 0.0) m_airDose = Gycm2; }
	double airDose(void) const { return m_airDose; } // Gycm2

	double getCalibrationValue(ProgressBar* = nullptr) const override;

	bool isValid(void) const override { return true; }
	bool validate(void) override {return true; }

protected:
	double m_photonEnergy = 100.0;
	double m_airDose = 1.0;
	std::uint64_t m_totalExposures = 10;
};


class IsotropicSource final : public Source
{
public:
	IsotropicSource();
	virtual ~IsotropicSource() = default;
	bool getExposure(Exposure& exposure, std::uint64_t exposureNumber) const override;
	double maxPhotonEnergyProduced() const override { 
		return m_maxPhotonEnergy; };
	void setTotalExposures(std::uint64_t nExposures) { 
		m_totalExposures = nExposures; };
	std::uint64_t totalExposures() const override { 
		return m_totalExposures; };
	
	double getCalibrationValue(ProgressBar* progress = nullptr) const override {
		return 1.0;
	};

	bool isValid() const override;
	bool validate() override;

	void setSpecter(const std::vector<double>& weights, const std::vector<double>& energies);
	void setCollimationAngles(double xRad, double yRad);
	std::pair<double, double> collimationAngles() const;
protected:

private:
	std::uint64_t m_totalExposures = 1;
	std::array<double, 2> m_collimationAngles = { 0,0 };
	std::shared_ptr<SpecterDistribution> m_specterDistribution = nullptr;
	double m_maxPhotonEnergy = 1.0;
};

class DXSource final : public Source
{
public:
	DXSource();
	virtual ~DXSource() = default;
	bool getExposure(Exposure& exposure, std::uint64_t i) const override;
	
	Tube& tube(void) { m_specterValid = false; return m_tube; };
	const Tube& tube(void) const { return m_tube; }

	double maxPhotonEnergyProduced() const override { return m_tube.voltage(); }

	std::uint64_t totalExposures(void) const override;
	void setTotalExposures(std::uint64_t exposures);

	void setCollimationAngles(const std::array<double, 2>& radians);
	const std::array<double, 2>& collimationAngles(void) const;

	void setCollimationAnglesDeg(const std::array<double, 2>& degrees);
	const std::array<double, 2> collimationAnglesDeg(void) const;

	void setFieldSize(const std::array<double, 2>& mm);
	const std::array<double, 2>& fieldSize(void) const;

	void setSourceDetectorDistance(double mm);
	double sourceDetectorDistance() const;

	void setSourceAngles(double primaryAngle, double secondaryAngle);
	void setSourceAngles(const std::array<double, 2>& angles);
	std::array<double, 2> sourceAngles() const;
	void setSourceAnglesDeg(double primaryAngle, double secondaryAngle);
	void setSourceAnglesDeg(const std::array<double, 2>& angles);
	std::array<double, 2> sourceAnglesDeg() const;

	void setTubeRotation(double angle);
	double tubeRotation() const;
	void setTubeRotationDeg(double angle);
	double tubeRotationDeg() const;

	void setDap(double Gycm2) { if (Gycm2 > 0.0) m_dap = Gycm2; }
	double dap(void) const { return m_dap; } // Gycm2

	double getCalibrationValue(ProgressBar* = nullptr) const override;
	
	const std::array<double, 3> tubePosition(void) const;

	bool isValid(void) const override { return m_specterValid; }
	bool validate(void) override { updateSpecterDistribution(); return m_specterValid; }

	void setModelHeelEffect(bool on) { m_modelHeelEffect = on; };
	bool modelHeelEffect() const { return m_modelHeelEffect; };

protected:
	void updateFieldSize(const std::array<double, 2>& collimationAngles);
	void updateCollimationAngles(const std::array<double, 2>& fieldSize);
	void updateSpecterDistribution();
	std::array<double, 6> zeroDirectionCosines() const;

private:
	double m_sdd = 1000.0;
	double m_dap = 1.0; // Gycm2
	std::array<double, 2> m_fieldSize;
	std::array<double, 2> m_collimationAngles;
	std::uint64_t m_totalExposures = 1000;
	Tube m_tube;
	double m_tubeRotationAngle = 0.0;
	bool m_specterValid = false;
	std::shared_ptr<SpecterDistribution> m_specterDistribution=nullptr;
	std::shared_ptr<HeelFilter> m_heelFilter = nullptr;
	bool m_modelHeelEffect = true;
};


class CTSource : public Source
{
public:
    CTSource();
	virtual ~CTSource() = default;
	virtual bool getExposure(Exposure& exposure, std::uint64_t i) const override = 0;

	Tube& tube(void) { m_specterValid = false; return m_tube; };
	const Tube& tube(void) const { return m_tube; }

	virtual double maxPhotonEnergyProduced() const override {return m_tube.voltage(); }

	void setBowTieFilter(std::shared_ptr<BowTieFilter> filter) { m_bowTieFilter = filter; }
	std::shared_ptr<BowTieFilter> bowTieFilter(void) { return m_bowTieFilter; }
	const std::shared_ptr<BowTieFilter> bowTieFilter() const { return m_bowTieFilter; }

	void setSourceDetectorDistance(double sdd);
    double sourceDetectorDistance(void) const;
    
	void setCollimation(double mmCollimation);
    double collimation(void) const;
    
	void setFieldOfView(double fov);
    double fieldOfView(void) const;
	
	void setGantryTiltAngle(double angle);
	double gantryTiltAngle() const;
	void setGantryTiltAngleDeg(double angle);
	double gantryTiltAngleDeg() const;

    void setStartAngle(double angle);
    double startAngle(void) const;
	void setStartAngleDeg(double angle);
	double startAngleDeg(void) const;
    
	void setExposureAngleStep(double);
    double exposureAngleStep(void) const;
	void setExposureAngleStepDeg(double);
	double exposureAngleStepDeg(void) const;
    
    virtual void setScanLenght(double scanLenght);
    double scanLenght(void) const;

	void setCtdiVol(double ctdivol) { if (ctdivol > 0.0) m_ctdivol = ctdivol; }
	double ctdiVol(void) const { return m_ctdivol; }

	bool useXCareFilter() const { return m_useXCareFilter; }
	void setUseXCareFilter(bool use) { m_useXCareFilter = use; }
	XCareFilter& xcareFilter() { return m_xcareFilter; }
	const XCareFilter& xcareFilter() const { return m_xcareFilter; }

	virtual std::uint64_t totalExposures(void) const override = 0;

	void setCtdiPhantomDiameter(std::uint64_t mm) { if (mm > 160) m_ctdiPhantomDiameter = mm; else m_ctdiPhantomDiameter = 160; }
	std::uint64_t ctdiPhantomDiameter(void) const { return m_ctdiPhantomDiameter; }

	virtual double getCalibrationValue(ProgressBar* = nullptr) const = 0;
	
	virtual std::uint64_t exposuresPerRotatition() const;
	
	bool isValid(void) const override { return m_specterValid; };
	virtual bool validate(void) override { updateSpecterDistribution(); return m_specterValid; };
	
	void setAecFilter(std::shared_ptr<AECFilter> filter) { m_aecFilter = filter; }
	std::shared_ptr<AECFilter> aecFilter(void) { return m_aecFilter; }

	virtual void updateFromWorld(const World& world) override { if (m_aecFilter) m_aecFilter->updateFromWorld(world); }

	void setModelHeelEffect(bool on) { m_modelHeelEffect = on; };
	bool modelHeelEffect() const { return m_modelHeelEffect; };

protected:
	template<typename T>
	static double ctCalibration(T& ctSource, ProgressBar* progress = nullptr);
	virtual void updateSpecterDistribution();

    double m_sdd;
    double m_collimation;
    double m_fov;
    double m_startAngle;
    double m_exposureAngleStep;
    double m_scanLenght;
	double m_ctdivol = 1.0;
	double m_gantryTiltAngle = 0.0;
	std::shared_ptr<AECFilter> m_aecFilter = nullptr;
	std::uint64_t m_ctdiPhantomDiameter = 320;
	std::shared_ptr<BowTieFilter> m_bowTieFilter=nullptr;
	XCareFilter m_xcareFilter;
	bool m_useXCareFilter = false;
	bool m_specterValid = false;
	Tube m_tube;
	std::shared_ptr<SpecterDistribution> m_specterDistribution=nullptr;
	std::shared_ptr<HeelFilter> m_heelFilter = nullptr;
	bool m_modelHeelEffect = true;
};

class CTSpiralSource final : public CTSource
{
public:
	CTSpiralSource();
	bool getExposure(Exposure& exposure, std::uint64_t i) const override;
	void setPitch(double pitch);
	double pitch(void) const;
	void setScanLenght(double scanLenght) override;
	std::uint64_t totalExposures(void) const override;
	double getCalibrationValue(ProgressBar* = nullptr) const override;
protected:
private:
	double m_pitch;
};

class CTAxialSource final : public CTSource
{
public:
	CTAxialSource();
	bool getExposure(Exposure& exposure, std::uint64_t i) const override;
	void setStep(double step);
	double step(void) const;
	void setScanLenght(double scanLenght) override;
	std::uint64_t totalExposures(void) const override;
	double getCalibrationValue(ProgressBar* progressBar) const override;
protected:
private:
	double m_step;
};

class CTDualSource final : public CTSource
{
public:

	//Source overrides
	CTDualSource();
	bool getExposure(Exposure& exposure, std::uint64_t i) const override;

	double tubeAmas() const { return m_tubeAmas; }
	double tubeBmas() const { return m_tubeBmas; }
	void setTubeAmas(double mas) { m_specterValid = false; m_tubeAmas = std::max(0.0, mas); }
	void setTubeBmas(double mas) { m_specterValid = false; m_tubeBmas = std::max(0.0, mas); }

	Tube& tubeB(void) { m_specterValid = false; return m_tubeB; };
	const Tube& tubeB(void) const { return m_tubeB; }

	double maxPhotonEnergyProduced()const override { return std::max(m_tube.voltage(), m_tubeB.voltage()); }

	std::uint64_t totalExposures(void) const override;
	std::uint64_t exposuresPerRotatition() const override;
	double getCalibrationValue(ProgressBar* = nullptr) const override;

	void setBowTieFilterB(std::shared_ptr<BowTieFilter> filter) { m_bowTieFilterB = filter; }
	std::shared_ptr<BowTieFilter> bowTieFilterB(void) { return m_bowTieFilterB; }
	const std::shared_ptr<BowTieFilter> bowTieFilterB() const { return m_bowTieFilterB; }

	void setSourceDetectorDistanceB(double sdd) { m_specterValid = false; m_sddB = std::abs(sdd); }
	double sourceDetectorDistanceB(void) const { return m_sddB; }	

	void setFieldOfViewB(double fov) { m_fovB = std::abs(fov); }
	double fieldOfViewB(void) const { return m_fovB; }

	void setPitch(double pitch);
	double pitch(void) const;
	void setScanLenght(double scanLenght) override;

	void setStartAngleB(double angle) { m_startAngleB = angle; }
	double startAngleB(void) const { return m_startAngleB; }
	void setStartAngleDegB(double angle);
	double startAngleDegB(void) const;

	bool validate(void) override { updateSpecterDistribution(); return m_specterValid; };
protected:
	void updateSpecterDistribution() override;
private:
	Tube m_tubeB;
	std::shared_ptr<SpecterDistribution> m_specterDistributionB = nullptr;
	double m_sddB;
	double m_fovB;
	double m_startAngleB;
	double m_pitch = 1.0;
	double m_tubeAmas = 100.0;
	double m_tubeBmas = 100.0;
	double m_tubeBweight = -1.0;
	double m_tubeAweight = -1.0;
	std::shared_ptr<BowTieFilter> m_bowTieFilterB = nullptr;
	std::shared_ptr<HeelFilter> m_heelFilterB = nullptr;
};
