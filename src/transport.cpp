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

#include "dxmc/transport.h"

#include <algorithm>
#include <mutex>
#include <atomic>
#include <execution>
#include <thread>
#include <memory>



namespace transport {

	constexpr double ELECTRON_REST_MASS = 510.9989461; // keV
	constexpr double PI_VAL = 3.14159265358979323846264338327950288;
	constexpr double PI_VAL2 = PI_VAL + PI_VAL;
	constexpr double RUSSIAN_RULETTE_PROBABILITY = 0.8;
	constexpr double RUSSIAN_RULETTE_ENERGY_THRESHOLD = 10.0; // keV
	constexpr double ENERGY_CUTOFF_THRESHOLD = 1.0; // keV
	constexpr double KEV_TO_MJ = 1.6021773e-13;

	std::mutex DOSE_MUTEX;
	std::mutex TALLY_MUTEX;
	std::mutex VARIANCE_MUTEX;

	constexpr double N_ERROR = 1.0e-9;

	
	/*
	//waiting for c++ 20 where we get atomic references. We then have lockless threadsafe update of values in dose array 
	inline void safeEnergyAdd(double& value, const double addValue)
	{
		std::atomic_ref<double> atomicValue(value);
		atomicValue.fetch_add(addValue);
	}
	*/

	inline void safeEnergyAdd(double& value, const double addValue)
	{
		std::scoped_lock guard(DOSE_MUTEX);
		//std::lock_guard<std::mutex> lock(DOSE_MUTEX);
		value += addValue;
	}

	template<typename T>
	inline void safeTallyAdd(T& value)
	{
		std::scoped_lock guard(TALLY_MUTEX);
		//std::lock_guard<std::mutex> lock(TALLY_MUTEX);
		value += static_cast<T>(1);
	}
	inline void safeVarianceAdd(double& value, const double addValue)
	{
		std::scoped_lock guard(VARIANCE_MUTEX);
		//std::lock_guard<std::mutex> lock(VARIANCE_MUTEX);
		value += value;
	}


	void rayleightScatterLivermore(Particle& particle, unsigned char materialIdx, const AttenuationLut& attLut, std::uint64_t seed[2], double& cosAngle)
	{
		// theta is scattering angle
		// see http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

		//finding qmax
		const double qmax = attLut.momentumTransferMax(particle.energy);
		const double amax = attLut.cumFormFactorSquared(materialIdx, qmax);

		do {
			const double r1 = randomUniform<double>(seed);

			const double aatq = amax * r1;

			const double q = attLut.momentumTransfer(static_cast<size_t>(materialIdx),aatq);

			cosAngle = attLut.cosAngle(particle.energy, q);

		} while ((0.5 + cosAngle * cosAngle * 0.5) < randomUniform<double>(seed));

		
		const double theta = std::acos(cosAngle);
		const double phi = randomUniform<double>(seed, 0.0, PI_VAL2);
		vectormath::peturb<double>(particle.dir, theta, phi);
	}
	

	double comptonScatterLivermore(Particle& particle, unsigned char materialIdx, const AttenuationLut& lut, std::uint64_t seed[2], double& cosAngle)
		// see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
		// and
		// https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
	{
		const double k = particle.energy / ELECTRON_REST_MASS;
		const double emin = 1.0 / (1.0 + 2.0 * k);
		const double gmax = 1.0 / emin + emin;
		double e;
		bool rejected;
		do {
			const double r1 = randomUniform<double>(seed);
			e = r1 + (1.0 - r1) * emin;
			cosAngle = 1.0 + 1.0 / k - 1.0 / (k * e);
			const double sinthetasqr = 1.0 - cosAngle * cosAngle;

			const double g = (1 / e + e - sinthetasqr) / gmax;
			
			const double q = lut.momentumTransferFromCos(particle.energy, cosAngle);
			const double scatterFactor = lut.comptonScatterFactor(static_cast<std::size_t>(materialIdx), q);
			const double r2 = randomUniform<double>(seed);
			rejected = r2 > g* scatterFactor;
			
		} while (rejected);

		
		const double theta = std::acos(cosAngle);
		const double phi = randomUniform<double>(seed, PI_VAL2);
		vectormath::peturb<double>(particle.dir, theta, phi);

		const double E0 = particle.energy;
		particle.energy *= e;

		return E0 * (1.0 - e);
	}


	double comptonScatter(Particle& particle, std::uint64_t seed[2], double& cosAngle)
		// see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
		// and
		// https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
	{
		const double k = particle.energy / ELECTRON_REST_MASS;
		const double emin = 1.0 / (1.0 + 2.0 * k);
		const double gmax = 1.0 / emin + emin;

		bool rejected;
		double sinthetasqr, e, t;
		do {
			const double r1 = randomUniform<double>(seed);
			e = r1 + (1.0 - r1) * emin;

			t = (1.0 - e) / (k * e);
			sinthetasqr = t * (2.0 - t);

			const double g = (1.0 / e + e - sinthetasqr) / gmax;

			const double r2 = randomUniform<double>(seed);
			rejected = r2 > g;
		} while (rejected);

		cosAngle = 1.0 - t;
		const double theta = std::acos(cosAngle);
		const double phi = randomUniform<double>(seed, PI_VAL2);
		vectormath::peturb<double>(particle.dir, theta, phi);

		const double E0 = particle.energy;
		particle.energy *= e;

		return E0 * (1.0 - e);
	}


	inline bool particleInsideWorld(const World& world, const Particle& particle)
	{
		const std::array<double, 6>& extent = world.matrixExtent();
		return (particle.pos[0] > extent[0] && particle.pos[0] < extent[1]) &&
			(particle.pos[1] > extent[2] && particle.pos[1] < extent[3]) &&
			(particle.pos[2] > extent[4] && particle.pos[2] < extent[5]);
	}

	inline std::size_t indexFromPosition(const double pos[3], const World& world)
	{
		//assumes particle is inside world
		std::size_t arraypos[3];
		const std::array<double, 6>& wpos = world.matrixExtent();
		const std::array<std::size_t, 3>& wdim = world.dimensions();
		const std::array<double, 3>& wspac = world.spacing();

		for (std::size_t i = 0; i < 3; i++)
			arraypos[i] = static_cast<std::size_t>((pos[i] - wpos[i * 2]) / wspac[i]);
		std::size_t idx = arraypos[2] * wdim[0] * wdim[1] + arraypos[1] * wdim[0] + arraypos[0];
		return idx;
	}

	bool computeInteractionsForced(const double eventProbability, const AttenuationLut& lutTable, Particle& p, const unsigned char matIdx, double& doseAdress, std::uint32_t& tallyAdress, double& varianceAdress, std::uint64_t seed[2], bool& updateMaxAttenuation)
	{
		const auto atts = lutTable.photoComptRayAttenuation(matIdx, p.energy);
		const double attPhoto = atts[0];
		const double attCompt = atts[1];
		const double attRayl = atts[2];
		const double attenuationTotal = attPhoto + attCompt + attRayl;

		const double weightCorrection = eventProbability * attPhoto / attenuationTotal;
		const double energyImparted = p.energy * p.weight * weightCorrection;
		safeEnergyAdd(doseAdress, energyImparted);
		safeTallyAdd(tallyAdress);
		safeVarianceAdd(varianceAdress, energyImparted * energyImparted);
		p.weight = p.weight * (1.0 - weightCorrection); // to prevent bias		

		const double r2 = randomUniform<double>(seed);
		if (r2 < eventProbability)
		{

			const double r3 = randomUniform(seed, attCompt + attRayl);

			if (r3 <= attCompt) // Compton event
			{
				double cosangle;
#ifdef DXMC_USE_LOWENERGY_COMPTON
				const double e = comptonScatterLivermore(p, matIdx, lutTable, seed, cosangle);
#else
				const double e = comptonScatter(p, seed, cosangle);
#endif
				const double energyImparted = e * p.weight;
				safeEnergyAdd(doseAdress, energyImparted);
				safeTallyAdd(tallyAdress);
				safeVarianceAdd(varianceAdress, energyImparted*energyImparted);
				updateMaxAttenuation = true;
				if (p.energy < ENERGY_CUTOFF_THRESHOLD)
				{
					safeEnergyAdd(doseAdress, p.energy * p.weight);
					return false;
				}
			}
			else // Rayleigh scatter event
			{
				double cosangle;
				rayleightScatterLivermore(p, matIdx, lutTable, seed, cosangle);
			}

		}
		return true;
	}
	bool computeInteractions(const AttenuationLut& lutTable, Particle& p, const unsigned char matIdx, double& doseAdress, std::uint32_t& tallyAdress, double& varianceAdress, std::uint64_t seed[2], bool& updateMaxAttenuation)
	{
		const auto atts = lutTable.photoComptRayAttenuation(matIdx, p.energy);
		const double attPhoto = atts[0];
		const double attCompt = atts[1];
		const double attRayl = atts[2];

		const double r3 = randomUniform(seed, attPhoto + attCompt + attRayl);
		if (r3 <= attPhoto) // Photoelectric event
		{
			const double energyImparted = p.energy * p.weight;
			safeEnergyAdd(doseAdress, energyImparted);
			safeTallyAdd(tallyAdress);
			safeVarianceAdd(varianceAdress, energyImparted * energyImparted);
			p.energy = 0.0;
			return false;
		}
		else if (r3 <= attPhoto + attCompt) // Compton event
		{
			double cosangle;
#ifdef DXMC_USE_LOWENERGY_COMPTON
			const double e = comptonScatterLivermore(p, matIdx, lutTable, seed, cosangle);
#else
			const double e = comptonScatter(p, seed, cosangle);
#endif
			const double energyImparted = e * p.weight;
			safeEnergyAdd(doseAdress, energyImparted);
			safeTallyAdd(tallyAdress);
			safeVarianceAdd(varianceAdress, energyImparted * energyImparted);
			updateMaxAttenuation = true;
			if (p.energy < ENERGY_CUTOFF_THRESHOLD)
			{
				safeEnergyAdd(doseAdress, p.energy * p.weight);
				return false;
			}
		}
		else // Rayleigh scatter event
		{
			double cosangle;
			rayleightScatterLivermore(p, matIdx, lutTable, seed, cosangle);
		}
		return true;
	}
	
	void sampleParticleSteps(const World& world, Particle& p, std::uint64_t seed[2], Result* result)
	{
		const AttenuationLut& lutTable = world.attenuationLut();
		const double* densityBuffer = world.densityBuffer();
		const unsigned char* materialBuffer = world.materialIndexBuffer();
		const std::uint8_t* measurementBuffer = world.measurementMapBuffer();
		double* energyImparted = result->dose.data();
		auto* tally = result->nEvents.data();
		double* variance = result->variance.data();

		double maxAttenuationInv;
		bool updateMaxAttenuation = true;
		bool continueSampling = true;
		bool ruletteCandidate = true;
		while (continueSampling)
		{
			if (updateMaxAttenuation)
			{
				maxAttenuationInv = 1.0 / lutTable.maxMassTotalAttenuation(p.energy);
				updateMaxAttenuation = false;
			}
			const double r1 = randomUniform<double>(seed);
			const double stepLenght = -std::log(r1) * maxAttenuationInv * 10.0; // cm -> mm
			for (std::size_t i = 0; i < 3; i++)
				p.pos[i] += p.dir[i] * stepLenght;

			if (particleInsideWorld(world, p))
			{
				const std::size_t bufferIdx = indexFromPosition(p.pos, world);
				const auto matIdx = materialBuffer[bufferIdx];
				const double density = densityBuffer[bufferIdx];
				const double attenuationTotal = lutTable.totalAttenuation(matIdx, p.energy) * density;
				const double eventProbability = attenuationTotal * maxAttenuationInv;

				if (measurementBuffer[bufferIdx] == 0) // naive sampling
				{
					const double r2 = randomUniform<double>(seed);
					if (r2 < eventProbability) // an event will happend
					{
						continueSampling = computeInteractions(lutTable, p, matIdx, energyImparted[bufferIdx], tally[bufferIdx], variance[bufferIdx], seed, updateMaxAttenuation);
					}
				}
				else // forced photoelectric effect 
				{
					continueSampling = computeInteractionsForced(eventProbability, lutTable, p, matIdx, energyImparted[bufferIdx], tally[bufferIdx], variance[bufferIdx], seed, updateMaxAttenuation);
				}

				if (continueSampling)
				{
					// russian rulette
					if ((p.energy < RUSSIAN_RULETTE_ENERGY_THRESHOLD) && ruletteCandidate)
					{
						ruletteCandidate = false;
						const double r4 = randomUniform<double>(seed);
						if (r4 < RUSSIAN_RULETTE_PROBABILITY) {
							continueSampling = false;
						}
						else {
							constexpr double factor = 1.0 / (1.0 - RUSSIAN_RULETTE_PROBABILITY);
							p.weight *= factor;
						}
					}
				}
			}
			else // particle not inside world
			{
				continueSampling = false;
			}
		}
	}

	bool transportParticleToWorld(const World& world, Particle& particle)
		//returns false if particle do not intersect world
	{
		bool isInside = particleInsideWorld(world, particle);
		if (isInside)
			return true;

		double amin = std::numeric_limits<double>::min();
		double amax = std::numeric_limits<double>::max();

		const std::array<double, 6>& extent = world.matrixExtent();
		for (std::size_t i = 0; i < 3; i++)
		{
			if (std::abs(particle.dir[i]) > N_ERROR)
			{
				double a0, an;
				a0 = (extent[i * 2] - particle.pos[i]) / particle.dir[i];
				an = (extent[i * 2 + 1] - particle.pos[i]) / particle.dir[i];
				amin = std::max(amin, std::min(a0, an));
				amax = std::min(amax, std::max(a0, an));
			}
		}
		if (amin < amax && amin > 0.0)
		{
			for (std::size_t i = 0; i < 3; i++)
			{
				particle.pos[i] += (amin + N_ERROR) * particle.dir[i]; // making sure particle is inside world
			}
			return true;
		}
		return false;
	}


	void transport(const World& world, const Exposure& exposure, std::uint64_t seed[2], Result* result)
	{
		Particle particle;
		const std::size_t nHistories = exposure.numberOfHistories();
		for (std::size_t i = 0; i < nHistories; ++i)
		{
			//Draw a particle
			exposure.sampleParticle(particle, seed);

			// Is particle intersecting with world
			const bool isInWorld = transportParticleToWorld(world, particle);
			if (isInWorld)
				sampleParticleSteps(world, particle, seed, result);
		}
	}
	
	void parallellRun(const World& w, const Source* source, Result* result,
		const std::uint64_t expBeg, const std::uint64_t expEnd, std::uint64_t nJobs, ProgressBar* progressbar)
	{
		const std::uint64_t len = expEnd - expBeg;
		if ((len <= 1) || (nJobs <= 1))
		{
			std::uint64_t seed[2];
			randomSeed(seed);  // get random seed from OS
			Exposure exposure;
			const auto& worldBasis = w.directionCosines();
			for (std::size_t i = expBeg; i < expEnd; i++)
			{
				if (progressbar)
					if (progressbar->cancel())
						return;
				source->getExposure(exposure, i);
				exposure.alignToDirectionCosines(worldBasis);
				transport(w, exposure, seed, result);
				if (progressbar)
					progressbar->exposureCompleted();
			}
			return;
		}
		
		const auto threadLen = len / nJobs;
		std::vector<std::thread> jobs;
		jobs.reserve(nJobs-1);
		std::uint64_t expBegThread = expBeg;
		std::uint64_t expEndThread = expBegThread + threadLen;
		for (std::size_t i = 0; i < nJobs-1; ++i)
		{
			jobs.push_back(std::thread(parallellRun, w, source, result, expBegThread, expEndThread, 1, progressbar));
			expBegThread = expEndThread;
			expEndThread = expBegThread + threadLen;
		}
		expEndThread = expEnd;
		parallellRun(w, source, result, expBegThread, expEndThread, 1, progressbar);
		for (std::size_t i = 0; i < nJobs-1; ++i)
			jobs[i].join();
	}

	void parallellRunCtdi(const CTDIPhantom& w, const CTSource* source, Result* result, 
		const std::uint64_t expBeg, const std::uint64_t expEnd, std::uint64_t nJobs, ProgressBar* progressbar)
	{
		const std::uint64_t len = expEnd - expBeg;

		if ((len == 1 ) || (nJobs <= 1))
		{
			std::uint64_t seed[2];
			randomSeed(seed);  // get random seed from OS
			Exposure exposure;
			const auto& worldBasis = w.directionCosines();
			for (std::size_t i = expBeg; i < expEnd; i++)
			{
				if (progressbar)
					if (progressbar->cancel())
						return;
				source->getExposure(exposure, i);
				auto pos = source->position();
				exposure.subtractPosition(pos); // aligning to center of phantom
				exposure.setPositionZ(0.0);
				exposure.alignToDirectionCosines(worldBasis);
				exposure.setBeamIntensityWeight(1.0);
				transport(w, exposure, seed, result);
				if (progressbar)
					progressbar->exposureCompleted();
			}
			return;
		}
		
		const auto threadLen = len / nJobs;
		std::vector<std::thread> jobs;
		jobs.reserve(nJobs-1);
		std::uint64_t expBegThread = expBeg;
		std::uint64_t expEndThread = expBegThread + threadLen;
		for (std::size_t i = 0; i < nJobs-1; ++i)
		{
			//jobs.push_back(std::async(std::launch::async, parallellRunCtdi, w, source, result, expBegThread, expEndThread, 1, progressbar));
			jobs.push_back(std::thread(parallellRunCtdi, w, source, result, expBegThread, expEndThread, 1, progressbar));
			expBegThread = expEndThread;
			expEndThread = expBegThread + threadLen;
		}
		expEndThread = expEnd;
		parallellRunCtdi(w, source, result, expBegThread, expEndThread, 1, progressbar);
		for (std::size_t i = 0; i < nJobs-1; ++i)
			jobs[i].join();
	}

	void computeResultVariance(Result& res)
	{
		//Computing variance by: Var[X] = E[X**2] - E[X]**2
		//and Var[A]+ Var[A] = 2Var[A]
		auto eBeg = res.dose.cbegin();
		auto eEnd = res.dose.cend();
		auto tBeg = res.nEvents.cbegin();
		auto vBeg = res.variance.begin();
		while (eBeg != eEnd)
		{
			const double nEv = static_cast<double>(*tBeg);
			*vBeg = (*vBeg / nEv - (*eBeg * *eBeg) / (nEv * nEv)) * nEv;
			++eBeg;
			++tBeg;
			++vBeg;
		}

	}

	void energyImpartedToDose(const World & world, Result& res, const double calibrationValue)
	{
		auto spacing = world.spacing();
		const double voxelVolume = spacing[0] * spacing[1] * spacing[2] / 1000.0; // cm3
		auto density = world.densityArray();

		std::transform(std::execution::par_unseq, res.dose.cbegin(), res.dose.cend(), density->cbegin(), res.dose.begin(),
			[=](double ei, double de)->double {
				const double voxelMass = de * voxelVolume * 0.001; //kg
				return de > 0.0 ? calibrationValue * KEV_TO_MJ * ei / voxelMass : 0.0;
		});
		std::transform(std::execution::par_unseq, res.variance.cbegin(), res.variance.cend(), density->cbegin(), res.variance.begin(),
			[=](double var, double de)->double {
				const double voxelMass = de * voxelVolume * 0.001; //kg
				return de > 0.0 ? calibrationValue * KEV_TO_MJ * var / voxelMass : 0.0;
			});
	}
	
	Result run(const World & world, Source* source, ProgressBar* progressbar, bool calculateDose)
	{
		Result result;
		result.dose.resize(world.size());
		result.nEvents.resize(world.size());
		result.variance.resize(world.size());
		std::fill(result.dose.begin(), result.dose.end(), 0.0);
		std::fill(result.nEvents.begin(), result.nEvents.end(), 0);
		std::fill(result.variance.begin(), result.variance.end(), 0.0);

		if (!source)
			return result;
	
		if (!world.isValid())
			return result;

		source->updateFromWorld(world);
		source->validate();
		if(!source->isValid())
			return result;

		const std::uint64_t totalExposures = source->totalExposures();
		
		const std::uint64_t nThreads = std::thread::hardware_concurrency();
		const std::uint64_t nJobs = std::max(nThreads, static_cast<std::uint64_t>(1));
		
		if (progressbar)
		{
			progressbar->setTotalExposures(totalExposures);
			progressbar->setDoseData(result.dose.data(), world.dimensions(), world.spacing());
		}

		const auto start = std::chrono::system_clock::now();
		parallellRun(world, source, &result, 0, totalExposures, nJobs, progressbar);
		result.simulationTime = std::chrono::system_clock::now() - start;

		if (progressbar)
		{
			progressbar->clearDoseData();
			if (progressbar->cancel())
			{
				std::fill(result.dose.begin(), result.dose.end(), 0.0);
				std::fill(result.nEvents.begin(), result.nEvents.end(), 0);
				std::fill(result.variance.begin(), result.variance.end(), 0.0);
				return result;
			}
		}

		//compute variance
		computeResultVariance(result);

		if (calculateDose)
		{
			double calibrationValue = source->getCalibrationValue(progressbar);
			//energy imparted to dose
			energyImpartedToDose(world, result, calibrationValue);
		}
		return result;
	}


	Result run(const CTDIPhantom & world, CTSource* source, ProgressBar* progressbar)
	{
		Result result;
		result.dose.resize(world.size());
		result.nEvents.resize(world.size());
		result.variance.resize(world.size());
		std::fill(result.dose.begin(), result.dose.end(), 0.0);
		std::fill(result.nEvents.begin(), result.nEvents.end(), 0);
		std::fill(result.variance.begin(), result.variance.end(), 0.0);

		if (!source)
			return result;

		if (!world.isValid())
			return result;
		
		if (!source->isValid())
			return result;

		const std::uint64_t totalExposures = source->exposuresPerRotatition();
		
		const std::uint64_t nThreads = std::thread::hardware_concurrency();
		const std::uint64_t nJobs = std::max(nThreads, static_cast<std::uint64_t>(1));
		
		if (progressbar)
		{
			progressbar->setTotalExposures(totalExposures, "CTDI Calibration");
			progressbar->setDoseData(result.dose.data(), world.dimensions(), world.spacing(), ProgressBar::Axis::Z);
		}

		const auto start = std::chrono::system_clock::now();
		parallellRunCtdi(world, source, &result, 0, totalExposures, nJobs, progressbar);
		result.simulationTime = std::chrono::system_clock::now() - start;

		//compute variance
		computeResultVariance(result);

		if(progressbar)
			progressbar->clearDoseData();
		
		energyImpartedToDose(world, result, 1.0);
		return result;
	}
}
