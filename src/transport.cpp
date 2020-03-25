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
#include <future>
#include <thread>
#include <memory>



namespace transport {

	constexpr double ELECTRON_REST_MASS = 510.9989461; // keV
	constexpr double PI_VAL = 3.14159265358979323846264338327950288;
	constexpr double PI_VAL2 = PI_VAL + PI_VAL;
	constexpr double RUSSIAN_RULETTE_PROBABILITY = 0.8;
	constexpr double RUSSIAN_RULETTE_ENERGY_THRESHOLD = 10.0; // keV
	constexpr double RUSSIAN_RULETTE_WEIGHT_THRESHOLD = 0.1; 
	constexpr double ENERGY_CUTOFF_THRESHOLD = 1.0; // keV
	constexpr double KEV_TO_MJ = 1.6021773e-13;

	std::mutex DOSE_MUTEX;
	std::mutex TALLY_MUTEX;

	constexpr double N_ERROR = 1.0e-9;

	/*template<typename T>
	void findNearestIndices(const T value, const std::vector<T>& vec, std::size_t& first, std::size_t& last)
	{
		auto beg = vec.begin();
		auto end = vec.end();
		auto upper = std::lower_bound(beg, end, value);

		if (upper == end)
			last = vec.size() - 1;
		else if (upper == beg)
			last = 2;
		else
			last = static_cast<std::size_t>(std::distance(beg, upper));
		first = last - 1;
		return;
	}

	template <typename T>
	inline T interp(T x, T x0, T x1, T y0, T y1)
	{
		return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
	}



	template <typename T>
	T interpolate(const std::vector<T>& xArr, const std::vector<T>& yArr, const T x)
	{
		auto first = xArr.begin();
		auto last = xArr.end();
		auto upper = std::lower_bound(first, last, x);

		if (upper == xArr.cend())
			return yArr.back();

		if (upper == xArr.cbegin())
			return yArr.front();

		const double x1 = *upper;

		auto yupper = yArr.begin();
		auto dist = std::distance(xArr.begin(), upper);
		std::advance(yupper, dist);

		const double y1 = *yupper;

		std::advance(yupper, -1);
		const double y0 = *yupper;

		std::advance(upper, -1);
		const double x0 = *upper;

		return interp(x, x0, x1, y0, y1);
	}*/

	/*
	//waiting for c++ 20 where we get atomic references. We then have lockless threadsafe update of values in dose array 
	inline void safeValueAdd(double& value, const double addValue)
	{
		std::atomic_ref<double> atomicValue(value);
		atomicValue.fetch_add(addValue);
	}
	*/

	inline void safeValueAdd(double& value, const double addValue)
	{
		std::lock_guard<std::mutex> lock(DOSE_MUTEX);
		value += addValue;
	}

	template<typename T>
	inline void safeTallyAdd(T& value)
	{
		std::lock_guard<std::mutex> lock(TALLY_MUTEX);
		value += (T)1;
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

	
	void sampleParticleSteps(const World& world, Particle& p, std::uint64_t seed[2], Result* result)
	{
		const AttenuationLut& lutTable = world.attenuationLut();
		const double* densityBuffer = world.densityBuffer();
		const unsigned char* materialBuffer = world.materialIndexBuffer();
		double* energyImparted = result->dose.data();
		auto tally = result->nEvents.data();

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
								
				const double attenuationTotal = lutTable.totalAttenuation(matIdx, p.energy) * densityBuffer[bufferIdx];

				const double r2 = randomUniform<double>(seed);
				if (r2 < (attenuationTotal * maxAttenuationInv)) // An event will happend
				{
					const auto atts = lutTable.photoComptRayAttenuation(matIdx, p.energy);
					const double attPhoto = atts[0];
					const double attCompt = atts[1];
					const double attRayl = atts[2];

					const double r3 = randomUniform(seed, attPhoto + attCompt + attRayl);
					if (r3 <= attPhoto) // Photoelectric event
					{
						safeValueAdd(energyImparted[bufferIdx], p.energy * p.weight);
						safeTallyAdd(tally[bufferIdx]);
						p.energy = 0.0;
						continueSampling = false;
					}
					else if (r3 <= attPhoto + attCompt) // Compton event
					{
						double cosangle;
#ifdef DXMC_USE_LOWENERGY_COMPTON
						const double e = comptonScatterLivermore(p, matIdx, world.attenuationLut(), seed, cosangle);
#else
						const double e = comptonScatter(p, seed, cosangle);
#endif
						safeValueAdd(energyImparted[bufferIdx], e * p.weight);
						safeTallyAdd(tally[bufferIdx]);
						updateMaxAttenuation = true;
					}
					else // Rayleigh scatter event
					{
						double cosangle;
						rayleightScatterLivermore(p, matIdx, world.attenuationLut(), seed, cosangle);
					}

					if (continueSampling)
					{
						if (p.energy < ENERGY_CUTOFF_THRESHOLD)
						{
							safeValueAdd(energyImparted[bufferIdx], p.energy * p.weight);
							continueSampling = false;
						}

						// russian rulette
						if ((p.energy < RUSSIAN_RULETTE_ENERGY_THRESHOLD || p.weight < RUSSIAN_RULETTE_WEIGHT_THRESHOLD) && ruletteCandidate)
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
	
	std::uint64_t parallellRun(const World& w, const Source* source, Result* result,
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
						return 0;
				source->getExposure(exposure, i);
				exposure.alignToDirectionCosines(worldBasis);
				transport(w, exposure, seed, result);
				if (progressbar)
					progressbar->exposureCompleted();
			}
			return source->historiesPerExposure() * len;
		}
		
		const auto threadLen = len / nJobs;
		std::vector<std::future<std::uint64_t>> jobs;
		jobs.reserve(nJobs);
		std::uint64_t expBegThread = expBeg;
		std::uint64_t expEndThread = expBegThread + threadLen;
		for (std::size_t i = 0; i < nJobs; ++i)
		{
			if (i == nJobs - 1)
				expEndThread = expEnd;
			jobs.push_back(std::async(std::launch::async, parallellRun, w, source, result, expBegThread, expEndThread, 1, progressbar));
			expBegThread = expEndThread;
			expEndThread = expBegThread + threadLen;
		}
		std::uint64_t sum = 0;
		for (std::size_t i = 0; i < nJobs; ++i)
			sum += jobs[i].get();
		return sum;
	}

	std::uint64_t parallellRunCtdi(const CTDIPhantom& w, const CTSource* source, Result* result, 
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
						return 0;
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
			return source->historiesPerExposure() * len;
		}
		
		const auto threadLen = len / nJobs;
		std::vector<std::future<std::uint64_t>> jobs;
		jobs.reserve(nJobs);
		std::uint64_t expBegThread = expBeg;
		std::uint64_t expEndThread = expBegThread + threadLen;
		for (std::size_t i = 0; i < nJobs; ++i)
		{
			if (i == nJobs - 1)
				expEndThread = expEnd;
			jobs.push_back(std::async(std::launch::async, parallellRunCtdi, w, source, result, expBegThread, expEndThread, 1, progressbar));
			expBegThread = expEndThread;
			expEndThread = expBegThread + threadLen;
		}
		std::uint64_t sum = 0;
		for (std::size_t i = 0; i < nJobs; ++i)
			sum += jobs[i].get();
		return sum;
	}


	void energyImpartedToDose(const World & world, std::vector<double>& energyImparted, const double calibrationValue)
	{
		auto spacing = world.spacing();
		const double voxelVolume = spacing[0] * spacing[1] * spacing[2] / 1000.0; // cm3
		auto density = world.densityArray();
		
		std::transform(std::execution::par_unseq, energyImparted.begin(), energyImparted.end(), density->begin(), energyImparted.begin(), 
			[=](double ei, double de)->double {
				const double voxelMass = de * voxelVolume * 0.001; //kg
				return de > 0.0 ? calibrationValue * KEV_TO_MJ * ei / voxelMass : 0.0;
		});
	}
	
	Result run(const World & world, Source* source, ProgressBar* progressbar, bool calculateDose)
	{
		Result result;
		result.dose.resize(world.size());
		result.nEvents.resize(world.size());
		std::fill(result.dose.begin(), result.dose.end(), 0.0);
		std::fill(result.nEvents.begin(), result.nEvents.end(), 0);

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

		auto nHistories = parallellRun(world, source, &result, 0, totalExposures, nJobs, progressbar);
		
		if (progressbar)
		{
			progressbar->clearDoseData();
			if (progressbar->cancel())
			{
				std::fill(result.dose.begin(), result.dose.end(), 0.0);
				std::fill(result.nEvents.begin(), result.nEvents.end(), 0);

				return result;
			}
		}
		if (calculateDose)
		{
			double calibrationValue = source->getCalibrationValue(progressbar);
			//energy imparted to dose
			energyImpartedToDose(world, result.dose, calibrationValue);
		}
		return result;
	}


	Result run(const CTDIPhantom & world, CTSource* source, ProgressBar* progressbar)
	{
		Result result;
		result.dose.resize(world.size());
		result.nEvents.resize(world.size());
		std::fill(result.dose.begin(), result.dose.end(), 0.0);
		std::fill(result.nEvents.begin(), result.nEvents.end(), 0);

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
		parallellRunCtdi(world, source, &result, 0, totalExposures + 1, nJobs, progressbar);
		if(progressbar)
			progressbar->clearDoseData();
		
		energyImpartedToDose(world, result.dose, 1.0);
		return result;
	}
}
