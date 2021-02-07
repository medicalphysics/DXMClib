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

#include "xraylib.h"

#include "dxmc/material.h"
#include "dxmc/source.h"
#include "dxmc/transport.h"

#include "testutils.h"

#include <array>
#include <chrono>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

using namespace dxmc;

template <typename T>
void normalize(std::vector<T>& v)
{
    const T sum = std::reduce(v.cbegin(), v.cend());

    std::transform(v.cbegin(), v.cend(), v.begin(), [=](auto val) { return val / sum; });
}

template <typename T>
class Test : public Transport<T> {
public:
    Test()
        : Transport<T>()
    {
    }

    template <int Correction = 0>
    bool testCompton(const T energy, const Material& mat)
    {
        auto& att = this->attenuationLut();
        std::vector<Material> mats;
        mats.push_back(mat);
        att.generate(mats, T { 1 }, energy);

        std::array<T, 3> pos = { 0, 0, 0 };
        std::array<T, 3> dir = { 1, 0, 0 };
        Particle<T> p { .pos = pos, .dir = dir };
        p.energy = energy;
        const T emin = std::max(T { 1 } / (1 + 2 * energy / ELECTRON_REST_MASS<T>()) - T { 0.1 }, T { 0 });

        const T emax = 1;
        const auto k = energy / ELECTRON_REST_MASS<T>();

        std::vector<T> earr(180, 0);

        for (std::size_t i = 0; i < earr.size(); ++i) {

            earr[i] = emin + (emax - emin) / (earr.size()) * i;
        }

        std::vector<T> res(earr.size(), 0);
        RandomState state;
        std::size_t nSamp = 3e6;
        for (std::size_t i = 0; i < nSamp; ++i) {
            auto psamp = p;
            this->comptonScatter<Correction>(psamp, 0, state);
            const T e = psamp.energy / p.energy;
            auto it = std::upper_bound(earr.cbegin(), earr.cend(), e);
            auto idx = std::distance(earr.cbegin(), it);
            ++res[idx - 1];
        }

        const auto estep = (earr[1] - earr[0]) / 2;
        std::transform(earr.cbegin(), earr.cend(), earr.begin(), [=](auto e) { return e + estep; });

        std::vector<T> meas(earr.size(), 0);

        for (std::size_t i = 0; i < meas.size(); ++i) {
            const auto e = earr[i];
            if (e > emin) {
                const auto angle = std::acos(1.0 + 1.0 / k - 1.0 / (e * k));
                xrl_error* error = nullptr;
                const auto val = DCS_Compt_CP(mat.name().c_str(), e * p.energy, angle, &error);
                if (val > 0)
                    meas[i] = val;
            }
        }

        normalize(meas);
        normalize(res);
        std::cout << "Compton differential scattering cross section\n";
        std::cout << "with initial energy " << energy << " keV in " << mat.prettyName() << " material\n";
        std::cout << "Sampling " << nSamp << " interactions with RMS differential cross section deviation of\n";
        const auto ms = std::transform_reduce(res.cbegin(), res.cend(), meas.cbegin(), T { 0 }, std::plus<>(), [](auto m, auto e) { return (m - e) * (m - e); });
        const auto rmsd = std::sqrt(ms / res.size());
        std::cout << rmsd << "\n";
        std::cout << "K2/K1, dxmclib, meas, diff, diff [%]\n";

        for (std::size_t i = 0; i < meas.size(); ++i) {
            std::cout << earr[i] << ", " << res[i] << ", " << meas[i] << ", " << res[i] - meas[i] << ", " << (res[i] - meas[i]) / meas[i] * 100 << "\n";
        }

        if (rmsd < 0.001)
            return true;
        return false;
    }
    template <int Correction = 0>
    bool testRayleight(const T energy, const Material& mat)
    {
        auto& att = this->attenuationLut();
        std::vector<Material> mats;
        mats.push_back(mat);
        att.generate(mats, T { 1 }, energy * 2);

        std::array<T, 3> pos = { 0, 0, 0 };
        std::array<T, 3> dir = { 1, 0, 0 };
        Particle<T> p { .pos = pos, .dir = dir };
        p.energy = energy;

        std::vector<T> cosAng(180, 0);

        for (std::size_t i = 0; i < cosAng.size(); ++i) {
            cosAng[i] = -1 + T { 2 } / (cosAng.size()) * i;
        }

        std::vector<T> res(cosAng.size(), 0);
        RandomState state;
        std::size_t nSamp = 3e6;
        for (std::size_t i = 0; i < nSamp; ++i) {
            auto psamp = p;
            this->rayleightScatter<Correction>(psamp, 0, state);
            const auto cos = vectormath::dot(p.dir.data(), psamp.dir.data());
            auto it = std::upper_bound(cosAng.cbegin(), cosAng.cend(), cos);
            auto idx = std::distance(cosAng.cbegin(), it);
            ++res[idx - 1];
        }

        const auto step = (cosAng[1] - cosAng[0]) / 2;
        std::transform(cosAng.cbegin(), cosAng.cend(), cosAng.begin(), [=](auto e) { return e + step; });

        std::vector<T> meas(cosAng.size(), 0);

        for (std::size_t i = 0; i < meas.size(); ++i) {
            const auto angle = std::acos(cosAng[i]);
            xrl_error* error = nullptr;
            const auto val = DCS_Rayl_CP(mat.name().c_str(), p.energy, angle, &error);
            if (val > 0)
                meas[i] = val;
        }

        normalize(meas);
        normalize(res);
        std::cout << "Rayleight differential scattering cross section\n";
        std::cout << "with initial energy " << energy << " keV in " << mat.prettyName() << " material\n";
        std::cout << "Sampling " << nSamp << " interactions with RMS differential cross section deviation of\n";
        const auto ms = std::transform_reduce(res.cbegin(), res.cend(), meas.cbegin(), T { 0 }, std::plus<>(), [](auto m, auto e) { return (m - e) * (m - e); });
        const auto rmsd = std::sqrt(ms / res.size());
        std::cout << rmsd << "\n";
        std::cout << "angle, dxmclib, meas, diff, diff [%]\n";

        for (std::size_t i = 0; i < meas.size(); ++i) {
            std::cout << std::acos(cosAng[i]) << ", " << res[i] << ", " << meas[i] << ", " << res[i] - meas[i] << ", " << (res[i] - meas[i]) / meas[i] * 100 << "\n";
        }

        if (rmsd < 0.001)
            return true;
        return false;
    }

    template <bool locks = true>
    void testSafeValueAddWorkerDose(Result<T>* res, T addValue, std::size_t N, std::size_t idx = 0) const
    {
        if constexpr (locks) {
            for (std::size_t i = 0; i < N; ++i) {
                this->safeValueAdd(res->dose[idx], addValue, res->locks[idx].dose);
            }
        } else {
            for (std::size_t i = 0; i < N; ++i) {
                this->safeValueAdd(res->dose[idx], addValue);
            }
        }
    }
    template <bool locks = true>
    void testSafeValueAddWorkerEvents(Result<T>* res, std::uint32_t addValue, std::size_t N, std::size_t idx = 0) const
    {
        if constexpr (locks) {
            for (std::size_t i = 0; i < N; ++i) {
                this->safeValueAdd(res->nEvents[idx], addValue, res->locks[idx].nEvents);
            }
        } else {
            for (std::size_t i = 0; i < N; ++i) {
                this->safeValueAdd(res->nEvents[idx], addValue);
            }
        }
    }

    bool testSafeValueAdd()
    {

        const std::size_t N = 1E6;

        Result<T> res(1);

        const T addValue { 1 };
        const std::uint64_t nThreads = std::thread::hardware_concurrency();
        const std::uint64_t nJobs = std::max(nThreads, static_cast<std::uint64_t>(8));
        std::vector<std::thread> jobs;
        jobs.reserve(nJobs);
        auto start = std::chrono::high_resolution_clock::now();
        for (std::size_t i = 0; i < nJobs; ++i) {
            jobs.emplace_back(&Test<T>::testSafeValueAddWorkerDose<false>, this, &res, addValue, N, 0);
        }
        for (auto& job : jobs) {
            job.join();
        }
        auto duration = std::chrono::high_resolution_clock::now() - start;
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(duration);

        const T expected = N * addValue * nJobs;
        const T value = res.dose[0];
        std::cout << "Test safe atomic add:\n";
        std::cout << "Time: " << time.count() / 1000.0;
        std::cout << "s\nAdded " << addValue << " " << N << " times in " << nJobs;
        std::cout << " threads.\nResult is expected to be " << expected << " and is ";
        std::cout << value << "\n";

        jobs.clear();
        jobs.reserve(nJobs);
        const std::uint32_t addValueEvents { 1 };
        for (std::size_t i = 0; i < nJobs; ++i) {
            jobs.emplace_back(&Test<T>::testSafeValueAddWorkerEvents<false>, this, &res, addValueEvents, N, 0);
        }
        for (auto& job : jobs) {
            job.join();
        }

        const std::uint32_t expectedEvents = N * addValueEvents * nJobs;
        const std::uint32_t valueEvents = res.nEvents[0];
        std::cout << "Test safe atomic add:\nAdded " << addValueEvents << " " << N << " times in " << nJobs;
        std::cout << " threads.\nResult is expected to be " << expectedEvents << " and is ";
        std::cout << valueEvents << "\n";

        return isEqual(expected, value) && isEqual(expectedEvents, valueEvents);
    }
};

template <typename T>
void testTransport()
{
    std::vector<Material> mats;
    mats.emplace_back((6));
    mats.emplace_back((8));
    mats.emplace_back((13));
    mats.emplace_back((82));
    mats.emplace_back(("H2O"));
    mats.back().setStandardDensity(1.0);

    const T Rayenergy = 10;
    const T Comenergy = 50;
    Test<T> t;
    t.testRayleight<1>(Rayenergy, mats[4]);
    t.testCompton<1>(Comenergy, mats[2]);
}

int main()
{
    testTransport<float>();

    return 0;
}