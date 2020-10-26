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

#include "dxmc/transport.h"

#include <array>
#include <functional>
#include <iostream>
#include <thread>
#include <vector>

using namespace dxmc;

constexpr double ERRF = 1e-6;

template <Floating T>
bool isEqual(T f1, T f2)
{
    return std::abs(f1 - f2) < ERRF;
}

template <typename T>
bool isEqual(T f1, T f2)
{
    return f1 == f2;
}

template <typename T>
class Test : public Transport<T> {
public:
    Test()
        : Transport<T>()
    {
    }
    bool testCompton(T energy, std::size_t Z)
    {

        auto& att = attenuationLut();
        std::vector<Material> materials;
        materials.emplace_back(Material(Z));
        att.generate(materials, T { 1 }, T { 60 });

        std::array<T, 3> pos = { 0, 0, 0 };
        std::array<T, 3> dir = { 1, 0, 0 };
        Particle<T> p { .pos = pos, .dir = dir };

        const T emin = T { 1 } / (1 + 2 * energy / ELECTRON_REST_MASS<T>());
        const T emax = 1;
        const auto k = energy / ELECTRON_REST_MASS<T>();

        std::array<T, 128> res;
        std::fill(res.begin(), res.end(), T { 0 });
        std::array<T, res.size()> earr;
        for (std::size_t i = 0; i < res.size(); ++i) {
            earr[i] = emin + ((emax - emin) / res.size()) * i;
        }

        RandomState state;

        for (std::size_t i = 0; i < 1e7; ++i) {
            p.energy = energy;
            T dummy;
            //comptonScatterLivermore(p, 0, state, dummy);
            comptonScatter(p, 0, state, dummy);
            const T e = p.energy / energy;
            auto it = std::upper_bound(earr.cbegin(), earr.cend(), e);
            auto idx = std::distance(earr.cbegin(), it);

            ++res[idx - 1];
        }

        std::array<T, res.size()> meas;
        std::fill(meas.begin(), meas.end(), T { 0 });
        for (std::size_t i = 0; i < meas.size(); ++i) {
            const auto estep = (earr[1] - earr[0]) / 2;
            const auto e = earr[i] + estep;

            const T a1 = std::log(1 + 2 * k);
            const T a2 = (1 - emin * emin) / 2;
            const T t = (1 - e) / (k * e);
            const T sin2 = t * (2 - t);
            const T cosAngle = 1 - t;
            const T gmax = 1 / emin + emin;
            const T p = (1 / e + e - sin2) / gmax;

            meas[i] = p;

            const auto q = att.momentumTransferFromCos(energy, cosAngle);
            const auto scatterFactor = att.comptonScatterFactor(0, q);
            meas[i] *= scatterFactor;
        }

        const auto meas_sum = std::accumulate(meas.cbegin(), meas.cend(), T { 0 });
        const auto res_sum = std::accumulate(res.cbegin(), res.cend(), T { 0 });
        std::transform(meas.cbegin(), meas.cend(), meas.begin(), [=](auto v) { return v / meas_sum; });
        std::transform(res.cbegin(), res.cend(), res.begin(), [=](auto v) { return v / res_sum; });

        std::cout << "Compton differential scattering cross section\n";
        std::cout << "with initial energy " << energy << " keV in material Z=" << Z << "\n";
        std::cout << "Sampling 1E7 interactions with RMS differential cross section deviation of\n";
        const auto ms = std::transform_reduce(res.cbegin(), res.cend(), meas.cbegin(), T { 0 }, std::plus<>(), [](auto m, auto e) { return (m - e) * (m - e); });
        const auto rmsd = std::sqrt(ms / res.size());
        std::cout << rmsd << "\n";
        /*
        std::cout << "e, DXMC, xraylib, scatterFactor dxmc, scatterFactor xraylib, momtrans dxmc, momtrans xraylib\n";
        for (std::size_t i = 0; i < res.size(); ++i) {
            const T e = emin + ((emax - emin) / res.size()) * i;
            const auto t = (1 / e - 1) / k;
            const auto sin2 = t * (2 - t);
            const auto cosAngle = 1 - t;
            const auto angle = std::acos(cosAngle) * RAD_TO_DEG<T>();

            std::cout << e << ", " << res[i] << ", ";
            std::cout << meas[i] << ", ";

            const auto q = att.momentumTransferFromCos(p.energy, cosAngle);
            const auto scatterFactor = att.comptonScatterFactor(0, q);
            std::cout << scatterFactor << ", " << SF_Compt(13, q, nullptr);
            std::cout << ", " << q << ", " << MomentTransf(p.energy, angle, nullptr) << "\n";
        }
        */
        if (rmsd < 0.001)
            return true;
        return false;
    }

    void testSafeValueAddWorkerDose(Result<T>* res, T addValue, std::size_t N, std::size_t idx = 0) const
    {
        for (std::size_t i = 0; i < N; ++i) {
            safeValueAdd(res->dose[idx], addValue, res->locks[idx].dose);
        }
    }
    void testSafeValueAddWorkerEvents(Result<T>* res, std::uint32_t addValue, std::size_t N, std::size_t idx = 0) const
    {
        for (std::size_t i = 0; i < N; ++i) {
            safeValueAdd(res->nEvents[idx], addValue, res->locks[idx].nEvents);
        }
    }

    bool testSafeValueAdd()
    {

        const std::size_t N = 1E6;

        Result<T> res(1);

        const T addValue { 1 };
        const std::uint64_t nThreads = std::thread::hardware_concurrency();
        const std::uint64_t nJobs = std::max(nThreads, static_cast<std::uint64_t>(2));
        std::vector<std::thread> jobs;
        jobs.reserve(nJobs);
        for (std::size_t i = 0; i < nJobs; ++i) {
            jobs.emplace_back(&Test<T>::testSafeValueAddWorkerDose, this, &res, addValue, N, 0);
        }
        for (auto& job : jobs) {
            job.join();
        }

        const T expected = N * addValue * nJobs;
        const T value = res.dose[0];
        std::cout << "Test safe atomic add:\n Added " << addValue << " " << N << " times in " << nJobs;
        std::cout << " threads.\nResult is expected to be " << expected << " and is ";
        std::cout << value << "\n";

        jobs.clear();
        jobs.reserve(nJobs);
        const std::uint32_t addValueEvents { 1 };
        for (std::size_t i = 0; i < nJobs; ++i) {
            jobs.emplace_back(&Test<T>::testSafeValueAddWorkerEvents, this, &res, addValueEvents, N, 0);
        }
        for (auto& job : jobs) {
            job.join();
        }

        const std::uint32_t expectedEvents = N * addValueEvents * nJobs;
        const std::uint32_t valueEvents = res.nEvents[0];
        std::cout << "Test safe atomic add:\n Added " << addValueEvents << " " << N << " times in " << nJobs;
        std::cout << " threads.\nResult is expected to be " << expectedEvents << " and is ";
        std::cout << valueEvents << "\n";

        return isEqual(expected, value) && isEqual(expectedEvents, valueEvents);
    }
};

template <typename T>
bool testTransport()
{
    Test<T> t;

    T energy = 54;
    bool success = t.testCompton(energy, 13);
    success = success && t.testSafeValueAdd();
    return success;
}

int main()
{
    bool success_d = testTransport<double>();
    bool success_f = testTransport<float>();
    if (success_d && success_f)
        return 1;
    return 0;
}