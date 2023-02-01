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

Copyright 2023 Erlend Andersen
*/

#include "dxmc/constants.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"

#include "xraylib.h"

#include <chrono>
#include <iostream>

template <typename T>
class Histogram {
public:
    Histogram(T start, T stop, std::size_t Nbins)
        : m_start(start)
        , m_stop(stop)
    {
        intensity.resize(Nbins);
        std::fill(intensity.begin(), intensity.end(), 0);
        edges.resize(Nbins + 1);
        m_step = (stop - start) / Nbins;
        for (std::size_t i = 0; i <= Nbins; ++i) {
            edges[i] = m_start + m_step * i;
        }
    }
    inline std::size_t idx(T val) const
    {
        const T frac = (val - m_start) / (m_stop - m_start);
        const auto idx_f = std::nextafter(frac * intensity.size(), T { 0 });
        const auto idx = static_cast<std::size_t>(idx_f);
        return idx;
    }

    void operator()(T value)
    {
        if (m_start < value && value < m_stop) {
            const auto i = idx(value);
            intensity[i] += 1;
        }
    }
    std::pair<std::vector<T>, std::vector<T>> plotData(bool normalize = true)
    {
        std::vector<T> x(intensity.size(), 0);
        std::vector<T> y(intensity.size(), 0);

        const auto sum = normalize ? std::reduce(intensity.cbegin(), intensity.cend(), T { 0 }) : T { 1 };
        for (std::size_t i = 0; i < intensity.size(); ++i) {
            y[i] = static_cast<T>(intensity[i]) / sum;
            x[i] = (edges[i] + edges[i + 1]) / 2;
        }
        return std::make_pair(x, y);
    }

private:
    std::vector<std::size_t> intensity;
    std::vector<T> edges;
    T m_start;
    T m_stop;
    T m_step = 0;
};

template <typename T, int type = 0>
bool testCoherent(std::size_t Z = 13, T energy = 5, bool print = false)
{
    bool success = true;
    constexpr std::size_t N = 1E6;
    Histogram hist(T { -1 }, T { 1 }, 90);

    auto material_opt = dxmc::Material2<T>::byZ(Z);
    auto material = material_opt.value();

    constexpr std::array<T, 3> dir = { 0, 0, 1 };
    dxmc::Particle<T> particle;
    particle.pos = { 0, 0, 0 };
    particle.dir = dir;
    particle.energy = energy;
    particle.weight = 1;

    dxmc::RandomState state;
    const auto tstart = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < N; ++i) {
        particle.dir = dir;
        dxmc::interactions::rayleightScatter<T, type>(particle, material, state);
        const T costheta = dxmc::vectormath::dot(particle.dir, dir);
        hist(costheta);
    }
    const auto tend = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart);

    const auto [x, y] = hist.plotData();
    std::vector<T> y_ana(x.size(), 0);
    std::vector<T> y_xraylib(x.size(), 0);
    // analytical
    if constexpr (type == 0) {
        std::transform(x.cbegin(), x.cend(), y_ana.begin(), [](const auto xv) {
            return (1 + xv * xv);
        });
        std::transform(x.cbegin(), x.cend(), y_xraylib.begin(), [](const auto xv) {
            const auto angle = std::acos(xv);
            return static_cast<T>(DCS_Thoms(angle, nullptr));
        });
    } else {
        std::transform(x.cbegin(), x.cend(), y_ana.begin(), [energy, &material](const auto xv) {
            const auto angle = std::acos(xv);
            const auto q = material.momentumTransfer(energy, angle);
            const auto formFactor = material.formFactor(q);
            return (1 + xv * xv) * formFactor * formFactor / 2;
        });
        std::transform(x.cbegin(), x.cend(), y_xraylib.begin(), [Z, energy](const auto xv) {
            const auto angle = std::acos(xv);
            return static_cast<T>(DCS_Rayl(Z, energy, angle, nullptr));
        });
    }
    const auto y_ana_sum = std::reduce(y_ana.cbegin(), y_ana.cend(), T { 0 });
    std::for_each(y_ana.begin(), y_ana.end(), [y_ana_sum](auto& yv) { yv /= y_ana_sum; });
    const auto y_xraylib_sum = std::reduce(y_xraylib.cbegin(), y_xraylib.cend(), T { 0 });
    std::for_each(y_xraylib.begin(), y_xraylib.end(), [y_xraylib_sum](auto& yv) { yv /= y_xraylib_sum; });

    for (std::size_t i = 0; i < y.size(); ++i) {
        const auto diff = std::abs(y[i] / y_ana[i] - 1);
        success = success && diff < T { 0.2 };
    }

    const auto corr = type == 0 ? "None" : "Livermore or IA";
    auto s_string = success ? "SUCCESS" : "FAILURE";

    std::cout << s_string;
    std::cout << " Rayleight with floating point size of " << sizeof(T) << " and ";
    std::cout << corr << " correction [time: " << time.count() << "ms]\n";
    if (print) {
        std::cout << "costheta, samples, analytical, xraylib\n";
        for (std::size_t i = 0; i < y.size(); ++i) {
            std::cout << x[i] << ", " << y[i] << ", " << y_ana[i] << ", " << y_xraylib[i] << std::endl;
        }
    }
    return success;
}

template <typename T, int type = 0>
bool testIncoherent(std::size_t Z = 13, T energy = 50, bool print = false)
{
    bool success = true;
    constexpr std::size_t N = 1E7;

    constexpr std::size_t sampleStep = 360;

    Histogram hist_angle(T { -1 }, T { 1 }, sampleStep);
    Histogram hist_energy(T { 0 }, T { 1 }, sampleStep);

    constexpr int Nshells = 12;

    auto material_opt = dxmc::Material2<T, Nshells>::byZ(Z);
    auto material = material_opt.value();

    auto atom = dxmc::AtomHandler<T>::Atom(Z);

    constexpr std::array<T, 3> dir = { 0, 0, 1 };
    dxmc::Particle<T> particle;
    particle.pos = { 0, 0, 0 };
    particle.dir = dir;
    particle.energy = energy;
    particle.weight = 1;

    std::size_t rejectedInteractions = 0;

    dxmc::RandomState state;
    const auto tstart = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < N; ++i) {
        particle.dir = dir;
        particle.energy = energy;
        const auto ei = dxmc::interactions::comptonScatter<T, type>(particle, material, state);

        hist_angle(dxmc::vectormath::dot(particle.dir, dir));
        hist_energy(particle.energy / energy);
    }
    const auto tend = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart);

    const auto [x_angle, y_angle] = hist_angle.plotData();
    const auto [x_energy, y_energy] = hist_energy.plotData();
    std::vector<T> y_energy_ana(x_energy.size(), 0);
    std::vector<T> y_angle_ana(x_angle.size(), 0);
    std::vector<T> y_xraylib(x_angle.size(), 0);

    if constexpr (type == 0) {
        std::transform(x_angle.cbegin(), x_angle.cend(), y_angle_ana.begin(), [energy](const auto cosangle) {
            const auto k = energy / dxmc::ELECTRON_REST_MASS<T>();
            const auto kc = k / (1 + k * (1 - cosangle));
            const auto e = kc / k;
            const auto sinangleSqr = 1 - cosangle * cosangle;
            const auto KN = e * e * (e + 1 / e - sinangleSqr);
            return KN;
        });
        std::transform(x_energy.cbegin(), x_energy.cend(), y_energy_ana.begin(), [energy](const auto e) {
            constexpr auto mec = dxmc::ELECTRON_REST_MASS<T>();
            const auto emin = mec / (mec + 2 * energy);
            if (e < emin)
                return T { 0 };
            const auto cosangle = 1 + mec / energy - mec / (energy * e);
            const auto sinangleSqr = 1 - cosangle * cosangle;
            const auto KN = (1 / e + e) * (1 - e * sinangleSqr / (1 + e * e));

            return KN;
        });
        std::transform(x_angle.cbegin(), x_angle.cend(), y_xraylib.begin(), [energy](const auto cosangle) {
            const auto theta = std::acos(cosangle);
            return static_cast<T>(DCS_KN(energy, theta, nullptr));
        });
    } else if constexpr (type > 0) {
        std::transform(x_angle.cbegin(), x_angle.cend(), y_angle_ana.begin(), [energy, &material](const auto cosangle) {
            const auto k = energy / dxmc::ELECTRON_REST_MASS<T>();
            const auto kc = k / (1 + k * (1 - cosangle));
            const auto e = kc / k;
            const auto sinangleSqr = 1 - cosangle * cosangle;
            const auto KN = e * e * (e + 1 / e - sinangleSqr);
            const auto q = material.momentumTransferCosAngle(energy, cosangle);
            return KN * material.scatterFactor(q);
        });
        std::transform(x_energy.cbegin(), x_energy.cend(), y_energy_ana.begin(), [energy, &material](const auto e) {
            constexpr auto mec = dxmc::ELECTRON_REST_MASS<T>();
            const auto emin = mec / (mec + 2 * energy);
            if (e < emin)
                return T { 0 };
            const auto cosangle = 1 + mec / energy - mec / (energy * e);
            const auto sinangleSqr = 1 - cosangle * cosangle;
            const auto KN = (1 / e + e) * (1 - e * sinangleSqr / (1 + e * e));
            const auto q = material.momentumTransferCosAngle(energy, cosangle);

            return KN * material.scatterFactor(q);
        });
        std::transform(x_angle.cbegin(), x_angle.cend(), y_xraylib.begin(), [Z, energy, &material](const auto cosangle) {
            const auto theta = std::acos(cosangle);
            const auto q = material.momentumTransferCosAngle(energy, cosangle);
            return static_cast<T>(DCS_Compt(Z, energy, theta, nullptr));
        });
    }

    auto normalize = [](std::vector<T>& v) -> void {
        const auto sum = std::reduce(v.cbegin(), v.cend(), T { 0 });
        if (sum > T { 0 })
            std::for_each(v.begin(), v.end(), [sum](auto& el) { el /= sum; });
    };

    normalize(y_angle_ana);
    normalize(y_energy_ana);
    normalize(y_xraylib);

    T max_diff = 0;
    for (std::size_t i = 0; i < y_angle.size() - 1; ++i) {
        const auto diff = std::abs(y_angle[i] / y_angle_ana[i] - 1);
        success = success && diff < T { 0.03 };
        max_diff = std::max(diff, max_diff);
    }

    std::string corr = "None";
    if (type == 1)
        corr = "Livermore";
    if (type > 1)
        corr = "IA";

    auto s_string = success ? "SUCCESS" : "FAILURE";

    std::cout << s_string;
    std::cout << " Compton with floating point size of " << sizeof(T) << " and ";
    std::cout << corr << " correction [time: " << time.count() << "ms]";
    std::cout << " Max error of " << max_diff << std::endl;
    if (print) {
        std::cout << "energyFraction, samples, analytical, angle, samples, analytical, xraylib\n";
        for (std::size_t i = 0; i < y_angle.size(); ++i) {
            std::cout << x_energy[i] << ", " << y_energy[i] << ", " << y_energy_ana[i] << ", ";
            std::cout << x_angle[i] << ", " << y_angle[i] << ", " << y_angle_ana[i] << ", " << y_xraylib[i] << std::endl;
        }
    }
    return success;
}
template <typename T>
bool testPhotoelectricEffectIA(std::size_t Z = 13, T energy = 50.0, bool print = false)
{
    bool success = true;
    constexpr std::size_t N = 1E6;
    constexpr std::size_t sampleStep = 360;

    constexpr int Nshells = 12;
    auto material_opt = dxmc::Material2<T, Nshells>::byZ(Z);
    auto material = material_opt.value();

    T min_energy = energy;
    for (const auto& s : material.shells()) {
        if (s.bindingEnergy < energy)
            min_energy = std::min(min_energy, energy - s.energyOfPhotonsPerInitVacancy);
    }

    Histogram hist_energy(min_energy - 1, energy + 1, sampleStep);

    auto atom = dxmc::AtomHandler<T>::Atom(Z);

    constexpr std::array<T, 3> dir = { 0, 0, 1 };
    dxmc::Particle<T> particle;
    particle.pos = { 0, 0, 0 };
    particle.dir = dir;
    particle.energy = energy;
    particle.weight = 1;

    dxmc::RandomState state;
    const auto tstart = std::chrono::high_resolution_clock::now();
    const T photoAtt = material.attenuationPhotoeletric(energy);

    for (std::size_t i = 0; i < N; ++i) {
        particle.dir = dir;
        particle.energy = energy;
        particle.weight = 1;
        const auto ei = dxmc::interactions::photoelectricEffectIA<T>(photoAtt, particle, material, state);
        hist_energy(ei);
    }
    const auto tend = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart);

    auto [x_energy, y_energy] = hist_energy.plotData();

    auto normalize = [](std::vector<T>& v) -> void {
        const auto sum = std::reduce(v.cbegin(), v.cend(), T { 0 });
        if (sum > T { 0 })
            std::for_each(v.begin(), v.end(), [sum](auto& el) { el /= sum; });
    };
    normalize(y_energy);

    auto s_string = success ? "SUCCESS" : "FAILURE";

    std::cout << s_string;
    std::cout << " PE with floating point size of " << sizeof(T) << " and ";
    std::cout << "[time:" << time.count() << "ms] ";
    std::cout << std::endl;
    if (print) {
        std::cout << "energy, samples\n";
        for (std::size_t i = 0; i < y_energy.size(); ++i) {
            std::cout << x_energy[i] << ", " << y_energy[i] << std::endl;
        }
    }
    return success;
}
/*
void printShell()
{
    std::cout << "Z, shell, bindingEnergy, photons per vacancy, PE10, PE30, PE50, PE70, PE100, PE150\n";

    constexpr std::array en = { 10., 30., 50., 70., 100., 150. };
    for (std::size_t Z = 1; Z < 90; ++Z) {

        const auto& atom = dxmc::AtomHandler<double>::Atom(Z);

        for (const auto& [shIdx, shell] : atom.shells) {

            if (shell.bindingEnergy > 1.0) {
                std::cout << Z << ", ";
                std::cout << shIdx << ", ";
                std::cout << shell.bindingEnergy << ", ";
                std::cout << shell.numberOfPhotonsPerInitVacancy << ", ";
                for (const auto& e : en) {
                    std::cout << dxmc::interpolate<double, false, false>(shell.photoel, e) / dxmc::interpolate(atom.photoel, e) << ", ";
                }
                std::cout << std::endl;
            }
        }
    }
}
*/
int main()
{
    // printShell();

    std::cout << "Testing interactions" << std::endl;
    bool success = true;
    success = success && testPhotoelectricEffectIA<double>(79, 100.0, true);

    success = success && testCoherent<double, 0>(13, 5.0, false);
    success = success && testCoherent<float, 0>(13, 5.0, false);
    success = success && testCoherent<double, 1>(13, 5.0, false);
    success = success && testCoherent<float, 1>(13, 5.0, false);

    success = success && testIncoherent<double, 0>(13, 50., false);
    success = success && testIncoherent<float, 0>(13, 50., false);
    success = success && testIncoherent<double, 1>(13, 50., false);
    success = success && testIncoherent<float, 1>(13, 50., false);
    success = success && testIncoherent<double, 2>(13, 50., false);
    success = success && testIncoherent<float, 2>(13, 50., false);

    success = success && testPhotoelectricEffectIA<double>(79, 100.0, false);
    success = success && testPhotoelectricEffectIA<float>(79, 100.0, false);

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
