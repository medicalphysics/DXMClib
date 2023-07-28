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

    auto material_opt = dxmc::Material<T>::byZ(Z);
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

    auto material_opt = dxmc::Material<T, Nshells>::byZ(Z);
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
template <typename T, int Nshells>
bool testPhotoelectricEffectIA(const std::map<std::size_t, T>& massFractions, const T energy = 50.0, bool print = false)
{
    bool success = true;
    constexpr std::size_t N = 1E6;

    const auto material = dxmc::Material<T, Nshells>::byWeight(massFractions).value();

    struct ShellStatistics {
        T prob = 0;
        T energy = 0;
        T bindingEnergy = 0;
        T weight = 0;
        T s_prob = 0;
        std::size_t count = 0;
    };

    std::vector<ShellStatistics> shell_stat;
    for (const auto& [Z, w] : massFractions) {
        const auto& atom = dxmc::AtomHandler<T>::Atom(Z);
        const auto pe = dxmc::interpolate(atom.photoel, energy);
        for (const auto& [Si, shell] : atom.shells) {
            const auto spe = dxmc::interpolate<T, false, false>(shell.photoel, energy);

            ShellStatistics s;
            s.prob = w * spe / pe;
            s.energy = shell.energyOfPhotonsPerInitVacancy;
            s.weight = shell.numberOfPhotonsPerInitVacancy;
            s.bindingEnergy = shell.bindingEnergy;
            shell_stat.push_back(s);
        }
    }
    // normalise probs:
    {
        const auto p_sum = std::reduce(shell_stat.cbegin(), shell_stat.cend(), T { 0 }, [](const T sum, const auto& shell) { return sum + shell.prob; });
        std::for_each(shell_stat.begin(), shell_stat.end(), [p_sum](auto& s) { s.prob /= p_sum; });
    }

    constexpr std::array<T, 3> dir = { 0, 0, 1 };
    dxmc::Particle<T> particle;
    particle.pos = { 0, 0, 0 };
    const auto photoAtt = material.attenuationPhotoeletric(energy);
    dxmc::RandomState state;
    const auto tstart = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < N; ++i) {
        particle.dir = dir;
        particle.energy = energy;
        particle.weight = 1;
        const auto ei = dxmc::interactions::photoelectricEffectIA<T>(photoAtt, particle, material, state);
        auto el = std::min_element(shell_stat.begin(), shell_stat.end(), [&particle](const auto left, const auto right) {
            return std::abs(particle.energy - left.energy) < std::abs(particle.energy - right.energy);
        });
        ++(el->count);
    }
    const auto tend = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart);

    std::for_each(shell_stat.begin(), shell_stat.end(), [N](auto& s) { s.s_prob = static_cast<T>(s.count) / N; });

    for (const auto& sh : shell_stat) {
        if (sh.energy > dxmc::MIN_ENERGY<T>())
            success = success && std::abs(sh.prob - sh.s_prob) < T { 0.1 };
    }

    auto s_string = success ? "SUCCESS" : "FAILURE";

    std::cout << s_string;
    std::cout << " PE with floating point size of " << sizeof(T) << " and ";
    std::cout << "[time:" << time.count() << "ms] ";
    std::cout << std::endl;
    if (print) {
        std::cout << "BindingEnergy, FL_energy, ana_prob , sample_prob\n";
        for (const auto& sh : shell_stat) {
            std::cout << sh.bindingEnergy << ", ";
            std::cout << sh.energy << ", ";
            std::cout << sh.prob << ", ";
            std::cout << sh.s_prob << std::endl;
        }
    }
    return success;
}

template <typename T>
bool testPhotoIAtemplate(const T energy = 50.0, bool print = false)
{
    bool success = true;
    std::map<std::size_t, T> mf1 {
        { 82, 1.0f },
        { 79, 1.0f },
        { 74, 5.0f }
    };
    success = success && testPhotoelectricEffectIA<T, 13>(mf1, energy, print);

    const auto& mf2 = dxmc::NISTMaterials<T>::Composition("Bone, Compact (ICRU)");
    success = success && testPhotoelectricEffectIA<T, 13>(mf2, energy, print);
    return success;
}

int main()
{
    // printShell();

    std::cout << "Testing interactions" << std::endl;
    bool success = true;

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

    success = success && testPhotoIAtemplate<double>(150.0, false);
    success = success && testPhotoIAtemplate<float>(150.0, false);

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
