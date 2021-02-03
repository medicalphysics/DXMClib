

#include "dxmc/attenuationinterpolator.h"
#include "dxmc/attenuationlut.h"
#include "xraylib.h"

#include <cassert>
#include <iostream>
#include <vector>

using namespace dxmc;

constexpr double ERRF = 1E-4;

template <Floating T>
bool isEqual(T a, T b)
{
    return std::abs(a - b) < ERRF;
}

template <typename T, typename U=float>
std::vector<U> normalize(const std::vector<T>& v)
{
    const U sum = std::reduce(v.cbegin(), v.cend());
    std::vector<U> r(v.size(), 0);
    std::transform(v.cbegin(), v.cend(), r.begin(), [=](auto val) { return val / sum; });
    return r;
}

template <typename T>
void testComptonScatterFactor(const std::vector<Material>& mats, const T minEnergy = 1, const T maxEnergy = 150)
{

    AttenuationLut<T> lut;
    lut.generate(mats, minEnergy, maxEnergy);

    const T min = 0;
    const T max = maxEnergy / 12.4;

    std::vector<T> q(200);

    std::cout << "Compton scatter factor test\n";
    std::cout << "q ";
    for (const auto& mat : mats) {
        std::cout << ", " << mat.prettyName() << " xraylib, " << mat.prettyName() << " dxmc";
    }
    std::cout << std::endl;
    for (std::size_t i = 0; i < q.size(); ++i) {
        q[i] = min + (max - min) / (q.size() - 1) * i;
        std::cout << q[i] << ", ";
        for (std::size_t j = 0; j < mats.size(); ++j) {
            std::cout << mats[j].getComptonNormalizedScatterFactor(q[i]) << ", ";
            std::cout << lut.comptonScatterFactor(j, q[i]) << ", ";
            assert(std::abs(mats[j].getComptonNormalizedScatterFactor(q[i]) - lut.comptonScatterFactor(j, q[i])) < 0.1);
        }
        std::cout << std::endl;
    }
}
template <typename T>
void testFormFactor(const std::vector<Material>& mats, const T minEnergy = 1, const T maxEnergy = 30)
{

    AttenuationLut<T> lut;
    lut.generate(mats, minEnergy, maxEnergy);

    const T min = 0;
    const T max = maxEnergy / 12.4;

    std::vector<T> q(50);
    for (std::size_t i = 0; i < q.size(); ++i) {
        q[i] = min + (max - min) / (q.size() - 1) * i;
    }
    std::cout << "Compton scatter factor test\n";
    std::cout << "q ";
    for (const auto& mat : mats) {
        std::cout << ", " << mat.prettyName() << " xraylib, " << mat.prettyName() << " dxmc";
    }
    std::cout << std::endl;
    std::vector<std::vector<T>> res;
    RandomState state;
    for (std::size_t j = 0; j < mats.size(); ++j) {
        std::vector<std::size_t> hist(q.size(), 0);
        for (std::size_t i = 0; i < 1E6; ++i) {
            const T q_val = lut.momentumTransferFromFormFactor(j, maxEnergy, state);
            const std::size_t idx = (q_val - min) * (q.size() - 1) / (max - min);
            ++hist[idx];
        }
        auto hist_n = normalize<std::size_t, float>(hist);

        std::vector<T> hist_a(q.size(), 0);
        for (std::size_t i = 0; i < q.size(); ++i) {

            hist_a[i] = mats[j].getRayleightFormFactorSquared(q[i]);
        }
        auto hist_an = normalize<float, float>(hist_a);
        res.push_back(hist_an);
        res.push_back(hist_n);
    }
    for (std::size_t i = 0; i < q.size(); ++i) {
        std::cout << q[i] << ", ";
        for (std::size_t j = 0; j < mats.size() * 2; ++j) {
            std::cout << res[j][i] << ", ";
        }
        std::cout << std::endl;
    }
}
template <typename T>
void testAttenuation()
{

    std::vector<Material> mats;
    mats.emplace_back((6));
    mats.emplace_back((8));
    mats.emplace_back((13));
    mats.emplace_back((82));
    mats.emplace_back(("H2O"));
    mats.back().setStandardDensity(1.0);

    //testComptonScatterFactor<T>(mats);
    testFormFactor<T>(mats);
}

template <Floating T>
bool testInterpolation()
{
    const T minEnergy = 1;
    const T maxEnergy = 150;

    std::vector<Material> materials;
    materials.emplace_back(("H2O"));
    materials.emplace_back(("H2O"));
    // materials.emplace_back((13));
    //materials.emplace_back((82));

    materials[0].setStandardDensity(1.0);
    materials[1].setStandardDensity(2.0);

    std::vector<T> dens;
    std::transform(materials.cbegin(), materials.cend(), std::back_inserter(dens), [](const auto& m) { return static_cast<T>(m.standardDensity()); });
    std::vector<std::uint8_t> mats(dens.size());
    std::iota(mats.begin(), mats.end(), 0);

    AttenuationLutInterpolator<T> attLut(materials, dens.cbegin(), dens.cend(), mats.cbegin(), T { 1 }, T { 150 });

    std::vector<T> energies;
    const std::size_t resolution = maxEnergy - minEnergy;
    for (std::size_t i = 0; i < resolution; ++i) {
        const T e = minEnergy + (i * (maxEnergy - minEnergy)) / (resolution - 1);
        energies.push_back(e);
    }

    for (const auto& mat : materials) {
        auto binding = mat.getBindingEnergies(minEnergy);
        for (const auto e : binding) {
            energies.push_back(e);
        }
    }
    std::sort(energies.begin(), energies.end());

    const auto material = 0;
    std::cout << "Energy, Total, Photo, Incoher, Coher, Total dxmc, Photo dxmc, Incoher dxmc, Coher dxmc, max total inv, " << materials[material].standardDensity() << "\n ";
    for (const auto e : energies) {
        std::cout << e << ", ";
        std::array<T, 3> xlib {
            static_cast<T>(CS_Photo_CP(materials[material].name().c_str(), e, nullptr)),
            static_cast<T>(CS_Compt_CP(materials[material].name().c_str(), e, nullptr)),
            static_cast<T>(CS_Rayl_CP(materials[material].name().c_str(), e, nullptr))
        };
        std::cout << std::reduce(xlib.cbegin(), xlib.cend()) << ", ";
        for (auto t : xlib) {
            std::cout << t << ", ";
        }
        const auto dx = attLut(material, e);
        const auto d_tot = std::reduce(dx.cbegin(), dx.cend());
        std::cout << d_tot << ", ";
        for (auto t : dx) {
            std::cout << t << ", ";
        }
        for (int i = 0; i < 3; ++i) {
            const auto valid = (dx[i] - xlib[i]) / xlib[i] * 1000;
            if (valid > 0.6)
                return false;
        }
        std::cout << attLut.maxAttenuationInverse(e) << ", ";
        std::cout << "\n";
    }

    return true;
}

int main(int argc, char* argv[])
{
    testAttenuation<float>();

    //  assert(testInterpolation<float>());

    return EXIT_FAILURE;
}
