

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

template <typename T, typename U = float>
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
    // compare scatter factor
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
    // compare form factor
    AttenuationLut<T> lut;
    lut.generate(mats, minEnergy, maxEnergy);

    const T min = 0;
    const T max = maxEnergy / 12.4;

    std::vector<T> q(50);
    for (std::size_t i = 0; i < q.size(); ++i) {
        q[i] = min + (max - min) / (q.size() - 1) * i;
    }
    std::cout << "Rayleight form factor test\n";
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
            if (idx >= hist.size())
                ++hist[idx - 1];
            else
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

template <typename T, int N = 3>
void testPhotoAttenuation(const std::vector<Material>& mats, const T minEnergy = 1, const T maxEnergy = 150)
{
    //Compare attenuation values
    std::vector<T> energy(139);

    for (std::size_t i = 0; i < energy.size(); ++i) {
        energy[i] = minEnergy + (maxEnergy - minEnergy) / (energy.size() - 1) * i;
    }

    AttenuationLut<T> lut;
    lut.generate(mats, minEnergy, maxEnergy, true);

    std::cout << "Attenuation test ";
    if (N == 0)
        std::cout << "Photoelectric\n";
    if (N == 1)
        std::cout << "Compton\n";
    if (N == 2)
        std::cout << "Rayleight\n";
    if (N > 2)
        std::cout << "Total\n";
    std::cout << "Energy";
    for (const auto& mat : mats) {
        std::cout << ", " << mat.prettyName() << " xraylib, " << mat.prettyName() << " dxmc";
    }

    std::cout << std::endl;
    for (const auto e : energy) {
        std::cout << e;
        std::size_t mIdx = 0;
        for (const auto& m : mats) {
            const T p = CS_Photo_CP(m.name().c_str(), e, nullptr);
            const T c = CS_Compt_CP(m.name().c_str(), e, nullptr);
            const T r = CS_Rayl_CP(m.name().c_str(), e, nullptr);
            const auto a = lut.photoComptRayAttenuation(mIdx, e);

            T v, mc;

            if (N == 0) {
                v = p;
            }
            if (N == 1) {
                v = c;
            }
            if (N == 2) {
                v = r;
            }
            if (N > 2) {
                v = c + p + r;
                mc = std::reduce(a.cbegin(), a.cend(), T { 0 });
            } else {
                mc = a[N];
            }
            std::cout << ", " << v;
            std::cout << ", " << mc;
            ++mIdx;
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
    testPhotoAttenuation<T, 0>(mats);
    testPhotoAttenuation<T, 1>(mats);
    testPhotoAttenuation<T, 2>(mats);
    testPhotoAttenuation<T, 3>(mats);
    testComptonScatterFactor<T>(mats);
    testFormFactor<T>(mats);
}

int main(int argc, char* argv[])
{
    testAttenuation<float>();

    //  assert(testInterpolation<float>());

    return EXIT_FAILURE;
}
