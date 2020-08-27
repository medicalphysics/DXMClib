

#include "dxmc/source.h"
#include "dxmc/transport.h"
#include "dxmc/tube.h"
#include "dxmc/world.h"

#include <cassert>
#include <chrono>
#include <execution>
#include <fstream>
#include <iostream>
#include <numeric>

constexpr double ERRF = 1e-4;
constexpr std::size_t histPerExposure = 1e6;
constexpr std::size_t nExposures = 16;

constexpr double PI = 3.14159265359;

// energy weighs pair for spectre
/*RQR-8
W/Al
100 kVp
11 deg anode angle
0% Ripple
Al filter thickness: 2.708 mm
Mean Energy: 50.6 keV
HVL: 3.950 mm Al
QVL: 9.840 mm Al
*/
const std::vector<double> TG195_100KV_raw({ 16.25, 1.423E-04, 16.75, 2.157E-04, 17.25, 3.102E-04, 17.75, 4.324E-04, 18.25, 5.840E-04, 18.75, 7.644E-04, 19.25, 9.784E-04, 19.75, 1.222E-03, 20.25, 1.491E-03, 20.75, 1.803E-03, 21.25, 2.129E-03, 21.75, 2.490E-03, 22.25, 2.863E-03, 22.75, 3.263E-03, 23.25, 3.658E-03, 23.75, 4.093E-03, 24.25, 4.504E-03, 24.75, 4.912E-03, 25.25, 5.347E-03, 25.75, 5.769E-03, 26.25, 6.168E-03, 26.75, 6.582E-03, 27.25, 6.965E-03, 27.75, 7.360E-03, 28.25, 7.710E-03, 28.75, 8.067E-03, 29.25, 8.368E-03, 29.75, 8.671E-03, 30.25, 8.975E-03, 30.75, 9.213E-03, 31.25, 9.476E-03, 31.75, 9.694E-03, 32.25, 9.903E-03, 32.75, 1.009E-02, 33.25, 1.025E-02, 33.75, 1.040E-02, 34.25, 1.053E-02, 34.75, 1.063E-02, 35.25, 1.073E-02, 35.75, 1.081E-02, 36.25, 1.087E-02, 36.75, 1.092E-02, 37.25, 1.096E-02, 37.75, 1.099E-02, 38.25, 1.100E-02, 38.75, 1.100E-02, 39.25, 1.099E-02, 39.75, 1.098E-02, 40.25, 1.095E-02, 40.75, 1.091E-02, 41.25, 1.086E-02, 41.75, 1.081E-02, 42.25, 1.076E-02, 42.75, 1.069E-02, 43.25, 1.063E-02, 43.75, 1.055E-02, 44.25, 1.048E-02, 44.75, 1.039E-02, 45.25, 1.031E-02, 45.75, 1.022E-02, 46.25, 1.012E-02, 46.75, 1.003E-02, 47.25, 9.933E-03, 47.75, 9.828E-03, 48.25, 9.732E-03, 48.75, 9.628E-03, 49.25, 9.516E-03, 49.75, 9.412E-03, 50.25, 9.302E-03, 50.75, 9.193E-03, 51.25, 9.084E-03, 51.75, 8.970E-03, 52.25, 8.862E-03, 52.75, 8.749E-03, 53.25, 8.637E-03, 53.75, 8.526E-03, 54.25, 8.409E-03, 54.75, 8.300E-03, 55.25, 8.185E-03, 55.75, 8.072E-03, 56.25, 7.959E-03, 56.75, 7.847E-03, 57.25, 7.737E-03, 57.75, 2.568E-02, 58.25, 7.513E-03, 58.75, 7.405E-03, 59.25, 3.920E-02, 59.75, 7.181E-03, 60.25, 7.071E-03, 60.75, 6.962E-03, 61.25, 6.854E-03, 61.75, 6.746E-03, 62.25, 6.640E-03, 62.75, 6.530E-03, 63.25, 6.425E-03, 63.75, 6.321E-03, 64.25, 6.214E-03, 64.75, 6.107E-03, 65.25, 6.006E-03, 65.75, 5.901E-03, 66.25, 5.797E-03, 66.75, 1.673E-02, 67.25, 5.592E-03, 67.75, 5.491E-03, 68.25, 5.390E-03, 68.75, 8.223E-03, 69.25, 5.055E-03, 69.75, 4.296E-03, 70.25, 4.236E-03, 70.75, 4.171E-03, 71.25, 4.110E-03, 71.75, 4.048E-03, 72.25, 3.982E-03, 72.75, 3.919E-03, 73.25, 3.852E-03, 73.75, 3.787E-03, 74.25, 3.719E-03, 74.75, 3.654E-03, 75.25, 3.585E-03, 75.75, 3.516E-03, 76.25, 3.449E-03, 76.75, 3.379E-03, 77.25, 3.308E-03, 77.75, 3.240E-03, 78.25, 3.169E-03, 78.75, 3.098E-03, 79.25, 3.026E-03, 79.75, 2.954E-03, 80.25, 2.882E-03, 80.75, 2.809E-03, 81.25, 2.736E-03, 81.75, 2.665E-03, 82.25, 2.592E-03, 82.75, 2.519E-03, 83.25, 2.445E-03, 83.75, 2.370E-03, 84.25, 2.296E-03, 84.75, 2.222E-03, 85.25, 2.148E-03, 85.75, 2.073E-03, 86.25, 1.999E-03, 86.75, 1.925E-03, 87.25, 1.850E-03, 87.75, 1.776E-03, 88.25, 1.700E-03, 88.75, 1.625E-03, 89.25, 1.550E-03, 89.75, 1.476E-03, 90.25, 1.400E-03, 90.75, 1.326E-03, 91.25, 1.251E-03, 91.75, 1.177E-03, 92.25, 1.101E-03, 92.75, 1.027E-03, 93.25, 9.529E-04, 93.75, 8.781E-04, 94.25, 8.041E-04, 94.75, 7.302E-04, 95.25, 6.559E-04, 95.75, 5.823E-04, 96.25, 5.089E-04, 96.75, 4.353E-04, 97.25, 3.623E-04, 97.75, 2.892E-04, 98.25, 2.166E-04, 98.75, 1.441E-04, 99.25, 7.193E-05, 99.75, 5.990E-06 });

/*RQR-9
W/Al
120 kVp
11 deg anode angle
0% Ripple
Al filter thickness: 2.861 mm
Mean Energy: 56.4 keV
HVL: 5.010 mm Al
*/
const std::vector<double> TG195_120KV_raw({ 16.75, 1.107E-04, 17.25, 1.625E-04, 17.75, 2.308E-04, 18.25, 3.172E-04, 18.75, 4.220E-04, 19.25, 5.486E-04, 19.75, 6.956E-04, 20.25, 8.610E-04, 20.75, 1.056E-03, 21.25, 1.264E-03, 21.75, 1.499E-03, 22.25, 1.748E-03, 22.75, 2.019E-03, 23.25, 2.293E-03, 23.75, 2.601E-03, 24.25, 2.900E-03, 24.75, 3.203E-03, 25.25, 3.531E-03, 25.75, 3.858E-03, 26.25, 4.176E-03, 26.75, 4.511E-03, 27.25, 4.830E-03, 27.75, 5.163E-03, 28.25, 5.469E-03, 28.75, 5.786E-03, 29.25, 6.065E-03, 29.75, 6.349E-03, 30.25, 6.638E-03, 30.75, 6.879E-03, 31.25, 7.143E-03, 31.75, 7.372E-03, 32.25, 7.597E-03, 32.75, 7.804E-03, 33.25, 7.994E-03, 33.75, 8.171E-03, 34.25, 8.339E-03, 34.75, 8.483E-03, 35.25, 8.622E-03, 35.75, 8.745E-03, 36.25, 8.849E-03, 36.75, 8.949E-03, 37.25, 9.031E-03, 37.75, 9.109E-03, 38.25, 9.170E-03, 38.75, 9.219E-03, 39.25, 9.264E-03, 39.75, 9.297E-03, 40.25, 9.319E-03, 40.75, 9.332E-03, 41.25, 9.333E-03, 41.75, 9.332E-03, 42.25, 9.327E-03, 42.75, 9.307E-03, 43.25, 9.292E-03, 43.75, 9.259E-03, 44.25, 9.229E-03, 44.75, 9.187E-03, 45.25, 9.149E-03, 45.75, 9.101E-03, 46.25, 9.044E-03, 46.75, 8.996E-03, 47.25, 8.937E-03, 47.75, 8.871E-03, 48.25, 8.813E-03, 48.75, 8.747E-03, 49.25, 8.672E-03, 49.75, 8.605E-03, 50.25, 8.530E-03, 50.75, 8.456E-03, 51.25, 8.381E-03, 51.75, 8.300E-03, 52.25, 8.226E-03, 52.75, 8.145E-03, 53.25, 8.065E-03, 53.75, 7.985E-03, 54.25, 7.899E-03, 54.75, 7.820E-03, 55.25, 7.736E-03, 55.75, 7.652E-03, 56.25, 7.568E-03, 56.75, 7.486E-03, 57.25, 7.403E-03, 57.75, 3.335E-02, 58.25, 7.236E-03, 58.75, 7.155E-03, 59.25, 5.339E-02, 59.75, 6.986E-03, 60.25, 6.903E-03, 60.75, 6.821E-03, 61.25, 6.739E-03, 61.75, 6.658E-03, 62.25, 6.578E-03, 62.75, 6.494E-03, 63.25, 6.415E-03, 63.75, 6.338E-03, 64.25, 6.256E-03, 64.75, 6.175E-03, 65.25, 6.100E-03, 65.75, 6.021E-03, 66.25, 5.942E-03, 66.75, 2.242E-02, 67.25, 5.788E-03, 67.75, 5.712E-03, 68.25, 5.637E-03, 68.75, 9.988E-03, 69.25, 5.257E-03, 69.75, 4.045E-03, 70.25, 4.019E-03, 70.75, 3.988E-03, 71.25, 3.960E-03, 71.75, 3.932E-03, 72.25, 3.900E-03, 72.75, 3.871E-03, 73.25, 3.838E-03, 73.75, 3.808E-03, 74.25, 3.774E-03, 74.75, 3.743E-03, 75.25, 3.709E-03, 75.75, 3.674E-03, 76.25, 3.641E-03, 76.75, 3.606E-03, 77.25, 3.570E-03, 77.75, 3.537E-03, 78.25, 3.500E-03, 78.75, 3.463E-03, 79.25, 3.426E-03, 79.75, 3.389E-03, 80.25, 3.351E-03, 80.75, 3.313E-03, 81.25, 3.274E-03, 81.75, 3.238E-03, 82.25, 3.200E-03, 82.75, 3.160E-03, 83.25, 3.121E-03, 83.75, 3.079E-03, 84.25, 3.039E-03, 84.75, 3.000E-03, 85.25, 2.959E-03, 85.75, 2.919E-03, 86.25, 2.878E-03, 86.75, 2.838E-03, 87.25, 2.797E-03, 87.75, 2.756E-03, 88.25, 2.712E-03, 88.75, 2.671E-03, 89.25, 2.629E-03, 89.75, 2.588E-03, 90.25, 2.544E-03, 90.75, 2.502E-03, 91.25, 2.460E-03, 91.75, 2.418E-03, 92.25, 2.374E-03, 92.75, 2.331E-03, 93.25, 2.289E-03, 93.75, 2.244E-03, 94.25, 2.202E-03, 94.75, 2.159E-03, 95.25, 2.115E-03, 95.75, 2.072E-03, 96.25, 2.029E-03, 96.75, 1.984E-03, 97.25, 1.941E-03, 97.75, 1.896E-03, 98.25, 1.853E-03, 98.75, 1.809E-03, 99.25, 1.765E-03, 99.75, 1.722E-03, 100.25, 1.677E-03, 100.75, 1.634E-03, 101.25, 1.589E-03, 101.75, 1.546E-03, 102.25, 1.501E-03, 102.75, 1.458E-03, 103.25, 1.414E-03, 103.75, 1.370E-03, 104.25, 1.326E-03, 104.75, 1.282E-03, 105.25, 1.238E-03, 105.75, 1.195E-03, 106.25, 1.151E-03, 106.75, 1.107E-03, 107.25, 1.063E-03, 107.75, 1.019E-03, 108.25, 9.761E-04, 108.75, 9.323E-04, 109.25, 8.893E-04, 109.75, 8.456E-04, 110.25, 8.027E-04, 110.75, 7.592E-04, 111.25, 7.158E-04, 111.75, 6.731E-04, 112.25, 6.300E-04, 112.75, 5.874E-04, 113.25, 5.445E-04, 113.75, 5.017E-04, 114.25, 4.594E-04, 114.75, 4.168E-04, 115.25, 3.747E-04, 115.75, 3.324E-04, 116.25, 2.903E-04, 116.75, 2.485E-04, 117.25, 2.067E-04, 117.75, 1.650E-04, 118.25, 1.236E-04, 118.75, 8.222E-05, 119.25, 4.102E-05, 119.75, 3.417E-06 });

std::pair<std::vector<double>, std::vector<double>> TG195_specter(const std::vector<double>& raw)
{
    std::pair<std::vector<double>, std::vector<double>> s;
    s.first.resize(raw.size() / 2);
    s.second.resize(raw.size() / 2);

    for (std::size_t i = 0; i < s.first.size(); ++i) {
        s.first[i] = raw[i * 2] - 0.25;
        //s.first[i] = raw[i * 2];
        s.second[i] = raw[i * 2 + 1];
    }
    return s;
}

std::pair<std::vector<double>, std::vector<double>> TG195_120KV()
{
    return TG195_specter(TG195_120KV_raw);
}

std::pair<std::vector<double>, std::vector<double>> TG195_100KV()
{
    return TG195_specter(TG195_100KV_raw);
}

bool test120Specter()
{

    std::vector<double> tg_energy(TG195_120KV_raw.size() / 2), tg_weight(TG195_120KV_raw.size() / 2);
    for (std::size_t i = 0; i < TG195_120KV_raw.size() / 2; ++i) {
        tg_energy[i] = TG195_120KV_raw[i * 2];
        tg_weight[i] = TG195_120KV_raw[i * 2 + 1];
    }

    const double tg_weightSum = std::reduce(tg_weight.cbegin(), tg_weight.cend());
    const double tg_meanEnergy = std::transform_reduce(tg_energy.cbegin(), tg_energy.cend(), tg_weight.cbegin(), 0.0, std::plus<>(), [=](double e, double w) { return w * e; });

    const auto s = TG195_120KV();
    auto energy = s.first;
    auto weight = s.second;

    SpecterDistribution dist(weight, energy);
    const double weightSum = std::reduce(weight.cbegin(), weight.cend());

    std::vector<double> sampl_specter(energy.size(), 0.0);
    const std::size_t n_samples = 1e7;
    double meanEnergy = 0.0;
    for (std::size_t i = 0; i < n_samples; ++i) {
        const auto j = dist.sampleIndex();
        sampl_specter[j] += 1;
        if (!(j < sampl_specter.size()))
            bool test = 0;
        meanEnergy += dist.sampleValue();
    }
    meanEnergy /= n_samples;
    const auto sampl_specter_sum = std::reduce(sampl_specter.cbegin(), sampl_specter.cend());
    const double sampl_meanEnergy = std::transform_reduce(tg_energy.cbegin(), tg_energy.cend(), sampl_specter.cbegin(), 0.0, std::plus<>(), [=](double e, double w) { return w * e / sampl_specter_sum; });

    std::cout << "Test specter TG195  120 kV W/Al\n";
    std::cout << "Sum weights: TG195 specter: " << tg_weightSum << " dxmc: " << weightSum << std::endl;
    std::cout << "Mean energy [keV]: TG195 specter: " << tg_meanEnergy << " dxmc sampling: " << sampl_meanEnergy << " dxmc: " << meanEnergy << std::endl;
    std::cout << std::endl;
    return true;
}

bool inline isEqual(double a, double b)
{
    return std::abs(a - b) < ERRF;
}

bool betw(std::size_t v, std::size_t min, std::size_t max)
{
    if (v >= min && v < max)
        return true;
    return false;
}
bool betw(double v, double min, double max)
{
    if (v >= min && v < max)
        return true;
    return false;
}
World generateTG195Case2World()
{

    Material air("Air, Dry (near sea level)");
    air.setStandardDensity(0.001205);
    Material soft("H62.9539171935344C12.9077870263354N1.16702581276482O22.7840718642933Na0.026328553360443P0.0390933975009805S0.0566470278101205Cl0.0341543557411274K0.0309747686593447");
    soft.setStandardDensity(1.03);

    std::array<double, 3> spacing = { 5, 5, 5 }; // mm
    std::array<std::size_t, 3> dim = { 78, 78, 400 };
    const auto size = dim[0] * dim[1] * dim[2];
    auto dens = std::make_shared<std::vector<double>>(size, air.standardDensity());
    auto idx = std::make_shared<std::vector<unsigned char>>(size, 0);
    // fill first slice with soft tissue
    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                if (z < 40) {
                    const auto i = x + y * dim[0] + z * dim[0] * dim[1];
                    dens->data()[i] = soft.standardDensity();
                    idx->data()[i] = static_cast<unsigned char>(1);

                    //center boxes
                    if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 5, 11))
                        idx->data()[i] = static_cast<unsigned char>(9 + 1);
                    if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 11, 17))
                        idx->data()[i] = static_cast<unsigned char>(8 + 1);
                    if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 17, 23))
                        idx->data()[i] = static_cast<unsigned char>(3 + 1);
                    if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 23, 29))
                        idx->data()[i] = static_cast<unsigned char>(7 + 1);
                    if (betw(x, 36, 42) && betw(y, 36, 42) && betw(z, 29, 35))
                        idx->data()[i] = static_cast<unsigned char>(6 + 1);

                    //periphery y
                    if (betw(x, 36, 42) && betw(y, 6, 12) && betw(z, 17, 23))
                        idx->data()[i] = static_cast<unsigned char>(1 + 1);
                    if (betw(x, 36, 42) && betw(y, 66, 72) && betw(z, 17, 23))
                        idx->data()[i] = static_cast<unsigned char>(5 + 1);
                    if (betw(x, 6, 12) && betw(y, 36, 42) && betw(z, 17, 23))
                        idx->data()[i] = static_cast<unsigned char>(2 + 1);
                    if (betw(x, 66, 72) && betw(y, 36, 42) && betw(z, 17, 23))
                        idx->data()[i] = static_cast<unsigned char>(4 + 1);
                }
            }

    World w;
    w.setDimensions(dim);
    w.setSpacing(spacing);
    w.setDensityArray(dens);
    w.addMaterialToMap(air);
    for (std::size_t i = 0; i < 10; ++i)
        w.addMaterialToMap(soft);
    w.setMaterialIndexArray(idx);
    return w;
}

bool TG195Case2AbsorbedEnergyMono()
{
    std::cout << "TG195 Case 2\n";
    std::cout << "Monochromatic source of 56.4 keV:\n";
    auto w = generateTG195Case2World();
    IsotropicSource src;
    std::vector<double> s({ 1.0 }), e({ 56.4 });
    w.setAttenuationLutMaxEnergy(*std::max_element(e.cbegin(), e.cend()));
    w.setAttenuationLutResolution(0.5);
    w.validate();
    assert(w.isValid());

    src.setSpecter(s, e);
    auto extent = w.matrixExtent();
    src.setPosition(0, 0, 1750 + extent[4]);
    const double halfAng = std::atan(390 * 0.5 / 1800.0);
    src.setCollimationAngles(halfAng * 2.0, halfAng * 2.0);
    std::array<double, 6> cosines = { -1, 0, 0, 0, 1, 0 };
    src.setDirectionCosines(cosines);
    src.setHistoriesPerExposure(histPerExposure);
    src.setTotalExposures(nExposures);
    src.validate();
    auto res = transport::run(w, &src, nullptr, false);
    auto& dose = res.dose;

    const double total_hist = static_cast<double>(src.totalExposures() * src.historiesPerExposure());
    double total_ev = 1000 * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [](double d, unsigned char i) -> double { return i > 0 ? d : 0.0; });
    total_ev /= total_hist;
    std::array<double, 9> subvol_ev;
    for (std::size_t i = 0; i < 9; ++i) {
        subvol_ev[i] = 1000 * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m) -> double { return m == i + 2 ? d : 0.0; });
        subvol_ev[i] /= total_hist;
    }

    const double sim_ev = 33171.4;
    std::array<double, 9> sim_subvol = { 27.01, 27.00, 36.67, 27.01, 27.01, 72.86, 53.35, 23.83, 14.60 };

    std::cout << "Total body dxmc: " << total_ev << ", TG195: " << sim_ev << ", difference aboslute: " << total_ev - sim_ev << ", difference [%]: " << (total_ev - sim_ev) / sim_ev * 100 << "\n";
    for (std::size_t i = 0; i < subvol_ev.size(); ++i)
        std::cout << "VOI " << i + 1 << " dxmc: " << subvol_ev[i] << ", TG195: " << sim_subvol[i] << ", difference absolute: " << subvol_ev[i] - sim_subvol[i] << ", difference [%]: " << (subvol_ev[i] - sim_subvol[i]) / sim_subvol[i] * 100 << "\n";
    std::cout << "\n";
    return true;
}

bool TG195Case2AbsorbedEnergy120()
{
    std::cout << "TG195 Case 2\n";
    std::cout << "Specter source of 120 kV:\n";

    auto w = generateTG195Case2World();
    IsotropicSource src;

    const auto specter = TG195_120KV();

    const auto max_energy = *(std::max_element(specter.first.cbegin(), specter.first.cend()));
    w.setAttenuationLutMaxEnergy(max_energy);
    //w.setAttenuationLutResolution(0.4);
    w.validate();
    assert(w.isValid());
    src.setSpecter(specter.second, specter.first);
    auto extent = w.matrixExtent();
    src.setPosition(0, 0, 1750 + extent[4]);
    const double halfAng = std::atan(390 * 0.5 / 1800.0);
    src.setCollimationAngles(halfAng * 2.0, halfAng * 2.0);
    std::array<double, 6> cosines = { -1, 0, 0, 0, 1, 0 };
    src.setDirectionCosines(cosines);
    src.setHistoriesPerExposure(histPerExposure);
    src.setTotalExposures(nExposures);
    src.validate();
    auto res = transport::run(w, &src, nullptr, false);
    auto& dose = res.dose;

    const double total_hist = static_cast<double>(src.totalExposures() * src.historiesPerExposure());
    double total_ev = 1000.0 * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [](double d, unsigned char i) -> double { return i > 0 ? d : 0.0; });
    total_ev /= total_hist;
    std::array<double, 9> subvol_ev;
    for (std::size_t i = 0; i < 9; ++i) {
        subvol_ev[i] = 1000 * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m) -> double { return m == i + 2 ? d : 0.0; });
        subvol_ev[i] /= total_hist;
    }

    const double sim_ev = 33125.98;
    std::array<double, 9> sim_subvol = { 24.97, 24.95, 33.52, 24.96, 24.97, 72.70, 49.99, 21.73, 13.48 };

    std::cout << "Total body dxmc: " << total_ev << ", TG195: " << sim_ev << ", difference aboslute: " << total_ev - sim_ev << ", difference [%]: " << (total_ev - sim_ev) / sim_ev * 100 << "\n";
    for (std::size_t i = 0; i < subvol_ev.size(); ++i)
        std::cout << "VOI " << i + 1 << " dxmc: " << subvol_ev[i] << ", TG195: " << sim_subvol[i] << ", difference absolute: " << subvol_ev[i] - sim_subvol[i] << ", difference [%]: " << (subvol_ev[i] - sim_subvol[i]) / sim_subvol[i] * 100 << "\n";
    std::cout << "\n";
    return true;
}

std::vector<std::size_t> circleIndices(const double center_x, const double center_y, const std::array<std::size_t, 3>& dim, const std::array<double, 3>& spacing, const double radii = 1.0)
{
    const double xo = dim[0] * spacing[0] * 0.5;
    const double yo = dim[1] * spacing[1] * 0.5;
    std::vector<std::size_t> ind;
    const double radii2 = radii * radii;
    for (std::size_t yi = 0; yi < dim[1]; ++yi)
        for (std::size_t xi = 0; xi < dim[0]; ++xi) {
            const double x = xi * spacing[0] - center_x - xo;
            const double y = yi * spacing[1] - center_y - yo;
            if ((x * x + y * y) <= radii2) {
                const auto i = xi + yi * dim[0];
                ind.push_back(i);
            }
        }
    return ind;
}

World generateTG195Case4World1()
{

    const std::array<std::size_t, 3> dim = { 323, 323, 600 };
    const std::array<double, 3> spacing = { 1, 1, 5 };
    const auto size = std::accumulate(dim.cbegin(), dim.cend(), (std::size_t)1, std::multiplies<std::size_t>());

    Material air("Air, Dry (near sea level)");
    air.setStandardDensity(0.001205);
    Material pmma("H53.2813989847746C33.3715774096566O13.3470236055689");
    pmma.setStandardDensity(1.19);

    auto mat = std::make_shared<std::vector<unsigned char>>(size, static_cast<unsigned char>(0));
    auto dens = std::make_shared<std::vector<double>>(size, air.standardDensity());
    //generate cylindar
    auto circ_ind = circleIndices(.0, .0, dim, spacing, 160.0);
    for (std::size_t z = 0; z < dim[2]; ++z)
        for (const auto ind : circ_ind) {
            const auto i = z * dim[0] * dim[1] + ind;
            (mat->data())[i] = static_cast<unsigned char>(1);
            (dens->data())[i] = pmma.standardDensity();
            if (betw((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5, -5, 5))
                mat->data()[i] = static_cast<unsigned char>(2);
            if (betw((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5, -15, -5))
                mat->data()[i] = static_cast<unsigned char>(3);
            if (betw((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5, -25, -15))
                mat->data()[i] = static_cast<unsigned char>(4);
            if (betw((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5, -35, -25))
                mat->data()[i] = static_cast<unsigned char>(5);
        }

    World w;
    w.setSpacing(spacing);
    w.setDimensions(dim);
    w.setDensityArray(dens);
    w.setMaterialIndexArray(mat);
    w.addMaterialToMap(air);
    w.addMaterialToMap(pmma);
    for (int i = 0; i < 4; ++i)
        w.addMaterialToMap(pmma);
    w.validate();
    return w;
}

World generateTG195Case4World2(bool force_interactions = false)
{

    const std::array<std::size_t, 3> dim = { 323, 323, 60 };
    const std::array<double, 3> spacing = { 1, 1, 50 };
    const auto size = std::accumulate(dim.cbegin(), dim.cend(), (std::size_t)1, std::multiplies<>());

    Material air("Air, Dry (near sea level)");
    air.setStandardDensity(0.001205);
    Material pmma("H53.2813989847746C33.3715774096566O13.3470236055689");
    pmma.setStandardDensity(1.19);

    auto mat = std::make_shared<std::vector<unsigned char>>(size, static_cast<unsigned char>(0));
    auto dens = std::make_shared<std::vector<double>>(size, air.standardDensity());
    auto meas = std::make_shared<std::vector<std::uint8_t>>(size, 0);
    //generate cylindar
    auto circ_ind = circleIndices(.0, .0, dim, spacing, 320. * 0.5);
    auto pmma2_ind = circleIndices(-150.0, .0, dim, spacing, 5.0);
    auto pmma1_ind = circleIndices(.0, .0, dim, spacing, 5.0);
    for (std::size_t z = 0; z < dim[2]; ++z) {
        for (const auto i : circ_ind) {
            const auto ind = z * dim[0] * dim[1] + i;
            mat->data()[ind] = static_cast<unsigned char>(1);
            dens->data()[ind] = pmma.standardDensity();
        }
        if (std::abs((z + 0.5) * spacing[2] - spacing[2] * dim[2] * 0.5) < 50) {
            const auto offset = z * dim[0] * dim[1];
            for (const auto i : pmma1_ind) {
                mat->data()[offset + i] = static_cast<unsigned char>(2);
                if (force_interactions)
                    meas->data()[offset + i] = static_cast<std::uint8_t>(1);
            }
            for (const auto i : pmma2_ind) {
                mat->data()[offset + i] = static_cast<unsigned char>(3);
                if (force_interactions)
                    meas->data()[offset + i] = static_cast<std::uint8_t>(1);
            }
        }
    }

    World w;
    w.setSpacing(spacing);
    w.setDimensions(dim);
    w.setDensityArray(dens);
    w.setMaterialIndexArray(mat);
    w.setMeasurementMapArray(meas);
    w.addMaterialToMap(air);
    w.addMaterialToMap(pmma);
    for (int i = 0; i < 2; ++i)
        w.addMaterialToMap(pmma);
    w.validate();
    return w;
}

bool TG195Case41AbsorbedEnergy()
{
    IsotropicSource src;
    src.setPosition(-600.0, 0., 0.);
    std::array<double, 6> cos = { 0, 1, 0, 0, 0, 1 };
    src.setDirectionCosines(cos);

    double direction[3];
    vectormath::cross(cos.data(), direction);

    std::cout << "TG195 Case 4.1:\nMonoenergetic specter of 56.4 kev\n";
    std::cout << "Collimation: 10 mm:\n";
    std::vector<double> s(1, 1.0), e(1, 56.4);
    src.setSpecter(s, e);
    src.setHistoriesPerExposure(histPerExposure);
    src.setTotalExposures(nExposures);
    src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(5. / 600.) * 2.);
    auto w = generateTG195Case4World1();
    w.setAttenuationLutMaxEnergy(*std::max_element(e.cbegin(), e.cend()));
    w.validate();

    auto res = transport::run(w, &src, nullptr, false);

    auto& dose = res.dose;
    const double total_hist = static_cast<double>(src.totalExposures() * src.historiesPerExposure());
    std::array<double, 4> voi_ev;
    for (std::size_t i = 0; i < 4; ++i) {
        voi_ev[i] = 1000.0 / total_hist * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m) -> double { return m == i + 2 ? d : 0.0; });
    }
    std::array<double, 4> sim_ev = { 11592.27, 2576.72, 1766.85, 1330.53 };
    for (std::size_t i = 0; i < voi_ev.size(); ++i)
        std::cout << "VOI " << i + 1 << " dxmc: " << voi_ev[i] << ", TG195: " << sim_ev[i] << ", difference absolute: " << voi_ev[i] - sim_ev[i] << ", difference [%]: " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";

    std::cout << "Collimation: 80 mm:\n";
    src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(40. / 600.) * 2.);

    res = transport::run(w, &src, nullptr, false);
    dose = res.dose;

    for (std::size_t i = 0; i < 4; ++i) {
        voi_ev[i] = 1000 / total_hist * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m) -> double { return m == i + 2 ? d : 0.0; });
    }

    sim_ev = { 3380.39, 3332.64, 3176.44, 2559.58 };

    for (std::size_t i = 0; i < voi_ev.size(); ++i)
        std::cout << "VOI " << i + 1 << " dxmc: " << voi_ev[i] << ", TG195: " << sim_ev[i] << ", difference absolute: " << voi_ev[i] - sim_ev[i] << ", difference [%]: " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";

    std::cout << "\nTG195 Case 4.1:\nSpecter of 120 kV W/Al\n";
    std::cout << "Collimation: 10 mm:\n";

    const auto specter = TG195_120KV();
    src.setSpecter(specter.second, specter.first);
    src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(5. / 600.) * 2.);

    w.setAttenuationLutMaxEnergy(*std::max_element(specter.first.cbegin(), specter.first.cend()));
    w.validate();
    res = transport::run(w, &src, nullptr, false);

    dose = res.dose;
    for (std::size_t i = 0; i < 4; ++i) {
        voi_ev[i] = 1000.0 / total_hist * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m) -> double { return m == i + 2 ? d : 0.0; });
    }
    sim_ev = { 13137.02, 2585.47, 1706.86, 1250.61 };

    for (std::size_t i = 0; i < voi_ev.size(); ++i)
        std::cout << "VOI " << i + 1 << " dxmc: " << voi_ev[i] << ", TG195: " << sim_ev[i] << ", difference absolute: " << voi_ev[i] - sim_ev[i] << ", difference [%]: " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";
    std::cout << "Collimation: 80 mm:\n";
    src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(40. / 600.) * 2.);

    res = transport::run(w, &src, nullptr, false);

    dose = res.dose;
    for (std::size_t i = 0; i < 4; ++i) {
        voi_ev[i] = 1000 / total_hist * std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](double d, unsigned char m) -> double { return m == i + 2 ? d : 0.0; });
    }
    sim_ev = { 3586.59, 3537.84, 3378.99, 2672.21 };

    for (std::size_t i = 0; i < voi_ev.size(); ++i)
        std::cout << "VOI " << i + 1 << " dxmc: " << voi_ev[i] << ", TG195: " << sim_ev[i] << ", difference absolute: " << voi_ev[i] - sim_ev[i] << ", difference [%]: " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";

    std::cout << "\n";

    return true;
}

struct TG19542Res {
    double doseC = 0.0;
    double doseP = 0.0;
    std::size_t tallyC = 0;
    std::size_t tallyP = 0;
    std::size_t nHist = 0;
    std::size_t angle = 0;
    double doseC_TG195 = 0.0;
    double doseP_TG195 = 0.0;
    double time_seconds = 0.0;

    template <typename D, typename M, typename T>
    TG19542Res(D dBeg, D dEnd, M mBeg, T tBeg, T tEnd, std::size_t nHistories, std::size_t ang, std::chrono::duration<double> time)
    {
        nHist = nHistories;
        angle = ang;
        time_seconds = time.count();
        doseC = std::transform_reduce(std::execution::par_unseq, dBeg, dEnd, mBeg, 0.0, std::plus<>(), [](auto d, auto ind) -> double { return ind == 2 ? d : 0; }) * 1000 / nHistories;
        doseP = std::transform_reduce(std::execution::par_unseq, dBeg, dEnd, mBeg, 0.0, std::plus<>(), [](auto d, auto ind) -> double { return ind == 3 ? d : 0; }) * 1000 / nHistories;
        tallyC = std::transform_reduce(std::execution::par_unseq, tBeg, tEnd, mBeg, (std::size_t)0, std::plus<>(), [](auto d, auto ind) -> std::size_t { return ind == 2 ? d : 0; });
        tallyP = std::transform_reduce(std::execution::par_unseq, tBeg, tEnd, mBeg, (std::size_t)0, std::plus<>(), [](auto d, auto ind) -> std::size_t { return ind == 3 ? d : 0; });
    }
};

bool TG195Case42AbsorbedEnergy()
{

    //TG195 results
    const std::vector<double> res42_tg_mcn({ 12.17, 12.11, 12.16, 12.14, 12.09, 12.10, 12.14, 12.12, 12.15, 12.15, 12.15, 12.16, 12.13, 12.15, 12.15, 12.14, 12.13, 12.15, 12.16, 12.15, 12.16, 12.16, 12.16, 12.16, 12.14, 12.16, 12.16, 12.14, 12.15, 12.10, 12.15, 12.12, 12.12, 12.14, 12.13, 12.14 });
    const std::vector<double> res42_tg_mpn({ 101.29, 99.80, 96.11, 89.97, 81.37, 70.54, 56.28, 38.13, 22.35, 12.17, 6.52, 3.64, 2.18, 1.41, 0.98, 0.75, 0.62, 0.55, 0.52, 0.55, 0.62, 0.75, 0.99, 1.41, 2.17, 3.63, 6.51, 12.20, 22.29, 38.15, 56.29, 70.64, 81.42, 89.82, 96.17, 99.92 });
    const std::vector<double> res42_tg_mcw({ 11.63, 11.63, 11.62, 11.62, 11.63, 11.60, 11.62, 11.62, 11.59, 11.63, 11.61, 11.61, 11.61, 11.62, 11.61, 11.62, 11.62, 11.63, 11.61, 11.63, 11.63, 11.62, 11.62, 11.62, 11.62, 11.62, 11.64, 11.62, 11.63, 11.65, 11.63, 11.63, 11.62, 11.62, 11.60, 11.61 });
    const std::vector<double> res42_tg_mpw({ 99.67, 98.30, 94.51, 88.49, 80.11, 69.26, 55.12, 37.24, 21.75, 11.77, 6.27, 3.46, 2.07, 1.34, 0.95, 0.72, 0.59, 0.52, 0.50, 0.52, 0.59, 0.72, 0.94, 1.35, 2.07, 3.47, 6.24, 11.80, 21.73, 37.28, 55.19, 69.26, 80.04, 88.40, 94.64, 98.33 });
    const std::vector<double> res42_tg_pcn({ 11.35, 11.40, 11.39, 11.40, 11.40, 11.37, 11.40, 11.38, 11.39, 11.41, 11.40, 11.38, 11.37, 11.36, 11.39, 11.38, 11.37, 11.37, 11.38, 11.36, 11.38, 11.38, 11.38, 11.37, 11.38, 11.38, 11.39, 11.38, 11.37, 11.38, 11.37, 11.40, 11.38, 11.37, 11.40, 11.40 });
    const std::vector<double> res42_tg_ppn({ 116.79, 115.31, 110.53, 103.15, 92.99, 79.72, 62.57, 40.92, 22.97, 12.17, 6.42, 3.57, 2.16, 1.42, 1.01, 0.78, 0.65, 0.58, 0.56, 0.58, 0.65, 0.78, 1.01, 1.41, 2.16, 3.58, 6.43, 12.17, 23.04, 40.93, 62.45, 79.76, 93.07, 103.27, 110.60, 115.27 });
    const std::vector<double> res42_tg_pcw({ 10.88, 10.92, 10.88, 10.90, 10.87, 10.90, 10.90, 10.88, 10.88, 10.89, 10.90, 10.88, 10.90, 10.89, 10.89, 10.89, 10.88, 10.89, 10.89, 10.89, 10.89, 10.90, 10.89, 10.89, 10.89, 10.90, 10.89, 10.91, 10.89, 10.88, 10.88, 10.88, 10.89, 10.89, 10.89, 10.90 });
    const std::vector<double> res42_tg_ppw({ 115.34, 113.76, 109.17, 101.71, 91.56, 78.39, 61.39, 40.09, 22.47, 11.78, 6.15, 3.42, 2.06, 1.35, 0.96, 0.74, 0.62, 0.55, 0.53, 0.55, 0.62, 0.74, 0.96, 1.35, 2.06, 3.42, 6.16, 11.79, 22.46, 40.14, 61.43, 78.33, 91.48, 101.61, 109.10, 113.84 });

    std::cout << "TG195 Case 4.2:\n";
    std::vector<double> specter_mono_weights(1, 1.0);
    std::vector<double> specter_mono_energies(1, 56.4);
    auto specter = TG195_120KV();
    std::vector<double> specter_poly_weights = specter.second;
    std::vector<double> specter_poly_energies = specter.first;

    IsotropicSource src;
    src.setHistoriesPerExposure(histPerExposure);
    src.setTotalExposures(nExposures);

    auto w = generateTG195Case4World2();
    w.setAttenuationLutMaxEnergy(120.0);
    w.validate();

    double rot_axis[3] = { 0, 0, 1 };

    std::vector<TG19542Res> results_mono_10;
    std::vector<TG19542Res> results_mono_80;
    std::vector<TG19542Res> results_poly_10;
    std::vector<TG19542Res> results_poly_80;
    //simulate 36 projections
    for (std::size_t i = 0; i < 360; i += 10) {
        std::cout << '\r' << "Prosessing angle (0-360) [" << static_cast<int>(i / 360.0 * 100.0) << " %] "
                  << " current angle: " << i;
        const auto nHistories = src.historiesPerExposure() * src.totalExposures();
        const double angle = i * PI / 180.0;
        std::array<double, 3> pos { -600, 0, 0 };
        std::array<double, 6> cos = { 0, 1, 0, 0, 0, 1 };
        vectormath::rotate(pos.data(), rot_axis, angle);
        src.setPosition(pos);
        vectormath::rotate(cos.data(), rot_axis, angle);
        vectormath::rotate(&cos[3], rot_axis, angle);
        src.setDirectionCosines(cos);

        //test_direction
        double dir[3];
        vectormath::cross(cos.data(), dir);

        src.setSpecter(specter_mono_weights, specter_mono_energies);
        //10mm collimation mono
        {
            src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(5. / 600.) * 2.);
            auto start = std::chrono::high_resolution_clock::now();
            auto res = transport::run(w, &src, nullptr, false);
            auto end = std::chrono::high_resolution_clock::now();
            results_mono_10.push_back(TG19542Res(res.dose.cbegin(), res.dose.cend(), w.materialIndexBuffer(), res.nEvents.cbegin(), res.nEvents.cend(), nHistories, i, end - start));
        }
        //80mm collimation mono
        {
            src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(40. / 600.) * 2.);
            auto start = std::chrono::high_resolution_clock::now();
            auto res = transport::run(w, &src, nullptr, false);
            auto end = std::chrono::high_resolution_clock::now();
            results_mono_80.push_back(TG19542Res(res.dose.cbegin(), res.dose.cend(), w.materialIndexBuffer(), res.nEvents.cbegin(), res.nEvents.cend(), nHistories, i, end - start));
        }

        src.setSpecter(specter_poly_weights, specter_poly_energies);
        //10mm collimation poly
        {
            src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(5. / 600.) * 2.);
            auto start = std::chrono::high_resolution_clock::now();
            auto res = transport::run(w, &src, nullptr, false);
            auto end = std::chrono::high_resolution_clock::now();
            results_poly_10.push_back(TG19542Res(res.dose.cbegin(), res.dose.cend(), w.materialIndexBuffer(), res.nEvents.cbegin(), res.nEvents.cend(), nHistories, i, end - start));
        }
        //80mm collimation poly
        {
            src.setCollimationAngles(std::atan(160. / 600.) * 2., std::atan(40. / 600.) * 2.);
            auto start = std::chrono::high_resolution_clock::now();
            auto res = transport::run(w, &src, nullptr, false);
            auto end = std::chrono::high_resolution_clock::now();
            results_poly_80.push_back(TG19542Res(res.dose.cbegin(), res.dose.cend(), w.materialIndexBuffer(), res.nEvents.cbegin(), res.nEvents.cend(), nHistories, i, end - start));
        }
    }
    std::cout << '\r';
    for (std::size_t i = 0; i < 36; ++i) {
        results_mono_10[i].doseC_TG195 = res42_tg_mcn[i];
        results_mono_10[i].doseP_TG195 = res42_tg_mpn[i];
        results_mono_80[i].doseC_TG195 = res42_tg_mcw[i];
        results_mono_80[i].doseP_TG195 = res42_tg_mpw[i];
        results_poly_10[i].doseC_TG195 = res42_tg_pcn[i];
        results_poly_10[i].doseP_TG195 = res42_tg_ppn[i];
        results_poly_80[i].doseC_TG195 = res42_tg_pcw[i];
        results_poly_80[i].doseP_TG195 = res42_tg_ppw[i];
    }
    std::cout << "Position, Angle [deg], Collimation, Specter, Dose [eV/hist], N events, TG195 result [ev/hist], Difference, Difference [%], Total time[s], Total histories \n";

    for (const auto& r : results_mono_10) {
        std::cout << "Center,     " << r.angle << ", 10, mono, " << r.doseC << ", " << r.tallyC << ", " << r.doseC_TG195 << ", " << r.doseC - r.doseC_TG195 << ", " << (r.doseC / r.doseC_TG195 - 1) * 100 << ", " << r.time_seconds << ", " << r.nHist << "\n";
        std::cout << "Pheriphery, " << r.angle << ", 10, mono, " << r.doseP << ", " << r.tallyP << ", " << r.doseP_TG195 << ", " << r.doseP - r.doseP_TG195 << ", " << (r.doseP / r.doseP_TG195 - 1) * 100 << ", " << r.time_seconds << ", " << r.nHist << "\n";
        "\n";
    }
    for (const auto& r : results_mono_80) {
        std::cout << "Center,     " << r.angle << ", 80, mono, " << r.doseC << ", " << r.tallyC << ", " << r.doseC_TG195 << ", " << r.doseC - r.doseC_TG195 << ", " << (r.doseC / r.doseC_TG195 - 1) * 100 << ", " << r.time_seconds << ", " << r.nHist << "\n";
        std::cout << "Pheriphery, " << r.angle << ", 80, mono, " << r.doseP << ", " << r.tallyP << ", " << r.doseP_TG195 << ", " << r.doseP - r.doseP_TG195 << ", " << (r.doseP / r.doseP_TG195 - 1) * 100 << ", " << r.time_seconds << ", " << r.nHist << "\n";
    }
    for (const auto& r : results_poly_10) {
        std::cout << "Center,     " << r.angle << ", 10, poly, " << r.doseC << ", " << r.tallyC << ", " << r.doseC_TG195 << ", " << r.doseC - r.doseC_TG195 << ", " << (r.doseC / r.doseC_TG195 - 1) * 100 << ", " << r.time_seconds << ", " << r.nHist << "\n";
        std::cout << "Pheriphery, " << r.angle << ", 10, poly, " << r.doseP << ", " << r.tallyP << ", " << r.doseP_TG195 << ", " << r.doseP - r.doseP_TG195 << ", " << (r.doseP / r.doseP_TG195 - 1) * 100 << ", " << r.time_seconds << ", " << r.nHist << "\n";
    }
    for (const auto& r : results_poly_80) {
        std::cout << "Center,     " << r.angle << ", 80, poly, " << r.doseC << ", " << r.tallyC << ", " << r.doseC_TG195 << ", " << r.doseC - r.doseC_TG195 << ", " << (r.doseC / r.doseC_TG195 - 1) * 100 << ", " << r.time_seconds << ", " << r.nHist << "\n";
        std::cout << "Pheriphery, " << r.angle << ", 80, poly, " << r.doseP << ", " << r.tallyP << ", " << r.doseP_TG195 << ", " << r.doseP - r.doseP_TG195 << ", " << (r.doseP / r.doseP_TG195 - 1) * 100 << ", " << r.time_seconds << ", " << r.nHist << "\n";
    }

    return true;
}

bool testAttenuation()
{
    const double energy = 56.4;
    Material m("Tissue, Soft (ICRP)");
    m.setStandardDensity(1.3);
    std::array<double, 3> spacing = { .1, .1, 1 };
    std::array<std::size_t, 3> dim = { 1, 1, 200 };
    const auto size = std::accumulate(dim.cbegin(), dim.cend(), (std::size_t)1, std::multiplies<>());
    auto dens = std::make_shared<std::vector<double>>(size, m.standardDensity());
    auto mat = std::make_shared<std::vector<unsigned char>>(size, static_cast<unsigned char>(0));
    auto meas = std::make_shared<std::vector<unsigned char>>(size, 0);
    World w;
    w.setDimensions(dim);
    w.setSpacing(spacing);
    w.setDensityArray(dens);
    w.setMaterialIndexArray(mat);
    w.setMeasurementMapArray(meas);
    w.addMaterialToMap(m);
    w.setAttenuationLutMaxEnergy(std::ceil(energy));
    w.validate();

    PencilSource pen;
    pen.setHistoriesPerExposure(histPerExposure);
    pen.setPhotonEnergy(energy);
    pen.setTotalExposures(nExposures);
    std::array<double, 6> cos = { 1, 0, 0, 0, 1, 0 };
    pen.setDirectionCosines(cos);
    pen.setPosition(0, 0, -400);

    const auto tot_hist = pen.historiesPerExposure() * pen.totalExposures();

    auto res_naive = transport::run(w, &pen, nullptr, false);
    std::fill(meas->begin(), meas->end(), 1);
    w.validate();
    auto res_force = transport::run(w, &pen, nullptr, false);

    std::vector<double> att(res_naive.dose.size());
    for (int i = 0; i < dim[2]; ++i)
        att[i] = std::exp(-(i + 1) * spacing[2] * 0.1 * m.standardDensity() * m.getTotalAttenuation(56.4));

    const auto att_max = *std::max_element(att.cbegin(), att.cend());
    const auto naive_max = *std::max_element(res_naive.dose.cbegin(), res_naive.dose.cend());
    const auto force_max = *std::max_element(res_force.dose.cbegin(), res_force.dose.cend());

    const double rms_naive = std::sqrt((1.0 / att.size()) * std::transform_reduce(res_naive.dose.cbegin(), res_naive.dose.cend(), att.cbegin(), 0.0, std::plus<>(), [=](double e, double a) -> double { return (e / naive_max - a / att_max) * (e / naive_max - a / att_max); }));
    const double rms_force = std::sqrt((1.0 / att.size()) * std::transform_reduce(res_force.dose.cbegin(), res_force.dose.cend(), att.cbegin(), 0.0, std::plus<>(), [=](double e, double a) -> double { return (e / force_max - a / att_max) * (e / force_max - a / att_max); }));

    std::cout << "Test attenuation for pencil beam in 1mm^2 tissue rod: \n";
    std::cout << "Monochromatic beam of " << energy << " kev. \n";
    std::cout << "RMS differense [%] from analytical attenuation; naive: " << rms_naive * 100.0 << ", forced: " << rms_force * 100.0 << "\n";
    /*std::cout << "Data: Position (midtpoint) [mm], naive dose, force dose, naive rel. dose, force rel.dose, analytical rel. attenuation, # Events naive, # Events force\n";
	for (int i = 0; i < dim[2]; ++i)
	{
		std::cout << i * spacing[2] + spacing[2] * .5 << ", " << res_naive.dose[i] << ", " << res_force.dose[i] << ", ";
		std::cout << res_naive.dose[i]/naive_max << ", " << res_force.dose[i]/force_max << ", ";
		std::cout << att[i] / att_max << ", " << res_naive.nEvents[i] << ", " << res_force.nEvents[i] << "\n";
	}
	*/
    if (rms_naive * 100 < 0.2 && rms_force * 100 < 0.2) {
        std::cout << "SUCCESS\n\n";
        return true;
    }
    std::cout << "FAILURE\n\n";
    return false;
}

int main(int argc, char* argv[])
{
    //test120Specter();
    auto success = true; //testAttenuation();
    success = success && TG195Case2AbsorbedEnergyMono();
    success = success && TG195Case2AbsorbedEnergy120();
    success = success && TG195Case41AbsorbedEnergy();
    success = success && TG195Case42AbsorbedEnergy();
    std::cout << "Press any key to exit";
    std::string dummy;
    std::cin >> dummy;
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}