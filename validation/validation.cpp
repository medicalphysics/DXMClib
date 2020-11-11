

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

using namespace dxmc;

constexpr double ERRF = 1e-4;
constexpr std::size_t histPerExposure = 1e7;
constexpr std::size_t nExposures = 4;

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

/*
RQR-M3
Mo/Mo
30 kVp
15 deg anode angle
0% Ripple
Mo filter thickness: 0.0386 mm
Mean Energy: 16.8 keV
HVL: 0.3431 mm Al
QVL: 0.7663 mm Al
*/
const std::vector<double> TG195_30KV_raw({ 7.25, 1.551E-04, 7.75, 4.691E-04, 8.25, 1.199E-03, 8.75, 2.405E-03, 9.25, 4.263E-03, 9.75, 6.797E-03, 10.25, 9.761E-03, 10.75, 1.314E-02, 11.25, 1.666E-02, 11.75, 2.013E-02, 12.25, 2.349E-02, 12.75, 2.666E-02, 13.25, 2.933E-02, 13.75, 3.167E-02, 14.25, 3.365E-02, 14.75, 3.534E-02, 15.25, 3.644E-02, 15.75, 3.741E-02, 16.25, 3.796E-02, 16.75, 3.823E-02, 17.25, 3.445E-01, 17.75, 3.770E-02, 18.25, 3.704E-02, 18.75, 3.639E-02, 19.25, 9.200E-02, 19.75, 2.178E-03, 20.25, 2.048E-03, 20.75, 2.043E-03, 21.25, 2.098E-03, 21.75, 2.193E-03, 22.25, 2.327E-03, 22.75, 2.471E-03, 23.25, 2.625E-03, 23.75, 2.770E-03, 24.25, 2.907E-03, 24.75, 3.000E-03, 25.25, 3.062E-03, 25.75, 3.058E-03, 26.25, 2.988E-03, 26.75, 2.823E-03, 27.25, 2.575E-03, 27.75, 2.233E-03, 28.25, 1.815E-03, 28.75, 1.290E-03, 29.25, 6.696E-04, 29.75, 4.086E-05 });

template <typename T>
std::pair<std::vector<T>, std::vector<T>> TG195_specter(const std::vector<double>& raw)
{
    std::pair<std::vector<T>, std::vector<T>> s;
    s.first.resize(raw.size() / 2);
    s.second.resize(raw.size() / 2);

    for (std::size_t i = 0; i < s.first.size(); ++i) {
        s.first[i] = raw[i * 2] - 0.25;
        //s.first[i] = raw[i * 2];
        s.second[i] = raw[i * 2 + 1];
    }
    return s;
}
template <typename T>
std::pair<std::vector<T>, std::vector<T>> TG195_120KV()
{
    return TG195_specter<T>(TG195_120KV_raw);
}
template <typename T>
std::pair<std::vector<T>, std::vector<T>> TG195_100KV()
{
    return TG195_specter<T>(TG195_100KV_raw);
}
template <typename T>
std::pair<std::vector<T>, std::vector<T>> TG195_30KV()
{
    return TG195_specter<T>(TG195_30KV_raw);
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

    const auto s = TG195_120KV<double>();
    auto energy = s.first;
    auto weight = s.second;

    SpecterDistribution<double> dist(weight, energy);
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
template <typename T>
bool inline isEqual(T a, T b)
{
    return std::abs(a - b) < ERRF;
}

template <typename T, typename U, typename W>
bool betw(T v, U min, W max)
{
    if (v >= min && v < max)
        return true;
    return false;
}

template <Floating T = double>
std::size_t indexFromPosition(const std::array<T, 3>& pos, const World<T>& world)
{
    //assumes particle is inside world
    std::size_t arraypos[3];
    const auto& wpos = world.matrixExtent();
    const auto& wdim = world.dimensions();
    const auto& wspac = world.spacing();

    for (std::size_t i = 0; i < 3; i++)
        arraypos[i] = static_cast<std::size_t>((pos[i] - wpos[i * 2]) / wspac[i]);
    const auto idx = arraypos[2] * wdim[0] * wdim[1] + arraypos[1] * wdim[0] + arraypos[0];
    return idx;
}

template <Floating T = double>
std::vector<T> getEVperHistory(const Result<T>& res, std::shared_ptr<std::vector<T>> density, const std::array<T, 3>& spacing, std::uint64_t histories)
{
    std::vector<T> ev(res.dose.size());
    const T voxelVolume = spacing[0] * spacing[1] * spacing[2] * T { 0.001 }; // mm->cm
    const T nHist = static_cast<T>(histories);

    std::transform(std::execution::par_unseq, res.dose.cbegin(), res.dose.cend(), density->cbegin(), ev.begin(), [=](auto e, auto d) -> T {
        const T voxelMass = d * voxelVolume * T { 0.001 }; // kg
        return 1000 * e * voxelMass / nHist; // returns eV/hist
    });
    return ev;
}

template <typename T>
World<T> generateTG195Case2World(bool forcedInteractions = false)
{
    //Material air("Air, Dry (near sea level)");
    Material air("C0.0150228136551869N78.439632744437O21.0780510531616Ar0.467293388746132");
    air.setStandardDensity(0.001205);
    Material soft("H62.9539171935344C12.9077870263354N1.16702581276482O22.7840718642933Na0.026328553360443P0.0390933975009805S0.0566470278101205Cl0.0341543557411274K0.0309747686593447");
    soft.setStandardDensity(1.03);

    std::array<T, 3> spacing = { 5, 5, 5 }; // mm
    std::array<std::size_t, 3> dim = { 78, 78, 360 };
    const auto size = dim[0] * dim[1] * dim[2];
    auto dens = std::make_shared<std::vector<T>>(size, static_cast<T>(air.standardDensity()));
    auto idx = std::make_shared<std::vector<std::uint8_t>>(size, 0);
    // fill first slice with soft tissue
    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const T zc = z * spacing[2] + spacing[2] / 2;
                if (betw(zc, 1550, 1550 + 200)) {
                    const T xc = x * spacing[0] + spacing[0] / 2;
                    const T yc = y * spacing[1] + spacing[1] / 2;

                    const auto i = x + y * dim[0] + z * dim[0] * dim[1];
                    dens->data()[i] = static_cast<T>(soft.standardDensity());
                    idx->data()[i] = static_cast<std::uint8_t>(1);

                    //center boxes
                    if (betw(xc, 180, 180 + 30) && betw(yc, 180, 180 + 30) && betw(zc, 1550 + 25 + 30 * 4, 1550 + 25 + 30 * 5))
                        idx->data()[i] = static_cast<std::uint8_t>(9 + 1);
                    if (betw(xc, 180, 180 + 30) && betw(yc, 180, 180 + 30) && betw(zc, 1550 + 25 + 30 * 3, 1550 + 25 + 30 * 4))
                        idx->data()[i] = static_cast<std::uint8_t>(8 + 1);
                    if (betw(xc, 180, 180 + 30) && betw(yc, 180, 180 + 30) && betw(zc, 1550 + 25 + 30 * 2, 1550 + 25 + 30 * 3))
                        idx->data()[i] = static_cast<std::uint8_t>(3 + 1);
                    if (betw(xc, 180, 180 + 30) && betw(yc, 180, 180 + 30) && betw(zc, 1550 + 25 + 30 * 1, 1550 + 25 + 30 * 2))
                        idx->data()[i] = static_cast<std::uint8_t>(7 + 1);
                    if (betw(xc, 180, 180 + 30) && betw(yc, 180, 180 + 30) && betw(zc, 1550 + 25 + 30 * 0, 1550 + 25 + 30 * 1))
                        idx->data()[i] = static_cast<std::uint8_t>(6 + 1);

                    //periphery y, x
                    if (betw(xc, 180, 180 + 30) && betw(yc, 30, 30 + 30) && betw(zc, 1550 + 25 + 30 * 2, 1550 + 25 + 30 * 3))
                        idx->data()[i] = static_cast<std::uint8_t>(1 + 1);
                    if (betw(xc, 180, 180 + 30) && betw(yc, 330, 330 + 30) && betw(zc, 1550 + 25 + 30 * 2, 1550 + 25 + 30 * 3))
                        idx->data()[i] = static_cast<std::uint8_t>(5 + 1);
                    if (betw(xc, 30, 30 + 30) && betw(yc, 180, 180 + 30) && betw(zc, 1550 + 25 + 30 * 2, 1550 + 25 + 30 * 3))
                        idx->data()[i] = static_cast<std::uint8_t>(2 + 1);
                    if (betw(xc, 330, 330 + 30) && betw(yc, 180, 180 + 30) && betw(zc, 1550 + 25 + 30 * 2, 1550 + 25 + 30 * 3))
                        idx->data()[i] = static_cast<std::uint8_t>(4 + 1);
                }
            }

    std::array<T, 3> origin = { 0, 0, 900 };

    World<T> w;
    w.setDimensions(dim);
    w.setSpacing(spacing);
    w.setOrigin(origin);

    /*for (std::size_t i = 1545; i < 1800;++i) {
        const T pz = i;
        std::array<T, 3> pos = { 0, 0, pz };
        auto ind = indexFromPosition(pos, w);
        std::cout << pz << ": " << static_cast<unsigned int>(idx->data()[ind]) << "\n";
    }

    for (int i = -390/2; i < 390/2; ++i) {
        const T px = i;
        std::array<T, 3> pos = { px, 0, 1550+100 };
        auto ind = indexFromPosition(pos, w);
        std::cout << px << ": " << static_cast<unsigned int>(idx->data()[ind]) << "\n";
    }*/

    w.setDensityArray(dens);
    w.addMaterialToMap(air);
    for (std::size_t i = 0; i < 10; ++i)
        w.addMaterialToMap(soft);
    w.setMaterialIndexArray(idx);

    if (forcedInteractions) {
        auto fmat = std::make_shared<std::vector<std::uint8_t>>(w.size(), 0);
        std::transform(
            std::execution::par_unseq, idx->cbegin(), idx->cend(), fmat->begin(), [](auto i) -> std::uint8_t { return i > 1 ? 1 : 0; });
        w.setMeasurementMapArray(fmat);
    }

    w.makeValid();
    return w;
}

template <typename T>
bool TG195Case2AbsorbedEnergy(bool specter = false, bool tomo = false, bool forcedInteractions = false)
{
    std::cout << "TG195 Case 2\n";
    if (forcedInteractions)
        std::cout << "Forcing interactions in subvolumes\n";
    auto w = generateTG195Case2World<T>();

    IsotropicSource<T> src;

    if (specter) {
        const auto specter = TG195_120KV<T>();
        src.setSpecter(specter.second, specter.first);
        std::cout << "Specter source of 120 kV specter\n";
    } else {
        std::vector<T> s({ 1.0 }), e({ 56.4 });
        src.setSpecter(s, e);
        std::cout << "Monochromatic source of 56.4 keV\n";
    }

    src.setHistoriesPerExposure(histPerExposure);
    src.setTotalExposures(nExposures);

    if (tomo) {
        std::array<T, 6> cosines = { 1, 0, 0, 0, 1, 0 };
        std::array<T, 3> rotaxis = { 1, 0, 0 };

        //calculating angles to align beam
        const T angle = DEG_TO_RAD<T>() * T { 15 };
        const T g = 1800;
        const T s = g * std::tan(angle);
        const T h2 = 390 / 2;
        const T a1 = std::atan(s / g - h2 / g);
        const T a2 = std::atan(s / g + h2 / g);
        const T a = (a2 - a1) / 2 + a1;
        const T angY = 2 * (a - a1);

        const T angX = 2 * std::atan(h2 / g);

        src.setCollimationAngles(angX, angY);
        vectormath::rotate(cosines.data(), rotaxis.data(), -a);
        vectormath::rotate(&cosines[3], rotaxis.data(), -a);
        src.setDirectionCosines(cosines);

        src.setPosition(0, -s, 0);
        auto exp = src.getExposure(0);
        std::cout << "Incident angle is 15 degrees\n";
    } else {
        src.setPosition(0, 0, 0);
        const T halfAng = std::atan((T { 390 } / 2) / T { 1800 });
        src.setCollimationAngles(halfAng * 2, halfAng * 2);
        std::array<T, 6> cosines = { 1, 0, 0, 0, 1, 0 };
        src.setDirectionCosines(cosines);
        std::cout << "Incident angle is 0 degrees\n";
    }

    src.validate();

    const auto total_hist = static_cast<T>(src.totalExposures() * src.historiesPerExposure());

    Transport<T> transport;
    auto res = transport(w, &src, nullptr, false);
    auto dose = getEVperHistory(res, w.densityArray(), w.spacing(), total_hist);

    std::array<T, 9> subvol_ev;
    std::array<std::uint64_t, 9> subvol_events;
    auto total_ev = std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), T { 0 }, std::plus<>(), [=](auto d, auto m) -> T { return m > 0 ? d : 0; });
    auto total_events = std::transform_reduce(std::execution::par_unseq, res.nEvents.cbegin(), res.nEvents.cend(), w.materialIndexArray()->begin(), 0, std::plus<>(), [=](auto d, auto m) { return m > 0 ? d : 0; });
    for (std::size_t i = 0; i < 9; ++i) {
        subvol_ev[i] = std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), T { 0 }, std::plus<>(), [=](auto d, auto m) -> T { return m == i + 2 ? d : 0; });
        subvol_events[i] = std::transform_reduce(std::execution::par_unseq, res.nEvents.cbegin(), res.nEvents.cend(), w.materialIndexArray()->begin(), 0, std::plus<>(), [=](auto d, auto m) { return m == i + 2 ? d : 0; });
    }

    T sim_ev;
    std::array<T, 9> sim_subvol;

    if (specter) {
        if (tomo) {
            sim_ev = 30923.13;
            sim_subvol = { 30.35, 23.52, 31.64, 23.52, 8.90, 70.53, 47.74, 20.31, 12.51 };
        } else {
            sim_ev = 33125.98;
            sim_subvol = { 24.97, 24.95, 33.52, 24.96, 24.97, 72.70, 49.99, 21.73, 13.48 };
        }
    } else {
        if (tomo) {
            sim_ev = 30883.83;
            sim_subvol = { 33.0807985, 25.475272, 34.62570725, 25.50542125, 9.79069025, 70.80499875, 51.0616275, 22.2764985, 13.54431025 };
        } else {
            sim_ev = 33171.4;
            sim_subvol = { 27.01, 27.00, 36.67, 27.01, 27.01, 72.86, 53.35, 23.83, 14.60 };
        }
    }
    const T simtime = std::chrono::duration_cast<std::chrono::seconds>(res.simulationTime).count();
    std::cout << "Simulation time " << simtime << " seconds";
    std::cout << " with " << simtime / total_hist << " seconds*CPU core per history\n";
    std::cout << "VOI, dxmc, dxmc nEvents, TG195, difference [eV/hist], difference [%]\n";

    std::cout << "Total body, " << total_ev << ", " << total_events << ", " << sim_ev << ", " << total_ev - sim_ev << ", " << (total_ev - sim_ev) / sim_ev * 100 << "\n";
    for (std::size_t i = 0; i < subvol_ev.size(); ++i)
        std::cout << "VOI " << i + 1 << ", " << subvol_ev[i] << ", " << subvol_events[i] << ", " << sim_subvol[i] << ", " << subvol_ev[i] - sim_subvol[i] << ", " << (subvol_ev[i] - sim_subvol[i]) / sim_subvol[i] * 100 << "\n";
    std::cout << "\n";
    return true;
}

template <typename T>
std::vector<std::size_t> circleIndices(const T center_x, const T center_y, const std::array<std::size_t, 3>& dim, const std::array<T, 3>& spacing, const T radii = 1.0)
{
    const T xo = dim[0] * spacing[0] / 2 + spacing[0] / 2;
    const T yo = dim[1] * spacing[1] / 2 + spacing[1] / 2;
    std::vector<std::size_t> ind;
    ind.reserve(dim[0] * dim[1]);
    const T radii2 = radii * radii;
    for (std::size_t yi = 0; yi < dim[1]; ++yi)
        for (std::size_t xi = 0; xi < dim[0]; ++xi) {
            const T x = xi * spacing[0] + spacing[0] / 2 - center_x - xo;
            const T y = yi * spacing[1] + spacing[1] / 2 - center_y - yo;
            if ((x * x + y * y) <= radii2) {
                const auto i = xi + yi * dim[0];
                ind.push_back(i);
            }
        }
    return ind;
}

template <typename T>
World<T> generateTG195Case3World()
{
    //std::array<T, 3> spacing = { 1, 1, 1 };
    //std::array<std::size_t, 3> dim = { 340, 300, 1320 };
    //std::array<T, 3> spacing = { 2, 2, 2 };
    //std::array<std::size_t, 3> dim = { 170, 150, 660};
    std::array<T, 3> spacing = { 1, 1, 2 };
    std::array<std::size_t, 3> dim = { 340, 300, 660 };
    World<T> w;
    w.setSpacing(spacing);
    w.setDimensions(dim);
    const auto size = w.size();

    //materials
    Material air("C0.0150228136551869N78.439632744437O21.0780510531616Ar0.467293388746132", "air");
    air.setStandardDensity(0.001205);
    Material pmma("H53.2813989847746C33.3715774096566O13.3470236055689", "pmma");
    pmma.setStandardDensity(1.19);
    Material breast("H61.9873215815672C25.2115870352038N0.812500094703561O11.959397097091P0.00826728025671175S0.00798629068764982K0.00655039238754296Ca0.00639022810261801", "breast");
    breast.setStandardDensity(0.952);
    Material breastSkin("H61.6819253067427C9.42172044694575N2.26874188115669O26.5008052432356P0.0359095410401931S0.034689042139856K0.02845211205691Ca0.027756426682265", "skin");
    breastSkin.setStandardDensity(1.09);
    Material water("H66.6220373399527O33.3779626600473", "water");
    water.setStandardDensity(1.0);

    std::vector<Material> materials({ air, water, pmma, breastSkin, breast });

    //arrays
    auto densArr = std::make_shared<std::vector<T>>(size, static_cast<T>(air.standardDensity()));
    auto matArr = std::make_shared<std::vector<std::uint8_t>>(size, 0);
    auto dens = densArr->data();
    auto mat = matArr->data();

    //making volume
    for (std::size_t z = 0; z < dim[2]; ++z) {
        const T zpos = z * spacing[2] + spacing[2] / 2 - dim[2] * spacing[2] / 2;
        for (std::size_t y = 0; y < dim[1]; ++y) {
            const T ypos = y * spacing[1] + spacing[1] / 2 - dim[1] * spacing[1] / 2;
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const T xpos = x * spacing[0] + spacing[0] / 2 - dim[0] * spacing[0] / 2;
                const auto idx = x + y * dim[0] + z * dim[0] * dim[1];

                //body
                if (betw(zpos, -110, 190) && betw(xpos, -170, 0) && betw(ypos, -150, 150)) {
                    dens[idx] = water.standardDensity();
                    mat[idx] = 1;
                }
                //pmma plates
                if (betw(zpos, 13, 15) || betw(zpos, 65, 67)) {
                    if (betw(xpos, 0, 140) && betw(ypos, -130, 130)) {
                        dens[idx] = pmma.standardDensity();
                        mat[idx] = 2;
                    }
                }

                //breast skin
                if (betw(zpos, 15, 65) && betw(xpos, 0, 100) && betw(ypos, -100, 100)) {
                    constexpr T radius_skin_sqr = T { 100 } * T { 100 };
                    const T skin_dist_sqr = xpos * xpos + ypos * ypos;
                    if (skin_dist_sqr < radius_skin_sqr) {
                        dens[idx] = breastSkin.standardDensity();
                        mat[idx] = 3;
                    }

                    //breast tissue
                    if (betw(zpos, 17, 63)) {
                        constexpr T radius_tissue_sqr = T { 98 } * T { 98 };
                        const T tissue_dist_sqr = xpos * xpos + ypos * ypos;
                        if (tissue_dist_sqr < radius_tissue_sqr) {
                            dens[idx] = breast.standardDensity();
                            mat[idx] = 4;
                        }

                        //making measurements VOIS, index start at 5
                        if (betw(zpos, 35, 45)) {
                            //central vois
                            //voi 1
                            if (betw(xpos, 40, 60) && betw(ypos, -60, -40)) {
                                mat[idx] = 4 + 1;
                            }
                            //voi 2
                            if (betw(xpos, 10, 30) && betw(ypos, -10, 10)) {
                                mat[idx] = 4 + 2;
                            }
                            //voi 3
                            if (betw(xpos, 40, 60) && betw(ypos, -10, 10)) {
                                mat[idx] = 4 + 3;
                            }
                            //voi 4
                            if (betw(xpos, 70, 90) && betw(ypos, -10, 10)) {
                                mat[idx] = 4 + 4;
                            }
                            //voi 5
                            if (betw(xpos, 40, 60) && betw(ypos, 40, 60)) {
                                mat[idx] = 4 + 5;
                            }
                        }
                        //voi 6
                        if (betw(zpos, 20, 30) && betw(xpos, 40, 60) && betw(ypos, -10, 10)) {
                            mat[idx] = 4 + 6;
                        }
                        //voi 7
                        if (betw(zpos, 50, 60) && betw(xpos, 40, 60) && betw(ypos, -10, 10)) {
                            mat[idx] = 4 + 7;
                        }
                    }
                }
            }
        }
    }
    w.setDensityArray(densArr);
    w.setMaterialIndexArray(matArr);

    auto meas = std::make_shared<std::vector<std::uint8_t>>(w.size());
    std::transform(std::execution::par_unseq, matArr->cbegin(), matArr->cend(), meas->begin(), [](auto m) -> std::uint8_t { return m > 4; });
    w.setMeasurementMapArray(meas);

    w.addMaterialToMap(air);
    w.addMaterialToMap(water);
    w.addMaterialToMap(pmma);
    w.addMaterialToMap(breastSkin);
    w.addMaterialToMap(breast);
    for (int i = 0; i < 7; ++i)
        w.addMaterialToMap(breast);

    w.makeValid();
    return w;
}

template <typename T>
bool TG195Case3AbsorbedEnergy(bool specter = false, bool tomo = false)
{
    std::cout << "TG195 Case 3\n";

    auto w = generateTG195Case3World<T>();

    IsotropicSource<T> src;

    if (specter) {
        const auto specter = TG195_30KV<T>();
        src.setSpecter(specter.second, specter.first);
        std::cout << "Specter source of 30 kV specter\n";
    } else {
        std::vector<T> s({ 1.0 }), e({ 16.8 });
        src.setSpecter(s, e);
        std::cout << "Monochromatic source of 16.8 keV\n";
    }

    src.setHistoriesPerExposure(histPerExposure);
    src.setTotalExposures(nExposures);

    if (tomo) {
        std::array<T, 6> cosines = { 1, 0, 0, 0, -1, 0 };
        std::array<T, 3> rotaxis = { 1, 0, 0 };

        //calculating angles to align beam
        const T angle = DEG_TO_RAD<T>() * T { 15 };
        const T g = 660;
        const T s = g * std::cos(angle);
        const T k = g * std::sin(angle);
        const T h2 = 260 / 2;
        const T a2 = std::atan((h2 + k) / s);
        const T a1 = std::atan((k - h2) / s);
        const T a = a2 / 2;

        const T angY = (a2 - a1) / 2;
        const T angX = std::atan(140 / (g * std::cos(angle)));

        src.setCollimationAngles(0, angX, -angY, angY);
        vectormath::rotate(cosines.data(), rotaxis.data(), -a);
        vectormath::rotate(&cosines[3], rotaxis.data(), -a);
        src.setDirectionCosines(cosines);

        src.setPosition(0, std::sin(angle) * g, std::cos(angle) * g);
        auto exp = src.getExposure(0);
        std::cout << "Incident angle is 15 degrees\n";
    } else {
        src.setPosition(0, 0, 660);
        const T halfAngX = std::atan(T { 140 } / T { 660 });
        const T halfAngY = std::atan(T { 130 } / T { 660 });

        src.setCollimationAngles(0, halfAngX, -halfAngY, halfAngY);
        std::array<T, 6> cosines = { 1, 0, 0, 0, -1, 0 };
        src.setDirectionCosines(cosines);
        std::cout << "Incident angle is 0 degrees\n";
    }

    src.validate();

    const auto total_hist = static_cast<T>(src.totalExposures() * src.historiesPerExposure());

    Transport<T> transport;
    transport.setBindingEnergyCorrection(true);
    transport.setLivermoreComptonModel(true);
    auto res = transport(w, &src, nullptr, false);
    auto dose = getEVperHistory(res, w.densityArray(), w.spacing(), total_hist);

    std::array<T, 7> subvol_ev;
    std::array<std::uint64_t, 7> subvol_events;
    auto total_ev = std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), T { 0 }, std::plus<>(), [=](auto d, auto m) -> T { return m > 3 ? d : 0; });
    auto total_events = std::transform_reduce(std::execution::par_unseq, res.nEvents.cbegin(), res.nEvents.cend(), w.materialIndexArray()->begin(), 0, std::plus<>(), [=](auto d, auto m) { return m > 3 ? d : 0; });
    for (std::size_t i = 5; i < 12; ++i) {
        subvol_ev[i - 5] = std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), T { 0 }, std::plus<>(), [=](auto d, auto m) -> T { return m == i ? d : 0; });
        subvol_events[i - 5] = std::transform_reduce(std::execution::par_unseq, res.nEvents.cbegin(), res.nEvents.cend(), w.materialIndexArray()->begin(), 0, std::plus<>(), [=](auto d, auto m) { return m == i ? d : 0; });
    }

    T sim_ev;
    std::array<T, 7> sim_subvol;

    if (specter) {
        if (tomo) {
            sim_ev = 4188.833;
            sim_subvol = { 14.390, 15.825, 15.972, 15.445, 17.171, 5.619, 49.022 };
        } else {
            sim_ev = 4293.433;
            sim_subvol = { 16.502, 16.658, 16.814, 16.249, 16.521, 6.041, 50.041 };
        }
    } else {
        if (tomo) {
            sim_ev = 4577.743;
            sim_subvol = { 15.217, 16.836, 16.943, 16.431, 18.370, 5.043, 54.974 };
        } else {
            sim_ev = 4697.333;
            sim_subvol = { 17.692, 18.070, 17.865, 17.262, 17.768, 5.417, 56.017 };
        }
    }
    const T simtime = std::chrono::duration_cast<std::chrono::seconds>(res.simulationTime).count();
    std::cout << "Simulation time " << simtime << " seconds";
    std::cout << " with " << simtime / total_hist << " seconds*CPU core per history\n";
    std::cout << "VOI, dxmc, dxmc nEvents, TG195, difference [eV/hist], difference [%]\n";

    std::cout << "Total body, " << total_ev << ", " << total_events << ", " << sim_ev << ", " << total_ev - sim_ev << ", " << (total_ev - sim_ev) / sim_ev * 100 << "\n";
    for (std::size_t i = 0; i < subvol_ev.size(); ++i)
        std::cout << "VOI " << i + 1 << ", " << subvol_ev[i] << ", " << subvol_events[i] << ", " << sim_subvol[i] << ", " << subvol_ev[i] - sim_subvol[i] << ", " << (subvol_ev[i] - sim_subvol[i]) / sim_subvol[i] * 100 << "\n";
    std::cout << "\n";
    return true;
}

template <typename T>
World<T> generateTG195Case4World1()
{

    const std::array<std::size_t, 3> dim = { 300, 300, 600 };
    const std::array<T, 3> spacing = { 4, 4, 5 };
    const auto size = std::accumulate(dim.cbegin(), dim.cend(), (std::size_t)1, std::multiplies<std::size_t>());

    //Material air("Air, Dry (near sea level)");
    Material air("C0.0150228136551869N78.439632744437O21.0780510531616Ar0.467293388746132");
    air.setStandardDensity(0.001205);
    Material pmma("H53.2813989847746C33.3715774096566O13.3470236055689");
    pmma.setStandardDensity(1.19);

    auto mat = std::make_shared<std::vector<std::uint8_t>>(size, static_cast<std::uint8_t>(0));
    auto dens = std::make_shared<std::vector<T>>(size, static_cast<T>(air.standardDensity()));
    //generate cylindar
    auto circ_ind = circleIndices(T { 0 }, T { 0 }, dim, spacing, T { 160 });
    for (std::size_t z = 0; z < dim[2]; ++z) {
        const auto zpos = z * spacing[2] - (dim[2] * spacing[2]) / 2;
        for (const auto ind : circ_ind) {
            const auto i = z * dim[0] * dim[1] + ind;
            (mat->data())[i] = static_cast<std::uint8_t>(1);
            (dens->data())[i] = pmma.standardDensity();
            if (betw(zpos + spacing[2] / 2, -5, 5))
                mat->data()[i] = static_cast<std::uint8_t>(2);
            if (betw(zpos + spacing[2] / 2, -15, -5))
                mat->data()[i] = static_cast<std::uint8_t>(3);
            if (betw(zpos + spacing[2] / 2, -25, -15))
                mat->data()[i] = static_cast<std::uint8_t>(4);
            if (betw(zpos + spacing[2] / 2, -35, -25))
                mat->data()[i] = static_cast<std::uint8_t>(5);
        }
    }
    World<T> w;
    w.setSpacing(spacing);
    w.setDimensions(dim);
    w.setDensityArray(dens);
    w.setMaterialIndexArray(mat);
    w.addMaterialToMap(air);
    w.addMaterialToMap(pmma);
    for (int i = 0; i < 4; ++i)
        w.addMaterialToMap(pmma);
    w.makeValid();
    return w;
}

template <typename T>
World<T> generateTG195Case4World2()
{
    const std::array<std::size_t, 3> dim = { 1201, 1201, 60 };
    const std::array<T, 3> spacing = { 1, 1, 50 };
    const auto size = std::accumulate(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());

    Material air("Air, Dry (near sea level)");
    air.setStandardDensity(0.001205);
    Material pmma("H53.2813989847746C33.3715774096566O13.3470236055689");
    pmma.setStandardDensity(1.19);

    auto mat = std::make_shared<std::vector<std::uint8_t>>(size, 0);
    auto dens = std::make_shared<std::vector<T>>(size, static_cast<T>(air.standardDensity()));
    auto meas = std::make_shared<std::vector<std::uint8_t>>(size, 0);
    //generate cylindar
    auto circ_ind = circleIndices(T { 0 }, T { 0 }, dim, spacing, T { 320 } * T { 0.5 });
    auto pmma2_ind = circleIndices(T { -150 }, T { 0 }, dim, spacing, T { 5 });
    auto pmma1_ind = circleIndices(T { 0 }, T { 0 }, dim, spacing, T { 5 });
    for (std::size_t z = 0; z < dim[2]; ++z) {
        const T zpos = z * spacing[2] - (spacing[2] * dim[2]) / 2;
        for (const auto i : circ_ind) {
            const auto ind = z * dim[0] * dim[1] + i;
            mat->data()[ind] = static_cast<std::uint8_t>(1);
            dens->data()[ind] = pmma.standardDensity();
        }
        if (betw(zpos + spacing[2] / 2, -50, 50)) {
            const auto offset = z * dim[0] * dim[1];
            for (const auto i : pmma1_ind) {
                mat->data()[offset + i] = static_cast<std::uint8_t>(2);
                meas->data()[offset + i] = static_cast<std::uint8_t>(1);
            }
            for (const auto i : pmma2_ind) {
                mat->data()[offset + i] = static_cast<std::uint8_t>(3);
                meas->data()[offset + i] = static_cast<std::uint8_t>(1);
            }
        }
    }

    World<T> w;
    w.setSpacing(spacing);
    w.setDimensions(dim);
    w.setDensityArray(dens);
    w.setMaterialIndexArray(mat);
    w.setMeasurementMapArray(meas);
    w.addMaterialToMap(air);
    w.addMaterialToMap(pmma);
    for (int i = 0; i < 2; ++i)
        w.addMaterialToMap(pmma);
    w.makeValid();
    return w;
}

template <typename T>
bool TG195Case41AbsorbedEnergy(bool specter = false, bool wide_collimation = false)
{
    std::cout << "TG195 Case 4.1:\n";
    IsotropicSource<T> src;
    src.setPosition(-600.0, 0., 0.);
    std::array<T, 6> cos = { 0, 1, 0, 0, 0, 1 };
    src.setDirectionCosines(cos);

    T direction[3];
    vectormath::cross(cos.data(), direction);
    if (specter) {
        std::cout << "Specter of 120 kV W/Al\n";
        const auto specter = TG195_120KV<T>();
        src.setSpecter(specter.second, specter.first);
    } else {
        std::cout << "Monoenergetic specter of 56.4 kev\n";
        std::vector<T> s(1, 1.0), e(1, 56.4);
        src.setSpecter(s, e);
    }
    src.setHistoriesPerExposure(histPerExposure);
    src.setTotalExposures(nExposures);

    if (wide_collimation) {
        std::cout << "Collimation: 80 mm:\n";
        src.setCollimationAngles(std::atan(T { 160 } / 600) * 2, std::atan(T { 40 } / 600) * 2);
    } else {
        std::cout << "Collimation: 10 mm:\n";
        src.setCollimationAngles(std::atan(T { 160 } / 600) * 2, std::atan(T { 5 } / 600) * 2);
    }
    auto w = generateTG195Case4World1<T>();

    Transport<T> transport;

    auto res = transport(w, &src, nullptr, false);
    const auto total_hist = static_cast<T>(src.totalExposures() * src.historiesPerExposure());

    const T simtime = std::chrono::duration_cast<std::chrono::seconds>(res.simulationTime).count();
    std::cout << "Simulation time " << simtime << " seconds";
    std::cout << " with " << simtime / total_hist << " seconds*CPU core per history\n";
    auto dose = getEVperHistory(res, w.densityArray(), w.spacing(), total_hist);

    std::array<T, 4> voi_ev;
    for (std::size_t i = 0; i < 4; ++i) {
        voi_ev[i] = std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](auto d, auto m) -> T { return m == i + 2 ? d : 0; });
    }

    std::array<std::uint64_t, 4> voi_nevent;
    for (std::size_t i = 0; i < 4; ++i) {
        voi_nevent[i] = std::transform_reduce(std::execution::par_unseq, res.nEvents.cbegin(), res.nEvents.cend(), w.materialIndexArray()->begin(), 0.0, std::plus<>(), [=](auto d, auto m) { return m == i + 2 ? d : 0; });
    }

    std::array<T, 4> sim_ev;
    if (specter) {
        if (wide_collimation) {
            sim_ev = { 3586.59, 3537.84, 3378.99, 2672.21 };
        } else {
            sim_ev = { 13137.02, 2585.47, 1706.86, 1250.61 };
        }
    } else {
        if (wide_collimation) {
            sim_ev = { 3380.39, 3332.64, 3176.44, 2559.58 };
        } else {
            sim_ev = { 11592.27, 2576.72, 1766.85, 1330.53 };
        }
    }

    std::cout << "VOI, dxmc, dxmc #events, TG195, difference [ev/hist], difference [%]\n";
    for (std::size_t i = 0; i < voi_ev.size(); ++i)
        std::cout << i + 1 << ", " << voi_ev[i] << ", " << voi_nevent[i] << ", " << sim_ev[i] << ", " << voi_ev[i] - sim_ev[i] << ", " << (voi_ev[i] - sim_ev[i]) / sim_ev[i] * 100 << "\n";
    return true;
}

template <typename T>
bool TG195Case42AbsorbedEnergy(bool specter = false, bool wide_collimation = false)
{

    std::array<T, 36> sim_ev_center, sim_ev_pher;
    if (specter) {
        if (wide_collimation) {
            sim_ev_center = { 10.878025, 10.9243, 10.884625, 10.89795, 10.87265, 10.902675, 10.8994, 10.880875, 10.875475, 10.8862, 10.895975, 10.88105, 10.8996, 10.886225, 10.8934, 10.8942, 10.879025, 10.8855, 10.894125, 10.8898, 10.8916, 10.895875, 10.889525, 10.889775, 10.89365, 10.901875, 10.894475, 10.906975, 10.888025, 10.877475, 10.883325, 10.875925, 10.8881, 10.886775, 10.88975, 10.900075 };
            sim_ev_pher = { 115.34325, 113.76275, 109.16925, 101.706, 91.562975, 78.39105, 61.388325, 40.08625, 22.471075, 11.781725, 6.14551, 3.42218, 2.05605, 1.35319, 0.96088275, 0.743808, 0.61922025, 0.55457575, 0.5309405, 0.55428325, 0.6219885, 0.74445025, 0.96480125, 1.3481875, 2.0611025, 3.4154, 6.15532, 11.7854, 22.461525, 40.13715, 61.42595, 78.328975, 91.481375, 101.61325, 109.10425, 113.8365 };
        } else {
            sim_ev_center = { 11.35475, 11.399475, 11.38755, 11.402175, 11.400225, 11.370625, 11.402625, 11.37715, 11.385375, 11.4096, 11.399825, 11.376675, 11.3698, 11.3613, 11.38865, 11.377925, 11.3692, 11.371525, 11.3807, 11.36, 11.37645, 11.379075, 11.379975, 11.368975, 11.377675, 11.38375, 11.3871, 11.3844, 11.37395, 11.3821, 11.371825, 11.395575, 11.379075, 11.372125, 11.4001, 11.39895 };
            sim_ev_pher = { 116.7915, 115.31075, 110.532, 103.1535, 92.9892, 79.7158, 62.57135, 40.92065, 22.966025, 12.169, 6.4216075, 3.5657225, 2.1609475, 1.418235, 1.014725, 0.780394, 0.64618, 0.5763575, 0.559597, 0.5776625, 0.648646, 0.7762825, 1.008508, 1.412125, 2.1584, 3.577945, 6.4296425, 12.1714, 23.039025, 40.926825, 62.451325, 79.759575, 93.074725, 103.269, 110.595, 115.27 };
        }
    } else {
        if (wide_collimation) {
            sim_ev_center = { 11.630625, 11.632925, 11.617175, 11.624825, 11.633075, 11.600225, 11.61655, 11.6235, 11.592875, 11.6258, 11.612, 11.612775, 11.608675, 11.623975, 11.611325, 11.6174, 11.6234, 11.627975, 11.60745, 11.632875, 11.628275, 11.6239, 11.61645, 11.617375, 11.621775, 11.6178, 11.6444, 11.61515, 11.626375, 11.64605, 11.63335, 11.628425, 11.622, 11.6198, 11.59835, 11.609925 };
            sim_ev_pher = { 99.665175, 98.300175, 94.509325, 88.48595, 80.1113, 69.261125, 55.124425, 37.2351, 21.754, 11.767, 6.269845, 3.460205, 2.073845, 1.3435025, 0.94542875, 0.71714775, 0.58643225, 0.52293525, 0.49996925, 0.5225535, 0.5875545, 0.719903, 0.94029425, 1.3461175, 2.07283, 3.4740625, 6.2371725, 11.79715, 21.73405, 37.28175, 55.1853, 69.2588, 80.036275, 88.397, 94.640025, 98.332825 };
        } else {
            sim_ev_center = { 12.168675, 12.11255, 12.163, 12.1358, 12.0881, 12.10275, 12.135925, 12.12425, 12.14915, 12.149575, 12.1482, 12.15955, 12.130775, 12.148475, 12.153925, 12.1383, 12.13475, 12.148675, 12.1576, 12.1452, 12.15665, 12.16225, 12.156625, 12.156875, 12.13905, 12.1561, 12.159325, 12.137025, 12.149725, 12.09575, 12.1471, 12.1231, 12.1211, 12.1368, 12.134825, 12.13955 };
            sim_ev_pher = { 101.29375, 99.802475, 96.108725, 89.966275, 81.369475, 70.5424, 56.2835, 38.1283, 22.348775, 12.166625, 6.515755, 3.643485, 2.180365, 1.41379, 0.984371, 0.7510235, 0.6165785, 0.54522075, 0.52163, 0.54515775, 0.61643625, 0.750186, 0.9914865, 1.4138375, 2.174915, 3.6299975, 6.5142425, 12.200725, 22.2865, 38.146275, 56.287325, 70.64165, 81.424625, 89.82015, 96.16995, 99.92435 };
        }
    }

    std::cout << "TG195 Case 4.2:\n";
    IsotropicSource<T> src;
    src.setHistoriesPerExposure(histPerExposure);
    src.setTotalExposures(nExposures);

    if (specter) {
        std::cout << "Specter of 120 kV W/Al\n";
        const auto specter = TG195_120KV<T>();
        src.setSpecter(specter.second, specter.first);
    } else {
        std::cout << "Monoenergetic specter of 56.4 kev\n";
        std::vector<T> s(1, T { 1 }), e(1, T { 56.4 });
        src.setSpecter(s, e);
    }
    if (wide_collimation) {
        std::cout << "Collimation: 80 mm:\n";
        src.setCollimationAngles(std::atan(T { 160 } / 600) * 2, std::atan(T { 40 } / 600) * 2);
    } else {
        std::cout << "Collimation: 10 mm:\n";
        src.setCollimationAngles(std::atan(T { 160 } / 600) * 2, std::atan(T { 5 } / 600) * 2);
    }

    auto w = generateTG195Case4World2<T>();

    Transport<T> transport;

    std::cout << "Angle, center dxmc [eV/hist], nEvents, pher dxmc [eV/hist], nEvents, center TG195 [eV/hist], pher TG195 [eV/hist], simtime [s], diff center[%], diff pher[%]\n";

    //simulate 36 projections
    for (std::size_t i = 0; i < 36; ++i) {
        std::cout << "Prosessing angle " << i * 10 << " [" << static_cast<int>(i / 36.0 * 100.0) << " %]";
        const auto nHistories = src.historiesPerExposure() * src.totalExposures();
        const T angle = (i * 10) * DEG_TO_RAD<T>();
        std::array<T, 3> rot_axis = { 0, 0, 1 };
        std::array<T, 3> pos { -600, 0, 0 };
        std::array<T, 6> cos = { 0, 1, 0, 0, 0, 1 };
        vectormath::rotate(pos.data(), rot_axis.data(), angle);
        src.setPosition(pos);
        vectormath::rotate(cos.data(), rot_axis.data(), angle);
        vectormath::rotate(&cos[3], rot_axis.data(), angle);
        src.setDirectionCosines(cos);

        T dir[3];
        vectormath::cross(cos.data(), dir);

        auto res = transport(w, &src, nullptr, false);
        const auto total_hist = static_cast<T>(src.totalExposures() * src.historiesPerExposure());

        const T simtime = std::chrono::duration_cast<std::chrono::milliseconds>(res.simulationTime).count();
        auto dose = getEVperHistory(res, w.densityArray(), w.spacing(), total_hist);
        auto matBegin = w.materialIndexArray()->cbegin();
        auto dose_pher = std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), matBegin, T { 0 }, std::plus<>(), [](auto d, auto m) -> T { return m == 3 ? d : 0; });
        auto dose_cent = std::transform_reduce(std::execution::par_unseq, dose.cbegin(), dose.cend(), matBegin, T { 0 }, std::plus<>(), [](auto d, auto m) -> T { return m == 2 ? d : 0; });
        auto nevents_pher = std::transform_reduce(std::execution::par_unseq, res.nEvents.cbegin(), res.nEvents.cend(), matBegin, 0, std::plus<>(), [](auto d, auto m) { return m == 3 ? d : 0; });
        auto nevents_cent = std::transform_reduce(std::execution::par_unseq, res.nEvents.cbegin(), res.nEvents.cend(), matBegin, 0, std::plus<>(), [](auto d, auto m) { return m == 2 ? d : 0; });

        std::cout << "\r" << i * 10 << ", " << dose_cent << ", " << nevents_cent << ", " << dose_pher << ", " << nevents_pher << ", ";
        std::cout << sim_ev_center[i] << ", " << sim_ev_pher[i] << ", " << simtime / 1000 << "s, ";
        std::cout << (dose_cent - sim_ev_center[i]) / sim_ev_center[i] * 100 << ", " << (dose_pher - sim_ev_pher[i]) / sim_ev_pher[i] * 100 << "\n";
    }

    return true;
}

template <typename T>
bool testAttenuation()
{

    Material pmma("H53.2813989847746C33.3715774096566O13.3470236055689");
    pmma.setStandardDensity(1.19);

    T nist_att = 1.970E-01;
    T dxmc_att = pmma.getTotalAttenuation(56.4);
    auto diff = (dxmc_att - nist_att) / nist_att * 100;

    const T energy = 56.4;
    Material m("Tissue, Soft (ICRP)");
    m.setStandardDensity(1.3);
    std::array<T, 3> spacing = { .1, .1, 1 };
    std::array<std::size_t, 3> dim = { 1, 1, 200 };
    const auto size = std::accumulate(dim.cbegin(), dim.cend(), (std::size_t)1, std::multiplies<>());
    const auto standardDensity = static_cast<T>(m.standardDensity());
    auto dens = std::make_shared<std::vector<T>>(size, standardDensity);
    auto mat = std::make_shared<std::vector<std::uint8_t>>(size, static_cast<std::uint8_t>(0));
    auto meas = std::make_shared<std::vector<std::uint8_t>>(size, 0);
    World<T> w;
    w.setDimensions(dim);
    w.setSpacing(spacing);
    w.setDensityArray(dens);
    w.setMaterialIndexArray(mat);
    w.setMeasurementMapArray(meas);
    w.addMaterialToMap(m);

    w.makeValid();

    PencilSource<T> pen;
    pen.setHistoriesPerExposure(histPerExposure);
    pen.setPhotonEnergy(energy);
    pen.setTotalExposures(nExposures);
    std::array<T, 6> cos = { 1, 0, 0, 0, 1, 0 };
    pen.setDirectionCosines(cos);
    pen.setPosition(0, 0, -400);

    const auto tot_hist = pen.historiesPerExposure() * pen.totalExposures();

    Transport<T> transport;
    auto res_naive = transport(w, &pen, nullptr, false);
    std::fill(meas->begin(), meas->end(), 1);
    w.makeValid();
    auto res_force = transport(w, &pen, nullptr, false);

    std::vector<T> att(res_naive.dose.size());
    for (int i = 0; i < dim[2]; ++i)
        att[i] = std::exp(-(i + 1) * spacing[2] * T { 0.1 } * m.standardDensity() * m.getTotalAttenuation(56.4));

    const auto att_max = *std::max_element(att.cbegin(), att.cend());
    const auto naive_max = *std::max_element(res_naive.dose.cbegin(), res_naive.dose.cend());
    const auto force_max = *std::max_element(res_force.dose.cbegin(), res_force.dose.cend());

    const T rms_naive = std::sqrt((1.0 / att.size()) * std::transform_reduce(res_naive.dose.cbegin(), res_naive.dose.cend(), att.cbegin(), T { 0.0 }, std::plus<>(), [=](auto e, auto a) -> T { return (e / naive_max - a / att_max) * (e / naive_max - a / att_max); }));
    const T rms_force = std::sqrt((1.0 / att.size()) * std::transform_reduce(res_force.dose.cbegin(), res_force.dose.cend(), att.cbegin(), T { 0.0 }, std::plus<>(), [=](auto e, auto a) -> T { return (e / force_max - a / att_max) * (e / force_max - a / att_max); }));

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
    if (rms_naive * 100 < T { 0.2 } && rms_force * 100 < T { 0.2 }) {
        std::cout << "SUCCESS\n\n";
        return true;
    }
    std::cout << "FAILURE\n\n";
    return false;
}

int main(int argc, char* argv[])
{
    //test120Specter();
    //testAttenuation<float>();

    auto success = true;

    success = success && TG195Case2AbsorbedEnergy<float>(false, false, false); // mono, tomo
    success = success && TG195Case2AbsorbedEnergy<float>(false, true, false);
    success = success && TG195Case2AbsorbedEnergy<float>(true, false, false);
    success = success && TG195Case2AbsorbedEnergy<float>(true, true, false);
    
    success = success && TG195Case3AbsorbedEnergy<float>(false, false);
    success = success && TG195Case3AbsorbedEnergy<float>(false, true);
    success = success && TG195Case3AbsorbedEnergy<float>(true, false);
    success = success && TG195Case3AbsorbedEnergy<float>(true, true);
    
    success = success && TG195Case41AbsorbedEnergy<float>(false, false);
    success = success && TG195Case41AbsorbedEnergy<float>(false, true);
    success = success && TG195Case41AbsorbedEnergy<float>(true, false);
    success = success && TG195Case41AbsorbedEnergy<float>(true, true);

    success = success && TG195Case42AbsorbedEnergy<float>(false, false);
    success = success && TG195Case42AbsorbedEnergy<float>(false, true);
    success = success && TG195Case42AbsorbedEnergy<float>(true, false);
    success = success && TG195Case42AbsorbedEnergy<float>(true, true);
    
    std::cout << "Press any key to exit";
    std::string dummy;
    std::cin >> dummy;
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}