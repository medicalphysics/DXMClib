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

#include "dxmc/beams/isotropicbeam.hpp"
#include "dxmc/beams/isotropicmonoenergybeam.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldboxgrid.hpp"

#include <fstream>
#include <iostream>

using namespace dxmc;

constexpr bool SAMPLE_RUN = true;

const auto N_THREADS = std::max(std::thread::hardware_concurrency(), 1u);
const std::uint64_t N_EXPOSURES = SAMPLE_RUN ? N_THREADS * 10 : N_THREADS * 100u;
const std::uint64_t N_HISTORIES = SAMPLE_RUN ? 10000 : 100000;

template <Floating T>
struct ResultKeys {
    std::string rCase;
    std::string volume;
    std::string specter;
    std::string modus;
    std::string model;
    bool forced = false;
    T result = 0;
    T result_std = 0;
    int precision = sizeof(T);
};

class ResultPrint {
private:
    std::ofstream m_myfile;

public:
    ResultPrint()
    {
        m_myfile.open("validationTable.txt", std::ios::out | std::ios::app);
    }
    ~ResultPrint()
    {
        m_myfile.close();
    }
    void header()
    {
        m_myfile << "Case, Volume, Specter, Model, Mode, Result, Stddev, Forced, Precision\n";
    }

    template <typename T>
    void operator()(const ResultKeys<T>& r, bool terminal = true)
    {
        print(r, terminal);
    }

    template <typename T>
    void print(const ResultKeys<T>& r, bool terminal = true)
    {
        const std::string forced = r.forced ? "1" : "0";

        m_myfile << r.rCase << ", ";
        m_myfile << r.volume << ", ";
        m_myfile << r.specter << ", ";
        m_myfile << r.model << ", ";
        m_myfile << r.modus << ", ";
        m_myfile << r.result << ", ";
        m_myfile << r.result_std << ", ";
        m_myfile << forced << ", ";
        m_myfile << r.precision << '\n';

        if (terminal) {
            std::cout << r.rCase << ", ";
            std::cout << r.volume << ", ";
            std::cout << r.specter << ", ";
            std::cout << r.model << ", ";
            std::cout << r.modus << ", ";
            std::cout << r.result << ", ";
            std::cout << r.result_std << ", ";
            std::cout << forced << ", ";
            std::cout << r.precision << '\n';
        }
    }
};

class Print {
private:
    std::ofstream m_myfile;

public:
    Print()
    {
        m_myfile.open("validation.txt", std::ios::out | std::ios::app);
    }
    ~Print()
    {
        m_myfile.close();
    }
    void date()
    {
        auto now = std::chrono::system_clock::now();
        auto in_time_t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&in_time_t), "%c %Z");
        m_myfile << "Time: " << ss.str() << std::endl;
    }
    template <typename T>
    void operator()(const T& msg)
    {
        std::cout << msg;
        m_myfile << msg;
    }
    template <typename T, typename... Args>
    void operator()(const T& msg, Args... args)
    {
        std::cout << msg;
        m_myfile << msg;
        this->operator()(args...);
    }
    template <typename T>
    void display(const T& msg)
    {
        std::cout << msg << "\r";
    }
};

template <typename T>
void saveBinaryArray(const std::vector<T>& data, const std::string& name)
{
    auto myfile = std::fstream(name, std::ios::out | std::ios::binary);
    const auto bytes = data.size() * sizeof(T);
    myfile.write((char*)&data[0], bytes);
    myfile.close();
}

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
std::vector<std::pair<T, T>> TG195_specter(const std::vector<double>& raw)
{
    std::vector<std::pair<T, T>> s;
    s.reserve(raw.size() / 2);
    for (std::size_t i = 0; i < raw.size(); i = i + 2) {
        s.push_back(std::make_pair(static_cast<T>(raw[i] - 0.25), static_cast<T>(raw[i + 1])));
    }
    return s;
}
template <typename T>
auto TG195_120KV()
{
    return TG195_specter<T>(TG195_120KV_raw);
}
template <typename T>
auto TG195_100KV()
{
    return TG195_specter<T>(TG195_100KV_raw);
}
template <typename T>
auto TG195_30KV()
{
    return TG195_specter<T>(TG195_30KV_raw);
}

template <typename T, int NShells = 5>
std::pair<T, Material2<T, NShells>> TG195_soft_tissue()
{
    std::map<std::size_t, T> Zs;
    Zs[1] = T { 10.5 };
    Zs[6] = T { 25.6 };
    Zs[7] = T { 2.7 };
    Zs[8] = T { 60.2 };
    Zs[11] = T { 0.1 };
    Zs[15] = T { 0.2 };
    Zs[16] = T { 0.3 };
    Zs[17] = T { 0.2 };
    Zs[19] = T { 0.2 };

    auto mat_cand = Material2<T, NShells>::byWeight(Zs).value();

    return std::make_pair(T { 1.03 }, mat_cand);
}

template <typename T, int NShells = 5>
std::pair<T, Material2<T, NShells>> TG195_pmma()
{
    std::map<std::size_t, T> pmma_w;
    pmma_w[1] = T { 53.28139898 };
    pmma_w[6] = T { 33.37157741 };
    pmma_w[8] = T { 13.34702361 };

    auto mat_cand = Material2<T, NShells>::byWeight(pmma_w).value();
    return std::make_pair(T { 1.19 }, mat_cand);
}

template <typename T, typename W, typename B>
auto runDispatcher(T& transport, W& world, const B& beam)
{
    dxmc::TransportProgress progress;

    bool running = true;
    std::thread job([&]() {
        transport(world, beam, &progress);
        running = false;
    });
    std::string message;
    while (running) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

template <Floating T, int LOWENERGYCORRECTION = 2>
bool TG195Case2AbsorbedEnergy(bool specter = false, bool tomo = false)
{
    constexpr int NShells = 5;
    using Box = WorldBoxGrid<T, NShells, LOWENERGYCORRECTION>;
    using Material = Material2<T, NShells>;

    auto [mat_dens, mat] = TG195_soft_tissue<T, NShells>();

    World2<T, Box> world;

    {
        const T box_halfside = T { 39 } / 2;
        const T box_height = 20;
        const T box_zbegin = 155;
        const std::array<T, 6> box_aabb = { -box_halfside, -box_halfside, box_zbegin, box_halfside, box_halfside, box_zbegin + box_height };
        Box box(box_aabb);
        box.setVoxelDimensions({ 78, 78, 40 });
        box.setMaterial(mat);
        box.setMaterialDensity(mat_dens);
        world.addItem(std::move(box));
    }

    world.build(T { 160 });
    if (specter) {
        IsotropicBeam<T> beam;
        const auto specter = TG195_120KV<T>();
        beam.setEnergySpecter(specter);
        beam.setNumberOfExposures(N_EXPOSURES);
        beam.setNumberOfParticlesPerExposure(N_HISTORIES);
        if (tomo) {
            constexpr T alpha = 15 * DEG_TO_RAD<T>();
            const T h = 180 * std::atan(alpha);
            beam.setPosition({ 0, -h, 0 });
            const T y_ang_min = std::atan((h - T { 39 } / 2) / 180);
            const T y_ang_max = std::atan((h + T { 39 } / 2) / 180);
            const T x_ang = std::atan(T { 39 } / (2 * 180));
            beam.setCollimationAngles(-x_ang, y_ang_min, x_ang, y_ang_max);
        } else {
            const T collangle = std::atan(T { 39 } / (2 * T { 180 }));
            beam.setCollimationAngles(-collangle, -collangle, collangle, collangle);
        }
        Transport<T> transport;
        auto time_elapsed = runDispatcher(transport, world, beam);
    } else {
        IsotropicMonoEnergyBeam<T> beam;
        beam.setEnergy(T { 56.4 });
        beam.setNumberOfExposures(N_EXPOSURES);
        beam.setNumberOfParticlesPerExposure(N_HISTORIES);

        if (tomo) {
            constexpr T alpha = 15 * DEG_TO_RAD<T>();
            const T h = 180 * std::atan(alpha);
            beam.setPosition({ 0, -h, 0 });
            std::array<T, 3> x_cos = { 1, 0, 0 };
            std::array<T, 3> y_cos = { 0, 1, 0 };
            y_cos = vectormath::rotate(y_cos, x_cos, 15 * DEG_TO_RAD<T>());
            const T beta = std::atan((h - T { 39 } / 2) / 180);
            const T y_ang_min = alpha - beta;
            const T phi = std::atan((h + T { 39 } / 2) / 180);
            const T y_ang_max = phi - beta - alpha;

            const T collangle = std::atan(T { 39 } / (2 * T { 180 }));
            beam.setCollimationAngles(-collangle, y_ang_min, collangle, y_ang_max);
            beam.setDirectionCosines(x_cos, y_cos);
        } else {
            const T collangle = std::atan(T { 39 } / (2 * T { 180 }));
            beam.setCollimationAngles(-collangle, -collangle, collangle, collangle);
        }

        Transport<T> transport;
        auto time_elapsed = runDispatcher(transport, world, beam);
    }

    const auto& boxes = world.getItems<Box>();
    const auto& box = boxes.front();

    const auto total_hist = static_cast<T>(N_EXPOSURES * N_HISTORIES);

    T total_ev = 0;
    T total_ev_var = 0;
    for (std::size_t i = 0; i < box.totalNumberOfVoxels(); ++i) {
        total_ev += box.dose(i).energyImparted();
        total_ev_var += box.dose(i).varianceEnergyImparted();
    }

    const auto ev_history = total_ev / (total_hist / 1000);
    const auto ev_history_var = total_ev_var / (total_hist / 1000);

    std::string model;
    if (LOWENERGYCORRECTION == 0)
        model = "None";
    if (LOWENERGYCORRECTION == 1)
        model = "Livermore";
    if (LOWENERGYCORRECTION == 2)
        model = "IA";

    ResultKeys<T> res;
    res.specter = specter ? "120kVp" : "56.4keV";
    res.rCase = "Case 2";
    res.volume = "Total body";
    res.model = model;
    res.result = ev_history;
    res.result_std = std::sqrt(ev_history_var);
    res.modus = tomo ? "tomosyntesis" : "radiography";

    double TG195_value;
    if (specter) {
        if (tomo) {
            TG195_value = 30923.13;
        } else {
            TG195_value = 33125.98;
        }
    } else {
        if (tomo) {
            TG195_value = 30883.83;
        } else {
            TG195_value = 33171.40;
        }
    }

    ResultPrint print;
    print(res, false);

    if (LOWENERGYCORRECTION == 0) {
        res.model = "TG195";
        res.result = TG195_value;
        res.result_std = 0;
        print(res, false);
    }

    std::cout << "TG195 Case 2 for " << res.modus << " orientation and " << res.specter << " photons\n";
    std::cout << "Whole body: " << ev_history << " eV/history, TG195: " << TG195_value << " eV/history";
    std::cout << ", differense: " << ev_history - TG195_value << " ev/history [" << (ev_history / TG195_value - 1) * 100 << "]%  \n";

    std::cout << "eV per history: " << ev_history << std::endl;
    std::cout << "eV per history std: " << std::sqrt(ev_history_var) << std::endl;

    return true;
}

template <typename T>
bool runAll()
{
    auto success = true;

    success = success && TG195Case2AbsorbedEnergy<T, 0>(false, false);
    success = success && TG195Case2AbsorbedEnergy<T, 0>(true, false);
    success = success && TG195Case2AbsorbedEnergy<T, 0>(false, true);
    success = success && TG195Case2AbsorbedEnergy<T, 0>(true, true);

    success = success && TG195Case2AbsorbedEnergy<T, 1>(false, false);
    success = success && TG195Case2AbsorbedEnergy<T, 1>(true, false);
    success = success && TG195Case2AbsorbedEnergy<T, 1>(false, true);
    success = success && TG195Case2AbsorbedEnergy<T, 1>(true, true);

    success = success && TG195Case2AbsorbedEnergy<T, 2>(false, false);
    success = success && TG195Case2AbsorbedEnergy<T, 2>(true, false);
    success = success && TG195Case2AbsorbedEnergy<T, 2>(false, true);
    success = success && TG195Case2AbsorbedEnergy<T, 2>(true, true);

    return success;
}

void printStart()
{
    Print print;
    print("TG195 validation cases\n");
    print.date();

    ResultPrint resPrint;
    resPrint.header();
}

int main(int argc, char* argv[])
{
    printStart();

    auto success = true;
    // success = success && selectOptions<double>();
    success = success && runAll<double>();

    success = success && runAll<float>();

    // std::cout << "Press any key to exit";
    // std::string dummy;
    // std::cin >> dummy;
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}