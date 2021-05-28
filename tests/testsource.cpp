#include "dxmc/source.hpp"

#include <cassert>
#include <chrono>
#include <future>
#include <iostream>
#include <thread>

using namespace dxmc;
constexpr double ERRF = 1e-4;

template <Floating T>
bool isEqual(T f1, T f2)
{
    return std::abs(f1 - f2) < ERRF;
}

template <typename T>
void initiateAll()
{
    PencilSource<T> pen;
    DXSource<T> dx;
    IsotropicSource<T> iso;
    CTAxialSource<T> ax;
    auto tube = ax.tube();

    CTTopogramSource<T> top;

    CTAxialDualSource<T> de_ax;

    CTSpiralSource<T> spiral;
    CTSpiralDualSource<T> de;
    CTAxialSource<T> from_spiral(spiral);
    CTAxialDualSource<T> from_spiral_de(de);
    CBCTSource<T> cbct;
}

template <typename T>
bool testDXSourceAngles(T pang, T sang, T tubeRotation)
{
    DXSource<T> src;

    std::array<T, 2> angles = { pang, sang };
    src.setTubeRotationDeg(tubeRotation);
    //std::cout << "angles set: " << angles[0] << ", " << angles[1];
    src.setSourceAnglesDeg(angles);
    auto anglesres = src.sourceAnglesDeg();
    //std::cout << " angles res: " << anglesres[0] << ", " << anglesres[1];
    //std::cout << " tube rot: " << tubeRotation << '\n';
    return isEqual(angles[0], anglesres[0]) && isEqual(angles[1], anglesres[1]);
}

template <typename T>
bool testDXSourceAnglesMany()
{
    std::cout << "Testing DX source angles: ";
    bool success = true;
    std::array<T, 7> angles = { -80, -60, -30, 0, 30, 60, 80 };
    auto tube_rot = angles;
    for (auto ap : angles)
        for (auto as : angles)
            for (auto tr : tube_rot) {
                success = success && testDXSourceAngles(ap, as, tr);
                assert(success);
            }

    assert(success);
    if (success)
        std::cout << "Success\n";
    else
        std::cout << "Failure\n";
    return success;
}

template <typename T>
bool testDXCalibration()
{
    std::cout << "Testing DX source calibration: ";
    DXSource<T> src;
    CTDIPhantom<T> world(320);

    auto dens = world.densityArray();
    auto mat = world.materialIndexArray();
    auto measArr = world.measurementMapArray();
    std::fill(dens->begin(), dens->end(), dens->at(0));
    std::fill(mat->begin(), mat->end(), 0);
    std::fill(measArr->begin(), measArr->end(), 1);

    const auto& dim = world.dimensions();
    const auto& spacing = world.spacing();

    auto idx = dim[0] / 2 + dim[1] / 2 * dim[0];

    src.setTotalExposures(4);
    src.setHistoriesPerExposure(500000);
    std::array<T, 6> cosines = { 1, 0, 0, 0, 1, 0 };
    src.setDirectionCosines(cosines);
    std::array<T, 2> fieldSize = { 10, 10 };
    src.setFieldSize(fieldSize);
    src.setDap(T { 1 });
    src.validate();

    Transport<T> transport;
    auto res = transport(world, &src);

    // dose at (0,0,0)
    const T measureArea = spacing[0] * spacing[1] / 10 / 10;
    const auto calcdose = src.dap() / measureArea;

    const auto measuredose = res.dose[idx];

    const auto diff = 100 * (calcdose - measuredose) / calcdose;

    if ((diff > -10) && (diff < 10)) {
        std::cout << "Success\n";
        return true;
    }
    std::cout << "Failure\n";
    return false;
}

template <typename T>
bool testCTForcedCalibration()
{
    std::cout << "Testing CT axial source calibration: \n";
    CTAxialSource<T> src;
    CTDIPhantom<T> world(320);

    //src.setPitch(0.5);
    src.setExposureAngleStepDeg(1);
    src.setHistoriesPerExposure(100000);

    using holePosition = typename CTDIPhantom<T>::HolePosition;
    std::array<holePosition, 5> position = { holePosition::Center, holePosition::West, holePosition::East, holePosition::South, holePosition::North };

    auto measurement = world.measurementMapArray();
    std::fill(measurement->begin(), measurement->end(), 0);

    Transport<T> transport;
    transport.setOutputMode(Transport<T>::OUTPUTMODE::EV_PER_HISTORY);
    auto res_n = transport(world, &src);

    for (auto pos : position) {
        const auto& holeIndices = world.holeIndices(pos);
        auto meas_ptr = measurement->data();
        for (const auto& idx : holeIndices)
            meas_ptr[idx] = 1;
    }

    auto res_m = transport(world, &src);

    struct Result {
        T mean_n = 0;
        T mean_m = 0;
        std::uint64_t events_n = 0;
        std::uint64_t events_m = 0;
    };

    std::array<Result, 5> measureDose;

    for (std::size_t i = 0; i < position.size(); ++i) {
        const auto& holeIndices = world.holeIndices(position[i]);
        for (const auto& idx : holeIndices) {
            measureDose[i].mean_n += res_n.dose[idx];
            measureDose[i].mean_m += res_m.dose[idx];
            measureDose[i].events_n += res_n.nEvents[idx];
            measureDose[i].events_m += res_m.nEvents[idx];
        }
        measureDose[i].mean_n /= static_cast<T>(holeIndices.size());
        measureDose[i].mean_m /= static_cast<T>(holeIndices.size());
    }

    bool success = true;

    for (auto m : measureDose) {
        std::cout << m.mean_n << ", " << m.mean_m << ", ";
        std::cout << m.events_n << ", " << m.events_m << "\n";

        success = success && std::abs(m.mean_n - m.mean_n) < T { 0.1 };
    }

    if (success)
        std::cout << "Success\n";
    else
        std::cout << "Failure\n";
    return success;
}

template <typename T>
bool testCTCalibration()
{
    std::cout << "Testing CT axial source calibration: ";
    CTAxialSource<T> src;
    CTDIPhantom<T> world(320);

    //src.setPitch(0.5);
    src.setExposureAngleStepDeg(1);
    src.setHistoriesPerExposure(10000);

    Transport<T> transport;
    auto res = transport(world, &src);

    using holePosition = typename CTDIPhantom<T>::HolePosition;
    std::array<holePosition, 5> position = { holePosition::Center, holePosition::West, holePosition::East, holePosition::South, holePosition::North };

    std::array<T, 5> measureDose;
    measureDose.fill(0);
    for (std::size_t i = 0; i < 5; ++i) {
        const auto& holeIndices = world.holeIndices(position[i]);
        for (const auto& idx : holeIndices)
            measureDose[i] += res.dose[idx];
        measureDose[i] /= static_cast<T>(holeIndices.size());
    }
    T pher = 0;
    for (int i = 1; i < 5; ++i)
        pher += measureDose[i];
    pher /= T { 4.0 };
    auto cent = measureDose[0];
    T ctdi = (pher * 2) / 3 + cent / 3;
    if (ctdi > 0.990 || ctdi < 1.01) {
        std::cout << "Success\n";
        return true;
    }
    std::cout << "Failure\n";
    return false;
}

template <typename T>
bool testTopogramCalibration()
{
    CTTopogramSource<T> top;
    top.setCollimation(2);
    top.getCalibrationValue(LOWENERGYCORRECTION::NONE);

    return true;
}
int main(int argc, char* argv[])
{
    std::cout << "Testing sources\n";
    initiateAll<float>();
    initiateAll<double>();
    bool success = true;
    /*success = success && testTopogramCalibration<float>();
    success = success && testDXCalibration<float>();
    success = success && testDXCalibration<double>();
    success = success && testCTCalibration<float>();
    success = success && testCTCalibration<double>();
    success = success && testCTForcedCalibration<float>();
    success = success && testDXSourceAnglesMany<float>();
    success = success && testDXSourceAnglesMany<double>();*/
    success = success && testCTForcedCalibration<float>();

    if (success)
        std::cout << "Test sources : Success\n";
    else
        std::cout << "Test sources : Failure\n";

    return !success;
}
