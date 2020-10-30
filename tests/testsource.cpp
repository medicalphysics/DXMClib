#include "dxmc/source.h"

#include <cassert>
#include <chrono>
#include <future>
#include <iostream>
#include <thread>

using namespace dxmc;

constexpr double RAD2DEG = 180.0 / 3.14159265359;
constexpr double DEG2RAD = 1.0 / RAD2DEG;
constexpr double ERRF = 1e-3;

template <typename T>
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
    CTAxialDualSource<T> de_ax;

    CTSpiralSource<T> spiral;
    CTSpiralDualSource<T> de;
    CTAxialSource<T> from_spiral = spiral;
    CTAxialDualSource<T> from_spiral_de = de;
}

template <typename T>
bool testDXSourceAngles(T pang, T sang, T tubeRotation)
{

    DXSource<T> src;

    std::array<T, 2> angles = { pang, sang };
    src.setTubeRotationDeg(tubeRotation);
    std::cout << "angles set: " << angles[0] << ", " << angles[1];
    src.setSourceAnglesDeg(angles);
    auto anglesres = src.sourceAnglesDeg();
    std::cout << " angles res: " << anglesres[0] << ", " << anglesres[1];
    std::cout << " tube rot: " << tubeRotation << '\n';
    return isEqual(angles[0], anglesres[0]) && isEqual(angles[1], anglesres[1]);
}

template <typename T>
bool testDXSourceAnglesMany()
{
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

    typedef CTDIPhantom<T>::HolePosition holePosition;
    std::array<holePosition, 5> position = { holePosition::Center, holePosition::West, holePosition::East, holePosition::South, holePosition::North };

    std::array<T, 5> measureDose;
    measureDose.fill(0);
    for (std::size_t i = 0; i < 5; ++i) {
        auto holeIndices = world.holeIndices(position[i]);
        for (auto idx : holeIndices)
            measureDose[i] += res.dose[idx];
        measureDose[i] /= static_cast<double>(holeIndices.size());
    }
    T pher = 0;
    for (int i = 1; i < 5; ++i)
        pher += measureDose[i];
    pher /= 4.0;
    auto cent = measureDose[0];
    T ctdi = pher * 2.0 / 3.0 + cent / 3.0;
    if (ctdi > 0.990 || ctdi < 1.01) {
        std::cout << "Success\n";
        return true;
    }
    std::cout << "Failure\n";
    return false;
}

int main(int argc, char* argv[])
{
    initiateAll<float>();
    initiateAll<double>();
    bool success = true;
    success = success && testCTCalibration<float>();
    success = success && testCTCalibration<double>();
    success = success && testDXSourceAnglesMany<float>();
    success = success && testDXSourceAnglesMany<double>();

    return !success;
}
