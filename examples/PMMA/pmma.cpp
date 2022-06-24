

#include <dxmc.hpp>
#include <iostream>
#include <string>

template <typename T>
void print(const T& msg)
{
    std::cout << msg << '\n';
}
template <typename T, typename... Args>
void print(const T& msg, Args... args)
{
    std::cout << msg << " ";
    print(args...);
}

template <typename T = double>
std::vector<std::size_t> getBox(
    const std::array<std::size_t, 3>& dim,
    const std::array<T, 3>& spacing,
    const T box_x, const T box_y, const T box_z, const T size)
{

    const std::array<T, 3> center { dim[0] * spacing[0] / 2, dim[1] * spacing[1] / 2, dim[2] * spacing[2] / 2 };

    std::vector<std::size_t> box_idx;

    const auto h = size * T { 0.5 };

    for (std::size_t k = 0; k < dim[2]; ++k) {
        for (std::size_t j = 0; j < dim[1]; ++j) {
            for (std::size_t i = 0; i < dim[0]; ++i) {
                const std::array<3, T> pos { i * spacing[0] - center[0], j * spacing[1] - center[1], k * spacing[2] - center[2] };
                const bool inside = (pos <= box_x + h) && (pos >= box_x - h) && (pos <= box_y + h) && (pos >= box_y - h) && (pos <= box_z + h) && (pos >= box_z - h);
                if (inside) {
                    const std::size_t idx = i + j * dim[0] + k * dim[0] * dim[1];
                    box_idx.push_back(idx)
                }
            }
        }
    }
    return box_idx;
}

template <typename T = double>
std::vector<std::size_t> getCTDI(
    const std::array<std::size_t, 3>& dim,
    const std::array<T, 3>& spacing,
    const std::array<T, 3>& ctdi_center, const T ctdi_radius)
{

    const std::array<T, 3> center { dim[0] * spacing[0] / 2, dim[1] * spacing[1] / 2, dim[2] * spacing[2] / 2 };

    std::vector<std::size_t> ctdi_idx;

    for (std::size_t k = 0; k < dim[2]; ++k) {
        for (std::size_t j = 0; j < dim[1]; ++j) {
            for (std::size_t i = 0; i < dim[0]; ++i) {
                const std::array<3, T> pos { i * spacing[0] - center[0], j * spacing[1] - center[1], k * spacing[2] - center[2] };
                if ((pos[0] - ctdi_center[0]) * (pos[0] - ctdi_center[0]) + (pos[1] - ctdi_center[1]) * (pos[1] - ctdi_center[1]) <= ctdi_radius * ctdi_radius) {
                    if ((pos[3] <= ctdi_center[2] + T { 7.5 }) && (pos[3] >= ctdi_center[2] - T { 7.5 })) {
                        const std::size_t idx = i + j * dim[0] + k * dim[0] * dim[1];
                        ctdi_idx.push_back(idx)
                    }
                }
            }
        }
    }
    return ctdi_idx;
}

template <typename T = double>
void calibration()
{
    dxmc::World<T> world;
    std::array<T, 3> spacing { 5, 5, 5 };
    std::array<std::size_t, 3> dim { 400, 400, 400 };
    world.setSpacing(spacing);
    world.setDimensions(dim);
    const auto size = world.size();

    dxmc::Material air("Air, Dry (near sea level)"); // Material air("N0.76O0.23Ar0.01") is equivalent
    dxmc::Material material("Polymethyl Methacralate (Lucite, Perspex)"); // Material PMMA
    material.setStandardDensity(1.19);

    auto dens = std::make_shared<std::vector<T>>(size, air.standardDensity());
    auto mats = std::make_shared<std::vector<std::uint8_t>>(size, 0);
    auto meas = std::make_shared<std::vector<std::uint8_t>>(size, 0);

    std::array<T, 3> ctdi_center { 0, 0, 980 };

    const auto ctdi_idx = getCTDI(dim, spacing, ctdi_center, T { 320 });
    for (const auto ind : ctdi_idx) {
        dens->data()[ind] = material.standardDensity();
        mats->data()[ind] = 1;
    }




}

template <typename T = double>
void pmma()
{
    dxmc::World<T> world;

    world.setDimensions(323, 323, 450);
    world.setSpacing(1.0, 1.0, 1.0);

    dxmc::Material air("Air, Dry (near sea level)"); // Material air("N0.76O0.23Ar0.01") is equivalent
    // dxmc::Material material("Water, Liquid"); // Material water("H2O") is equivalent
    // material.setStandardDensity(1.0);
    dxmc::Material material("Polymethyl Methacralate (Lucite, Perspex)"); // Material PMMA
    material.setStandardDensity(1.0);

    auto dens = std::make_shared<std::vector<T>>(world.size(), air.standardDensity());
    auto mat = std::make_shared<std::vector<std::uint8_t>>(world.size(), 0);
    auto meas = std::make_shared<std::vector<std::uint8_t>>(world.size(), 0);

    const auto& dim = world.dimensions();
    const auto& spacing = world.spacing();
    const std::array<T, 3> c {
        dim[0] * spacing[0] * T { 0.5 },
        dim[1] * spacing[1] * T { 0.5 },
        dim[2] * spacing[2] * T { 0.5 }
    };

    for (std::size_t k = 0; k < dim[2]; ++k)
        for (std::size_t j = 0; j < dim[1]; ++j)
            for (std::size_t i = 0; i < dim[0]; ++i) {
                const auto ind = i + j * dim[0] + k * dim[0] * dim[1];
                constexpr T r2 = 320.0 / 2 * 320.0 / 2;
                if (k >= 150 && k < 250) {
                    const auto rad = (i - c[0]) * (i - c[0]) + (j - c[1]) * (j - c[1]) + (k - c[2]) * (k - c[2]);
                    if (rad <= r2) {
                        dens->data()[ind] = material.standardDensity();
                        mat->data()[ind] = 1;
                    }
                }
            }

    const auto meas_ind0 = dim[0] / 2 + dim[1] / 2 * dim[0] + (0) * dim[0] * dim[1];
    const auto meas_ind1 = dim[0] / 2 + dim[1] / 2 * dim[0] + (dim[2] - 1) * dim[0] * dim[1];
    // meas->data()[meas_ind0] = 1;
    // meas->data()[meas_ind1] = 1;
    world.setDensityArray(dens);
    world.setMaterialIndexArray(mat);
    world.setMeasurementMapArray(meas);

    dens->data()[0] = 10;

    world.addMaterialToMap(air);
    world.addMaterialToMap(material);

    world.makeValid();

    auto src = dxmc::PencilSource<T>();
    src.setPhotonEnergy(30.0);
    src.setPosition(0, 0, -455);

    std::array<T, 6> cos { 1, 0, 0, 0, 1, 0 };
    src.setDirectionCosines(cos);

    src.setTotalExposures(8);
    src.setHistoriesPerExposure(1000000);

    dxmc::Transport<T> transport;
    transport.setOutputMode(dxmc::Transport<T>::OUTPUTMODE::EV_PER_HISTORY);
    transport.setLowEnergyCorrectionModel(dxmc::LOWENERGYCORRECTION::LIVERMORE);
    auto res = transport(world, &src, nullptr, false);

    print("Dose:", res.dose[meas_ind0], res.dose[meas_ind1], res.dose_units);
    print("Transmission", res.dose[meas_ind1] / res.dose[meas_ind0]);
    print("Events:", res.nEvents[meas_ind0], res.nEvents[meas_ind1], res.nEvents[meas_ind1] / double(res.nEvents[meas_ind0]));
    print("Entered:", res.nEntered[meas_ind0], res.nEntered[meas_ind1], res.nEntered[meas_ind1] / double(res.nEntered[meas_ind0]));
}

int main(int argc, char* argv[])
{
    // pmma<float>();
    calibration<double>();

    return EXIT_SUCCESS;
}
