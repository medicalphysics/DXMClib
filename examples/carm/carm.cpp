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

#include "dxmc/beams/dxbeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/transportprogress.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/enclosedroom.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include "phantomreader.hpp"

#include <iostream>
#include <vector>

dxmc::AAVoxelGrid<double, 5, 1, 0> testPhantom()
{
    auto d = ICRP110PhantomReader::readFemalePhantom("AF.dat", "AF_media.dat", "AF_organs.dat");

    dxmc::AAVoxelGrid<double, 5, 1, 0> phantom;
    using Material = dxmc::Material<double, 5>;
    std::vector<Material> materials;
    for (auto& w : d.mediaComposition()) {
        auto mat_cand = Material::byWeight(w);
        if (mat_cand)
            materials.push_back(mat_cand.value());
        else
            throw std::runtime_error("error");
    }
    phantom.setData(d.dimensions(), d.densityData(), d.mediaData(), materials);
    phantom.setSpacing(d.spacing());
    return phantom;
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
        std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << std::flush << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

int main()
{

    using CTDIPhantom = dxmc::CTDIPhantom<double, 5, 1>;
    using Mesh = dxmc::TriangulatedMesh<double, 5, 1>;
    using Sphere = dxmc::WorldSphere<double, 5, 1>;
    using VGrid = dxmc::AAVoxelGrid<double, 5, 1, 0>;
    using Room = dxmc::EnclosedRoom<double, 5, 1>;
    using World = dxmc::World<double, CTDIPhantom, Mesh, Sphere, VGrid, Room>;
    using Viz = dxmc::VisualizeWorld<double>;

    World world {};
    world.reserveNumberOfItems(5);
    auto& carm = world.addItem<Mesh>({ "carm.stl" });
    auto& table = world.addItem<Mesh>({ "table.stl" });
    table.translate({ -30, 0, 0 });

    auto& room = world.addItem<Room>();
    room.setInnerRoomAABB({ -250, -150, -120, 150, 150, 120 });
    room.setWallThickness(2);
    const auto lead = dxmc::Material<double, 5>::byZ(82).value();
    const auto lead_atom = dxmc::AtomHandler<double>::Atom(82);
    const auto lead_dens = dxmc::AtomHandler<double>::Atom(82).standardDensity;
    room.setMaterial(lead, lead_dens * 0.2 / 2.0);

    /*auto& ctdi = world.addItem<CTDIPhantom>({});
    ctdi.translate({ 16, 0, 9 });
    table.translate({ 0, 0, -17 });
    ctdi.translate({ 0, 0, -17 });
    auto ctdi_aabb = ctdi.AABB();
*/

    auto& phantom = world.addItem(testPhantom());
    phantom.rollAxis(2, 0);
    phantom.rollAxis(2, 1);
    phantom.flipAxis(2);
    auto table_aabb = table.AABB();
    auto phantom_aabb = phantom.AABB();
    phantom.translate({ -40, 0, table_aabb[5] - phantom_aabb[2] });

    auto& doctor = world.addItem(testPhantom());
    // doctor.rollAxis(2, 0);
    // doctor.rollAxis(2, 1);
    doctor.flipAxis(1);
    auto doctor_aabb = doctor.AABB();
    doctor.translate({ -40, -40, -doctor_aabb[2] - 120 });

    world.build();

    // adding beam
    using Beam = dxmc::DXBeam<double>;
    const std::array<double, 3> source_pos = { 0, 0, -70 };
    Beam beam(source_pos);
    beam.setBeamSize(12, 12, 114);
    beam.setNumberOfExposures(200);
    beam.setNumberOfParticlesPerExposure(1000000);
    beam.setDAPvalue(1.0);

    dxmc::Transport transport;
    runDispatcher(transport, world, beam);

    double max_doctor_dose = 0;
    for (const auto& d : doctor.getDoseScores()) {
        max_doctor_dose = std::max(max_doctor_dose, d.dose());
    }

    constexpr double dosenorm = 0.0000005;
    std::cout << "Max doctor dose " << max_doctor_dose << " mGy, dose norm: " << dosenorm << std::endl;

    Viz viz(world);
    viz.addColorByValueItem(&doctor);
    viz.addColorByValueItem(&phantom);
    viz.setColorByValueMinMax(0, dosenorm);
    auto buffer = viz.createBuffer<double>(2048, 2048);
    viz.addLineProp(beam, 114, .2);

    viz.setDistance(400);
    viz.setAzimuthalAngleDeg(60);
    std::vector<double> angles;
    for (std::size_t i = 0; i < 12; ++i)
        angles.push_back(i * 30);

    for (auto a : angles) {
        viz.setPolarAngleDeg(a);
        viz.suggestFOV();
        viz.generate(world, buffer);
        std::string name = "test" + std::to_string(int(a)) + ".png";
        viz.savePNG(name, buffer);
        std::cout << "Rendertime " << buffer.renderTime.count() << " ms"
                  << "(" << 1000.0 / buffer.renderTime.count() << " fps)" << std::endl;
    }

    viz.setAzimuthalAngleDeg(120);
    for (auto a : angles) {
        viz.setPolarAngleDeg(a);
        viz.suggestFOV();
        viz.generate(world, buffer);
        std::string name = "test_low" + std::to_string(int(a)) + ".png";
        viz.savePNG(name, buffer);
        std::cout << "Rendertime " << buffer.renderTime.count() << " ms"
                  << "(" << 1000.0 / buffer.renderTime.count() << " fps)" << std::endl;
    }
    return 0;
}