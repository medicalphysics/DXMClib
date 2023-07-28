

#include "dxmc/beams/isotropicmonoenergybeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <typename T>
void writeImage(const std::vector<T>& buffer, const std::string& name)
{
    std::ofstream file;
    file.open(name, std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
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
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

template <typename T>
std::vector<dxmc::Triangle<T>> getPyramid()
{
    std::vector<std::array<T, 3>> p;
    constexpr T d = 30;
    p.push_back({ 1, 1, 0 }); // 0
    p.push_back({ 1, -1, 0 }); // 1
    p.push_back({ -1, -1, 0 }); // 2
    p.push_back({ -1, 1, 0 }); // 3
    p.push_back({ 0, 0, 1 });
    for (auto& i : p)
        for (auto& j : i)
            j *= d;

    std::vector<dxmc::Triangle<T>> t;
    t.push_back({ p[0], p[1], p[4] });
    t.push_back({ p[1], p[2], p[4] });
    t.push_back({ p[2], p[3], p[4] });
    t.push_back({ p[3], p[0], p[4] });

    // underside
    t.push_back({ p[0], p[3], p[2] });
    t.push_back({ p[2], p[1], p[0] });

    return t;
}

template <typename T>
std::vector<dxmc::Triangle<T>> getBox(T scale = 1)
{
    std::vector<std::array<T, 3>> p;
    p.push_back({ 1, 1, 1 }); // 0
    p.push_back({ 1, 1, -1 }); // 1
    p.push_back({ 1, -1, 1 }); // 2
    p.push_back({ -1, 1, 1 }); // 3
    p.push_back({ -1, -1, 1 }); // 4
    p.push_back({ -1, 1, -1 }); // 5
    p.push_back({ 1, -1, -1 }); // 6
    p.push_back({ -1, -1, -1 }); // 7
    for (auto& i : p)
        for (auto& j : i)
            j *= scale;

    std::vector<dxmc::Triangle<T>> t;
    t.push_back({ p[0], p[3], p[4] });
    t.push_back({ p[0], p[4], p[2] });
    t.push_back({ p[6], p[2], p[4] });
    t.push_back({ p[6], p[4], p[7] });
    t.push_back({ p[7], p[4], p[3] });
    t.push_back({ p[7], p[3], p[5] });
    t.push_back({ p[5], p[1], p[6] });
    t.push_back({ p[5], p[6], p[7] });
    t.push_back({ p[1], p[0], p[2] });
    t.push_back({ p[1], p[2], p[6] });
    t.push_back({ p[5], p[3], p[0] });
    t.push_back({ p[5], p[0], p[1] });
    return t;
}

template <typename T, std::size_t N = 5, int L = 2>
void testMeshVisualization()
{

    using Mesh = dxmc::TriangulatedMesh<T, 5, 2>;
    using World = dxmc::World2<T, Mesh>;
    using Material = dxmc::Material<T, 5>;

    World world;

    const auto waterComp = dxmc::NISTMaterials<T>::Composition("Water, Liquid");
    auto water = Material::byWeight(waterComp).value();

    int option;
    option = 0; // Box
    option = 1; // Triangle
    option = 2; // Bunny
    option = 3; // bunny_low
    //option = 4; // duck
    

    if (option == 0) {
        const auto triangles = getBox<T>();
        Mesh mesh(triangles);
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    } else if (option == 1) {
        const auto triangles = getPyramid<T>();
        Mesh mesh(triangles);
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    } else if (option == 2) {
        Mesh mesh("bunny.stl");
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    } else if (option == 3) {
        Mesh mesh("bunny_low.stl");
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    } else if (option == 4) {
        Mesh mesh("duck.stl");
        mesh.setMaterial(water, 1);
        world.addItem(std::move(mesh));
    }

    world.build(T { 0 });

    dxmc::VisualizeWorld<T> viz(world);

    int height = 1024;
    int width = 1024;
    std::vector<T> buffer(height * width * 4, T { 1 });

    for (std::size_t i = 0; i < 12; ++i) {
        viz.setDistance(500);
        viz.setPolarAngle(std::numbers::pi_v<T> / 3.0);
        viz.setAzimuthalAngle(std::numbers::pi_v<T> * i / 6.0);
        // viz.setCameraPosition({ -60, -30, -10 });
        viz.suggestFOV();
        viz.generate(world, buffer, width, height);
        std::string name = "color_" + std::to_string(i) + ".png";
        viz.savePNG(name, buffer, width, height);
        //writeImage(buffer, name);
    }
}

template <dxmc::Floating T, bool MESH = true>
T testScoring()
{
    using Mesh = dxmc::TriangulatedMesh<T, 5, 2>;
    using Box = dxmc::WorldBox<T, 5, 2>;
    using World = dxmc::World2<T, Mesh, Box>;
    using Material = dxmc::Material<T, 5>;

    World world;
    const auto waterComp = dxmc::NISTMaterials<T>::Composition("Water, Liquid");
    auto water = Material::byWeight(waterComp).value();

    if constexpr (MESH) {

        auto tri = getBox<T>();
        auto& mesh = world.addItem<Mesh>({ tri });
        mesh.setMaterial(water, T { 1 });
    } else {
        auto& box = world.addItem<Box>({ T { 1 } });
        box.setMaterial(water, T { 1 });
    }
    world.build();

    using Beam = dxmc::IsotropicMonoEnergyBeam<T>;
    Beam beam({ 0, 0, -1000 }, { 1, 0, 0, 0, 1, 0 }, 60);
    beam.setNumberOfExposures(8);
    beam.setNumberOfParticlesPerExposure(1E5);

    dxmc::Transport transport;

    runDispatcher(transport, world, beam);

    const auto* item = world.getItemPointers()[0];
    const auto& dose = item->dose();

    return dose.energyImparted();
}

int main(int argc, char* argv[])
{
    std::cout << "Testing tetrahedal mesh\n";
    testMeshVisualization<double>();

    std::cout << "Testing dose scoring of mesh\n";
    bool success = true;
    const auto dmesh = testScoring<double, true>();
    const auto dbox = testScoring<double, false>();
    const auto ddiff = (dmesh / dbox - 1) * 100;
    success = success && std::abs(ddiff) < 0.1;
    if (success)
        std::cout << "SUCCESS: ";
    else
        std::cout << "FAILURE: ";
    std::cout << "Dose mesh: " << dmesh << ", dose box: ";
    std::cout << dbox << ", difference[%] " << ddiff << std::endl;

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
