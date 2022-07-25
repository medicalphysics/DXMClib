

#include "dxmc/world/kdtree.hpp"
#include "dxmc/world/triangulatedmesh.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
    dxmc::STLReader<double> reader;

    auto mesh = reader("untitled.stl");
    std::cout << "Message: " << reader.message() << std::endl;

    const auto& vertices = mesh.getVertices();
    const auto& faceIdx = mesh.getFaceIndices();

    dxmc::KDTree<double> kdtree(vertices, faceIdx);


    dxmc::Particle<double> p;
    p.pos = { .0, .0, -5 };
    p.dir = { 0, 0, 1 };

    const auto& aabb = mesh.calculateAABB();
    const auto& vert = mesh.getVertices();

    auto res = kdtree.intersect(p, vert, aabb);
    if (res)
        std::cout << *res; 


    return EXIT_SUCCESS;
}
