

#include "dxmc/world/kdtree.hpp"
#include "dxmc/world/triangulatedmesh.hpp"


#include <iostream>

int main(int argc, char* argv[])
{
    dxmc::STLReader<double> reader;

    auto mesh = reader("bunny.stl");
    std::cout << "Message: " << reader.message() << std::endl;

    

    
    

    dxmc::Particle<double> p;
    p.pos = { 0,0,0 };
    p.dir = { 1, 0, 0 };
    dxmc::vectormath::normalize(p.dir);

    const auto& aabb_p = mesh.calculateAABB();
    std::array<double, 3> dist;
    for (std::size_t i = 0; i < 3; ++i)
        dist[i] = -(aabb_p[i] + aabb_p[i + 3] * 0.5);
    mesh.translate(dist);
    const auto& aabb = mesh.calculateAABB();

    auto triangles = mesh.getTriangles();
    dxmc::KDTree<double, dxmc::Triangle<double>> kdtree(triangles);
    std::cout << "KDTree depth : " << kdtree.depth() << std::endl;


    auto res = kdtree.intersect(p,  aabb);
    if (res)
        std::cout << "Intersect KDTree: " <<  *res<< std::endl; 

    //brute
    double bh_min = std::numeric_limits<double>::max();
    double bh_max = std::numeric_limits<double>::lowest();
    for (const auto& t : triangles) {
        auto h = t.intersect<1>(p);
        if (h) {
            bh_min = std::min(*h, bh_min);
            bh_max = std::max(*h, bh_max);
        }
    }
    std::cout << "Brute force hits: " << bh_min << ", " << bh_max << std::endl;




    return EXIT_SUCCESS;
}
