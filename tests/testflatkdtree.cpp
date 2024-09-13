

#include "dxmc/world/kdtreeflat.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <cstdint>
#include <iostream>
#include <variant>
#include <vector>

int main()
{
    std::vector<std::variant<dxmc::WorldSphere<>>> items;
    const double r = 4;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            items.push_back(dxmc::WorldSphere<> { r, { 3 * r * i, 3 * r * j, 0 } });

    std::vector<std::variant<dxmc::WorldSphere<>>*> itemptrs;
    for (auto& item : items)
        itemptrs.push_back(&item);

    dxmc::KDTreeFlat<dxmc::WorldSphere<>> tree;
    tree.setData(itemptrs);

    const auto aabb = tree.calculateAABB();

    dxmc::Particle p;
    p.pos = { 0, -100, 0 };
    p.dir = { 0, 1, 0 };
    auto t = tree.intersect(p, aabb);
    auto tv = tree.intersectVisualization(p, aabb);

    return 0;
}
