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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.hpp"
#include "dxmc/world/baseobject.hpp"
#include "dxmc/world/kdtree.hpp"

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

namespace dxmc {

template <Floating T>
class ItemCollection : public BaseObject<T> {
public:
    ItemCollection()
        : BaseObject<T>()
    {
    }

    template <typename U>
    requires std::is_base_of<BaseObject<T>, U>::value void addItem(std::shared_ptr<U> item)
    {
        if (item) {
            // making sure we add an unique object
            auto item_base = std::static_pointer_cast<BaseObject<T>>(item);
            if (std::find(m_items.begin(), m_items.end(), item_base) == m_items.end()) {
                m_items.push_back(std::static_pointer_cast<BaseObject<T>>(item));

                std::vector<BaseObject<T>*> raw_values;
                raw_values.reserve(m_items.size());
                for (auto& i : m_items)
                    raw_values.push_back(i.get());
                m_kdtree = KDTree(raw_values);
                this->m_aabb = m_kdtree.AABB();
            }
        }
    }
    void clearItems() { m_items.clear(); }
    void translate(const std::array<T, 3>& dist) override
    {
        m_kdtree.translate(dist);
        this->m_aabb = m_kdtree.AABB();
    }
    std::array<T, 3> center() const override
    {
        std::array<T, 3> res { 0, 0, 0 };
        for (const auto& i : m_items) {
            const auto c = i->center();
            res = vectormath::add(res, c);
        }
        for (auto& n : res)
            n /= m_items.size();
        return res;
    }
    std::optional<T> intersect(const Particle<T>& p) const override
    {
        return m_kdtree.intersect(p, this->m_aabb);
    }
    std::optional<T> intersect(const Particle<T>& p, const std::array<T, 2>& tbox) const override
    {
        return m_kdtree.intersect(p, this->m_aabb, tbox);
    }
    T transport(Particle<T>& p, RandomState& state) override
    {
        return 0;
    }

protected:
private:
    std::vector<std::shared_ptr<BaseObject<T>>> m_items;
    KDTree<BaseObject<T>*> m_kdtree;
};

}