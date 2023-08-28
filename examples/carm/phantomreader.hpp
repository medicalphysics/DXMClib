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

#include <array>
#include <map>
#include <string>
#include <vector>

#pragma once

class ICRP110PhantomReader {
public:
    ICRP110PhantomReader() = default;

    std::vector<std::uint8_t> media() const;
    std::vector<std::uint8_t> organs() const;
    std::array<double, 3> spacing() const { return m_spacing; }
    std::array<std::size_t, 3> dimensions() const { return m_dim; }

protected:
    static ICRP110PhantomReader readFemalePhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path);
    static ICRP110PhantomReader readMalePhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path);

private:
    std::vector<std::uint8_t> m_organ_data;
    std::array<double, 3> m_spacing;
    std::array<std::size_t, 3> m_dim;
    std::vector<std::uint8_t> m_organ_media;
    std::vector<std::string> m_organ_name;
    std::vector<double> m_organ_density;
    std::vector<std::map<std::size_t, double>> m_media_composition;
};