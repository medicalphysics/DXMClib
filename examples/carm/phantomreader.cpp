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

#include "phantomreader.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <numeric>
#include <optional>
#include <sstream>

// trim from start (in place)
inline void ltrim(std::string& s)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
        return !std::isspace(ch);
    }));
}

// trim from end (in place)
inline void rtrim(std::string& s)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
        return !std::isspace(ch);
    }).base(),
        s.end());
}

// trim from both ends (in place)
inline void trim(std::string& s)
{
    rtrim(s);
    ltrim(s);
}

std::vector<std::uint8_t> readASCIIData(const std::string& path)
{
    std::ifstream t;
    std::vector<std::uint8_t> buffer;
    t.open(path); // open input file
    if (t.is_open()) {
        t.seekg(0, std::ios::end); // go to the end
        auto length = t.tellg(); // report location (this is the length)
        buffer.resize(length);
        t.seekg(0, std::ios::beg); // go back to the beginning
        t.read(reinterpret_cast<char*>(&buffer[0]), length); // read the whole file into the buffer
        t.close();
    }
    return buffer;
}

struct Organ {
    std::string name;
    double density = 0;
    std::uint8_t ID = 0;
    std::uint8_t mediaID = 0;
};

std::vector<Organ> readASCIIOrgans(const std::string& organ_path)
{
    std::ifstream t(organ_path);
    std::stringstream buffer;
    if (t.is_open()) {
        buffer << t.rdbuf();
        t.close();
    } else {
        return;
    }

    auto to_uchar = [](std::string_view s) -> std::optional<std::uint8_t> {
        std::uint8_t value {};
        if (std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc {})
            return value;
        else
            return std::nullopt;
    };
    auto to_double = [](std::string_view s) -> std::optional<double> {
        double value {};
        if (std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc {})
            return value;
        else
            return std::nullopt;
    };

    std::vector<Organ> organs;

    for (std::string line; std::getline(buffer, line);) {
        if (line.size() >= 67) {
            Organ organ;
            // ID
            auto IDs = line.substr(0, 6);
            auto ID = to_uchar(IDs);
            if (ID)
                organ.ID = ID.value();
            auto name = line.substr(6, 55);
            trim(name);
            organ.name = name;
            auto TNRs = line.substr(55, 58);
            ltrim(TNRs);
            auto TNR = to_uchar(TNRs);
            if (TNR)
                organ.mediaID = TNR.value();
            auto denss = line.substr(61, 66);
            auto dens = to_double(denss);
            if (dens)
                organ.density = dens.value();

            organs.push_back(organ);
        }
    }

    return organs;
}

ICRP110PhantomReader ICRP110PhantomReader::readFemalePhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path)
{
    ICRP110PhantomReader data;
    data.m_dim = { 299, 137, 346 };
    data.m_spacing = { 1.775, 1.775, 4.84 };
    const auto size = std::reduce(data.m_dim.cbegin(), data.m_dim.cend(), 1, std::multiplies<>());
    data.m_organ_data = readASCIIData(phantom_path);
}