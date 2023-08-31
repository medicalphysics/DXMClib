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
#include <array>
#include <charconv>
#include <execution>
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

auto to_uchar(std::string_view s) -> std::optional<std::uint8_t>
{
    std::uint8_t value {};
    if (std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc {})
        return value;
    else
        return std::nullopt;
};
auto to_double(std::string_view s) -> std::optional<double>
{
    double value {};
    if (std::from_chars(s.data(), s.data() + s.size(), value).ec == std::errc {})
        return value;
    else
        return std::nullopt;
};

std::vector<std::uint8_t> readASCIIData(const std::string& path)
{
    std::ifstream t;
    std::vector<std::uint8_t> buffer;
    t.open(path); // open input file
    if (t.is_open()) {
        t.seekg(0, std::ios::end); // go to the end
        auto length = t.tellg(); // report location (this is the length)

        t.seekg(0, std::ios::beg); // go back to the beginning
        std::string test(length, '0');
        t.read(&test[0], length);
        t.close();

        test.erase(std::remove_if(std::execution::par_unseq, test.begin(), test.end(), [](unsigned char x) { return x == '\t' || x == '\n'; }), test.end());

        const auto buffer_size = test.size() / 4;

        buffer.resize(buffer_size, 255);
        for (std::size_t i = 0; i < buffer.size(); ++i) {
            char* start = test.data() + 4 * i;
            char* end = test.data() + 4 * (i + 1);

            while (std::isspace(*start))
                ++start;
            std::from_chars(start, end, buffer[i]);
        }
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
    std::vector<Organ> organs;
    std::ifstream t(organ_path);
    std::stringstream buffer;
    if (t.is_open()) {
        buffer << t.rdbuf();
        t.close();
    } else {
        return organs;
    }

    std::size_t lineNr = 0;

    for (std::string line; std::getline(buffer, line);) {
        lineNr++;
        if (line.size() >= 66 && lineNr > 4) {
            Organ organ;
            // ID
            auto IDs = line.substr(0, 6);
            auto ID = to_uchar(IDs);
            if (ID)
                organ.ID = ID.value();
            auto name = line.substr(6, 49);
            trim(name);
            organ.name = name;
            auto TNRs = line.substr(55, 3);
            ltrim(TNRs);
            auto TNR = to_uchar(TNRs);
            if (TNR)
                organ.mediaID = TNR.value();
            auto denss = line.substr(61, 5);
            auto dens = to_double(denss);
            if (dens)
                organ.density = dens.value();

            organs.push_back(organ);
        }
    }
    organs.push_back({ .name = "Air", .density = 0.001, .ID = 0, .mediaID = 0 });
    
    return organs;
}

struct Media {
    std::uint8_t ID = 0;
    std::map<std::size_t, double> composition;
    std::string name;
};

std::vector<Media> readASCIIMedia(const std::string& media_path)
{
    std::vector<Media> medias;
    std::ifstream t(media_path);
    std::stringstream buffer;
    if (t.is_open()) {
        buffer << t.rdbuf();
        t.close();
    } else {
        return medias;
    }

    std::size_t lineNr = 0;

    for (std::string line; std::getline(buffer, line);) {
        lineNr++;
        if (line.size() >= 154 && lineNr > 3) {
            Media media;
            // ID
            auto IDs = line.substr(0, 6);
            auto ID = to_uchar(IDs);
            if (ID)
                media.ID = ID.value();
            auto name = line.substr(6, 72);
            trim(name);
            media.name = name;

            // composition
            const std::array Zs = { 1, 6, 7, 8, 11, 12, 15, 16, 17, 19, 20, 26, 53 };
            for (std::size_t i = 0; i < Zs.size(); ++i) {
                constexpr std::size_t offset = 76;
                constexpr std::size_t step = 6;
                auto s = line.substr(offset + step * i, step);
                trim(s);
                auto w = to_double(s);
                if (w)
                    if (w.value() > 0.0)
                        media.composition[Zs[i]] = w.value();
            }
            medias.push_back(media);
        }
    }

    medias.push_back({ .ID = 0, .composition = { { 7, 80.0 }, { 8, 20.0 } }, .name = "Air" });    
    return medias;
}

bool sanitizeIDs(std::vector<std::uint8_t>& organ_data, std::vector<Organ>& organs, std::vector<Media>& media)
{
    // finding organ IDs in organ data;
    auto organIDs = organ_data; // copy of vector
    std::sort(std::execution::par_unseq, organIDs.begin(), organIDs.end());
    organIDs.erase(std::unique(std::execution::par_unseq, organIDs.begin(), organIDs.end()), organIDs.end());

    // do we have all organs, if not add air?
    {
        std::uint8_t teller = 0;
        while (teller < organs.size())


        bool missingOrgan = false;
        for (auto& id : organIDs) {
            auto pos = std::find_if(organs.cbegin(), organs.cend(), [id](const auto& o) { return o.ID == id; });
            missingOrgan = missingOrgan || pos == organs.cend();
            if (missingOrgan)
                bool test = false;
        }
    }
    if (missingOrgan)
        return false;
    //removing organs not in volume
    for (std::size_t i = 0; i < organs.size();) {
        const auto& o = organs[i];
        auto pos = std::find(organIDs.cbegin(), organIDs.cend(), o.ID);
        if (pos == organIDs.cend())
            organs.erase(organs.cbegin() + i);
        else
            ++i;
    }
    std::sort(organs.begin(), organs.end(), [](const auto& lh, const auto& rh) { return lh.ID < rh.ID; });
    for (std::uint8_t i = 0; i < organs.size(); ++i)
    {

    }


    std::replace(organ_data.begin(), organ_data.end(), )




    for (const auto oID : organIDs) {
        auto pos = std::find_if(organs.begin(), organs.end(), [oID](const auto& o) { return o.ID == oID; });
        if (pos != organs.cend())
            organs.erase(pos);
    }
    return true;
    // making organ IDs consecutive

    /*
    const auto maxOrganId = std::max_element(organIDs.cbegin(), organIDs.cend());
    std::vector<std::size_t> lin_map_organ(*maxOrganId + 1, 255);
    for (std::size_t i = 0; i < organIDs.size(); ++i) {
        lin_map_organ[organIDs[i]] = i;
    }
    for (auto& o : organs)
        o.ID = lin_map_organ[o.ID];
    std::transform(std::execution::par_unseq, organ_data.cbegin(), organ_data.cend(), organ_data.begin(), [lin_map_organ](const auto oIdx) {
        return lin_map_organ[oIdx];
    });
    std::sort(organs.begin(), organs.end(), [](const auto& lh, const auto& rh) { return lh.ID < rh.ID; });
    organs.erase(std::find_if(organs.cbegin(), organs.cend(), [](const auto& o) { return o.ID == 255; }), organs.end());

    // making media IDs consecutive
    std::vector<std::uint8_t> mediaIDs(media.size());
    std::transform(media.cbegin(), media.cend(), mediaIDs.begin(), [](const auto& o) { return o.ID; });
    // do we need all media?
    for (std::size_t i = 0; i < mediaIDs.size(); ++i) {
        bool found = false;
        if (mediaIDs[i] == 0)
            found = true;
        else
            for (const auto& o : organs)
                found = found || o.mediaID == mediaIDs[i];
        if (!found)
            mediaIDs[i] = 255;
    }
    std::sort(mediaIDs.begin(), mediaIDs.end());
    mediaIDs.erase(std::unique(mediaIDs.begin(), mediaIDs.end()), mediaIDs.end());
    mediaIDs.erase(std::find(mediaIDs.begin(), mediaIDs.end(), 255), mediaIDs.end());
    const auto maxMediaId = std::max_element(mediaIDs.cbegin(), mediaIDs.cend());
    std::vector<std::size_t> lin_map_media(*maxMediaId + 1, 0);
    for (std::size_t i = 0; i < mediaIDs.size(); ++i) {
        lin_map_media[mediaIDs[i]] = i;
    }

    for (auto& m : media)
        m.ID = lin_map_media[m.ID];
    for (auto& o : organs)
        o.mediaID = lin_map_media[o.mediaID];

    std::vector<std::uint8_t> media_data(organ_data.size(), 0);
    std::vector<std::uint8_t> organ_media_map(organs.size());
    for (std::size_t i = 0; i < organs.size(); ++i) {
        organ_media_map[i] = organs[i].mediaID;
    }
    std::transform(std::execution::par_unseq, organ_data.cbegin(), organ_data.cend(), media_data.begin(), [organ_media_map](const auto oIdx) {
        return organ_media_map[oIdx];
    });

    return;
    */
}

ICRP110PhantomReader ICRP110PhantomReader::readFemalePhantom(const std::string& phantom_path, const std::string& media_path, const std::string& organ_path)
{
    ICRP110PhantomReader data;
    data.m_dim = { 299, 137, 346 };
    data.m_spacing = { 1.775, 1.775, 4.84 };
    const auto size = std::reduce(data.m_dim.cbegin(), data.m_dim.cend(), 1, std::multiplies<>());
    data.m_organ_data = readASCIIData(phantom_path);
    // if (size != data.m_organ_data.size())
    //     return data;

    auto organs = readASCIIOrgans(organ_path);
    auto media = readASCIIMedia(media_path);
    sanitizeIDs(data.m_organ_data, organs, media);
    return data;
}