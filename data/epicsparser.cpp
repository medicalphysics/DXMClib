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

#include "epicsparser.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

struct DataSegment {
    std::uint8_t Z = 0;
    std::uint8_t Yi = 0;
    std::uint8_t Yo = 0;
    std::uint8_t Iflag = 0;
    std::uint8_t C = 0;
    std::uint8_t I = 0;
    std::uint8_t S = 0;
    std::uint8_t X1 = 0;
    std::uint8_t dataDim = 0;
    double AW = 0.0;
    double A = 0.0;
    std::vector<double> data;

    void clear()
    {
        Z = 0;
        Yi = 0;
        Yo = 0;
        Iflag = 0;
        C = 0;
        I = 0;
        S = 0;
        X1 = 0;
        dataDim = 0;
        AW = 0.0;
        A = 0.0;
        data.clear();
    }
};

void processFirstHeaderLine(const std::string& line, DataSegment& segment)
{
    segment.Z = std::stoi(line.substr(0, 3));
    segment.A = std::stod(line.substr(3, 3));
    segment.Yi = std::stoi(line.substr(7, 2));
    segment.Yo = std::stoi(line.substr(10, 2));
    segment.Iflag = std::stoi(line.substr(31, 1));
    segment.AW = std::stod(line.substr(13, 10));
}
void processSecondHeaderLine(const std::string& line, DataSegment& segment)
{
    segment.C = std::stoi(line.substr(0, 2));
    segment.I = std::stoi(line.substr(2, 3));
    segment.S = std::stoi(line.substr(5, 3));
    segment.X1 = static_cast<std::uint8_t>(std::stod(line.substr(21, 10)));
}

EPICSparser::EPICSparser(const std::string& path)
{
    read(path);
}

std::vector<double> split(const std::string& s)
{
    constexpr std::size_t sublen = 16;
    const auto n_data = s.size() / 16;
    std::vector<double> data(n_data);
    for (std::size_t i = 0; i < n_data; ++i) {
        data[i] = std::stod(s.substr(i * sublen, sublen));
    }
    return data;
}

void processSegments(const std::vector<DataSegment>& segments, std::map<std::uint8_t, AtomicElement>& elements)
{
    std::vector<int> shells;
    for (const auto& seg : segments) {
        if (!elements.contains(seg.Z)) { // adding uniqe element
            elements.emplace(seg.Z, seg.Z);
            //elements[seg.Z] = AtomicElement(seg.Z);
            elements[seg.Z].setAtomicWeight(seg.AW);
        }
        if (seg.Yi == 7) { // incoming photon
            if (seg.C == 71) { // coherent
                if (seg.X1 == 0) { // whole atom
                    if (seg.I == 0) {
                        elements[seg.Z].setCoherentData(seg.data);
                    }
                }
            }
            if (seg.C == 72) { // incoherent
                if (seg.X1 == 0) { // whole atom
                    if (seg.I == 0) {
                        elements[seg.Z].setIncoherentData(seg.data);
                    }
                }
            }
            if (seg.C == 73) { // photoelectric
                if (seg.X1 == 0) { // whole atom
                    if (seg.I == 0) {
                        elements[seg.Z].setPhotoelectricData(seg.data);
                    }
                } else {
                    if (seg.I == 0 && seg.Z==16)
                        shells.push_back(seg.X1);
                }
            }
        }
    }
    auto stop = true;
}

void EPICSparser::read(const std::string& path)
{
    std::ifstream stream;
    stream.open(std::string(path));

    std::vector<DataSegment> segments;

    if (stream.is_open()) {
        std::string line;
        std::uint8_t headerline = 1;
        DataSegment segment;
        std::size_t linenumber = 0;
        while (std::getline(stream, line)) {
            if (headerline == 0) {
                if (line.size() >= endIdx()) {
                    if (line.back() == '1') { // end of segment
                        segments.push_back(segment);
                        segment.clear();
                        headerline = 1;
                    }
                } else { // read data line
                    const auto data = split(line);
                    if (segment.dataDim == 0) { // find data dimenionality
                        segment.dataDim = data.size();
                    }
                    if (segment.dataDim != data.size()) {
                        bool test = false;
                    }
                    for (const auto d : data) {
                        segment.data.push_back(d);
                    }
                }
            } else if (headerline == 2) {
                processSecondHeaderLine(line, segment);
                headerline = 0;

            } else if (headerline == 1) {
                processFirstHeaderLine(line, segment);
                headerline++;
            }
            ++linenumber;
        }
    }

    processSegments(segments, m_elements);
}
