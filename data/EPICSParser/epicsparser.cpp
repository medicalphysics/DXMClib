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
#include "dxmc/material/atomserializer.hpp"

#include <algorithm>
#include <charconv>
#include <fstream>
#include <iostream>
#include <optional>
#include <sstream>
#include <vector>

struct DataSegment {
    int Z = 0;
    int Yi = 0;
    int Yo = 0;
    // int Iflag = 0;
    int C = 0;
    int I = 0;
    int S = 0;
    int X1 = 0;
    int dataDim = 0;
    double AW = 0.0;
    double A = 0.0;
    std::vector<double> data;

    void clear()
    {
        Z = 0;
        Yi = 0;
        Yo = 0;
        // Iflag = 0;
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
    segment.AW = std::stod(line.substr(13, 10));
    // segment.Iflag = std::stoi(line.substr(31, 1));
}
void processSecondHeaderLine(const std::string& line, DataSegment& segment)
{
    segment.C = std::stoi(line.substr(0, 2));
    segment.I = std::stoi(line.substr(2, 3));
    segment.S = std::stoi(line.substr(5, 3));
    segment.X1 = static_cast<int>(std::stod(line.substr(21, 10)));
}

EPICSparser::EPICSparser(const std::string& path)
{
    read(path);
}

EPICSparser::EPICSparser(std::vector<char>& data)
{
    deserializeElements(data);
}

std::vector<double> split( std::string& s)
{
    // For some reason the exponent indicated is 'D' not 'E' for EPICS 2023 data (sigh).
    std::replace(s.begin(), s.end(), 'D', 'E');

    constexpr std::size_t sublen = 16;
    const auto n_data = s.size() / 16;
    std::vector<double> data(n_data);
    for (std::size_t i = 0; i < n_data; ++i) {
        data[i] = std::stod(s.substr(i * sublen, sublen));
    }
    return data;
}

void processSegments(const std::vector<DataSegment>& segments, std::map<std::uint64_t, AtomicElementHandler>& elements)
{
    for (const auto& seg : segments) {
        if (!elements.contains(seg.Z)) { // adding uniqe element
            elements.emplace(seg.Z, seg.Z);
            elements[seg.Z].setAtomicWeight(seg.AW);
        }
        if (seg.Yi == 7) { // incoming photon
            if (seg.C == 71) { // coherent
                if (seg.X1 == 0) { // whole atom
                    if (seg.I == 0) { // integrated cross section
                        elements[seg.Z].setCoherentData(seg.data);
                    }
                }
            }
            if (seg.C == 72) { // incoherent
                if (seg.X1 == 0) { // whole atom
                    if (seg.I == 0) { // integrated cross section
                        elements[seg.Z].setIncoherentData(seg.data);
                    }
                    if (seg.Yo == 7) {
                        if (seg.I == 10) { // avg energy of scattered photon
                            elements[seg.Z].setIncoherentAvgEnergyScatteredPhoton(seg.data);
                        }
                    }
                }
            }
            if (seg.C == 73) { // photoelectric
                if (seg.X1 == 0) { // whole atom
                    if (seg.I == 0) { // integrated cross section
                        elements[seg.Z].setPhotoelectricData(seg.data);
                    }
                } else { // subshell
                    if (seg.I == 0) { // integrated cross section
                        elements[seg.Z].setShellPhotoelectricData(seg.X1, seg.data);
                    }
                }
            }
            if (seg.C == 93) { // coherent and incoherent data (not cross section)
                if (seg.X1 == 0) { // whole atom
                    if (seg.I == 941) { // form factor
                        elements[seg.Z].setFormFactor(seg.data);
                    }
                    /*
                    if (seg.I == 943) { // imaginary anomalous scattering factor
                        elements[seg.Z].setImaginaryAnomalousSF(seg.data);
                    }
                    if (seg.I == 944) { // imaginary anomalous scattering factor
                        elements[seg.Z].setRealAnomalousSF(seg.data);
                    }
                    */
                    if (seg.I == 942) { // incoherent scatter function
                        elements[seg.Z].setIncoherentSF(seg.data);
                    }
                }
            }
        }
        if (seg.Yi == 0) { // other atom/shell parameters
            if (seg.C == 91) { // subshell parameters
                if (seg.I == 912) { // Number of electrons
                    elements[seg.Z].setShellNumberOfElectrons(seg.data);
                } else if (seg.I == 913) { // Binding energy
                    elements[seg.Z].setShellBindingEnergy(seg.data);
                } else if (seg.I == 914) { // electron kinetic energy
                    elements[seg.Z].setShellKineticEnergy(seg.data);
                }
            }
            if (seg.C == 92) { // transition data
                if (seg.Yo == 7) { // photon data
                    if (seg.I == 933) { // number of particles per vacancy
                        elements[seg.Z].setShellNumberOfPhotonsPerInitVacancy(seg.data);
                    }
                    if (seg.I == 934) { // avg energy of particles per vacancy
                        elements[seg.Z].setShellEnergyOfPhotonsPerInitVacancy(seg.data);
                    }
                }
            }
        }
    }
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
        while (std::getline(stream, line)) {
            if (headerline == 0) {
                if (line.size() >= endIdx()) {
                    if (line.back() == '1' || line.ends_with("1\r")) { // end of segment
                        segments.push_back(segment);
                        segment.clear();
                        headerline = 1;
                    }
                } else { // read data line
                    const auto data = split(line);
                    if (segment.dataDim == 0) { // find data dimenionality
                        segment.dataDim = data.size();
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
        }
    }
    processSegments(segments, m_elements);
}

void EPICSparser::readHartreeFockProfiles(const std::string& path)
{
    std::ifstream stream;
    stream.open(std::string(path));

    if (stream.is_open()) {
        std::string line;
        std::size_t linenumber = 0;
        while (std::getline(stream, line)) {
            if (linenumber > 0) { // we dont read the header
                auto start = line.data();
                auto end = start + line.size();
                int Z = 0;
                int shell = 0;
                double J = 0.0;
                auto [ptrZ, errZ] = std::from_chars(start, end, Z);
                if (errZ == std::errc {}) {
                    auto [ptrS, errS] = std::from_chars(ptrZ + 1, end, shell);
                    if (errS == std::errc {}) {
                        auto [ptrJ, errJ] = std::from_chars(ptrS + 1, end, J);
                        if (errJ == std::errc {}) {
                            m_elements[Z].setShellHartreeFockProfile_0(shell, J);
                        }
                    }
                }
            }
            linenumber++;
        }
    }
}
void EPICSparser::readStandardDensities(const std::string& path)
{
    std::ifstream stream;
    stream.open(std::string(path));

    if (stream.is_open()) {
        std::string line;
        std::size_t linenumber = 0;
        while (std::getline(stream, line)) {
            if (linenumber > 0) { // we dont read the header
                auto start = line.data();
                auto end = start + line.size();
                int Z = 0;
                double dens = 0;

                auto [ptrZ, errZ] = std::from_chars(start, end, Z);
                if (errZ == std::errc {}) {
                    auto [ptrD, errD] = std::from_chars(ptrZ + 1, end, dens);
                    if (errD == std::errc {}) {
                        m_elements[Z].setStandardDensity(dens);
                    }
                }
            }
            linenumber++;
        }
    }
}

std::vector<char> EPICSparser::serializeElements() const
{
    std::map<std::uint64_t, dxmc::AtomicElement<double>> map;
    for (const auto& [key, element] : m_elements) {
        map[key] = element.atom();
    }
    return dxmc::AtomSerializer::serializeAtoms(map);
}

void EPICSparser::deserializeElements(std::vector<char>& data)
{
    auto map = dxmc::AtomSerializer::deserializeAtoms(data);
    m_elements.clear();
    for (const auto& [key, atom] : map) {
        m_elements[key] = AtomicElementHandler(atom);
    }
}
