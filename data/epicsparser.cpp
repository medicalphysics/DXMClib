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
                    if (seg.I == 943) { // imaginary anomalous scattering factor
                        elements[seg.Z].setImaginaryAnomalousSF(seg.data);
                    }
                    if (seg.I == 944) { // imaginary anomalous scattering factor
                        elements[seg.Z].setRealAnomalousSF(seg.data);
                    }
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

std::string writePairVector(const std::vector<std::pair<double, double>>& vec)
{
    std::stringstream ss;
    for (const auto& [first, second] : vec) {

        ss << "{" << first << ", " << second << "}"
           << ",";
    }
    return "{" + ss.str() + "}";
}

std::string writePairVector(const std::string& name, const std::vector<std::pair<double, double>>& vec)
{
    std::string start = "auto " + name + " =";
    return start + writePairVector(vec) + ";\n";
}

std::string writeShellsMap(const std::string& name, const std::map<std::uint8_t, AtomicShell>& shells)
{
    std::string start = "std::map<std::uint8_t, Shell> " + name + ";\n";
    std::string end = "\n";

    std::stringstream ss;
    for (const auto& [shell_id, obj] : shells) {
        ss << name << "[" << static_cast<int>(shell_id) << "]= {\n";
        ss << ".shell=" << static_cast<int>(shell_id) << ",\n";
        ss << ".bindingEnergy=" << obj.bindingEnergy() << ",\n";
        ss << ".numberOfElectrons=" << obj.numberOfElectrons() << ",\n";
        ss << ".hartreeFockOrbital_0=" << obj.hartreeFockOrbital_0() << ",\n";
        ss << ".numberOfPhotonsPerInitVacancy=" << obj.numberOfPhotonsPerInitVacancy() << ",\n";
        ss << ".energyOfPhotonsPerInitVacancy=" << obj.energyOfPhotonsPerInitVacancy() << ",\n";
        ss << ".photo=" << writePairVector(obj.photoelectricData());
        ss << "\n};";
    }
    return start + ss.str() + end;
}

std::string writeElementsMap(const std::string& name, const std::map<std::uint8_t, AtomicElement>& elements)
{
    std::string start = "std::map<std::uint8_t, Atom> " + name + "; \n";
    std::string end = "\n";
    std::stringstream ss;
    for (const auto& [Z, obj] : elements) {
        std::string exp = name + "[" + std::to_string(static_cast<int>(Z)) + "]";
        ss << "{\n";
        ss << exp << "= {\n ";
        ss << ".Z=" << static_cast<int>(Z) << ",\n";
        ss << ".atomicWeight=" << obj.atomicWeight() << ",\n";
        ss << "};\n";
        ss << exp << ".photo=" << writePairVector(obj.photoelectricData()) << ";\n";
        ss << exp << ".incoherent=" << writePairVector(obj.incoherentData()) << ";\n";
        ss << exp << ".coherent=" << writePairVector(obj.coherentData()) << ";\n";
        ss << exp << ".formFactor=" << writePairVector(obj.formFactor()) << ";\n";
        ss << exp << ".scatterFactor=" << writePairVector(obj.incoherentSF()) << ";\n";
        const std::string shellname = "shells";
        ss << writeShellsMap(shellname, obj.shells());
        ss << exp << ".shells =" << shellname << ";\n";
        ss << "}\n";
    }
    return start + ss.str() + end;
}

bool EPICSparser::writeMaterialHeaderFile(const std::string& filename) const
{
    if (m_elements.size() == 0)
        return false;
    std::ofstream fstream(filename);
    if (!fstream.is_open())
        return false;

    auto start = R"V0G0N(
//This file is autogenerated by dxmclib
#pragma once
#include <vector>
#include <map>
#include <utility>
namespace dxmclib{
namespace material{
struct Shell{
    std::uint8_t shell = 0;
    double bindingEnergy = 0;
    double numberOfElectrons = 0; 
    double hartreeFockOrbital_0 = 0;
    double numberOfPhotonsPerInitVacancy = 0;
    double energyOfPhotonsPerInitVacancy = 0;
    std::vector<std::pair<double, double>> photo;
};
struct Atom{
    std::uint8_t Z=0;
    double atomicWeight=0;
    std::vector<std::pair<double, double>> photo;
    std::vector<std::pair<double, double>> incoherent;
    std::vector<std::pair<double, double>> coherent;
    std::vector<std::pair<double, double>> formFactor;
    std::vector<std::pair<double, double>> scatterFactor;
    std::map<std::uint8_t, Shell> shells;
};

static std::map<std::uint8_t, Atom> generateAtoms() {
)V0G0N";

    auto end = R"V0G0N(return atoms;}}})V0G0N";

    fstream << start;

    fstream << std::endl;

    const auto& photo = m_elements.at(16).photoelectricData();

    fstream << writeElementsMap("atoms", m_elements);

    fstream << end;

    fstream.close();
    return true;
}
