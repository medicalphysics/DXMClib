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
#include "dxmc/material/atomicelement.hpp"
#include "dxmc/material/atomserializer.hpp"

#include <execution>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace dxmc {

template <Floating T>
class AtomHandler {

public:
    static const AtomicElement<T>& Atom(std::uint64_t Z)
    {
        auto& instance = Instance();
        if (instance.m_elements.contains(Z)) {
            return instance.m_elements.at(Z);
        }
        return instance.m_dummyElement;
    }

    static const std::map<std::uint64_t, AtomicElement<T>>& allAtoms()
    {
        const auto& instance = Instance();
        return instance.m_elements;
    }

    AtomHandler(const AtomHandler&) = delete;
    void operator=(const AtomHandler&) = delete;

protected:
    static const AtomHandler& Instance()
    {
        static AtomHandler instance;
        return instance;
    }
    AtomHandler()
    {
        // reading data
        std::string datapath;
        if (std::filesystem::exists("physicslists.bin")) {
            datapath = "physicslists.bin";
        } else {
            const std::string EPICSdataBuildPath { DXMCLIB_PHYSICSLISTSPATH };
            if (std::filesystem::exists(EPICSdataBuildPath)) {
                datapath = EPICSdataBuildPath;
            }
        }

        // reading buffer
        std::ifstream buffer_file(datapath, std::ios::binary);
        if (buffer_file.good()) {
            std::vector<char> data(std::istreambuf_iterator<char>(buffer_file), {});
            if constexpr (std::is_same<T, double>::value) {
                m_elements = AtomSerializer::deserializeAtoms(data);
            } else {
                auto elements = AtomSerializer::deserializeAtoms(data);
                for (const auto& [key, element] : elements) {
                    m_elements[key] = typecastAtom(element);
                }
            }
        }
    }

    static std::vector<std::pair<T, T>> typecastPairVector(const std::vector<std::pair<double, double>>& r)
    {
        std::vector<std::pair<T, T>> l(r.size());
        std::transform(std::execution::par_unseq, r.cbegin(), r.cend(), l.begin(), [](const auto& p) {
            return std::make_pair(static_cast<T>(p.first), static_cast<T>(p.second));
        });
        return l;
    }
    static AtomicShell<T> typecastShell(const AtomicShell<double>& r)
    {
        AtomicShell<T> l;
        l.shell = r.shell;
        l.numberOfElectrons = static_cast<T>(r.numberOfElectrons);
        l.bindingEnergy = static_cast<T>(r.bindingEnergy);
        l.kineticEnergy = static_cast<T>(r.kineticEnergy);
        l.HartreeFockOrbital_0 = static_cast<T>(r.HartreeFockOrbital_0);
        l.numberOfPhotonsPerInitVacancy = static_cast<T>(r.numberOfPhotonsPerInitVacancy);
        l.energyOfPhotonsPerInitVacancy = static_cast<T>(r.energyOfPhotonsPerInitVacancy);
        l.photoel = typecastPairVector(r.photoel);
        return l;
    }
    static AtomicElement<T> typecastAtom(const AtomicElement<double>& r)
    {
        AtomicElement<T> l;
        l.Z = r.Z;
        l.atomicWeight = static_cast<T>(r.atomicWeight);
        l.standardDensity = static_cast<T>(r.standardDensity);
        l.coherent = typecastPairVector(r.coherent);
        l.formFactor = typecastPairVector(r.formFactor);
        l.incoherentSF = typecastPairVector(r.incoherentSF);
        l.incoherent = typecastPairVector(r.incoherent);
        l.photoel = typecastPairVector(r.photoel);
        l.incoherentMeanScatterEnergy = typecastPairVector(r.incoherentMeanScatterEnergy);
        for (const auto& [key, shell] : r.shells) {
            l.shells[key] = typecastShell(shell);
        }
        return l;
    }

private:
    AtomicElement<T> m_dummyElement;
    std::map<std::uint64_t, AtomicElement<T>> m_elements;
};

}