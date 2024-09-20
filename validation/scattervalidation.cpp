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

Copyright 2024 Erlend Andersen
*/

#include "dxmc/interactions.hpp"

#include <algorithm>
#include <fstream>
#include <iostream>

struct Histogram {
    double start = 0;
    double stop = 1;
    int N = 100;
    std::vector<std::size_t> intensity;
    Histogram(int n = 100, double start_val = 0, double stop_val = 1)
    {
        N = std::clamp(n, 2, 1000);
        intensity.resize(N);
        start = std::min(start_val, stop_val);
        stop = std::max(start_val, stop_val);
    }
    void operator()(double v)
    {
        double fidx = ((v - start) / (stop - start)) * (N - 1);
        if (fidx > 0 && fidx < N - 1) {
            auto idx = static_cast<int>(fidx);
            intensity[idx]++;
        }
    }

    std::vector<double> getBins(bool midpoint = true) const
    {
        if (midpoint) {
            std::vector<double> bins(N);
            const auto step = (stop - start) / N;
            for (std::size_t i = 0; i < N; ++i) {
                bins[i] = start + (i + 0.5) * step;
            }
            return bins;
        } else {
            std::vector<double> bins(N + 1);
            const auto step = (stop - start) / N;
            for (std::size_t i = 0; i < N + 1; ++i) {
                bins[i] = start + i * step;
            }
            return bins;
        }
    }
    std::vector<double> getValues() const
    {
        std::vector<double> vals(N);
        const auto sum = static_cast<double>(std::accumulate(intensity.cbegin(), intensity.cend(), std::size_t { 0 }));
        std::transform(intensity.cbegin(), intensity.cend(), vals.begin(), [=](auto i) { return i / sum; });
        return vals;
    }
    void save(const std::string& fname) const
    {
        std::ofstream myfile;
        myfile.open(fname, std::ios::out);
        auto bins = getBins();
        auto vals = getValues();
        for (std::size_t i = 0; i < vals.size(); ++i) {
            myfile << bins[i] << "," << vals[i] << '\n';
        }
        myfile.close();
    }
};

class ResultPrint {
public:
    ResultPrint()
    {
        m_myfile.open("validationScatterTable.txt", std::ios::out);
    }
    ResultPrint(const std::string& fname)
    {
        m_myfile.open(fname, std::ios::out);
    }

    ~ResultPrint()
    {
        m_myfile.close();
    }
    void header()
    {
        m_myfile << "Model,InteractionType,x,y\n";
    }
    void operator()(const std::string& model, const std::string& type, double x, double y)
    {
        m_myfile << model << ",";
        m_myfile << type << ",";
        m_myfile << x << ",";
        m_myfile << y << '\n';
    }
    void operator()(const Histogram& hist, const std::string& model, const std::string& type)
    {
        const auto b = hist.getBins();
        const auto v = hist.getValues();
        for (std::size_t i = 0; i < v.size(); ++i) {
            m_myfile << model << ',';
            m_myfile << type << ',';
            m_myfile << b[i] << ',';
            m_myfile << v[i] << '\n';
        }
    }

private:
    std::ofstream m_myfile;
};

template <int Model = 1>
Histogram comptonScatterAngle(const dxmc::Material<5>& material, double energy, std::size_t N = 1e6)
{
    dxmc::RandomState state;
    Histogram h(180, 0, 180);
    for (std::size_t i = 0; i < N; ++i) {
        dxmc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        dxmc::interactions::comptonScatter<5, Model>(p, material, state);
        const auto angle = dxmc::vectormath::angleBetween({ 0, 0, 1 }, p.dir) * dxmc::RAD_TO_DEG<>();
        h(angle);
    }
    return h;
}

template <int Model = 1>
Histogram comptonScatterEnergy(const dxmc::Material<5>& material, double energy, std::size_t N = 1e6)
{
    dxmc::RandomState state;
    Histogram h(100, .9 / (1 + energy / dxmc::ELECTRON_REST_MASS() * 2), 1);
    for (std::size_t i = 0; i < N; ++i) {
        dxmc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        dxmc::interactions::comptonScatter<5, Model>(p, material, state);
        h(p.energy / energy);
    }
    return h;
}

template <int Model = 1>
Histogram rayleightScatterAngle(const dxmc::Material<5>& material, double energy, std::size_t N = 1e6)
{
    dxmc::RandomState state;
    Histogram h(180, 0, 180);
    for (std::size_t i = 0; i < N; ++i) {
        dxmc::Particle p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 }, .energy = energy, .weight = 1 };
        dxmc::interactions::rayleightScatter<5, Model>(p, material, state);
        const auto angle = dxmc::vectormath::angleBetween({ 0, 0, 1 }, p.dir) * dxmc::RAD_TO_DEG<>();
        h(angle);
    }
    return h;
}

template <int I = 0>
void saveHist(ResultPrint& p, const auto& material, double energy)
{
    auto h_ang = comptonScatterAngle<I>(material, energy, 1000000);
    auto h_en = comptonScatterEnergy<I>(material, energy, 1000000);
    auto hr_ang = rayleightScatterAngle<I>(material, energy, 1000000);
    std::string model = "None";
    if (I == 1)
        model = "Livermore";
    if (I == 2)
        model = "IA";

    p(h_ang, model, "ComptonAngle");
    p(h_en, model, "ComptonEnergy");
    p(hr_ang, model, "RayleighAngle");
}

int main()
{
    auto material = dxmc::Material<5>::byNistName("Water, Liquid").value();

    ResultPrint p;
    p.header();
    saveHist<0>(p, material, 60.0);
    saveHist<1>(p, material, 60.0);
    saveHist<2>(p, material, 60.0);

    return 0;
}