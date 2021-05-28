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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.hpp"

#include <array>
#include <atomic>
#include <chrono>
#include <execution>
#include <iomanip>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace dxmc {
template <Floating T = double>
struct DoseProgressImageData {
    std::array<std::size_t, 2> dimensions = { 0, 0 };
    std::array<T, 2> spacing = { 0, 0 };
    std::vector<std::uint8_t> image;
};

template <Floating T>
class ProgressBar {
public:
    enum class Axis { X,
        Y,
        Z };
    ProgressBar() {}
    ProgressBar(std::uint64_t totalExposures) { setTotalExposures(totalExposures); }
    void setTotalExposures(std::uint64_t totalExposures)
    {
        m_totalExposures = totalExposures;
        m_currentExposures = 0;
        m_startTime = std::chrono::system_clock::now();
    } //not thread safe

    void setPrefixMessage(const std::string& msg) { m_message = msg; } // not threadsafe

    void exposureCompleted() // threadsafe
    {
        m_currentExposures++;
        auto duration = std::chrono::system_clock::now() - m_startTime;
        auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
        m_secondsElapsed.exchange(static_cast<T>(seconds));
    }
    std::string getETA() const
    {
        const auto cExp = m_currentExposures.load();
        if (cExp == 0)
            return "ETA: estimating... [0%]";
        const auto secondsElapsed = m_secondsElapsed.load();
        const auto totExp = m_totalExposures.load();
        const auto secondsRemaining = secondsElapsed / cExp * (totExp - cExp);
        const auto percent = (T { 100 } * cExp) / totExp;
        return makePrettyTime(secondsRemaining, percent);
    }
    void setCancel(bool cancel) // threadsafe
    {
        m_cancel.exchange(cancel);
    }
    bool cancel(void) const // threadsafe
    {
        return m_cancel.load();
    }
    void setPlaneNormal(Axis planeNormal) { m_doseAxis = planeNormal; }
    void setDoseData(const T* doseData, const std::array<std::size_t, 3>& doseDimensions, const std::array<T, 3>& doseSpacing)
    {
        std::scoped_lock guard(m_doseMutex);
        m_doseData = doseData;
        m_doseDimensions = doseDimensions;
        m_doseSpacing = doseSpacing;
    }
    void clearDoseData()
    {
        std::scoped_lock guard(m_doseMutex);
        m_doseData = nullptr;
        std::fill(m_doseSpacing.begin(), m_doseSpacing.end(), T { 0.0 });
        std::fill(m_doseDimensions.begin(), m_doseDimensions.end(), 0);
    }

    std::shared_ptr<DoseProgressImageData<T>> computeDoseProgressImage()
    {
        std::scoped_lock guard(m_doseMutex);
        if (m_doseAxis == Axis::Y)
            return computeDoseProgressImageY();
        else if (m_doseAxis == Axis::Z)
            return computeDoseProgressImageZ();
        return computeDoseProgressImageX();
    }

protected:
    std::string makePrettyTime(T seconds, T percent) const
    {
        std::stringstream ss;
        ss << std::fixed << std::setprecision(0);
        ss << m_message << "ETA: about ";
        if (seconds > 120)
            ss << seconds / 60 << " minutes";
        else
            ss << seconds << " seconds";
        ss << " [" << percent << "%]";
        return ss.str();
    }

    std::shared_ptr<DoseProgressImageData<T>> computeDoseProgressImageX()
    {
        if (!m_doseData)
            return nullptr;

        auto doseProgressImage = std::make_shared<DoseProgressImageData<T>>();
        m_doseMipBuffer.resize(m_doseDimensions[1] * m_doseDimensions[2]);
        std::fill(m_doseMipBuffer.begin(), m_doseMipBuffer.end(), 0.0);
        doseProgressImage->image.resize(m_doseDimensions[1] * m_doseDimensions[2], 0);
        doseProgressImage->dimensions[0] = m_doseDimensions[1];
        doseProgressImage->dimensions[1] = m_doseDimensions[2];
        doseProgressImage->spacing[0] = m_doseSpacing[1];
        doseProgressImage->spacing[1] = m_doseSpacing[2];

        //doing mip over X axis
        T global_max = 0.0;
        for (std::size_t i = 0; i < m_doseDimensions[1] * m_doseDimensions[2]; ++i) {
            const auto beg = m_doseData + i * m_doseDimensions[0];
            const auto end = beg + m_doseDimensions[0];
            const auto max_val = *std::max_element(beg, end);
            m_doseMipBuffer[i] = max_val;
            global_max = std::max(global_max, max_val);
        }

        const T scalefactor = T { 255.0 } / global_max;

        std::transform(std::execution::par_unseq, m_doseMipBuffer.cbegin(), m_doseMipBuffer.cend(), doseProgressImage->image.begin(),
            [=](const T el) -> std::uint8_t { return static_cast<std::uint8_t>(el * scalefactor); });

        return doseProgressImage;
    }

    std::shared_ptr<DoseProgressImageData<T>> computeDoseProgressImageY()
    {
        if (!m_doseData)
            return nullptr;

        auto doseProgressImage = std::make_shared<DoseProgressImageData<T>>();
        m_doseMipBuffer.resize(m_doseDimensions[0] * m_doseDimensions[2]);
        std::fill(m_doseMipBuffer.begin(), m_doseMipBuffer.end(), 0.0);
        doseProgressImage->image.resize(m_doseDimensions[0] * m_doseDimensions[2], 0);
        doseProgressImage->dimensions[0] = m_doseDimensions[0];
        doseProgressImage->dimensions[1] = m_doseDimensions[2];
        doseProgressImage->spacing[0] = m_doseSpacing[0];
        doseProgressImage->spacing[1] = m_doseSpacing[2];

        //doing mip over Y axis
        T global_max = 0.0;
        for (std::size_t k = 0; k < m_doseDimensions[2]; ++k)
            for (std::size_t j = 0; j < m_doseDimensions[1]; ++j)
                for (std::size_t i = 0; i < m_doseDimensions[0]; ++i) {
                    const auto dp_idx = i + m_doseDimensions[0] * k;
                    const auto d_idx = i + m_doseDimensions[0] * j + m_doseDimensions[0] * m_doseDimensions[1] * k;
                    const auto max_val = std::max(m_doseData[d_idx], m_doseMipBuffer[dp_idx]);
                    m_doseMipBuffer[dp_idx] = max_val;
                    global_max = std::max(max_val, global_max);
                }
        const T scalefactor = 255.0 / global_max;

        std::transform(std::execution::par_unseq, m_doseMipBuffer.cbegin(), m_doseMipBuffer.cend(), doseProgressImage->image.begin(),
            [=](const T el) -> std::uint8_t { return static_cast<std::uint8_t>(el * scalefactor); });

        return doseProgressImage;
    }

    std::shared_ptr<DoseProgressImageData<T>> computeDoseProgressImageZ()
    {
        if (!m_doseData)
            return nullptr;

        auto doseProgressImage = std::make_shared<DoseProgressImageData<T>>();
        m_doseMipBuffer.resize(m_doseDimensions[0] * m_doseDimensions[1]);
        std::fill(m_doseMipBuffer.begin(), m_doseMipBuffer.end(), 0.0);
        doseProgressImage->image.resize(m_doseDimensions[0] * m_doseDimensions[1], 0);
        doseProgressImage->dimensions[0] = m_doseDimensions[0];
        doseProgressImage->dimensions[1] = m_doseDimensions[1];
        doseProgressImage->spacing[0] = m_doseSpacing[0];
        doseProgressImage->spacing[1] = m_doseSpacing[1];

        //doing mip over Z axis
        T global_max = 0.0;
        for (std::size_t k = 0; k < m_doseDimensions[2]; ++k)
            for (std::size_t j = 0; j < m_doseDimensions[1]; ++j)
                for (std::size_t i = 0; i < m_doseDimensions[0]; ++i) {
                    const auto dp_idx = i + m_doseDimensions[0] * j;
                    const auto d_idx = i + m_doseDimensions[0] * j + m_doseDimensions[0] * m_doseDimensions[1] * k;
                    const auto max_val = std::max(m_doseData[d_idx], m_doseMipBuffer[dp_idx]);
                    m_doseMipBuffer[dp_idx] = max_val;
                    global_max = std::max(global_max, max_val);
                }
        const T scalefactor = 255.0 / global_max;

        std::transform(std::execution::par_unseq, m_doseMipBuffer.cbegin(), m_doseMipBuffer.cend(), doseProgressImage->image.begin(),
            [=](const T el) -> std::uint8_t { return static_cast<std::uint8_t>(el * scalefactor); });
        return doseProgressImage;
    }

private:
    std::atomic<std::uint64_t> m_totalExposures = 0;
    std::atomic<std::uint64_t> m_currentExposures = 0;
    std::chrono::system_clock::time_point m_startTime;
    std::atomic<T> m_secondsElapsed;
    std::string m_message;
    std::atomic<bool> m_cancel = false;
    std::mutex m_doseMutex;
    const T* m_doseData = nullptr;
    std::vector<T> m_doseMipBuffer;
    std::array<std::size_t, 3> m_doseDimensions = { 0, 0, 0 };
    std::array<T, 3> m_doseSpacing = { 1, 1, 1 };
    Axis m_doseAxis = Axis::Y;
};
}