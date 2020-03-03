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
#include <atomic>
#include <chrono>
#include <string>
#include <sstream>
#include <iomanip> 
#include <mutex>
#include <execution>
#include <array>
#include <vector>
#include <memory>

struct DoseProgressImageData
{
	std::array<std::size_t, 2> dimensions = { 0,0 };
	std::vector<unsigned char> image;
};


class ProgressBar
{
public:
	ProgressBar() {};
	ProgressBar(std::uint64_t totalExposures) { setTotalExposures(totalExposures); }
	void setTotalExposures(std::uint64_t totalExposures, const std::string& message = "")
	{
		m_totalExposures = totalExposures;
		m_currentExposures = 0;
		m_message = message;
		m_startTime = std::chrono::system_clock::now();
	} //not thread safe
	void exposureCompleted() // threadsafe
	{
		m_currentExposures++;
		auto duration = std::chrono::system_clock::now() - m_startTime;
		auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();
		m_secondsElapsed.exchange(static_cast<double>(seconds));
	}
	std::string getETA() const
	{
		auto secondsRemaining = m_secondsElapsed.load() / m_currentExposures.load() * (m_totalExposures.load() - m_currentExposures.load());
		return makePrettyTime(secondsRemaining);
	}
	void setCancel(bool cancel) // threadsafe
	{
		m_cancel.exchange(cancel);
	}
	bool cancel(void) const // threadsafe
	{
		return m_cancel.load();
	}

	
	void setDoseData(double* doseData, const std::array<std::size_t, 3>& doseDimensions)
	{
		const std::lock_guard<std::mutex> lock(m_doseMutex);
		m_doseData = doseData;
		m_doseDimensions = doseDimensions;
	}

	
	std::unique_ptr<DoseProgressImageData> computeDoseProgressImage()
	{
		const std::lock_guard<std::mutex> lock(m_doseMutex);
		
		std::unique_ptr<DoseProgressImageData> doseProgressImage = nullptr;
		if (!m_doseData)
		{
			return doseProgressImage;
		}
		doseProgressImage = std::make_unique<DoseProgressImageData>();
		std::vector<double> doseProgressData(m_doseDimensions[0] * m_doseDimensions[2], 0.0);
		doseProgressImage->image.resize(m_doseDimensions[0] * m_doseDimensions[2], 0);
		doseProgressImage->dimensions[0] = m_doseDimensions[0];
		doseProgressImage->dimensions[1] = m_doseDimensions[2];

		//doing mip over Y axis
		
		for (std::size_t k = 0; k < m_doseDimensions[2]; ++k)
			for (std::size_t j = 0; j < m_doseDimensions[1]; ++j)
				for (std::size_t i = 0; i < m_doseDimensions[0]; ++i)
				{
					const auto dp_idx = i + m_doseDimensions[0] * k;
					const auto d_idx = i + m_doseDimensions[0] * j + m_doseDimensions[0] * m_doseDimensions[1] * k;
					doseProgressData[dp_idx] = std::max(m_doseData[d_idx], doseProgressData[dp_idx]);
				}
		const auto max_element = std::max_element(std::execution::par_unseq, doseProgressData.begin(), doseProgressData.end());
		const double scalefactor = 255.0 / *max_element;

		std::transform(std::execution::par_unseq, doseProgressData.cbegin(), doseProgressData.cend(), doseProgressImage->image.begin(),
			[=](const double el)->unsigned char {return static_cast<unsigned char>(el * scalefactor); }
		);

		return doseProgressImage; 
	}


protected:
	std::string makePrettyTime(double seconds) const {
		std::stringstream ss;
		ss << std::fixed << std::setprecision(0);
		ss << m_message << " ETA: about ";
		if (seconds > 120.0)
			ss << seconds / 60.0 << " minutes";
		else
			ss << seconds << " seconds";
		return ss.str();
	}

private:
	std::atomic<std::uint64_t> m_totalExposures = 0;
	std::atomic<std::uint64_t> m_currentExposures = 0;
	std::chrono::system_clock::time_point m_startTime;
	std::atomic<double> m_secondsElapsed;
	std::string m_message;
	std::atomic<bool> m_cancel = false;
	std::mutex m_doseMutex;
	double* m_doseData = nullptr;
	std::array<std::size_t, 3> m_doseDimensions;

};