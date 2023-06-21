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

#include "lodepng/dxmclodepngwrapper.hpp"
#include "lodepng.h"

namespace dxmc {
namespace dxmclodepng {
    std::vector<std::uint8_t> encodePNG(const std::vector<std::uint8_t>& image, std::size_t width, size_t height)
    {

        std::vector<unsigned char> png;

        if (image.size() == width * height * 4) {
            auto error = lodepng::encode(png, image, width, height);
        }
        return png;
    }
};
}