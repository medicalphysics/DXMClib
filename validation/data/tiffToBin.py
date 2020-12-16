"""
This file is part of DXMClib.

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

Copyright 2020 Erlend Andersen
"""


import numpy as np
from skimage import io

#Path to tiff file from report TG195 case 5 describing the phantom
TIFF_FILE = r"./TG 195 Case 5 Voxelized Volume.tif"
TIFF_FILE=r"C:\Users\erlend\OneDrive\OpenDXMCvalidering\RPT_195ElectronicResources\05-Case 5 Voxelized Volume\TG 195 Case 5 Voxelized Volume.tif"

#Procedure to convert the tiff file to binary file that can be easily read by c++
def convert_tiff(path):
    im = io.imread(path)
    # import pdb;pdb.set_trace()
    flat = im.flatten().astype(np.uint8)
    flat.tofile("case5world.bin")
    return 
   
if __name__ == '__main__':
    convert_tiff(TIFF_FILE)
