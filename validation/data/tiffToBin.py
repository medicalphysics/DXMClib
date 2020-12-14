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
from matplotlib import pylab as plt

TIFF_FILE = r"C:\Users\erlend\OneDrive\OpenDXMCvalidering\RPT_195ElectronicResources\05-Case 5 Voxelized Volume\TG 195 Case 5 Voxelized Volume.tif"


def convert_tiff(path):
    im = io.imread(path)
    # c = np.moveaxis(im, 0, -1)
    plt.imshow(im[:,:,150])
    plt.show()
    # import pdb;pdb.set_trace()
    flat = im.flatten().astype(np.uint8)
    flat.tofile("mat.bin")
    print(im.shape)
    print(flat.min(), flat.max())
    
    
    materials = ["Voxel Value,Material,Density (g/cm3),H,C,N,O,Na,Mg,P,S,Cl,Ar,K,Ca,I,Total,Reference",
        "0,Air,0.001205,0,0.0124,75.5268,23.1781,0,0,0,0,0,1.2827,0,0,0,100,****",
        "1,Cushion/Foam,0.075,7.8,64.7,8.4,19.1,0,0,0,0,0,0,0,0,0,100,*****",
        "2,Carbon fiber,1.2,0,100,0,0,0,0,0,0,0,0,0,0,0,100,",
        "3,Soft tissue,1.03,10.5,25.6,2.7,60.2,0.1,0,0.2,0.3,0.2,0,0.2,0,0,100,**",
        "4,Heart,1.05,10.4,13.9,2.9,71.8,0.1,0,0.2,0.2,0.2,0,0.3,0,0,100,**",
        "5,Lung,0.26,10.3,10.5,3.1,74.9,0.2,0,0.2,0.3,0.3,0,0.2,0,0,100,**",
        "6,Liver,1.06,10.2,13.9,3.0,71.6,0.2,0,0.3,0.3,0.2,0,0.3,0,0,100,**",
        "7,Gallbladder,1.03,10.5,25.6,2.7,60.2,0.1,0,0.2,0.3,0.2,0,0.2,0,0,100,*",
        "8,Spleen,1.06,10.3,11.3,3.2,74.1,0.1,0,0.3,0.2,0.2,0,0.3,0,0,100,**",
        "9,Stomach,1.03,10.6,11.5,2.2,75.1,0.1,0,0.1,0.1,0.2,0,0.1,0,0,100,**",
        "10,Large Intestine,1.03,10.6,11.5,2.2,75.1,0.1,0,0.1,0.1,0.2,0,0.1,0,0,100,**",
        "11,Pancreas,1.04,10.6,16.9,2.2,69.4,0.2,0,0.2,0.1,0.2,0,0.2,0,0,100,**",
        "12,Adrenal,1.03,10.5,25.6,2.7,60.2,0.1,0,0.2,0.3,0.2,0,0.2,0,0,100,*",
        "13,Thyroid,1.05,10.4,11.9,2.4,74.5,0.2,0,0.1,0.1,0.2,0,0.1,0,0.1,100,**",
        "14,Thymus,1.03,10.5,25.6,2.7,60.2,0.1,0,0.2,0.3,0.2,0,0.2,0,0,100,*",
        "15,Small Intestine,1.03,10.6,11.5,2.2,75.1,0.1,0,0.1,0.1,0.2,0,0.1,0,0,100,**",
        "16,Esophagus,1.03,10.6,11.5,2.2,75.1,0.1,0,0.1,0.1,0.2,0,0.1,0,0,100,**",
        "17,Skin,1.09,10.0,20.4,4.2,64.5,0.2,0,0.1,0.2,0.3,0,0.1,0,0,100,**",
        "18,Breast,0.93,11.2,61.9,1.7,25.1,0,0,0.025,0.025,0,0,0.025,0.025,0,100,***",
        "19,Cortical Bone,1.92,3.4,15.5,4.2,43.5,0.1,0.2,10.3,0.3,0,0,0,22.5,0,100,**"]
    weights={'H':1.01,
             'C':12.01,
             'N':14.01,
             'O':16.0,
             'Na':23.0,
             'Mg':24.32,
             'P':30.98,
             'S':32.07,
             'Cl':35.46,
             'Ar':39.94,
             'K':39.1,
             'Ca':40.08,
             'I':126.91}

    header = materials.pop(0).split(',')
    lines =list()
    dens = list()
    for m in materials:
        mat = m.split(',')
        line = mat[0]+";"+mat[1]+";"
        atoms = header[3:]
        mass_frac = str()
        for ind, value in enumerate(mat[3:]):
            atom = atoms[ind]
            if atom in weights:
                weight = weights[atom]
                mass_dens = float(mat[ind+3])
                if mass_dens > 0.0:
                    mass_frac+="{}{}".format(atom, mass_dens*weight)
                
        line += mass_frac +'\n'
        lines.append(line)
        # import pdb;pdb.set_trace()
        dens.append((int(mat[0]), float(mat[2])))
        
    with open("mat.txt", 'w') as f:
        f.writelines(lines)
    
    densarr = np.zeros_like(flat, dtype=np.float32)
    
    for ind, d in dens:
        idx = flat == ind
        densarr[idx] = d
    densarr.tofile("dens.bin")


if __name__ == '__main__':
    convert_tiff(TIFF_FILE)
