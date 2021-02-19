# This file is part of DXMClib.

# DXMClib is free software : you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# DXMClib is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

# Copyright 2021 Erlend Andersen


import seaborn as sns
import pandas as pd
from matplotlib import pylab as plt
import numpy as np

def readData():    
    converters={"Result":float, "Stddev":float}
    dt = pd.read_csv("validationTable.txt", sep=", ", engine='python', converters=converters, error_bad_lines=False)
    
    return dt


def plotCase2(dt_full):
    dt=dt_full[dt_full['Case']=="Case 2"]
    dt=dt[dt['Volume'] != 'Total body']
    
    sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt )
    plt.savefig("Case2.png", dpi=900)
    plt.show()
    

def plotCase3(dt_full):
    dt=dt_full[dt_full['Case']=="Case 3"]
    dt=dt[dt['Volume'] != 'Total body']
    
    sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt )
    plt.savefig("Case3.png", dpi=900)
    plt.show()
    
def plotCase41(dt_full):
    dt=dt_full[dt_full['Case']=="Case 4.1"]        
    sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt )
    plt.savefig("Case41.png", dpi=900)
    plt.show()
    
def plotCase42(dt_full):
    dt=dt_full[dt_full['Case']=="Case 4.2"]
    
    dtp_ind=['Cent' in e for e in dt['Mode']]
    dtp = dt[dtp_ind]
    sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dtp )
    plt.savefig("Case42_Cent.png", dpi=900)
    plt.show()
    
    dtp_ind=['Pher' in e for e in dt['Mode']]
    dtp = dt[dtp_ind]
    sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dtp )
    plt.savefig("Case42_Pher.png", dpi=900)
    plt.show()
    
def plotCase5(dt_full):
    dt=dt_full[dt_full['Case']=="Case 5"]        
    vols = set([m for m in dt['Volume']])
    for m in vols:
        dtm = dt[dt['Volume']==m]
        sns.catplot(x="Mode", y="Result",hue='Model', row=None, col='Specter', data=dtm )
        plt.savefig("Case5_{}.png".format(m), dpi=900)
        plt.show()

    


if __name__=='__main__':
    dt = readData()
    plotCase2(dt)
    plotCase3(dt)
    plotCase41(dt)
    plotCase42(dt)
    plotCase5(dt)

