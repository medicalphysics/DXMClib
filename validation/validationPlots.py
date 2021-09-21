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
    nNaN= dt.isnull().sum().sum()
    if nNaN > 0:
        print("Warning: Found {} NaN values or missing data, dropping rows".format(nNaN))
        return dt.dropna()
    return dt

def fix_axis(fg, rotate_labels=True):
    fg.set(ylim=(0, None))
    fg.set_axis_labels(y_var="Energy [eV/history]")
    if rotate_labels:
        for ax in fg.axes_dict.values():             
            for tick in ax.get_xticklabels():
                tick.set_rotation(70)
    fg.tight_layout()

def plotCase2(dt_full):
    dt=dt_full[dt_full['Case']=="Case 2"]    
    dt_vol=dt[dt['Volume'] != 'Total body']
    for forced in [0, 1]:
        dt_f=dt_vol[dt_vol['Forced']==forced]
        g = sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt_f )
        fix_axis(g)
        plt.savefig("Case2_volume{}.png".format("_forced" if forced==1 else ""), dpi=900)
        plt.show()
        plt.clf()

    dt_tot = dt[dt['Volume'] == "Total body"]
    g = sns.catplot(x="Forced", y="Result",hue='Model', row='Specter', col='Mode', data=dt_tot )
    fix_axis(g, False)
    plt.savefig("Case2_total_body.png", dpi=900)
    plt.show()
    plt.clf()
 
def plotCase3(dt_full):
    dt=dt_full[dt_full['Case']=="Case 3"]    
    dt_vol=dt[dt['Volume'] != 'Total body']
    for forced in [0, 1]:
        dt_f=dt_vol[dt_vol['Forced']==forced]
        g = sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt_f )
        plt.savefig("Case3_volume{}.png".format("_forced" if forced==1 else ""), dpi=900)
        fix_axis(g)       
        plt.show()
        plt.clf()

    dt_tot = dt[dt['Volume'] == "Total body"]
    g = sns.catplot(x="Forced", y="Result",hue='Model', row='Specter', col='Mode', data=dt_tot )
    fix_axis(g, False)
    plt.savefig("Case3_total_body.png", dpi=900)
    plt.show()
    plt.clf()

def plotCase41(dt_full):
    dt=dt_full[dt_full['Case']=="Case 4.1"]        
    g = sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt )
    fix_axis(g, False)
    plt.savefig("Case41.png", dpi=900)
    plt.show()
    
def plotCase42(dt_full):
    dt_42=dt_full[dt_full['Case']=="Case 4.2"]
    for forced in [1, 0]:
        dt = dt_42[dt_42['Forced'] == forced]
        dtp_ind=['Cent' in e for e in dt['Mode']]
        dtp = dt[dtp_ind]
        g = sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dtp )
        fix_axis(g)
        plt.savefig("Case42_Cent{}.png".format("_forced" if forced==1 else ""), dpi=900)
        plt.show()
    
        dtp_ind=['Pher' in e for e in dt['Mode']]
        dtp = dt[dtp_ind]
        g = sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dtp )
        fix_axis(g)
        plt.savefig("Case42_Pher{}.png".format("_forced" if forced==1 else ""), dpi=900)
        plt.show()
    
def plotCase5(dt_full):
    dt=dt_full[dt_full['Case']=="Case 5"]        
    vols = set([m for m in dt['Volume']])
    for m in vols:
        dtm = dt[dt['Volume']==m]
        g = sns.catplot(x="Mode", y="Result",hue='Model', row=None, col='Specter', data=dtm )
        fix_axis(g)
        plt.savefig("Case5_{}.png".format(m), dpi=900)
        plt.show()

if __name__=='__main__':
    dt = readData()
    plotCase2(dt)
    plotCase3(dt)
    plotCase41(dt)
    plotCase42(dt)
    plotCase5(dt)
