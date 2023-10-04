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
import os

def readData():    
    converters={"Result":float, "Stddev":float}
    dt = pd.read_csv("validationTable.txt", sep=", ", engine='python', converters=converters)
    nNaN= dt.isnull().sum().sum()
    if nNaN > 0:
        print("Warning: Found {} NaN values or missing data, dropping rows".format(nNaN))
        return dt.dropna()
    return dt

def fix_axis(fg, rotate_labels=True, ylabel="Energy [eV/history]"):
    fg.set(ylim=(0, None))
    fg.set_axis_labels(y_var=ylabel)
    if rotate_labels:
        for ax in fg.axes_dict.values():             
            for tick in ax.get_xticklabels():
                tick.set_rotation(70)
    fg.tight_layout()

def plotCase2(dt_full, show=False):
    dt=dt_full[dt_full['Case']=="Case 2"]    
    dt_vol=dt[dt['Volume'] != 'Total body']
   
    g = sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt_vol )
    fix_axis(g)
    plt.savefig("plots/Case2_volume.png", dpi=900)
    if show:
        plt.show()
    plt.clf()

    dt_tot = dt[dt['Volume'] == "Total body"]
    g = sns.catplot(x="Specter", y="Result",hue='Model', col='Mode',  data=dt_tot )
    fix_axis(g, False)
    plt.savefig("plots/Case2_total_body.png", dpi=900)
    if show:
        plt.show()
    plt.clf()
 
def plotCase3(dt_full, show=False):
    dt=dt_full[dt_full['Case']=="Case 3"]    
    dt_vol=dt[dt['Volume'] != 'Total body']
   
    g = sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt_vol )
    fix_axis(g)
    plt.savefig("plots/Case3_volume.png", dpi=900)
    if show:
        plt.show()
    plt.clf()

    dt_tot = dt[dt['Volume'] == "Total body"]
    g = sns.catplot(x="Specter", y="Result",hue='Model', col='Mode',  data=dt_tot )
    fix_axis(g, False)
    plt.savefig("plots/Case3_total_body.png", dpi=900)
    if show:
        plt.show()
    plt.clf()

def plotCase41(dt_full, show=False):
    dt=dt_full[dt_full['Case']=="Case 4.1"]        
    g = sns.catplot(x="Volume", y="Result",hue='Model', row='Specter', col='Mode', data=dt )
    fix_axis(g, False)
    plt.savefig("plots/Case41.png", dpi=900)
    if show:
        plt.show()
    
def plotCase42(dt_full, show=False):
    dt=dt_full[dt_full['Case']=="Case 4.2"]
    
    dt['Volume [angle]'] = [float(d) for d in dt['Volume']]
    dtp_ind=['Cent' in e for e in dt['Mode']]
    dtp = dt[dtp_ind]
    g = sns.relplot(x="Volume [angle]", y="Result",hue='Model', row='Specter', col='Mode', data=dtp )
    fix_axis(g)
    plt.savefig("plots/Case42_Cent.png", dpi=900)
    if show:
        plt.show()
    
    dtp_ind=['Pher' in e for e in dt['Mode']]
    dtp = dt[dtp_ind]
    g = sns.relplot(x="Volume [angle]", y="Result",hue='Model', row='Specter', col='Mode', data=dtp )
    fix_axis(g)
    plt.savefig("plots/Case42_Pher.png", dpi=900)
    if show:
        plt.show()
    
def plotCase5(dt_full, show=False):
    dt=dt_full[dt_full['Case']=="Case 5"]        
    vols = set([m for m in dt['Volume']])
    for m in vols:
        dtm = dt[dt['Volume']==m]
        g = sns.catplot(x="Mode", y="Result",hue='Model', row=None, col='Specter', data=dtm )
        fix_axis(g)
        plt.savefig("plots/Case5_{}.png".format(m), dpi=900)
        if show:
            plt.show()

def plotRuntimes(dt, show=False):
    df = dt[dt['Model'] != 'TG195']    
    for cas in set(df['Case']):        
        dff = df[df['Case'] == cas]        
        g = sns.catplot(x="Mode", y="SimulationTime",col='Model', row='Specter', hue='Precision', data=dff, kind='bar')
        fix_axis(g, ylabel='Simulation time [ms]')
        plt.savefig("plots/Runtimes_{}.png".format(cas), dpi=900)
        if show:
            plt.show()
        plt.clf()

if __name__=='__main__':
    dt = readData()
    try:
        os.mkdir("plots")
    except Exception:
        pass
    
    plotCase2(dt)
    plotCase3(dt)
    plotCase41(dt)
    plotCase42(dt)
    plotCase5(dt)
    plotRuntimes(dt)
