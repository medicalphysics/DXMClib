# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 23:48:17 2022

@author: ander
"""

import pandas as pd
import seaborn as sns
from matplotlib import pylab as plt

def testAtom():
    path = r"atomTestData.csv"    
    df = pd.read_csv(path,sep=',')
    
    xmin = df['e'].min()
    xmax = df['e'].max()
    
    ## Edit these to view only a small energy segment
    # xmin = 5.5
    # xmax = 6.5
    
    for t in ['scatterfactor', 'formfactor','photoelectric', 'coherent', 'incoherent', 'total',]:
        idx = df['type'] == t
        dff = df[idx]
        _ = plt.figure(figsize=(9, 3), dpi=300)
        plt.title(t)
        ax = plt.subplot(131)
        hue_att = ['dxmc', 'lin', 'xlib']        
        sns.lineplot(data=dff, x='e', y='att', hue='kind', hue_order=hue_att, ax=ax)
        plt.yscale('log')
        plt.xscale('log')
        plt.xlim(xmin, xmax)
        ax = plt.subplot(132)
        sns.lineplot(data=dff, x='e', y='att', hue='kind', hue_order=['diff',], ax=ax)
        plt.xscale('log')
        plt.xlim(xmin, xmax)
        ax = plt.subplot(133)
        sns.lineplot(data=dff, x='e', y='att', hue='kind', hue_order=['diffp',], ax=ax)
        plt.xscale('log')
        plt.xlim(xmin, xmax)   
        
        plt.tight_layout()
        plt.show()
        
    sns.relplot(data=df, x='e', y='att', hue='type', col='kind', hue_order=['photoelectric', 'coherent', 'incoherent', 'total'], 
                col_order=['dxmc', 'lin', 'xlib'], kind='line')
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim(xmin,xmax)    
    plt.tight_layout()
    plt.show()    


if __name__=='__main__':    
    testAtom()
        
    