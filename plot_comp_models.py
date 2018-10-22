# -*- coding: utf-8 -*-
"""
Created on Mon May 21 08:45:37 2018

@author: dori
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from glob import glob
import matplotlib.pyplot as plt

models = ['HW', 'KC']
extension = '.png'

pprojD = []
pmassproj = []
pmassD = []
pvD = []
names = []


def powerLaw(x,a,b):
    return a*x**b

def linearLaw(x,a,b):
    return a+b*x

def fitPower(X,Y):
    popt, pcov = curve_fit(powerLaw, X, Y)
    return popt, pcov

def fitLinearLog(X,Y):
    popt, pcov = curve_fit(linearLaw, np.log10(X), np.log10(Y))
    popt[0] = 10.0**popt[0]
    return popt, pcov

datafiles = sorted(glob('./*HW.csv'))
badcol='Unnamed: 0'

for csvfiles in sorted(datafiles):
    csvfiles = csvfiles[:-6] # cut model and extension
    fig, ax = plt.subplots(1,2,figsize=(10,5))
    for i, model in enumerate(models):
        csvfile = csvfiles + model + '.csv'
        print(csvfile)
        data = pd.read_csv(csvfile)
        if badcol in data.columns:
            data.drop(badcol,axis=1,inplace=True)
        data.sort_values(by='Jdmax',inplace=True)
    
        ax[i].scatter(1e3*data['Jdmax'],data['vDproj'],linewidth=0.2)
        xlim = [0,1e3*max(data.Jdmax)]
        ylim = [0,max(data.vDproj)]
        ax[i].set_xlim(xlim)
        ax[i].set_ylim(ylim)
        xdata = np.linspace(xlim[0],xlim[1],100)
        popt, pcov = fitLinearLog(data['Jdmax'], data['vDproj'])
        v = powerLaw(xdata*1e-3,*popt)
        ax[i].plot(xdata,v,c='g',linewidth=2,label='%1.3fx$^{%1.1f}$' % tuple(popt))
        pvD.append(popt)
        if max(data['vDproj'])>3.5:
            for j in range(len(ax)):
                ax[j].set_ylim([0.0,7.0])
        else:
            for j in range(len(ax)):
                ax[j].set_ylim([0.0,4.0])
        ax[i].legend(loc=4)
        ax[i].set_ylabel('term speed  [m/s]')
        ax[i].set_xlabel('D$_{max}$ 3D  [mm]')
        ax[i].set_title('vel model '+model)
        ax[i].grid()
    fig.suptitle(csvfiles,fontweight='heavy')
    fig.tight_layout()
    plt.savefig(csvfiles+'_comp_models'+extension)
    names.append(csvfile[:-4])

