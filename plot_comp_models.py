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
plt.close('all')
models = ['HW', 'KC', 'B89', 'B92']
extension = '.png'

csvs = []
pprojD = []
pmassproj = []
pmassD = []
pvD = []
names = []


def powerLaw(x, a, b):
    return a*x**b

def linearLaw(x, a, b):
    return a+b*x

def atlasLaw(x, a, b, c):
    return a-b*np.exp(-c*x)

def fitPower(X, Y):
    popt, pcov = curve_fit(powerLaw, X, Y)
    return popt, pcov

def fitLinearLog(X, Y):
    popt, pcov = curve_fit(linearLaw, np.log10(X), np.log10(Y))
    popt[0] = 10.0**popt[0]
    return popt, pcov

def fitAtlas(X, Y):
    popt, pcov = curve_fit(atlasLaw, X, Y, p0=[1.5, 1.4, 1000], bounds=([0.5, 0.5, 10], [8.5, 2.5, 10000]))
    return popt, pcov

datafiles = sorted(glob('./*B92.csv'))
badcol='Unnamed: 0'

indexes = [d[2:-7] for d in datafiles]
coeffs = ['a','b','c','d','e']
columns = [b+a for a in models for b in coeffs]
data_coeffs = pd.DataFrame(index=indexes, columns=columns)

for csvfiles in sorted(datafiles):
    csvfiles = csvfiles[:-7] # cut model and extension
    particles = csvfiles[2:]
    fig, ax = plt.subplots(1,len(models),figsize=(18,5))
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
        #popt, pcov = fitPower(data['Jdmax'], data['vDproj'])
        v = powerLaw(xdata*1e-3,*popt)
        ax[i].plot(xdata,v,c='k',linewidth=2,label='%1.3fx$^{%1.1f}$' % tuple(popt))
        data_coeffs.loc[particles, 'a'+model]=popt[0]
        data_coeffs.loc[particles, 'b'+model]=popt[1]
        #pvD.append(popt)
        #popt, pcov = fitAtlas(data['Jdmax'], data['vDproj'])
        Deq = np.cbrt(6.0*(data['mass']/917.0)/np.pi)
        popt, pcov = fitAtlas(Deq, data['vDproj'])
        #v = atlasLaw(xdata*1e-3, *popt)
        v = atlasLaw(Deq, *popt)
        data_coeffs.loc[particles, 'c'+model]=popt[0]
        data_coeffs.loc[particles, 'd'+model]=popt[1]
        data_coeffs.loc[particles, 'e'+model]=popt[2]
        data['vcalc'] = v
        dat = data.sort_values(by='Jdmax')
        #ax[i].plot(xdata, v, c='r',linewidth=2,label='%1.3f-%1.3f exp(-%1.1fx)$' % tuple(popt))
        ax[i].plot(dat['Jdmax']*1.0e3, dat['vcalc'], c='r',linewidth=2,label='%1.3f-%1.3f exp(-%1.1fx)' % tuple(popt))
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
    DF = pd.DataFrame()
data_coeffs.to_csv('velocity_fits_coefficients.csv')