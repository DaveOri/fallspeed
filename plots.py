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

datafiles = glob('./*.csv')
badcol='Unnamed: 0'

pprojD = []
pmassproj = []
pmassD = []
pvD = []
names = []

for csvfile in sorted(datafiles):
    data = pd.read_csv(csvfile)
    if badcol in data.columns:
        data.drop(badcol,axis=1,inplace=True)
    data.sort_values(by='Jdmax',inplace=True)

    def powerLaw(x,a,b):
        return a*x**b
    
    def linearLaw(x,a,b):
        return a+b*x
    
    fig, ax = plt.subplots(2,3,figsize=(10,7))
    #ax[0,0].scatter(1e3*data['Jdmax'],1e3*data['Ddmax'])
    ax[0,0].scatter(1e3*data['Jdmax'],1e3*data['Dprojcirc'],linewidth=0.2)
    ax[0,0].set_xlim([0,35])
    ax[0,0].set_ylim([0,35])
    xylim = [0,1e3*max(max(data.Jdmax),max(data.Dprojcirc))]
    ax[0,0].plot(xylim, xylim, ls="-", c='r')
    ax[0,0].set_xlim(xylim)
    ax[0,0].set_ylim(xylim)
    ax[0,0].set_ylabel('D circ proj  [mm]')
    ax[0,0].set_xlabel('D$_{max}$ 3D  [mm]')
    
    ax[0,1].scatter(1e3*data['Jdmax'],1e6*data['projarea'],linewidth=0.2)
    xlim = [0,1e3*max(data.Jdmax)]
    ylim = [0,1e6*max(data.projarea)]
    xdata = np.linspace(xlim[0],xlim[1],100)
    popt, pcov = curve_fit(powerLaw, data['Jdmax'], data['projarea'])
    ax[0,1].plot(xdata,1e6*powerLaw(xdata*1e-3,*popt),c='g',linewidth=2,label='%1.3fx$^{%1.1f}$' % tuple(popt))
    pprojD.append(popt)
    ax[0,1].legend(loc=4)
    ax[0,1].set_xlim(xlim)
    ax[0,1].set_ylim(ylim)
    ax[0,1].set_ylabel('proj area  [mm$^2$]')
    ax[0,1].set_xlabel('D$_{max}$ 3D  [mm]')
    
    ax[0,2].scatter(1e6*data['projarea'],1e6*data['mass'],linewidth=0.2)
    xlim = [0,1e6*max(data.projarea)]
    ylim = [0,1e6*max(data.mass)]
    ax[0,2].set_xlim(xlim)
    ax[0,2].set_ylim(ylim)
    xdata = np.linspace(xlim[0],xlim[1],100)
    popt, pcov = curve_fit(linearLaw, data['projarea'], data['mass'])
    ax[0,2].plot(xdata,1e6*linearLaw(xdata*1e-6,*popt),c='g',linewidth=2,label='%1.1g + %1.1f*x' % tuple(popt))
    pmassproj.append(popt)
    ax[0,2].legend(loc=4)
    ax[0,2].set_xlabel('proj area  [mm$^2$]')
    ax[0,2].set_ylabel('mass     [mg]')
    
    ax[1,0].scatter(1e3*data['Jdmax'],1e6*data['mass'],linewidth=0.2)
    xlim = [1e3*min(data.Jdmax),1e3*max(data.Jdmax)]
    ylim = [1e6*min(data.mass),1e6*max(data.mass)]
    ax[1,0].set_xlim(xlim)
    ax[1,0].set_ylim(ylim)
    xdata = np.linspace(xlim[0],xlim[1],100)
    popt, pcov = curve_fit(powerLaw, data['Jdmax'], data['mass'])
    ax[1,0].plot(xdata,1e6*powerLaw(xdata*1e-3,*popt),c='g',linewidth=2,label='%1.1gx$^{%1.1f}$' % tuple(popt))
    pmassD.append(popt)    
    ax[1,0].legend(loc=4)
    ax[1,0].set_yscale('log')
    ax[1,0].set_xscale('log')
    ax[1,0].set_ylabel('mass   [mg]')
    ax[1,0].set_xlabel('D$_{max}$ 3D  [mm]')
    
    ax[1,1].scatter(1e3*data['Jdmax'],data['vDproj'],linewidth=0.2)
    xlim = [0,1e3*max(data.Jdmax)]
    ylim = [0,max(data.vDproj)]
    ax[1,1].set_xlim(xlim)
    ax[1,1].set_ylim(ylim)
    xdata = np.linspace(xlim[0],xlim[1],100)
    popt, pcov = curve_fit(powerLaw, data['Jdmax'], data['vDproj'])
    ax[1,1].plot(xdata,powerLaw(xdata*1e-3,*popt),c='g',linewidth=2,label='%1.3fx$^{%1.1f}$' % tuple(popt))
    pvD.append(popt)    
    ax[1,1].legend(loc=4)
    ax[1,1].set_ylabel('term speed  [m/s]')
    ax[1,1].set_xlabel('D$_{max}$ 3D  [mm]')
    
    ax[1,2].scatter(data['vDproj'],data['vDJ'],linewidth=0.2)
    ax[1,2].set_ylabel('speed D$_{max}$  [m/s]')
    ax[1,2].set_xlabel('speed D$_{proj}$  [m/s]')
    xylim = [0,max(max(data.vDJ),max(data.vDproj))]
    ax[1,2].plot(xylim, xylim, ls="-", c='r')
    ax[1,2].set_xlim(xylim)
    ax[1,2].set_ylim(xylim)
    fig.suptitle('s'+csvfile[2:-4],fontweight='heavy')
    fig.tight_layout()
    fig.savefig('s'+csvfile[2:-4]+'.png')
    names.append(csvfile[2:-4])
#%%
xdata = np.linspace(0,35,60)
plt.figure(figsize=(8,6))
for i,j in zip(pvD,names):
    if 'multaneous' in j:
        plt.plot(xdata,powerLaw(xdata*1e-3,*i),label=j[-3:])
    else:
        plt.plot(xdata,powerLaw(xdata*1e-3,*i),'h-',label=j[-3:])
plt.legend()
plt.title('simultaneous-lines         subsequent-markers')
plt.savefig('all_fitted_vD.png')