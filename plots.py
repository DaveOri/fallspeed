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

model = 'HW'
model = 'KC'
extension = '_'+model+'.png'

datafiles = glob('./*'+model+'.csv')
badcol='Unnamed: 0'

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

for csvfile in sorted(datafiles):
    data = pd.read_csv(csvfile)
    if badcol in data.columns:
        data.drop(badcol,axis=1,inplace=True)
    data.sort_values(by='Jdmax',inplace=True)

   
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
#    popt, pcov = curve_fit(powerLaw, data['Jdmax'], data['mass'])
    popt, pcov = fitLinearLog(data['Jdmax'], data['mass'])
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
#    popt, pcov = curve_fit(powerLaw, data['Jdmax'], data['vDproj'])
    popt, pcov = fitLinearLog(data['Jdmax'], data['vDproj'])
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
    fig.suptitle(csvfile[2:-4],fontweight='heavy')
    fig.tight_layout()
    fig.savefig(csvfile[2:-4]+extension)
    names.append(csvfile[:-4])
#%%
xdata = np.logspace(np.log10(0.001),np.log10(35),120)
plt.figure(figsize=(8,6))
for i,j in zip(pvD,names):
    if 'multaneous' in j:
        plt.plot(xdata,powerLaw(xdata*1e-3,*i),label=j[-5:-2])
        print(j)
    elif 'only' in j:
        plt.plot(xdata,powerLaw(xdata*1e-3,*i),'.',label='rime C')
        print(j)
    elif 'needle' in j:
        plt.plot(xdata,powerLaw(xdata*1e-3,*i),'--',label=j[0:11])
        print(j)
    else:
        plt.plot(xdata,powerLaw(xdata*1e-3,*i),'.-',label=j[-5:-2])
        print(j)
plt.ylim([0,8])
plt.ylabel('terminal fall velocity   [m/s]')
plt.xlabel('Maximum Dimension     [mm]')
plt.legend(ncol=3)
plt.grid()
plt.title('simultaneous-lines         subsequent-markers '+model)
plt.savefig('all_fitted_vD'+extension)

xdata = np.linspace(0,35,60)
plt.figure(figsize=(8,6))
for i,j in zip(pvD,names):
    if 'multaneous' in j:
        #plt.plot(xdata,powerLaw(xdata*1e-3,*i),label=j[-3:])
        print(j)
    elif 'only' in j:
        #plt.plot(xdata,powerLaw(xdata*1e-3,*i),'.',label='rime C')
        print(j)
    elif 'needle' in j:
        plt.plot(xdata,powerLaw(xdata*1e-3,*i),'--',label=j[2:16])
        print(j)
    else:
        #plt.plot(xdata,powerLaw(xdata*1e-3,*i),'h-',label=j[-3:])
        print(j)
plt.ylim([0,3])
plt.legend()
plt.grid()
plt.ylabel('terminal fall velocity   [m/s]')
plt.xlabel('Maximum Dimension     [mm]')
plt.title('unrimed needle aggregates '+model)
plt.savefig('all_needles_vD'+extension)
