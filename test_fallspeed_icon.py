#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 18:15:36 2019

@author: dori
"""

import numpy as np
from khvorostyanov_2001 import drops, spheres
from corr_power_law import corrPowLaw
import matplotlib.pyplot as plt
plt.close('all')
rw = 1000.0

hail = {'am':392.33, 'bm':3.0, 'av':106.33, 'bv':0.5, 'Dm':1.87e-4, 'Dx':1.1e-2, 'mod':'sphere'}
graupel = {'am':500.86, 'bm':3.18, 'av':406.67, 'bv':0.85, 'Dm':2.11e-4, 'Dx':1.3e-2, 'mod':'sphere'}
rain = {'am':rw*np.pi/6, 'bm':3.0, 'av':494.74, 'bv':0.70311, 'Dm':8.0e-5, 'Dx':0.8e-2, 'mod':'drop'}
cloud = {'am':rw*np.pi/6, 'bm':3.0, 'av':24388657.6, 'bv':2.0, 'Dm':2.0e-6, 'Dx':8.0e-5, 'mod':'drop'}

def compare(hydro):
  diams = np.linspace(hydro['Dm'], hydro['Dx'], 1000)
  vpowlaw = corrPowLaw(diams, a=hydro['av'], b=hydro['bv'],
                       rho_air=1.2041, T=293.15)
  
  if hydro['mod'] == 'sphere':
    vk = spheres(diams, rho_part=6.0*hydro['am']*diams**(hydro['bm']-3)/np.pi,
                 rho_air=1.2041, nu_air=1.516e-5)
  elif hydro['mod'] == 'drop':
    vk = drops(diams, rho_air=1.2041, nu_air=1.516e-5)
  
  return diams, vpowlaw, vk

fg, ax = plt.subplots(1,1)
f, axs = plt.subplots(2,2)
d, vpl, vk = compare(hail)
axs[0,0].plot([min(vk), max(vk)], [min(vk), max(vk)], c='r', label='1:1')
axs[0,0].scatter(vpl, vk, s=1, label='data')
axs[0,0].grid()
axs[0,0].set_ylabel('Khvorostianov sphere  [m/s]')
axs[0,0].set_xlabel('Power Law icon  [m/s]')
axs[0,0].set_title('Hail')
axs[0,0].legend()

ax.plot(d, vpl, c='C1', label='hail')
ax.plot(d, vk, '--', c='C1')

d, vpl, vk = compare(graupel)
axs[0,1].plot([min(vk), max(vk)], [min(vk), max(vk)], c='r', label='1:1')
axs[0,1].scatter(vpl, vk, s=1, label='data')
axs[0,1].grid()
axs[0,1].set_ylabel('Khvorostianov sphere  [m/s]')
axs[0,1].set_xlabel('Power Law icon  [m/s]')
axs[0,1].set_title('Graupel')
axs[0,1].legend()

ax.plot(d, vpl, c='C2', label='graupel')
ax.plot(d, vk, '--', c='C2')

d, vpl, vk = compare(rain)
axs[1,0].plot([min(vk), max(vk)], [min(vk), max(vk)], c='r', label='1:1')
axs[1,0].scatter(vpl, vk, s=1, label='data')
axs[1,0].grid()
axs[1,0].set_ylabel('Khvorostianov drop  [m/s]')
axs[1,0].set_xlabel('Power Law icon  [m/s]')
axs[1,0].set_title('Rain')
axs[1,0].legend()

ax.plot(d, vpl, c='C3', label='rain')
ax.plot(d, vk, '--', c='C3')

d, vpl, vk = compare(cloud)
axs[1,1].plot([min(vk), max(vk)], [min(vk), max(vk)], c='r', label='1:1')
axs[1,1].scatter(vpl, vk, s=1, label='data')
axs[1,1].grid()
axs[1,1].set_ylabel('Khvorostianov drop  [m/s]')
axs[1,1].set_xlabel('Power Law icon  [m/s]')
axs[1,1].set_title('Cloud')
axs[1,1].legend()

ax.plot(d,vpl, c='C4', label='cloud')
ax.plot(d,vk, '--', c='C4')

ax.legend()
ax.grid()
ax.set_xlabel('diameter   [m]')
ax.set_ylabel('fallspeed    [m/s]')
ax.set_title('___ powLaw --- Khvorostianov')
#ax.set_xscale('log')
#ax.set_yscale('log')

