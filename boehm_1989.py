#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 10:44:35 2019

@author: dori
"""

import numpy as np

def dia2vel(diam, rho_air, nu_air, mass, area):
  grav = 9.81
  rho_ice = 917.0
  q = area / (np.pi/4.0 * diam**2)
  eta_air = nu_air*rho_air # dynamic viscosity
  
  #alpha = np.array(as_ratio) #1.0
  #X_0 = 2.8e6
  X = 8.0*mass*grav*rho_air/(np.pi*(eta_air**2)*q**0.25)
  
  #k = np.minimum(np.maximum(0.82+0.18*alpha,np.ones_like(alpha)*0.85),0.37+0.63/alpha,1.33/(np.maximum(np.log(alpha),np.ones_like(alpha)*0.0)+1.19)) #k is 1 for alpha=1
  #gama_big = np.maximum(np.ones_like(alpha)*1.0, np.minimum(np.ones_like(alpha)*1.98,3.76-8.41*alpha+9.18*alpha**2-3.53*alpha**3)) #1 for alpha=1
  #C_DP = np.maximum(0.292*k*gama_big,0.492-0.2/np.sqrt(alpha)) #0.292 for alpha=1
  #C_DP = np.maximum(1.0,q*(1.46*q-0.46))*C_DP #0.292 for alpha=1
  #C_DP_prim = C_DP*(1.0+1.6*(X/X_0)**2)/(1.0+(X/X_0)**2) #0.292 for small particles; larger for bigger particles 
  #beta = np.sqrt(1.0+C_DP_prim/6.0/k*np.sqrt(X/C_DP_prim))-1
  #N_Re0 = 6.0*k/C_DP_prim*beta**2
  #C_DO = 4.5*k**2*np.maximum(alpha,np.ones_like(alpha)*1.0)
  #gama_small = (C_DO - C_DP)/4.0/C_DP
  #N_Re  = N_Re0*(1.0 + (2.0*beta*np.exp(-beta*gama_small))/((2.0+beta)*(1.0+beta)) )
  Re = 8.5*((1.0+0.1519*X**0.5)**0.5-1.0)**2
  vterm_bohm = Re*eta_air/diam/rho_air
  return vterm_bohm