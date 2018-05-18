# -*- coding: utf-8 -*-
"""
Created on Fri May 18 16:20:35 2018

@author: dori
"""
import numpy as np

def dia2vel(diaSpec_SI, rho_air_SI, nu_SI, mass, area, k=0.5):
#   !in
#      !nDia: no of diameters
#      !diaSpec_SI = diameter spectrum [m]
#      !rho_air_SI = density of air [kg/m³]
#      !nu_SI = kinematic viscosity of air [m²/s]
#      !mass = mass of the particle [SI]
#      !area = cross section area [SI]
#      !out
#      !velSpec: velocity spectrum [m/s]
#
#      !Heymsfield, A. J. & Westbrook, C. D. Advances in the Estimation of Ice Particle Fall Speeds
#      !Using Laboratory and Field Measurements. Journal of the Atmospheric Sciences 67, 2469–2482 (2010).
#
#      ! equations 9-11
#
#
#      ! Modules used:
#      use settings, only: verbose
#      use kinds
#      use constants
#      use report_module
#      implicit none
#      
#      integer, intent(in) :: nDia
#      real(kind=dbl), intent(in), dimension(ndia)::diaSpec_SI, mass,area
#      real(kind=dbl), intent(in) :: rho_air_SI, nu_SI, k
#      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
#      real(kind=dbl), dimension(ndia)::area_proj, Xstar, Re
#      real(kind=dbl) :: delta_0, C_0, eta
#
#
#      integer(kind=long), intent(out) :: errorstatus
#      integer(kind=long) :: err = 0
#      character(len=33) :: nameOfRoutine = 'dia2vel_heymsfield10_particles'

# no check for baundaries due to different possible particle types
#     k = 0.5   #defined in the paper
    delta_0 = 8.0
    C_0 = 0.35
    g = 9.81
    
    area_proj = area/((np.pi/4.)*diaSpec_SI**2)
    eta = nu_SI * rho_air_SI #!now dynamic viscosity

    Xstar = 8.0*rho_air_SI*mass*g/(np.pi*area_proj**(1.0-k)*eta**2)# !eq 9
    Re=0.25*delta_0**2*((1.0+((4.0*Xstar**0.5)/(delta_0**2.0*C_0**0.5)))**0.5 - 1 )**2 #!eq10
     
    velSpec = eta*Re/(rho_air_SI*diaSpec_SI)
    return velSpec