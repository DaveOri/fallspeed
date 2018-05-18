# -*- coding: utf-8 -*-
"""
Created on Fri May 18 16:20:35 2018

@author: dori
"""

subroutine dia2vel_heymsfield10_particles &
    (errorstatus,&     	! out
    nDia,& 		!in
    diaSpec_SI,&	!in
    rho_air_SI,&	!in
    nu_SI,&		!in
    mass,&              !in
    area,&
    k,&                 !in
    velSpec)		!out

      !in
      !nDia: no of diameters
      !diaSpec_SI = diameter spectrum [m]
      !rho_air_SI = density of air [kg/m³]
      !nu_SI = kinematic viscosity of air [m²/s]
      !mass = mass of the particle [SI]
      !area = cross section area [SI]
      !out
      !velSpec: velocity spectrum [m/s]

      !Heymsfield, A. J. & Westbrook, C. D. Advances in the Estimation of Ice Particle Fall Speeds
      !Using Laboratory and Field Measurements. Journal of the Atmospheric Sciences 67, 2469–2482 (2010).

      ! equations 9-11


      ! Modules used:
      use settings, only: verbose
      use kinds
      use constants
      use report_module
      implicit none

      integer, intent(in) :: nDia
      real(kind=dbl), intent(in), dimension(ndia)::diaSpec_SI, mass,area
      real(kind=dbl), intent(in) :: rho_air_SI, nu_SI, k
      real(kind=dbl), dimension(ndia), intent(out) :: velSpec
      real(kind=dbl), dimension(ndia)::area_proj, Xstar, Re
      real(kind=dbl) :: delta_0, C_0, eta


      integer(kind=long), intent(out) :: errorstatus
      integer(kind=long) :: err = 0
      character(len=33) :: nameOfRoutine = 'dia2vel_heymsfield10_particles'

      if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

      ! no check for baundaries due to different possible particle types
      err = success

!       k = 0.5d0 !defined in the paper
      delta_0 = 8.0d0
      C_0 = 0.35d0

      area_proj = area/((pi/4.d0)*diaSpec_SI**2)
      eta = nu_SI * rho_air_SI !now dynamic viscosity

      Xstar = 8.d0 * rho_air_SI * mass * g / (pi * area_proj**(1.d0-k) * eta**2) !eq 9
      Re = delta_0**2/4.d0 * ((1.d0+((4.d0*sqrt(Xstar))/(delta_0**2*sqrt(C_0))))**0.5 - 1 )**2 !eq10
     
      velSpec = eta * Re /(rho_air_SI*diaSpec_SI)

      errorstatus = err
      if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

      return

  end subroutine dia2vel_heymsfield10_particles