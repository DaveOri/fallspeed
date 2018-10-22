#
# terminal fall velocity of ice from Khvorostyanov and Curry (2005/2002)
# optional parameter to switch to smooth surface parameters
# 
# Note: Rassmussen and Heymsfield JAS 44, 2754pp, 1987 interpolate to Cd 0.6 at 
# Re > 2.5e4 on wet particles. This is not done here to avoid a transition to pure water drops
#
import numpy as np

def dia2vel(diam, rho_air, nu_air, mass, area, smooth=False ):
  #mtot = sp%m_r + sp%m_i + sp%m_w 
  grav = 9.807
  rho_ice = 917.0

  # Best number eq. (2.4b) with buoyancy
  Vb = mass/rho_ice
  Fb = rho_air * Vb * grav
  eta_air = nu_air*rho_air # dynamic viscosity
  Xbest = 2. * np.abs(mass*grav-Fb) * rho_air * diam**2 / (area * eta_air**2)
  if( smooth ):
    Cd  = X2Cd_kc05smooth(Xbest)
  else:
    Cd  = X2Cd_kc05rough(Xbest)
  return np.sqrt( 2*np.abs(mass*grav - Fb)/(rho_air * area * Cd))

def X2Cd_kc05rough(Xbest):
  do_i = 5.83
  co_i = 0.6
  Ct = 1.6
  X0_i = 0.35714285714285714285e-6 # 1.0/2.8e6
  # derived constants
  c1 = 4.0 / ( do_i**2 * np.sqrt(co_i))
  c2 = 0.25 * do_i**2

  # Re-X eq. (2.5)
  bracket = np.sqrt(1.0 + c1 * np.sqrt(Xbest)) - 1.0
  # turbulent Reynold's number, eq (3.3)
  psi = (1.0+(Xbest*X0_i)**2) / (1.0+Ct*(Xbest*X0_i)**2)
  Re  = c2*bracket**2 # * np.sqrt(psi) # TODO remove psi in Re?
  # eq. (2.1) from KC05 with (3.2)
  return co_i * (1.0 + do_i / np.sqrt(Re))**2 / psi


def X2Cd_kc05smooth(Xbest):
  do_i = 9.06
  co_i = 0.292
  Ct = 1.6
  X0_i = 1.0/6.7e6
  c1 = 4.0/(do_i**2 * np.sqrt(co_i))
  c2 = 0.25 * do_i**2

  # Re-X eq. (2.5)
  bracket = np.sqrt(1.0 + c1 * np.sqrt(Xbest)) - 1.0
  # turbulent Reynold's number, eq (3.3)
  psi = (1+(Xbest*X0_i)**2) / (1+Ct*(Xbest*X0_i)**2)
  Re  = c2*bracket**2 #* np.sqrt(psi) # TODO remove psi in Re?
  # eq. (2.1) from KC05 with (3.2)
  return co_i * (1. + do_i/np.sqrt(Re))**2 / psi


#   subroutine dia2vel_khvorostyanov01_particles &
#     (errorstatus,&      # out
#     nDia,&    #in
#     diaSpec_SI,&  #in
#     rho_air_SI,&  #in
#     nu_SI,&   #in
#     mass_SI,& #in
#     area_SI,& #in
#     velSpec)    #out

#       #in
#       #nDia: no of diameters
#       #diaSpec_SI = diameter spectrum [m]
#       #rho_air_SI = density of air [kg/m³]
#       #nu_SI = kinematic viscosity of air [m²/s]
#       #mass_size_a_SI,mass_size_b parameters of mass size relation m = a*D_max^b [SI]
#       #area_size_a_SI,area_size_b parameters of mass area relation A = a*D_max^b [SI]
#       #out
#       #velSpec: velocity spectrum [m/s]

#       # Khvorostyanov, V. I. & Curry, J. A. Terminal Velocities of Droplets and Crystals:
#       # Power Laws with Continuous Parameters over the Size Spectrum.
#       # Journal of the Atmospheric Sciences 59, 1872–1884 (2002).
#       # equation 3.3


#       # Modules used:
#       use settings, only: verbose
#       use kinds
#       use constants
#       use report_module
#       implicit none

#       integer, intent(in) :: nDia
#       real(kind=dbl), intent(in), dimension(ndia)::diaSpec_SI, mass_SI, area_SI
#       real(kind=dbl), intent(in) :: rho_air_SI, nu_SI
#       real(kind=dbl), dimension(ndia), intent(out) :: velSpec
#       real(kind=dbl) :: rho_air, g_cp, my
#       real(kind=dbl) :: c1, delta0
#       real(kind=dbl), dimension(ndia):: X, bRe, aRe#, Av, Bv
#       real(kind=dbl), dimension(ndia)::diaSpec, mass, area

#       integer(kind=long), intent(out) :: errorstatus
#       integer(kind=long) :: err = 0
#       character(len=33) :: nameOfRoutine = 'dia2vel_khvorostyanov01_particles'


#       if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

#       # no check for boundaries due to different possible particle types
#       err = success

#       #variables to cgs...
#       mass = mass_SI * 1000.d0
#       area = area_SI * 100.d0**2

#       g_cp = g*100.d0 #cm/s
#       rho_air = rho_air_SI/1000.d0 #g/cm³
#       diaSpec = 100.d0*diaSpec_SI #cm
#       my = nu_SI*10000.d0 #cm² /s
#       c1 = 0.0902d0
#       delta0 = 9.06d0

# #       X = 2*mass_size_a*g_cp*diaSpec**(mass_size_b + 2.d0 - area_size_b)/&
# #       (area_size_a*rho_air*my**2) #eq. 17 of Mitchell et al 1996
#       X = 2 * mass * g_cp  * diaSpec**2 /(area * rho_air* my**2)#eq. 3 of Mitchell et al 1996


#       bRe = 0.5d0*c1*X**0.5d0*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**(-1) * &
#       (1 + c1*X**0.5d0)**(-0.5d0) #2.12

#       aRe = (delta0**2/4.d0)*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**2 / X**bRe #2.13


# #       Av = aRe * my**(1.d0-2.d0*bRe) * ((2.d0*mass_size_a*g_cp)/(rho_air*area_size_a))**bRe #2.24
# #       Bv = (bRe * (mass_size_b-area_size_b+2.d0)) - 1.d0 #2.25
# # 
# #       velSpec = Av * diaSpec**Bv #2.23

#       velSpec = aRe *  my**(1.d0-2.d0*bRe) * ( (2.d0 * mass * g_cp) / (area * rho_air) )**bRe * &
#                     diaSpec**(2d0*bRe - 1d0)#2.20a

#       velSpec = velSpec/100.d0 #CGS to SI, now m/s

#       errorstatus = err
#       if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

#       return

#   end subroutine dia2vel_khvorostyanov01_particles


#   subroutine dia2vel_khvorostyanov01_spheres &
#     (errorstatus,&      # out
#     nDia,&        # in
#     diaSpec,&       # in
#     rho_air,&       # in
#     my,&        # in
#     rho_particle,&      # in
#     velSpec)    # out
    
#       #in
#       #nDia: no of diameters
#       #diaSpec = diameter spectrum [m]
#       #rho_air = density of air [kg/m³]
#       #rho_particle = density of particle [kg/m³]
#       #my = kinematic viscosity of air [m²/s]

#       #
#       #out
#       #velSpec: velocity spectrum [m/s]

#       # Khvorostyanov, V. I. & Curry, J. A. Terminal Velocities of Droplets and Crystals:
#       # Power Laws with Continuous Parameters over the Size Spectrum.
#       # Journal of the Atmospheric Sciences 59, 1872–1884 (2002).
#       # equation 3.3


#       use settings, only: verbose
#       use kinds
#       use constants
#       use report_module
#       implicit none

#       integer, intent(in) :: nDia
#       real(kind=dbl), intent(in), dimension(ndia)::diaSpec, rho_particle
#       real(kind=dbl), intent(in) :: rho_air, my
#       real(kind=dbl), dimension(ndia), intent(out) :: velSpec
#       real(kind=dbl) :: rho_air_cp, g_cp, my_cp
#       real(kind=dbl) :: c1, delta0
#       real(kind=dbl), dimension(ndia):: vB, A, X, bRe, aRe
#       real(kind=dbl), dimension(ndia)::diaSpec_cp, rho_particle_cp

#       integer(kind=long), intent(out) :: errorstatus
#       integer(kind=long) :: err = 0
#       character(len=31) :: nameOfRoutine = 'dia2vel_khvorostyanov01_particles'

#       if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

#       # no check for baundaries due to different possible particle types
#       err = success

#       #variables to cgs...
#       rho_particle_cp = rho_particle/1000.d0 #g/cm³
#       g_cp = g*100.d0 #cm/s
#       rho_air_cp = rho_air/1000.d0 #g/cm³
#       diaSpec_cp = 100.d0*diaSpec #cm
#       my_cp = my*10000.d0 #cm² /s
#       c1 = 0.0902d0
#       delta0 = 9.06d0

#       vB = 4d0/3d0 * pi * (diaSpec_cp*0.5d0)**3
#       A = pi * (diaSpec_cp*0.5d0)**2
#       X = 2*vB*(rho_particle_cp-rho_air_cp)*g_cp*diaSpec_cp**2/&
#       (A*rho_air_cp*my_cp**2) #eq 2.7

#       bRe = 0.5d0*c1*X**0.5d0*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**(-1) * &
#       (1 + c1*X**0.5d0)**(-0.5d0) #eq. 2.12

#       aRe = (delta0**2/4.d0)*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**2 / X**bRe #2.13

#       velSpec = aRe*my_cp**(1-2*bRe) * (2*vB*g_cp/A * ((rho_particle_cp-rho_air_cp)/rho_air_cp))**bRe *&
#       diaSpec_cp**(2*bRe-1)             #2.20b

#       velSpec = velSpec/100.d0 #CGS to SI, now m/s

#       errorstatus = err
#       if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

#       return
#   end subroutine dia2vel_khvorostyanov01_spheres


#   subroutine dia2vel_khvorostyanov01_drops &
#     (errorstatus,&      # out
#     nDia,&        # in
#     diaSpec,&       # in
#     rho_air,&       # in
#     my,&        # in
#     velSpec)    # out
    
#       #in
#       #nDia: no of diameters
#       #diaSpec = diameter spectrum [m]
#       #rho_air density of air [kg/m³]
#       #my = kinematic viscosity of air [m²/s]

#       #out
#       #velSpec: velocity spectrum [m/s]

#       # Khvorostyanov, V. I. & Curry, J. A. Terminal Velocities of Droplets and Crystals:
#       # Power Laws with Continuous Parameters over the Size Spectrum.
#       # Journal of the Atmospheric Sciences 59, 1872–1884 (2002).
#       # equation 2.20c

#       # defined for diameters from 0 to 8.5 mm, non spherical shape is corrected

#       use settings, only: verbose
#       use kinds
#       use constants
#       use report_module
#       implicit none

#       integer, intent(in) :: nDia
#       real(kind=dbl), intent(in), dimension(ndia)::diaSpec
#       real(kind=dbl), intent(in) :: rho_air, my
#       real(kind=dbl), dimension(ndia), intent(out) :: velSpec
#       real(kind=dbl) :: rho_air_cp, g_cp, rho_water_cp, my_cp
#       real(kind=dbl) :: c1, delta0, lam
#       real(kind=dbl), dimension(ndia):: vB, A, X, bRe, aRe, xi
#       real(kind=dbl), dimension(ndia)::diaSpec_cp

#       integer(kind=long), intent(out) :: errorstatus
#       integer(kind=long) :: err = 0
#       character(len=80) :: msg
#       character(len=29) :: nameOfRoutine = 'dia2vel_khvorostyanov01_drops'

#       if (verbose >= 2) call report(info,'Start of ', nameOfRoutine)

#       #check for boundaries (including 1% numerical tolerance)
#       if (MAXVAL(diaSpec) > 8.5d-3*1.01d0) then
#     print*, "Largest diameter for dia2vel_khvorostyanov01_spheres is 8.5 mm, got",&
#     MAXVAL(diaSpec), "[m]"
#     errorstatus = fatal
#     msg = 'Diameter out of specs#'
#     call report(errorstatus, msg, nameOfRoutine)
#     return
#       else
#     err = success
#       end if
#       #variables to cgs...
#       rho_water_cp = rho_water/1000.d0 #g/cm³
#       g_cp = g*100.d0 #cm/s
#       rho_air_cp = rho_air/1000.d0 #g/cm³
#       diaSpec_cp = 100.d0*diaSpec #cm
#       my_cp = my*10000.d0 #cm² /s
#       c1 = 0.0902d0
#       delta0 = 9.06d0


#       lam = 0.47
#       xi = exp(-diaSpec_cp/lam)+(1.d0-exp(-diaSpec_cp/lam))*(1.d0/(1.d0+(diaSpec_cp/lam))) #eq. 3.4

#       vB = 4d0/3d0 * pi * (diaSpec_cp*0.5d0)**3 * xi #correct for non-spherity
#       A = pi * (diaSpec_cp*0.5d0)**2 #no correction neccessary for A

#       X = 2*vB*(rho_water_cp-rho_air_cp)*g_cp*diaSpec_cp**2/&
#       (A*rho_air_cp*my_cp**2) #eq. 2.7

#       bRe = 0.5d0*c1*X**0.5d0*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**(-1) * &
#       (1 + c1*X**0.5d0)**(-0.5d0) #eq. 2.12

#       aRe = (delta0**2/4.d0)*((1.d0+c1*X**0.5d0)**0.5d0 - 1.d0)**2 / X**bRe #eq. 2.13



#       velSpec = aRe*my_cp**(1-2*bRe) * (4.d0/3.d0 * g_cp * xi * ((rho_water_cp-rho_air_cp)/rho_air_cp))**bRe *&
#       diaSpec_cp**(3*bRe-1)  #eq. 2.20c

#       velSpec = velSpec/100.d0 #CGS to SI, now m/s

#       errorstatus = err
#       if (verbose >= 2) call report(info,'End of ', nameOfRoutine)

#       return

#   end subroutine dia2vel_khvorostyanov01_drops


def fall_velocity_KC(area, mass, D_max, T, P):
    """The Khvorostyanov-Curry fall velocity.

    Args:
        area: Projected area [m^2].
        mass: Particle mass [kg].
        D_max: Particle maximum dimension [m].
        T: Ambient temperature [K].
        P: Ambient pressure [Pa].

    Returns:
        The fall velocity [m/s].
    """
    C0 = 0.6 # p. 4348
    delta0 = 5.83 # p. 4348
    C1 = 4.0/(delta0**2 * np.sqrt(C0)) # appendix
    Ct = 1.6 # p. 4345
    k = 2 # p. 4345
    X0 = 2.8e6 # p. 4345

    eta = air_dynamic_viscosity(T)
    rho_air = air_density(T, P)

    g = 9.807

    X = 2.0 * rho_air * mass * g * D_max**2 / (area * eta**2) # 2.4b
    X_sqrt = np.sqrt(X)

    h = np.sqrt(1 + C1*X_sqrt)
    b_Re = C1*X_sqrt / (2 * (h-1) * h) # 2.8
    a_Re = (delta0/4.0) * (h-1)**2 / X**b_Re # 2.7
    Re = a_Re * X**b_Re # 2.6

    # 2.2
    Cd = C0*(1+delta0/np.sqrt(Re))**2
    
    # 3.2
    X_rel = (X/X0)**k
    psi = (1.0 + X_rel)/(1.0 + Ct*X_rel)
    # 3.1
    Cd /= psi
    print(X_rel, Cd)
    # 2.1
    return np.sqrt(2*mass*g / (rho_air*area*Cd)) 
    

def air_kinematic_viscosity(T, P):
    """The kinematic viscosity of air.

    Args:
        T: Ambient temperature [K].
        P: Ambient pressure [Pa].

    Returns:
        The kinematic viscosity [m^2/s].
    """
    rho = air_density(T, P)
    mu = air_dynamic_viscosity(T)
    return mu/rho


def air_dynamic_viscosity(T):
    """The kinematic viscosity of air.

    Args:
        T: Ambient temperature [K].

    Returns:
        The kinematic viscosity [Pa/s].
    """
    mu0 = 1.716e-5
    T0 = 273.15
    C = 111.0
    return mu0 * ((T0+C)/(T+C)) * (T/T0)**1.5


def air_density(T, P):
    """The density of air.

    Args:
        T: Ambient temperature [K].
        P: Ambient pressure [Pa].

    Returns:
        The kinematic viscosity [Pa/s].
    """
    R = 28704e-2 # gas constant for air
    return P / (T*R)
