import numpy as np

def drops(diam, rho_air=1.2041, nu_air=1.516e-5):
    D = diam*100.0
    g = 9.81e2
    v = nu_air*1.0e4
    cl = 0.0902
    d0 = 9.06
    lam = 0.47
    rw = 1.0
    ra = rho_air*1.0e-3

    xi = np.exp(-D/lam)+(1.-np.exp(-D/lam))*(1./(1.+(D/lam)))
    vB = np.pi*D*D*D*xi/6. #correct for non-spherity
    A = np.pi*D*D*0.25 #no correction neccessary for A
    X = 2.0*vB*(rw-ra)*g*D*D/(A*ra*v*v)
    sX = np.sqrt(X)
    bRe = 0.5*cl*sX/((np.sqrt(1.+cl*sX)-1.)*np.sqrt(1+cl*sX))
    aRe = (d0*d0*0.25)*(np.sqrt(1.+cl*sX)-1.)**2 / X**bRe
    
    velSpec = aRe*v**(1.-2.*bRe) * (4./3.*g*xi*(rw-ra)/ra)**bRe *D**(3.*bRe-1.)
    return velSpec*0.01


def spheres(diam, rho_part, rho_air=1.2041, nu_air=1.516e-5):
    D = diam*100.0
    g = 9.81e2
    v = nu_air*1.0e4
    cl = 0.0902
    d0 = 9.06
    ra = rho_air*1.0e-3
    rp = rho_part*1.0e-3

    vB = np.pi*D*D*D/6.
    A = np.pi*D*D*0.25
    X = 2*vB*(rp-ra)*g*D**2/(A*ra*v**2)

    bRe = 0.5*cl*X**0.5/((1.+cl*X**0.5)**0.5 - 1.) * (1 + cl*X**0.5)**(-0.5) 
    aRe = (d0*d0*0.25)*((1.+cl*X**0.5)**0.5 - 1.)**2 / X**bRe

    velSpec = aRe*v**(1.-2.*bRe) * (2.*vB*g/A * (rp-ra)/ra)**bRe *D**(2.*bRe-1.)
    return velSpec*0.01