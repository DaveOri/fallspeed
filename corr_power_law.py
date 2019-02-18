import numpy as np

def corrPowLaw(diam, a, b, rho_air=1.2041, T=293.15):
    r0 = 1.2038631624242195
    if rho_air > r0:
      rho_air = r0 - 0.001
    Y = 0.43*np.log10(r0/rho_air)-0.4*(np.log10(r0/rho_air))**2.5
    return a*diam**b*10.**Y*(1.+(0.0023*(1.1-(rho_air/r0))*(293.15-T)))