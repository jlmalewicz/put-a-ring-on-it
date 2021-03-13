# Contains constants & functions common in GR

import numpy as np

G = 6.6743e-11 #Nm^2/kg^2




def Schwarzschild_metric(M, R, theta):
    # theta = sph. polar angle
    # geometrized units
    g = np.zeroes([4,4])
    g[0,0] = -(1-2*M/r)
    g[1,1] =  (1-2*M)**(-1)
    g[2,2] =  r**2
    g[3,3] =  r**2 * np.sin(theta)**2
    return g

def convert_sph_to_car(r, theta, phi):
    x = r * np.sin(phi) * np.cos(theta)
    y = r * np.sin(phi) * np.sin(theta)
    z = r * np.cos(phi)
    return x, y, z

def convert_car_to_sph(x, y, z):
    r     = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan(y/x)
    phi   = np.arctant(np.sqrt(x**2 + y**2)/z)

def convert_pol_to_car(rho, theta, z):
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    return x, y, z

def convert_car_to_pol(x, y, z):
    rho   = np.sqrt(x**2 + y**2) 
    theta = np.arctan(y/x)
