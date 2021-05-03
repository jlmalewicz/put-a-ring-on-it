from gaussxw import gaussxwab
from scipy import constants

def integral(r1,b):
    def h(z,r1,b):
        de = z/(1-z)+r1
        return 2/((z+r1*(1-z))**2) * (1/(b**2) - 1/(de**2) + (2*M*G)/(c**2*de**3))**(-1/2)
        
    lower_lim = 0
    up_lim = 1
    
    x1, w1 = gaussxwab(1000,lower_lim,up_lim)
    integral1 = 0

    for l in range(1000):
        integral1 += w1[l]*h(x1[l],r1,b)
        
    return integral1
