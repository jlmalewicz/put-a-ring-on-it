# integral
Let $r = \frac{z}{1-z}+r_1 \rightarrow dr = \frac{dz}{(1-z)^2}$. Then, the integral becomes:
$$\Delta\phi=2\int_{r_1}^{\infty}\frac{dr}{r^2}\left(\frac{1}{b^2}-W_{eff}(r)\right)^{-1/2}=2\int_{0}^{1}\frac{1}{\left(\frac{z}{1-z}+r_1\right)^2}\frac{dz}{(1-z)^2}\left(\frac{1}{b^2}-\frac{1}{\left(\frac{z}{1-z}+r_1\right)^2}+\frac{2MG}{c^2\left(\frac{z}{1-z}+r_1\right)^3}\right)^{-1/2}$$
$$\Delta\phi=2\int_{0}^{1}\frac{dz}{\left[z+r_1(1-z)\right]^2}\left(\frac{1}{b^2}-\frac{1}{\left(\frac{z}{1-z}+r_1\right)^2}+\frac{2MG}{c^2\left(\frac{z}{1-z}+r_1\right)^3}\right)^{-1/2}$$

from gaussxw import gaussxwab
from scipy import constants

b =  
r1 = 

G = constants.G
c = constants.c
solar_mass = 1.98847e30
M = solar_mass*1e3 
R = 2*G*M/(c**2)

def h(z,r1,b):
    de = z/(1-z)+r1
    return 2/((z+r1*(1-z))**2) * (1/(b**2) - 1/(de**2) + (2*M*G)/(c**2)/(de**3))**(-1/2)
        
lower_lim = 0
up_lim = 1
    
x1, w1 = gaussxwab(1000,lower_lim,up_lim)
integral = 0

for l in range(1000):
    integral += w1[l]*h(x1[l],r1,b)
integral
