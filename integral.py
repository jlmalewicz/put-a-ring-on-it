# integral
Let $r = \frac{z}{1-z}+r_1 \rightarrow dr = \frac{dz}{(1-z)^2}$. Then, the integral becomes:
$$\Delta\phi=2\int_{r_1}^{\infty}\frac{dr}{r^2}\left(\frac{1}{b^2}-W_{eff}(r)\right)^{-1/2}=2\int_{0}^{1}\frac{1}{\left(\frac{z}{1-z}+r_1\right)^2}\frac{dz}{(1-z)^2}\left(\frac{1}{b^2}-\frac{1}{\left(\frac{z}{1-z}+r_1\right)^2}+\frac{2M}{\left(\frac{z}{1-z}+r_1\right)^3}\right)^{-1/2}$$
$$\Delta\phi=2\int_{0}^{1}\frac{dz}{\left[z+r_1(1-z)\right]^2}\left(\frac{1}{b^2}-\frac{1}{\left(\frac{z}{1-z}+r_1\right)^2}+\frac{2M}{\left(\frac{z}{1-z}+r_1\right)^3}\right)^{-1/2}$$

from gaussxw import gaussxwab
from scipy import constants

b = 
M = 
r1 = 

def h(z):
    return 2/((z+r1*(1-z))**2) * (1/b**2 - 1/((z/(1-z)+r1)**2) + (2*M)/((z/(1-z)+r1)**3))**(-1/2)
  
# calculate integral using  the routines in the file gaussxb.py on Canvas
lower_lim = 0
up_lim = 1
    
x1, w1 = gaussxwab(50,lower_lim,up_lim)
integral = 0

for l in range(50):
    integral += w1[l]*h(x1[l])
print(integral)
