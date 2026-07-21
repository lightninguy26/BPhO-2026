# plotting probability density of a hydrogenic atom against radius

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from scipy.special import eval_genlaguerre
from scipy.constants import electron_mass,Planck,elementary_charge
from math import factorial,sqrt

# defining parameters & constants
m_e=electron_mass
e=elementary_charge
h=Planck
e_0=8.85*1e-12 # vacuum permittivity

n=1
l=0
Z=1
# Quantum numbers & nuclear charge

a_0_m=(e_0*h**2)/(pi*m_e*e**2)
a_0_a=a_0_m*1e10
a=a_0_a/Z # mu is approximately electron mass

r=np.linspace(0.01,4,1000) # radii to plot across in angstroms
x=(2*r)/(a*n)

def radial(r):
    first_term=sqrt((factorial(n-l-1))/
                    (2*n*factorial(n+l)))
    second_term=pow(2/(a*n),3/2)
    third_term=pow(x,l) * np.exp(-x/2) * eval_genlaguerre(n-l-1,2*l+1,x)
    return (first_term*second_term*third_term)**2 # multiply by r^2 to convert probability density to radial probability


plt.plot(r,radial(r))
plt.title(f'Hydrogenic Atom/Ion: Z={Z}, n={n}, l={l}')
plt.xlabel("Radius / $\mathrm{\AA}$")
plt.ylabel("Radial Probability density / $\mathrm{\AA}$ $^{-3}$")

plt.show()



