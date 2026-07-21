 # plotting probability density of a hydrogenic atom against radius

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from scipy.special import eval_genlaguerre,lpmv
from scipy.constants import electron_mass,Planck,elementary_charge
from math import factorial,sqrt

# defining parameters & constants
m_e=electron_mass
e=elementary_charge
h=Planck
e_0=8.85*1e-12 # vacuum permittivity

n=2
l=0
Z=1
m=0
# Quantum numbers & nuclear charge

a_0_m=(e_0*h**2)/(pi*m_e*e**2)
a_0_a=a_0_m*1e10
a=a_0_a/Z # mu is approximately electron mass

boundary=2 # in angstroms: 1/2 the length of the square slice to plot


def radial(r):
    x=(2*r)/(a*n)
    first_term=sqrt((factorial(n-l-1))/
                    (2*n*factorial(n+l)))
    second_term=pow(2/(a*n),3/2)
    third_term=pow(x,l) * np.exp(-x/2) * eval_genlaguerre(n-l-1,2*l+1,x)
    return (first_term*second_term*third_term) 

# second part of hydrogenic wavefunction
def angular(theta,phi):
    if m>0:
      return spherical(theta,phi,m) + spherical(theta,phi,-m)
    elif m<0:
        return spherical(theta,phi,-m) - spherical(theta,phi,m)
    else:
        return spherical(theta,phi,0)

def spherical(theta,phi,m):
    first_term= pow(-1,m)
    second_term= sqrt( ((2*l+1) * factorial(l-m)) / ((4*pi) * factorial(l+m)) )
    third_term = lpmv(m,l,np.cos(theta)) * np.exp(1j*m*phi)

    return first_term*second_term*third_term

def wavefunction(x,y,z):
    r=np.sqrt(x**2+y**2+z**2)
    r_safe=np.where(r==0,1e-12,r)
    # account for when r is 0
    phi=np.arctan2(y,x)
    theta=np.arccos(z/r_safe)

    psi=radial(r)*spherical(theta,phi,m)
    density=np.abs(psi)**2
    return density # square wavefunction to get probabilities

def checkvalue():
    # retrieve the name of the plane slice e.g. z=0
    for name, value in {"x": x, "y": y, "z": z}.items():
        if not isinstance(value, np.ndarray):
            return name

y=np.linspace(-boundary,boundary,2000)
x=np.linspace(-boundary,boundary,2000) 
z=0 # coordinates to plot across
X,Y=np.meshgrid(x,y)


densities=wavefunction(X,Y,z)
densities/=np.max(densities)

mask= densities>0.005*np.max(densities)
# boundary condition for axes cutoffs

x_interesting=X[mask]
y_interesting=Y[mask]

xmin,xmax = x_interesting.min(),x_interesting.max()
ymin,ymax = y_interesting.min(),y_interesting.max()
padding=0.1*max(xmax-xmin,ymax-ymin)



plt.figure(figsize=(6,5))
plt.imshow(
    densities,
    extent=[X.min(),X.max(),Y.min(),Y.max()],
    origin="lower",
    cmap='jet',
    aspect='equal'
)
plt.colorbar(label=r"Normalised $|\psi|^2$")
plt.xlabel(r"$x$ / $\mathrm{\AA}$")
plt.ylabel(r"$y$ / $\mathrm{\AA}$")
plt.title(fr"${checkvalue()}={0}$ slice, $n={n},\ l={l},\ m={m},\ Z={Z}$ protons ")

#plt.xlim(xmin-padding,xmax+padding)
#plt.ylim(ymin-padding,ymax+padding)
# automatically chosen 'interesting' regions of graph to display

plt.show()


    

