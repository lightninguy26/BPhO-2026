 # plotting probability density of a hydrogenic atom against radius

import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
from numpy import pi
from scipy.special import eval_genlaguerre,lpmv
from scipy.constants import electron_mass,Planck,elementary_charge
from math import factorial,sqrt

pio.renderers.default = "browser"

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
opacity=0.15
cutoff_probability=0.01

densities_list=[]
Z_plane_list=[]

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


y=np.linspace(-boundary,boundary,100)
x=np.linspace(-boundary,boundary,100) 
z=np.linspace(-boundary,boundary,30) # coordinates to plot across
X,Y=np.meshgrid(x,y)

fig=go.Figure()

for Z_value in z:
 
    densities=wavefunction(X,Y,Z_value)
    densities_list.append(densities)
    Z_plane_list.append(Z_value)


    #mask= densities>0.005*np.max(densities)
    # boundary condition for axes cutoffs

    #x_interesting=X[mask]
    #y_interesting=Y[mask]

    #xmin,xmax = x_interesting.min(),x_interesting.max()
    #ymin,ymax = y_interesting.min(),y_interesting.max()
    #padding=0.1*max(xmax-xmin,ymax-ymin)

densities_list=np.array(densities_list)
densities_list/=np.max(densities_list)
# normalise densities


for i,d in enumerate(densities_list):

    z=Z_plane_list[i]
    Z_plane=np.full_like(X,z)

    mask = d>cutoff_probability
    Z_visible=np.where(mask,Z_plane,np.nan)
    # only include 'interesting' regions
    fig.add_trace(go.Surface(
        x=X,
        y=Y,
        z=Z_visible,
        surfacecolor=d,
        opacity=opacity,
        colorscale='Jet',
        cmin=0,
        cmax=1,
        colorbar=dict(
            title="Normalised |ψ|²"
        ),
        showscale= i==0
    ))
# add each z slice to the 3D plot

fig.update_layout(
    title="Hydrogenic orbital probability density",
    scene=dict(
        xaxis_title="x / Å",
        yaxis_title="y / Å",
        zaxis_title="z / Å"
    )
)

fig.show()






    

