# plotting probability density of a hydrogenic atom against radius

import plotly.graph_objects as go
import plotly.io as pio
import numpy as np
from numpy import pi
from scipy.special import eval_genlaguerre,lpmv,sph_harm_y
from scipy.constants import electron_mass,Planck,elementary_charge
from math import factorial,sqrt

pio.renderers.default = "browser"

# defining parameters & constants
m_e=electron_mass
e=elementary_charge
h=Planck
e_0=8.85*1e-12 # vacuum permittivity


Z=1
n=2
orbital='sp²(-1/2,-d√3/2,0) = (1/√3)ψ₂s + √(2/3)(-½ψ₂px - √3/2 ψ₂py)'
# Quantum numbers for hybridised orbitals & nuclear charge

a_0_m=(e_0*h**2)/(pi*m_e*e**2)
a_0_a=a_0_m*1e10
a=a_0_a/Z # mu is approximately electron mass

boundary=20 # in angstroms: 1/2 the length of the square slice to plot
resolution=150
opacity=0.3
isomin=0.15 # threshold for 'interesting' points



def radial(r,l):
    x=(2*r)/(a*n)
    first_term=sqrt((factorial(n-l-1))/
                    (2*n*factorial(n+l)))
    second_term=pow(2/(a*n),3/2)
    third_term=pow(x,l) * np.exp(-x/2) * eval_genlaguerre(n-l-1,2*l+1,x)
    return (first_term*second_term*third_term)

# second part of hydrogenic wavefunction
def angular(theta,phi,m,l):
    if m==0:
        return spherical(theta,phi,0,l) # pz
    else:
        mu=abs(m)
        Y_pos = sph_harm_y(l, mu, theta, phi)
        Y_neg = sph_harm_y(l, -mu, theta, phi)
        if m>0:
            Y_real = (Y_neg + (-1)**mu * Y_pos)/sqrt(2)
            return np.real(Y_real) # px
        else:
            Y_real = (Y_neg - (-1)**mu * Y_pos) / (1j * sqrt(2))
            return np.real(Y_real) # py
    # use the library function for the more complex case when m is not 0
        

def spherical(theta,phi,m,l):
    first_term= pow(-1,m)
    second_term= sqrt( ((2*l+1) * factorial(l-m)) / ((4*pi) * factorial(l+m)) )
    third_term = lpmv(m,l,np.cos(theta)) * np.exp(1j*m*phi)

    return first_term*second_term*third_term

def wavefunction(x,y,z,m,l):
    r=np.sqrt(x**2+y**2+z**2)
    r_safe=np.where(r==0,1e-12,r)
    # account for when r is 0
    phi=np.arctan2(y,x)
    theta=np.arccos(z/r_safe)

    psi=radial(r,l=l)*angular(theta,phi,m,l=l)
    
    return psi # do not square yet when calculating hybridised orbitals

def checkvalue():
    # retrieve the name of the plane slice e.g. z=0
    for name, value in {"x": x, "y": y, "z": z}.items():
        if not isinstance(value, np.ndarray):
            return name


y=np.linspace(-boundary,boundary,resolution)
x=np.linspace(-boundary,boundary,resolution) 
z=np.linspace(-boundary,boundary,resolution) # coordinates to plot across
X,Y,Z_values=np.meshgrid(x,y,z)

fig=go.Figure()


combined_psi = (wavefunction(X,Y,Z_values,l=0,m=0)/sqrt(3) + (-wavefunction(X,Y,Z_values,l=1,m=1)/2 - (sqrt(3)/2)* wavefunction(X,Y,Z_values,l=1,m=-1) )*sqrt(2/3))
# normalized hybrid orbital

densities = np.abs(combined_psi)**2

densities/=np.max(densities)
# normalise probability densities
mask = densities>isomin
# boundary condition for axes cutoffs
x_interesting=X[mask]
y_interesting=Y[mask] 
z_interesting=Z_values[mask]
# apply filter and then calculate cutoffs with padding
xmin,xmax = x_interesting.min(),x_interesting.max()
ymin,ymax = y_interesting.min(),y_interesting.max()
zmin, zmax = z_interesting.min(), z_interesting.max()
padding=0.1*max(xmax-xmin,ymax-ymin,zmax-zmin)

fig.add_trace(go.Isosurface(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z_values.flatten(),
    value=densities.flatten(),
    isomin=isomin,
    isomax=0.6,
    surface_count=3,
    opacity=opacity,
    colorscale='Jet'
))
# create the iso surface plot and define density boundaries

def set_plot(n,orbital):
    fig.update_layout(
            title=f"Hydrogenic orbital probability density |ψ|²: n = {n}, Z = {Z} protons, hybridised orbital= {orbital}",
            scene=dict(
                xaxis=dict(title="x / Å",range=[xmin-padding,xmax+padding]),
                yaxis=dict(title="y / Å",range=[ymin-padding,ymax+padding]),
                zaxis=dict(title="z / Å",range=[zmin-padding,zmax+padding]),
                aspectmode='cube'
            )
        )

set_plot(n=n,orbital=orbital)

fig.show()

