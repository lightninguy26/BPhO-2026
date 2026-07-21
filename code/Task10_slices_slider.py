from matplotlib.widgets import Slider
 # plotting probability density of a hydrogenic atom against radius

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from scipy.special import eval_genlaguerre,lpmv,sph_harm_y
from scipy.constants import electron_mass,Planck,elementary_charge
from math import factorial,sqrt

# defining parameters & constants
m_e=electron_mass
e=elementary_charge
h=Planck
e_0=8.85*1e-12 # vacuum permittivity

n=2
l=1
Z=1
m=0
# Quantum numbers & nuclear charge

a_0_m=(e_0*h**2)/(pi*m_e*e**2)
a_0_a=a_0_m*1e10
a=a_0_a/Z # mu is approximately electron mass

resolution=600
boundary=4 # in angstroms: 1/2 the length of the square slice to plot
fig,ax=plt.subplots(figsize=(6,5))
image = None

def radial(r):
    x=(2*r)/(a*n)
    first_term=sqrt((factorial(n-l-1))/
                    (2*n*factorial(n+l)))
    second_term=pow(2/(a*n),3/2)
    third_term=pow(x,l) * np.exp(-x/2) * eval_genlaguerre(n-l-1,2*l+1,x)
    return (first_term*second_term*third_term) 

# second part of hydrogenic wavefunction
def angular(theta,phi,m,l):
    if m==0:
        return spherical(theta,phi,0,l=l)
    else:
        Y_pos = sph_harm_y(l, abs(m), theta, phi)
        Y_neg = sph_harm_y(l, -abs(m), theta, phi)
        if m>0:
            return Y_pos+Y_neg
        else:
            return (Y_pos-Y_neg)/1j
    # use the library function for the more complex case when m is not 0

def spherical(theta,phi,m,l):
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

    psi=radial(r)*angular(theta,phi,m,l=l)
    density=np.abs(psi)**2
    return density # square wavefunction to get probabilities

def checkvalue():
    # retrieve the name of the plane slice e.g. z=0
    for name, value in {"x": x, "y": y, "z": z}.items():
        if not isinstance(value, np.ndarray):
            return name


z=np.linspace(-boundary,boundary,resolution)
x=np.linspace(-boundary,boundary,resolution) 
y=0 # coordinates to plot across
X,Y=np.meshgrid(x,z)

def replot():
    global l
    global m
    l=n-1 if l>n-1 else l
    if m<-l:
        m=-l
    elif m>l:
        m=l
    # check if quantum values are valid and clamp them if not
    densities=wavefunction(X,Y,z)
    densities/=np.max(densities)

    mask= densities>0.005*np.max(densities)
    # boundary condition for axes cutoffs

    x_interesting=X[mask]
    y_interesting=Y[mask]

    xmin,xmax = x_interesting.min(),x_interesting.max()
    ymin,ymax = y_interesting.min(),y_interesting.max()
    padding=0.1*max(xmax-xmin,ymax-ymin)
    
    #fig.colorbar(image,ax=ax,label=r"Normalised $|\psi|^2$")
    image.set_data(densities)
    
    ax.set_title(fr"${checkvalue()}={0}$ slice, $n={n},\ l={l},\ m={m},\ Z={Z}$ protons ")

    ax.set_xlim(xmin-padding,xmax+padding)
    ax.set_ylim(ymin-padding,ymax+padding) 
    # automatically chosen 'interesting' regions of graph to display
    fig.canvas.draw_idle()


densities=wavefunction(X,Y,z)
densities/=np.max(densities)

mask= densities>0.005*np.max(densities)
    # boundary condition for axes cutoffs

x_interesting=X[mask]
y_interesting=Y[mask]

xmin,xmax = x_interesting.min(),x_interesting.max()
ymin,ymax = y_interesting.min(),y_interesting.max()
padding=0.1*max(xmax-xmin,ymax-ymin)
    
ax.set_xlabel(r"$x$ / $\mathrm{\AA}$")
ax.set_ylabel(r"$z$ / $\mathrm{\AA}$")
ax.set_title(fr"${checkvalue()}={0}$ slice, $n={n},\ l={l},\ m={m},\ Z={Z}$ protons ")

plt.subplots_adjust(bottom=0.25)

ax.set_xlim(xmin-padding,xmax+padding)
ax.set_ylim(ymin-padding,ymax+padding) 

image = ax.imshow(
        densities,
        extent=[X.min(),X.max(),Y.min(),Y.max()],
        origin="lower",
        cmap='jet',
        aspect='equal'
    )
fig.colorbar(image,ax=ax,label=r"Normalised $|\psi|^2$")

def adjust_n(N):
    global n
    n=int(N)
    replot()

def adjust_z(charge):
    global Z
    global a
    Z=int(charge)
    a=a_0_a/Z
    replot()

def adjust_m(M):
    global m
    m=int(M)
    replot()

def adjust_l(L):
    global l
    l=int(L)
    replot()

ax_n_slider = plt.axes([0.3,0.12,0.4,0.02],facecolor='teal')
n_slider= Slider(ax_n_slider,"n",valmin=1,valmax=5,valinit=2,valstep=1)
n_slider.on_changed(adjust_n)

ax_z_slider = plt.axes([0.3,0.09,0.4,0.02],facecolor='magenta')
z_slider= Slider(ax_z_slider,"Z",valmin=1,valmax=10,valinit=1,valstep=1)
z_slider.on_changed(adjust_z)

ax_l_slider = plt.axes([0.3,0.06,0.4,0.02],facecolor='green')
l_slider= Slider(ax_l_slider,"l",valmin=0,valmax=4,valinit=1,valstep=1)
l_slider.on_changed(adjust_l)

ax_m_slider = plt.axes([0.3,0.03,0.4,0.02],facecolor='red')
m_slider= Slider(ax_m_slider,"m",valmin=-4,valmax=4,valinit=0,valstep=1)
m_slider.on_changed(adjust_m)

plt.show()


    

