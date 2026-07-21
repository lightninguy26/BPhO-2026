# plotting rings on a phosphor screen due to electron diffraction

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.widgets import Slider
import numpy as np
from math import sqrt,asin,sin

# setting constants and parameters
h=6.626e-34
e=1.6e-19
V=2*1e3 # in V and to be adjusted with slider 
d=0.123*1e-9 #in m and to be toggled with slider 
r=65*1e-3 #in m and is the radius of the screen
m = 9.1e-31 # mass of an electron in kg

fig,ax=plt.subplots()

def plot_ring(X:float):
     # plot ring halos to simulate real screen
            for lw, a in [(10, 0.04), (6, 0.07), (3, 0.10)]:
                glow = Circle(
                    (0, 0),
                    X,
                    fill=False,
                    linewidth=lw,
                    alpha=a
                )
                ax.add_patch(glow)

            # sharp ring, drawn on top
            ring = Circle(
                (0, 0),
                X,
                fill=False,
                linewidth=1.5,
                alpha=0.85,
                color='green',
                label="Diffraction maximum"
            )
            ax.add_patch(ring)

def show_plot(v,d:float):

    screen=Circle((0,0),r,fill=False,linewidth=2,label="Screen boundary")
    ax.add_patch(screen)

    ax.set_xlim(-r,r)
    ax.set_ylim(-r,r)
    ax.set_aspect("equal",adjustable="box")
    # plot settings
    ax.set_xlabel("screen x / m")
    ax.set_ylabel("screen y / m")
    ax.set_title(f"Rings on a phosphor screen (V={v*1e-3:.1f}kV, d={d*1e9:.3f}nm)")

    ax.grid(alpha=0.2,linestyle='-',linewidth=0.8)
    ax.legend(fontsize=6)
    plt.subplots_adjust(bottom=0.25)


def check_sin_term(n:float,w:float,d:float)->bool:
    return (n*w)/(2*d)<=1

def update_curve_voltage(v):
    # update function
    v_V=v*1e3
    ax.clear()
    n=1
    while True:
        wavelength=h/sqrt(2*m*e*v_V) #calculate wavelength of electron

        if check_sin_term(n,wavelength,d):
            # sine term is less than one so ring exists
            theta= asin((n*wavelength)/(2*d))
            phi=theta*2
            x=r*sin(2*phi)
            plot_ring(x)

            n+=1
            # increase order until maximum not possible
        else:
            break

    show_plot(v_V,d)
    V=v_V
    plt.draw()

def update_curve_d(D):
    # update function
    D_m = D*1e-9
    ax.clear()
    n=1 
    while True:
        wavelength=h/sqrt(2*m*e*V) #calculate wavelength of electron

        if check_sin_term(n,wavelength,D_m):
            # sine term is less than one so ring exists
            theta= asin((n*wavelength)/(2*D_m))
            phi=theta*2
            x=r*sin(2*phi)
            plot_ring(x)

            n+=1
            # increase order until maximum not possible
        else:
            break

    show_plot(V,D_m)
    d=D_m
    plt.draw()
    

# plot for the first time round
n=1
while True:
    wavelength=h/sqrt(2*m*e*V) #calculate wavelength of electron

    if check_sin_term(n,wavelength,d):
        # sine term is less than one so ring exists
        theta= asin((n*wavelength)/(2*d))
        phi=theta*2
        x=r*sin(2*phi)
        plot_ring(x)

        n+=1
        # increase order until maximum not possible
    else:
        break

show_plot(V,d)

ax_voltage_slider = plt.axes([0.3,0.1,0.4,0.05],facecolor='teal')
voltage_slider= Slider(ax_voltage_slider,"Accelerating Voltage/kV",valmin=1,valmax=5,valinit=2,valstep=0.1)
voltage_slider.on_changed(update_curve_voltage)

ax_d_slider = plt.axes([0.4,0.02,0.3,0.05],facecolor='magenta')
d_slider= Slider(ax_d_slider,"Atomic spacing/nm",valmin=0.123,valmax=0.213,valinit=0.123,valstep=0.09)
d_slider.on_changed(update_curve_d)

plt.show()




