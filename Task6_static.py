# plotting rings on a phosphor screen due to electron diffraction

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from math import sqrt,asin,sin

# setting constants and parameters
h=6.626e-34
e=1.6e-19
V=2*1e3 # in V and to be adjusted with slider later
d1,d2=0.123*1e-9,0.213*1e-9 #in m and to be toggled with slider later
r=65*1e-3 #in m and is the radius of the screen
m = 9.1e-31 # mass of an electron in kg

fig,ax=plt.subplots()
screen=Circle((0,0),r,fill=False,linewidth=2,label="Screen boundary")
ax.add_patch(screen)

def check_sin_term(n:float,w:float,d:float)->bool:
    return (n*w)/(2*d)<=1

for d,clr in [(d1,'green'),(d2,'red')]:      
    n=1
    while True:
        wavelength=h/sqrt(2*m*e*V) #calculate wavelength of electron

        if check_sin_term(n,wavelength,d):
            # sine term is less than one so ring exists
            theta= asin((n*wavelength)/(2*d))
            phi=theta*2
            x=r*sin(2*phi)

            # plot ring halos to simulate real screen
            for lw, a in [(5, 0.04), (3, 0.07), (1.5, 0.10)]:
                glow = Circle(
                    (0, 0),
                    x,
                    fill=False,
                    linewidth=lw,
                    alpha=a
                )
                ax.add_patch(glow)

            # sharp ring, drawn on top
            ring = Circle(
                (0, 0),
                x,
                fill=False,
                linewidth=0.75,
                alpha=0.85,
                color=clr,
                label=f"Diffraction maxima: spacing={d}nm"
            )
            ax.add_patch(ring)
            
            n+=1
        else:
            break


ax.set_xlim(-r,r)
ax.set_ylim(-r,r)
ax.set_aspect("equal",adjustable="box")
# plot settings
ax.set_xlabel("screen x / m")
ax.set_ylabel("screen y / m")
ax.set_title("Rings on a phosphor screen (V=2kV, d=0.123nm & 0.213nm)")


ax.grid(alpha=0.2,linestyle='-',linewidth=0.8)

ax.legend(fontsize=3,loc='upper right')

plt.show()




