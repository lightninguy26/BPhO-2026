# plotting rings on a phosphor screen due to electron diffraction

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
from math import sqrt,asin,sin

# setting constants and parameters
h=6.626e-34
e=1.6e-19
V=[1e3,2e3,3e3,4e3,5e3] # in V and to be adjusted with slider later
d1,d2=0.123*1e-9,0.213*1e-9 #in m and to be toggled with slider later
r=65*1e-3 #in m and is the radius of the screen
m = 9.1e-31 # mass of an electron in kg
n=1

l1=[]
l2=[]

for d,l in [(d1,l1),(d2,l2)]:
    for v in V:      
        # collect data to plot inverse root voltage against sin of theta
        wavelength=h/sqrt(2*m*e*v) #calculate wavelength of electron
        theta= asin((n*wavelength)/(2*d))
        phi=theta*2
        x=r*sin(2*phi)

        sint= sin(theta)
        ivst=1/sqrt(v)

        l.append((sint,ivst))

print("--- 0.123nm LIST--- :")
print(l1)
print()
print("--- 0.213nm LIST--- :")
print(l2)
        
        
      






