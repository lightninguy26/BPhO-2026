# plotting einstein's model of the heat capacity C of solids
import matplotlib.pyplot as plt
import numpy as np
from math import cbrt

# parameters and constants

debye_temps = np.array([170,343.5,420,428,470,645,2230])
# values for gold, copper, titanium,  aluminium, iron, silicon and carbon
chemical_names=['Au','Cu','Ti','Al','Fe','Si','C']
einstein_temps = debye_temps * cbrt(np.pi/6)
# einstein temps in terms of debyes

R=8.314


temps = np.linspace(0.1,1000,1000)

# core plotting loop
for i in range(len(einstein_temps)):
    t=einstein_temps[i]
    # calculate x values
    x_values = t/temps
    
    exp_term = np.exp(x_values)
    # apply einstein's formula to calculate values of C
    second_term = (x_values**2 * exp_term) / ((exp_term -1)**2)
    C_values= 3*R*second_term

    plt.plot(temps,C_values,label=chemical_names[i])

plt.title('Einstein model of solid molar heat capacity')
plt.xlabel('Temperature/K')
plt.ylabel('Molar heat capacity / Jmol$^{-1}$ K$^{-1}$')
plt.axhline(3*R,linestyle='--',label='Dulong-Petit limit')

plt.legend()
plt.show()


    
