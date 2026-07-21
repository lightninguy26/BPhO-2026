# plotting einstein's model of the heat capacity C of solids
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np

R=8.314

fig,ax = plt.subplots(figsize=(9,6))

temps = np.linspace(0.1,1000,1000)

# core plotting loop
def update_curve(t):
    ax.clear()
    # calculate x values
    x_values = t/temps
    
    exp_term = np.exp(x_values)
    # apply einstein's formula to calculate values of C
    second_term = (x_values**2 * exp_term) / ((exp_term -1)**2)
    C_values= 3*R*second_term

    ax.set_title('Einstein model of solid molar heat capacity')
    ax.set_xlabel('Temperature/K')
    ax.set_ylabel('Molar heat capacity / Jmol$^{-1}$ K$^{-1}$')
    ax.axhline(3*R,linestyle='--',label='Dulong-Petit limit')
    ax.plot(temps,C_values)
    plt.draw()

# handle plotting
ax.set_title('Einstein model of solid molar heat capacity')
ax.set_xlabel('Temperature/K')
ax.set_ylabel('Molar heat capacity / Jmol$^{-1}$ K$^{-1}$')
ax.axhline(3*R,linestyle='--',label='Dulong-Petit limit')
#ax.legend()
plt.subplots_adjust(bottom=0.25)

#initialize a slider
ax_slider = plt.axes([0.2,0.1,0.7,0.05],facecolor='teal')
slider= Slider(ax_slider,"Einstein Temperature",valmin=50,valmax=2000,valinit=300,valstep=10)
slider.on_changed(update_curve)


plt.show()


    