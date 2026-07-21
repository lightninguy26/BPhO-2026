# plotting the Planck spectrum for various temperatures and analysing it to extract the Physics

import matplotlib.pyplot as plt
import numpy as np

# parameters
h=6.63*10**-34
c=2.998*10**8
kb=1.381*10**-23
dx=1/6000

area_dict={}
max_wavelength_dict={}
# temperatures in kelvin
temperatures = np.array([3500,4000,4500,5000,5500])
# range of wavelengths in nm to plot
wavelengths = np.linspace(0.1,3000.0,num=int(1/dx))
wavelengths_m = wavelengths * 1e-9

# create function
def planck(w:np.ndarray,t:float)->np.ndarray:
    # apply formula
    term_1 = (2*h*c**2)/w**5
    term_2_denominator = np.exp( (h*c) / (w*kb*t) )-1
    term_2 = 1/term_2_denominator
    # handle term-by-term to reduce clutter 
    return term_1*term_2

for t in temperatures:
    # pass wavelengths into function and get array of irradiances back out
    irradiances=planck(wavelengths_m,t)

    # integration using trapezium method: multiplication by pi to convert spectral radiance to total emitted power
    intensity = np.pi * np.trapezoid(irradiances,wavelengths_m)
    
    # append to area dict
    area_dict[t]=intensity
    # log max wavelength for later analysis of wien's law
    max_irradiance = max(irradiances)
    index = np.argmax(irradiances)
    max_wavelength_dict[t]= wavelengths_m[index]

    plt.plot(wavelengths,irradiances,label=f'{t}K')

# display max wavelength and area dictionaries 
print('MAX WAVELENGTHS:')
print(max_wavelength_dict)
print('AREAS:')
print(area_dict)

plt.title('Spectral Radiance vs. Wavelength',
          fontsize=20,
          family='Arial',
          fontweight='bold')

plt.xlabel('Wavelength/nm')
plt.ylabel('Solar irradiance/Wm^-3')

plt.legend()
plt.grid(linewidth=0.5)

plt.show()
