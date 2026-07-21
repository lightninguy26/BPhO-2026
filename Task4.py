# plotting stopping voltage vs. frequency for various metals 

import matplotlib.pyplot as plt
import numpy as np

# constants/parameters
workfunctions={
    'Ag':4.3,
    'Ca':2.9,
    'Au':5.1,
    'Cu':4.7,
    'Sn':4.4,
    'Cs':2.1,
    'Pt':5.6,
    'Ni':4.6,
    'Na':2.4
}
frequencies=np.linspace(0,2,num=300)
# frequencies in 10^15 Hz
e=1.602*10**-19
h=6.626*10**-34

for label,W in workfunctions.items():

    W*=10**-19
    f=frequencies * 10**15
    # apply the formula
    stopping_voltages = (h/e) * f - (W/e)

    plt.plot(f,stopping_voltages,label=label,linewidth=0.75)

plt.ylim(bottom=0)
plt.title('Photoelectric effect')
# vertical lines representing colours in the visible light spectrum
plt.vlines(x=[0.44*10**15,0.52*10**15,0.565*10**15,0.635*10**15],ymin=[0,0,0,0],ymax=[9,9,9,9],colors=['red','yellow','green','blue'],linestyles='--',lw=0.6)
plt.annotate('Visible spectrum',xy=(5.5e14,5.5),xytext=(-5e13,6),arrowprops={})

plt.text(1.2e15,0.5,r'$f_0 = \frac{W}{h}$',fontweight='bold')
plt.ylabel('Stopping voltage/volts')
plt.xlabel('Frequency/10$^{15}$Hz')
plt.legend()

plt.show()


