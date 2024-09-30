# -*- coding: utf-8 -*-
"""
Created on Wed May 24 14:36:46 2023

@author: Kiana
"""

# -*- coding: utf-8 -*-
"""
Created on Fri May 19 15:54:34 2023

@author: Kiana
"""


import numpy as np
from numpy import pi, r_
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.ticker as ticker


font = {'family' : 'Arial',
        'size'   : 14}

plt.rc('font', **font)
size = 14
plt.rcParams['axes.grid'] = False
fig, ax = plt.subplots(1, 1)
fig.set_figheight(6)
fig.set_figwidth(8) 
ax.yaxis.get_ticklocs(minor=True) 
ax.minorticks_on()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(which = 'minor',bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

plt.axvline(x=520, color='gray', linestyle='--')


data = np.genfromtxt('Emission_data.txt', unpack=True).T
plt.plot(data[0,1:],data[7,1:],'o-', markersize = '8', alpha=.4,linewidth = '3',color = 'Darkblue',label = '100 mM - Adduct')
plt.plot(data[0,1:],data[4,1:],'o-', markersize = '8', alpha=.4,linewidth = '3',color = 'Darkred',label = '100 mM - Reduced')
plt.plot(data[0,1:],data[1,1:],'o-', markersize = '8' , alpha=.4,linewidth = '3',color = 'Darkorange',label = '100 mM - Oxidized')
plt.xlabel('$\mathregular{Wavelength \, \, (nm)}$', labelpad=6,family='Arial')
plt.ylabel('$\mathregular{Emmission  \,\, (Excited \, \, at \, \, 405 \,\, nm)}$', labelpad=6,family='Arial')


a = data[0,1:]
b = data[7,1:]
c = [a,b]

#plt.plot(a,b)

with open("100 mM - Adduct - 405.txt", "w") as file:
    for x in zip(*c):
        file.write("{0}\t{1}\n".format(*x))
        
a = data[0,1:]
b = data[4,1:]
c = [a,b]

#plt.plot(a,b)

with open("100 mM - Reduced - 405.txt", "w") as file:
    for x in zip(*c):
        file.write("{0}\t{1}\n".format(*x))        

a = data[0,1:]
b = data[1,1:]
c = [a,b]

#plt.plot(a,b)

with open("100 mM - Oxidized - 405.txt", "w") as file:
    for x in zip(*c):
        file.write("{0}\t{1}\n".format(*x))      









plt.xlim(437,800)
plt.ylim(0,3100)
plt.legend()
fig.savefig('Figure4c_300dpi.jpg', dpi=600)
plt.show()