# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 18:34:09 2023

@author: Kiana
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 18:04:26 2023

@author: Kiana
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 18:00:30 2023

@author: Kiana
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 17:47:48 2023

@author: Kiana
"""

# -*- coding: utf-8 -*-
"""
Created on Wed May 24 15:42:05 2023

@author: Kiana
"""

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
size = 20


plt.rcParams['axes.grid'] = False
fig = plt.figure()
fig, ax = plt.subplots(1,1, figsize=(8,6))
ax.yaxis.get_ticklocs(minor=True) 
ax.minorticks_on()
ax.tick_params(bottom=True, top=True, left=True, right=True)
ax.tick_params(which = 'minor',bottom=True, top=True, left=True, right=True)
ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"



data = np.genfromtxt('ph14.txt', unpack=True).T
plt.plot(data[0,1:],data[5,1:],'o-', markersize = '8', alpha=.5,linewidth = '2',color = 'deeppink',label = 'pH 14 with dry ice')
plt.plot(data[0,1:],data[4,1:],'o-', markersize = '8', alpha=.5,linewidth = '2',color = 'blue',label = 'pH 14 ')

plt.xlabel('$\mathregular{Wavelength \, \, (nm)}$', labelpad=6,family='Arial')
plt.ylabel('$\mathregular{Emmission \, \, Spectrum \,\, (Excited \, \, at \, \, 510 \,\, nm)}$', labelpad=6,family='Arial')
plt.xlim(550,800)
plt.ylim(0,500)
plt.legend()
plt.show()

fig.savefig('ph14510.tif')
plt.show()