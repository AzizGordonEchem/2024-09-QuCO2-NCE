# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 18:22:58 2023

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
        'size'   : 16}

plt.rc('font', **font)
size = 14
plt.figsize=(18,10)
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



data = np.genfromtxt('22peaq510.txt', unpack=True).T
plt.plot(data[0,1:],data[3,1:],'o-', markersize = '8', alpha=.4,linewidth = '2',color = 'Darkblue',label = '50 mM - Adduct')
plt.plot(data[0,1:],data[1,1:],'o-', markersize = '8', alpha=.4,linewidth = '2',color = 'DarkRed',label = '50 mM - Reduced')
plt.plot(data[0,1:],data[2,1:],'o-', markersize = '8' , alpha=.4,linewidth = '2',color = 'DarkOrange',label = '50 mM - Oxidized')
plt.xlabel('$\mathregular{Wavelength \, \, (nm)}$', labelpad=6,family='Arial')
plt.ylabel('$\mathregular{Emmission \, \, Spectrum \,\, (Excited \, \, at \, \, 510 \,\, nm)}$', labelpad=6,family='Arial')
plt.xlim(550,800)
plt.ylim(0,4000)

plt.legend()
fig.savefig('22peaq510.tif')
plt.show()