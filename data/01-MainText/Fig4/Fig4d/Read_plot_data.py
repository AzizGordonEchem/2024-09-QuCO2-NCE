# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 18:55:41 2024

@author: amini
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

# Set font properties
plt.rc('font', family='Arial', size=14)
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.grid'] = False

# Create the plot figure and axis
fig, ax = plt.subplots(figsize=(8, 6))

# Enable minor ticks and customize tick parameters
ax.minorticks_on()
ax.tick_params(which='both', direction='in', bottom=True, top=True, left=True, right=True)

# Draw a vertical line at x=600 nm
plt.axvline(x=600, color='gray', linestyle='--')

# Load the data from the file
data = np.genfromtxt('Emission_data.txt', unpack=True).T

# Plot the data
plt.plot(data[0, 1:], data[5, 1:], 'o-', markersize=8, alpha=0.4, linewidth=3, color='darkblue', label='100 mM - Adduct')
plt.plot(data[0, 1:], data[3, 1:], 'o-', markersize=8, alpha=0.4, linewidth=3, color='darkred', label='100 mM - Reduced')
plt.plot(data[0, 1:], data[1, 1:], 'o-', markersize=8, alpha=0.4, linewidth=3, color='darkorange', label='100 mM - Oxidized')



a = data[0,1:]
b = data[5, 1:]
c = [a,b]

#plt.plot(a,b)

with open("100 mM - Adduct - 475.txt", "w") as file:
    for x in zip(*c):
        file.write("{0}\t{1}\n".format(*x))
        
a = data[0,1:]
b = data[3, 1:]
c = [a,b]

#plt.plot(a,b)

with open("100 mM - Reduced - 475.txt", "w") as file:
    for x in zip(*c):
        file.write("{0}\t{1}\n".format(*x))        

a = data[0,1:]
b = data[1, 1:]
c = [a,b]

#plt.plot(a,b)

with open("100 mM - Oxidized - 475.txt", "w") as file:
    for x in zip(*c):
        file.write("{0}\t{1}\n".format(*x))      







# Set the axis labels
plt.xlabel('$\mathregular{Wavelength \, \, (nm)}$', labelpad=6)
plt.ylabel('$\mathregular{Emission \,\, (Excited \, \, at \, \, 475 \,\, nm)}$', labelpad=6)

# Set the limits for the x and y axes
plt.xlim(510, 800)
plt.ylim(0, 2000)

# Add the legend and save the figure
plt.legend()
fig.savefig('Figure4d_300dpi.jpg', dpi=600)

# Display the plot
plt.show()