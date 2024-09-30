# -*- coding: utf-8 -*-


import matplotlib.pyplot as plt
import numpy as np

# Configure plot settings
plt.rc('font', family='Arial', size=16)
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams['axes.grid'] = False

# Set figure dimensions
fig, ax1 = plt.subplots(1, 1)
fig.set_figheight(6)
fig.set_figwidth(8)

# Enable minor ticks
ax1.minorticks_on()
ax1.tick_params(bottom=True, top=True, left=True, right=True, which='both')
ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

# Define equilibrium constant values
K1_values = np.logspace(4, 24, num=120)

# Initialize lists to store concentration data
dissociation_constant = []
adduct_conc = []
adduct2_conc = []
bicarbonate_conc = []
carbonate_conc = []
co2aq = []
sumsum = []

# Loop over K1 values and load corresponding data
for K1 in K1_values:
    filename = f"concentrationofcarbontt_{K1:.6e}.txt"
    
    # Load the data from the file
    try:
        a_k1 = np.loadtxt(filename, skiprows=1, unpack=True)
    except FileNotFoundError:
        print(f"File {filename} not found. Skipping this K1 value.")
        continue
    
    # Append data to corresponding lists
    dissociation_constant.append(a_k1[0])
    co2aq.append(a_k1[1])
    bicarbonate_conc.append(a_k1[2])
    carbonate_conc.append(a_k1[3])
    adduct_conc.append(a_k1[4])
    adduct2_conc.append(2 * a_k1[4])
    sumsum.append(a_k1[1] + a_k1[2] + a_k1[3] + 2 * a_k1[4])

# Plot settings
line_width = 3
alpha_value = 1

# Plot the data on a semi-logarithmic plot
ax1.semilogx(dissociation_constant, co2aq, linewidth=line_width, color='green', alpha=alpha_value, label='$\mathregular{CO_2(aq)}$')
ax1.semilogx(dissociation_constant, bicarbonate_conc, linewidth=line_width, color='black', alpha=alpha_value, label='$\mathregular{HCO_3^-}$')
ax1.semilogx(dissociation_constant, carbonate_conc, linewidth=line_width, color='purple', alpha=alpha_value, label='$\mathregular{CO_3^{2-}}$')
ax1.semilogx(dissociation_constant, adduct_conc, linewidth=line_width, color='blue', alpha=alpha_value, label='$\mathregular{AQ(CO_2)_2^{2-}}$')

# Annotations for each curve
plt.text(1e19, 0.11, '$\mathregular{[ AQ(CO_2)_2^{2-} ]}$', color='blue', alpha=alpha_value)
plt.text(1e7, 0.175, '$\mathregular{[HCO_3^-}]$', color='black', alpha=alpha_value)
plt.text(1e20, 0.0125, '$\mathregular{[CO_2(aq)}]$', color='green', alpha=alpha_value)
plt.text(1e7, 0.0125, '$\mathregular{[CO_3^{2-}}]$', color='purple', alpha=alpha_value)

# Set axis limits
ax1.set_xlim(1e4, 1e24)

# Set labels for the axes
plt.xlabel('Adduct formation equilibrium constant $\mathregular{(1/M^2)}$')
plt.ylabel('Concentration (M)')

# Highlight specific regions on the plot with shaded areas
plt.axvspan(10**4, 10**12, facecolor='gray', alpha=0.15)
plt.axvspan(10**12, 10**17, facecolor='red', alpha=0.15)
plt.axvspan(10**17, 10**24, facecolor='blue', alpha=0.1)

# Show the plot
plt.show()

# Save the figure to a file
fig.savefig('Figure1_conc__equilibriumconstant0p5bar.tif')

