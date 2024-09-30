# -*- coding: utf-8 -*-


# Import necessary libraries
import os
import dill
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
from matplotlib.legend_handler import HandlerBase

# Set font parameters
font = {'family': 'Arial', 'size': 16}
plt.rc('font', **font)

# Set figure and axis properties
plt.rcParams['axes.grid'] = False
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

# Create the plot
fig, ax1 = plt.subplots(1, 1, figsize=(8, 6))
ax1.minorticks_on()
ax1.tick_params(which='both', bottom=True, top=True, left=True, right=True)
ax1.tick_params(which='minor', bottom=True, top=True, left=True, right=True)
ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

# Load data from files
folder = ''
E, I, I_sims = [], [], []

# Iterate through files in the folder
for file in os.listdir(folder or '.'):
    if file.endswith('data_E.txt'):
        with open(os.path.join(folder, file), 'rb') as f:
            E.append(dill.load(f))
    elif file.endswith('data_I.txt'):
        with open(os.path.join(folder, file), 'rb') as f:
            I.append(dill.load(f))
    elif file.endswith('sim_I.txt'):
        with open(os.path.join(folder, file), 'rb') as f:
            I_sims.append(dill.load(f))

# Set scan rates and colors for the plot
scan_rates = [0.1, 0.5, 0.9]
colors = cm.tab20b(np.linspace(0, 2000, 25000))

# Plot experimental and simulation data
for i in range(3):
    adjusted_E = np.asarray(E[i]) - 0.1  # Adjust potential data
    plt.plot(adjusted_E, I[i], linewidth=2, color=colors[9 - 4 * i], label='Experiment')
    plt.plot(adjusted_E, I_sims[i], '--', linewidth=2, color=colors[9 - 4 * i], label='Simulation')

# Custom legend handler for experiment vs simulation
class AnyObjectHandler(HandlerBase):
    def create_artists(self, legend, orig_handle, x0, y0, width, height, fontsize, trans):
        l1 = plt.Line2D([x0, y0 + width], [0.7 * height, 0.7 * height], linestyle='--', color=orig_handle[0])
        l2 = plt.Line2D([x0, y0 + width], [0.3 * height, 0.3 * height], linestyle='-', color=orig_handle[0])
        return [l1, l2]

# Add legend for scan rates
plt.legend([(colors[9], colors[9]), (colors[5], colors[5]), (colors[1], colors[1])],
           [f'{scan_rates[0]} V/s', f'{scan_rates[1]} V/s', f'{scan_rates[2]} V/s'],
           handler_map={tuple: AnyObjectHandler()}, title='Experimental (-) and Simulated (--)')

# Set axis labels
plt.ylabel('$\mathregular{Current \, \, (mA)}$', labelpad=2)
plt.xlabel('$\mathregular{Potential \,\, Versus \, \, SHE \, \, (V)}$', labelpad=0.2)

# Define helper functions to annotate min and max points
def annot_max(x, y, ax=None):
    xmax = x[np.argmax(y)]
    ymax = y.max()
    text = "x={:.3f}, y={:.3f}".format(xmax, ymax)
    if not ax:
        ax = plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops = dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data', textcoords="axes fraction", arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    # ax.annotate(text, xy=(xmax, ymax), xytext=(xmax + 0.7, ymax + 0.7), **kw)

def annot_min(x, y, ax=None):
    xmin = x[np.argmin(y)]
    ymin = y.min()
    text = "x={:.3f}, y={:.3f}".format(xmin, ymin)
    if not ax:
        ax = plt.gca()
    bbox_props = dict(boxstyle="square,pad=0.3", fc="w", ec="k", lw=0.72)
    arrowprops = dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=60")
    kw = dict(xycoords='data', textcoords="axes fraction", arrowprops=arrowprops, bbox=bbox_props, ha="right", va="top")
    # ax.annotate(text, xy=(xmin, ymin), xytext=(xmin + 0.9, ymin + 0.7), **kw)

# Set axis limits and show the plot
plt.xlim(-0.95, 0.75)
plt.show()

# Save the figure
fig.savefig('figure2d-pdf.pdf', dpi=600)