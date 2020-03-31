"""
Created on Wed Mar 25 11:18:32 2020

Parse data from Master JSON output and plot

Usage: python plot_master_sim.py json_file png_file

@author: david
"""
import json
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

if len(sys.argv) < 2:
    print("Please specify JSON input file")
    sys.exit()
if len(sys.argv) == 2:
    json_file = sys.argv[1]
    png_file = 'SIR-curves.png'
if len(sys.argv) > 2:
    json_file = sys.argv[1]    
    png_file = sys.argv[2]

# Opening JSON file 
with open(json_file) as json_file: 
    data = json.load(json_file)

plotE = False
plotI = False
plotIl = False
plotIh = False

t = np.array(data['t'])
S = np.array(data['S'])
if "E" in data:
    E = np.array(data['E'])
    plotE = True
if "I" in data:
    I = np.array(data['I'])
    plotI = True
if "Il" in data:
    Il = np.array(data['Il'])
    plotIl = True
if "Ih" in data:
    Ih = np.array(data['Ih'])
    plotIh = True
    
# Set up plot
sns.set()
sns.set_context("talk")
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

# Plot dynamics
ax.plot(t, S, label="Susceptible")
if plotE:
    ax.plot(t, E, label="Exposed")
if plotI:
    ax.plot(t, I, label="Infected")
if plotIl:
    ax.plot(t, Il, label="Infected low")
if plotIh:
    ax.plot(t, Ih, label="Infected high")
ax.set_xlabel('Time')
ax.set_ylabel('Count')
ax.legend()
fig.tight_layout()
fig.savefig(png_file, dpi=200)