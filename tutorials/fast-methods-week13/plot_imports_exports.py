#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 14:36:00 2022

@author: david
"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# Split df into an imports and exports df
df = pd.read_csv('anc_loc_changes.csv')
imports_df = df[df['Destination'] == 'North Carolina']
exports_df = df[df['Origin'] == 'North Carolina']

# Plot imports into NC
sns.set_theme(style="darkgrid")
fig, ax = plt.subplots(figsize=(8, 5))
imports_by_state = imports_df['Origin'].value_counts()
sns.barplot(imports_by_state.index, imports_by_state.values)
ax.set_ylabel('Imports into NC')
plt.xticks(rotation=90) # rotate xtick labels so readable
fig.tight_layout()
fig.savefig('nc_imports_per_state.png', dpi=200)

# Plot exports into NC
fig, ax = plt.subplots(figsize=(8, 5))
exports_by_state = exports_df['Destination'].value_counts()
sns.barplot(exports_by_state.index, exports_by_state.values)
ax.set_ylabel('Exports into NC')
plt.xticks(rotation=90) # rotate xtick labels so readable
fig.tight_layout()
fig.savefig('nc_exports_per_state.png', dpi=200)