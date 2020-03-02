"""

Plot Bayesian Skyline plots using matplotlib and seaborn

@author: david
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def convert_dates(df):
    
    """
    Influenza seasons begins at week 40 according to CDC
    So for a given season:
        Weeks >40 fall in previous year
        Weeks <40 fall in next year
    
    WARNING: This will break for seasons prior to 2000
    """
    
    dates = []
    for index, row in df.iterrows():            
        season = row['SEASON']
        week = row['WEEK']
        years = season.split('-')
        if week >= 40:
            d = float(years[0]) + (week/52)
        else:
            d = 2000 + float(years[1]) + (week/52)
        dates.append(d)
        
    return dates


"Read in data from Skyline analysis"
tsv_file = 'nc_flu_skyline_data.tsv'
df = pd.read_table(tsv_file, sep="\t")

"Set up subplots"
fig, axs = plt.subplots(2, 1)

"Plot skyline in first subplot"
sns.set()
sns.lineplot(x="Time", y="Median", data=df, ax=axs[0])
axs[0].set_ylabel('Effective size', fontsize=12)
axs[0].set_xlabel('') # no label
axs[0].fill_between(df['Time'], df['Lower'], df['Upper'], alpha=.3) # alpha val sets transparency
date_lim = axs[0].get_xlim()

"Read in ILI data from CDC"
csv_file = 'cdc_state_ILI_reports.csv'
df = pd.read_table(csv_file, sep=",")
df = df[df['SUB AREA'] == 'North Carolina'] # get entries for North Carolina
df['DATES'] = convert_dates(df) # convert flu season weeks to decimal years

"Plot ILI data in second subplot"
sns.lineplot(x="DATES", y="PERCENT P&I", data=df, color='darkorange', ax=axs[1])
axs[1].set_xlabel('Time')
axs[1].set_xlim(date_lim) # set xlim so subplots have same x-axis

"Save image"
fig.set_size_inches(6, 6)
plt.show()
img_file = 'nc_flu_skyline_vs_ILI.png'
fig.savefig(img_file, dpi=200)