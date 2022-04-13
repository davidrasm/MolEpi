"""
Created on Wed Aug  4 12:59:55 2021

Extract dates from GISAID metadata file and write to tsv file

@author: david
"""

import pandas as pd    

# Set meta data file and output file
meta_file = "metadata_NC-10k_USA-10k.tsv"
dates_file = "dates.tsv" # output file

# Read in meta data and put in pandas dataframe
df = pd.read_csv(meta_file,sep="\t")

# Extract only Accession ID and Collection date
df = df[['Accession ID','Collection date']]

# Write to tab-delimited csv file
df.to_csv(dates_file,sep="\t",index=False)


