"""
Simple script for retreiving FASTA files from GenBank using NCBI's Entrez
Created on Wed May  8 09:59:01 2019
@author: david
"""

import pandas as pd
from Bio import SeqIO
from Bio import Entrez
    
def rename_seqs(records,df):
    
    "Iterate through list of records, renaming each as we go"
    new_records = []
    for record in records:
        accession = record.name # get accession number
        entry = df.index[df['Sequence Accession'] == accession].tolist() # get entry in df by accession number
        "Check to make sure entry exists and is unique"
        if len(entry) == 0:
            print("Entry not found for accession: " + accession)
        elif len(entry) > 1:
            print("More than one entry found for accession: " + accession)
        date = df.at[entry[0],'Collection Date'] # get date string
        record.id = accession + '_' + date # create new record id
        record.description = '' # set description blank
        new_records.append(record)

    return new_records
        
"Set ouput Fasta file name"
fasta_out = "influenzaA_H3N2_NC_2010-2019.fasta"

"We'll use Pandas to import our list of sequences" 
df = pd.read_csv('influenzaA_H3N2_NC_2010-2019.tsv','\t')
seq_list = df['Sequence Accession'] # get the seq accessions from the dataframe
seq_str = ",".join(seq_list) # a comma seperated list of all sequence accessions

"We'll use Biopython's Entrez inferface to download the sequences from NCBI GenBank"
Entrez.email = "your_email@important_univ.edu" # Always tell NCBI who you are
handle = Entrez.efetch(db="nucleotide", id=seq_str, rettype="gb", retmode="text")
records = SeqIO.parse(handle, "gb")

"Simplify seq names and add sampling dates"
records = rename_seqs(records,df)

"Print new records IDs"
for record in records:
    print(record.id)
    print(repr(record.seq))

SeqIO.write(records, fasta_out, "fasta")