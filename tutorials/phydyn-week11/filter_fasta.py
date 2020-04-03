"""
Filter sequences in fasta file
Created on Wed May  8 09:59:01 2019
@author: david
"""

from Bio import SeqIO
from Bio import Entrez
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

def read_accessions(fp):
    "Reads list from file"
    with open(fp) as acc_lines:
        return [line.strip() for line in acc_lines]

fasta_file = 'gisaid_cov2020_sequences.fasta'
records = SeqIO.parse(fasta_file, "fasta")

match = 'USA/WA'
#match = 'Wuhan'

sprob = 0.5 # sample inclusion prob

records_out = [];
for record in records:
    if match in record.name:
        print(record.name)
        if random.random() < sprob:
            records_out.append(record)
    
#    print(record.id)
#    print(record.name)
#    print(record.description)
#    print(repr(record.seq))

SeqIO.write(records_out, "WAhalf_cov2020_sequences.fasta", "fasta")