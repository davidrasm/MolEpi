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

fasta_file = 'WAhalf_cov2020_aligned.fasta'
records = SeqIO.parse(fasta_file, "fasta")

records_out = [];
for record in records:
	record.id = record.id.replace('|', '_') + '_Il'
	record.id = record.id.replace('2020_EPI_ISL_', '')
	print(record.id)
	record.name = ''
	record.description = ''
	records_out.append(record)

SeqIO.write(records_out, "WAhalf_cov2020_relabeled.fasta", "fasta")