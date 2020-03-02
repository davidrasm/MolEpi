"""
Created on Tue Feb 25 14:21:08 2020

# Convert between phylo tree formats
# Writes BEAST-style NEXUS tree files

@author: david
"""

from Bio import Phylo
import baltic as bt
import re

# '1-2936','2937-4936', '4937-6870','6871-8473','8474-19706','19707-20428'
in_tree = 'sars-like-CoVs-sub_19707-20428.tre'
out_tree = 'sars-like-CoVs-sub_19707-20428.nexus'

# Basic conversion with Phylo module doesn't work here because we need BEAST-style NEXUS files
#Phylo.convert(in_tree, 'newick', out_tree, 'nexus')

myTree=bt.loadNewick(in_tree, absoluteTime=False)
#myTree.setAbsoluteTime(2019.00) # need to set this to time of last sampled tip
myTree.traverse_tree() ## required to set heights
myTree.treeStats() ## report stats about tree
names = []
for idx,k in enumerate(myTree.Objects): ## iterate over a flat list of branches
    if k.branchType=='leaf':
        curr_name = k.numName
        names.append(curr_name)
        
#print(names)

date_str = '_2020.00'

# Write taxa names
nex=open(out_tree,"w")
nex.write("#NEXUS\n")
nex.write("Begin taxa;\n")
nex.write("\tDimensions ntax=" + str(len(names)) + ";\n")	
nex.write("\t\tTaxlabels\n")
for n in names:
    nex.write("\t\t\t" + n + date_str + "\n")    
nex.write("\t\t\t;\n")
nex.write("End;\n")		

# Write translation 	
nex.write("Begin trees;\n")	
nex.write("\tTranslate\n")	
for idx,n in enumerate(names):
    nex.write("\t\t" + str(idx+1) + ' ' + n + date_str + "\n")
nex.write(";\n")

# Write tree
with open(in_tree, 'r') as file:
    tree_str = file.read().replace('\n', '')
for idx,n in enumerate(names):
    tree_str = re.sub(n, str(idx+1), tree_str)    
nex.write("tree TREE1 = " + tree_str + "\n")
nex.write("End;\n")