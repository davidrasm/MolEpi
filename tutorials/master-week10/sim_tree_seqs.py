"""
Created on Fri Mar 27 08:30:47 2020

Simulate sequences along a phylogenetic tree using pyvolve

@author: david
"""
import baltic as bt
import pyvolve
import sys

def relabel_tips(treeFile,outTree):

    myTree=bt.loadNewick(treeFile, absoluteTime=False)
    #myTree.setAbsoluteTime(2019.00) # need to set this to time of last sampled tip
    myTree.traverse_tree() ## required to set heights
    myTree.treeStats() ## report stats about tree
    
    for k in myTree.Objects: ## iterate over a flat list of branches  
        if k.branchType=='leaf':
            curr_name = k.numName
            if 'type' in k.traits:
                ttype = k.traits['type']
            else:
                ttype = 'Unknown' #Assign tip unknown type if not given
            if 'time' in k.traits:
                time = k.traits['time']
            else:
                time = 'Unknown'
            k.numName = curr_name + '_' + str(time) + '_' + ttype
            print(k.numName)
            
    newick = myTree.toString(traits=[],numName=True,nexus=False)
    newick = newick.replace('\'','')
    
    tfile = open(outTree, "w")
    tfile.write(newick)
    tfile.close()


if len(sys.argv) < 2:
    print("Please specify input NEXUS file")
    sys.exit()
if len(sys.argv) == 2:
    treeFile = sys.argv[1]
    stem = treeFile.split(".")[0]
    relabledTreeFile = stem + '_relabeled.tre'
    fastaFile = stem + '.fasta'

"Relabel tips in tree so they contain tip types and times"
relabel_tips(treeFile,relabledTreeFile)

"User defined params"
mut_rate = 0.005 # rate per day
freqs = [0.25, 0.25, 0.25, 0.25] # equilibrium nuc freqs
seq_length = 2000 # length of sim sequences
kappa = 2.75 # kappa parameter in HKY model

"Read in phylogeny along which Pyvolve should simulate"
my_tree = pyvolve.read_tree(file = relabledTreeFile, scale_tree = mut_rate) #Scale_tree sets absolute mutation rate
#pyvolve.print_tree(my_tree) # Print the parsed phylogeny

"HKY model with kappa"
nuc_model = pyvolve.Model( "nucleotide", {"kappa":kappa, "state_freqs":freqs})

"Define a Partition object which evolves set # of positions according to my_model"
my_partition = pyvolve.Partition(models = nuc_model, size = seq_length)

"Define an Evolver instance to evolve a single partition"
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree) 

"Evolve sequences with custom file names"
my_evolver(ratefile = "sim_ratefile.txt", infofile = "sim_infofile.txt", seqfile = fastaFile)



