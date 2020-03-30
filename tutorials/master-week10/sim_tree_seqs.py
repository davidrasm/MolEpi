"""
Created on Fri Mar 27 08:30:47 2020

Simulate sequences along a phylogenetic tree using pyvolve

@author: david
"""
import baltic as bt
import pyvolve

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

treeFile = 'SEI2RTree_sim2_output.nexus'
outTree = 'SEI2RTree_sim2_relabeled.tre'
simTree = outTree
fastaFile = 'SEI2RTree_sim2_seqs.fasta'

relabel_tips(treeFile,outTree)

"User defined params"
mut_rate = 0.005 # rate per day
freqs = [0.25, 0.25, 0.25, 0.25]
seq_length = 1000
kappa = 2.75

"Read in phylogeny along which Pyvolve should simulate"
"Scale_tree sets absolute mutation rate"
my_tree = pyvolve.read_tree(file = simTree, scale_tree = mut_rate)
#pyvolve.print_tree(my_tree) # Print the parsed phylogeny

"Specify nucleotide substitution model with custom rates"
#custom_mu = {"AC":0.5, "AG":0.25, "AT":1.23, "CG":0.55, "CT":1.22, "GT":0.47}
#nuc_model = pyvolve.Model( "nucleotide", {"mu":custom_mu, "state_freqs":freqs} )

"Or just use an HKY model with kappa"
nuc_model = pyvolve.Model( "nucleotide", {"kappa":kappa, "state_freqs":freqs})

"Define a Partition object which evolves set # of positions according to my_model"
my_partition = pyvolve.Partition(models = nuc_model, size = seq_length)
#my_partition = pyvolve.Partition(models = nuc_model, root_sequence = "GATAGAAC") # Or with a root seq

"Define an Evolver instance to evolve a single partition"
my_evolver = pyvolve.Evolver(partitions = my_partition, tree = my_tree) 

"Evolve sequences with custom file names"
my_evolver(ratefile = "sim_ratefile.txt", infofile = "sim_infofile.txt", seqfile = fastaFile)



