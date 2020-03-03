"""
Created on Tue Mar  3 07:44:11 2020

Simulate sequences along a phylogenetic tree using pyvolve

@author: david
"""
import pyvolve

"User defined params"
mut_rate = 0.005
freqs = [0.25, 0.25, 0.25, 0.25]
seq_length = 1000
kappa = 2.75

"Read in phylogeny along which Pyvolve should simulate"
"Scale_tree sets absolute mutation rate"
my_tree = pyvolve.read_tree(file = "AMR-sim.tre", scale_tree = mut_rate)
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
my_evolver(ratefile = "AMR_ratefile.txt", infofile = "AMR_infofile.txt", seqfile = "AMR-seqsim.fasta" )