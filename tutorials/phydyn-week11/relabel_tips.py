"""
Created on Fri Mar 27 08:30:47 2020

Relabel tips

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
            curr_name = curr_name.replace('|', '_') + '_Il'
            curr_name = curr_name.replace('2020_EPI_ISL_', '')
            k.numName = curr_name
            print(k.numName)
            
    newick = myTree.toString(traits=[],numName=True,nexus=False)
    newick = newick.replace('\'','')
    
    tfile = open(outTree, "w")
    tfile.write(newick)
    tfile.close()


# if len(sys.argv) < 2:
#     print("Please specify input NEXUS file")
#     sys.exit()
# if len(sys.argv) == 2:
#     treeFile = sys.argv[1]
#     stem = treeFile.split(".")[0]
#     relabledTreeFile = stem + '_relabeled.tre'
#     fastaFile = stem + '.fasta'

"Relabel tips in tree so they contain tip types and times"
treeFile = 'WAhalf_cov2020_tree10001.tre'
relabledTreeFile = 'WAhalf_cov2020_tree10001_relabeled.tre'
relabel_tips(treeFile,relabledTreeFile)





