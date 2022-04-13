#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 10:41:23 2022

@author: david
"""

import dendropy
import numpy as np
import pandas as pd

def reconstruct_MP(tree,feature_dic):
    
    """
        Reconstruct ancestral states using Sankoff's max parsimony algorithm (see Felsenstein p15-18)
        Parsimony scores are computed assuming all transitions have a cost of one.
        Tip/ancestral features are input and output as dictionaries
        Internal nodes are labeled as n<X> where X is an int determined by the position of the node in a pre-order traversal
    """
    
    state_set = set([feature_dic[node.taxon.label.replace(' ','_')] for node in tree.leaf_node_iter()])
    states = len(state_set)
    
    "Add a state2int dict map so this works with non-integer data types"
    state2int = {state:index for index, state in enumerate(state_set)}
    int2state = {index:state for index, state in enumerate(state_set)}
    
    "Post-order traversal to compute costs in terms of required state transitions"
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            costs = [np.inf]*states
            tip_state = state2int[feature_dic[node.taxon.label.replace(' ','_')]] # feature_dic[node.name]
            costs[tip_state] = 0
            node.costs = costs
        else:
            costs = [0]*states
            for i in range(states):
                child_costs = []
                for child in node.child_node_iter():
                    temp_costs = [0]*states
                    for j in range(states):
                        temp_costs[j] = child.costs[j]
                        if i != j:
                            temp_costs[j] = temp_costs[j] + 1 # add cost for transitioning between i and j
                    child_costs.append(temp_costs)
                costs[i] = sum([min(c) for c in child_costs])
            node.costs = costs
    
    "Pre-order traversal to select anc states based on least cost parsimony score"
    anc_dic = {}
    internal_cntr = 0 # counter for internal nodes encountered
    for node in tree.preorder_node_iter():
        costs = node.costs
        if node == tree.seed_node: # i.e. root
            root_state = costs.index(min(costs)) # or np.argmin(node.costs)
            node.state = root_state
            anc_dic['root'] = root_state
        else:
            parent_state = node.parent_node.state
            least_cost = min(costs)
            least_cost_state = costs.index(least_cost)
            if parent_state == least_cost_state:
                anc_state = parent_state
            else:
                parent_st_cost = costs[parent_state]
                if parent_st_cost < (least_cost+1):
                    anc_state = parent_state # if parent state costs less than transitioning to least cost state
                else:
                    anc_state = least_cost_state
            node.state = anc_state
            if node.is_leaf():
                anc_dic[node.taxon.label.replace(' ','_')] = anc_state
            else:
                name = 'n' + str(internal_cntr)
                anc_dic[name] = anc_state
                internal_cntr += 1
    
    "Covert from integers back to original states"
    for k,v in anc_dic.items():
        anc_dic[k] = int2state[v]
        
    for node in tree.preorder_node_iter():
        node.state = int2state[node.state]

    "Should return tree and anc_dic"
    return tree

def ancestral_changes(tree):
    
    """
        Iterate over tree, identifying all branches in which the parent and child state differs
    """
    
    changes = 0
    times = []
    origins = []
    destinations = [] 
    tree.calc_node_ages(ultrametricity_precision=False)
    for node in tree.postorder_node_iter():
        if node != tree.seed_node: # i.e. if not root
            parent_state = node.parent_node.state
            if node.state != parent_state:
                changes += 1
                times.append(node.age)
                origins.append(parent_state)
                destinations.append(node.state)
                
    print("Total number of state changes: " + str(changes))
    
    # Place events in dataframe
    df = pd.DataFrame({'EventTime' : times})
    df['Origin'] = origins
    df['Destination'] = destinations
    
    return df

#Read in meta data
meta_file = "metadata_NC-10k_USA-10k.tsv"
meta_df = pd.read_csv(meta_file,sep="\t")

#Get state-level location from location column
meta_df['State'] = meta_df['Location'].str.split('/').apply(lambda x: x[2].strip() if (len(x) > 2) else 'None') # direct way
state_counts = meta_df['State'].value_counts()

# Load in tree using DendroPy 
tree_file = 'hcov_NC10k_USA10k_dated.tre'
taxa = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(file=open(tree_file, 'r'), schema="newick", rooting="default-rooted", taxon_namespace=taxa) 

#Make a dict mapping tip names (key) to tip states (value)
tip_state_dic = {}
for index, row in meta_df.iterrows():
    tip_state_dic[row['Accession ID']] = row['State']

# Run ancestral location reconstruction      
tree = reconstruct_MP(tree,tip_state_dic)

# Find changes in state along tree
changes_df = ancestral_changes(tree)
changes_df.to_csv('anc_loc_changes.csv',index=False)

