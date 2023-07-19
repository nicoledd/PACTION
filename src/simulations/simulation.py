#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 2021

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import itertools
import math
import numpy as np
import random
from itertools import combinations

import networkx as nx
import itertools

from write import writeFullTree, writeFullTreeDOT, writeSnvTree, writeSnvTreeDOT, writeCnaTree, writeCnaTreeDOT, writeFullClones, writeSnvClones, writeCnaClones, writeNoisySnvClones, writeNoisyCnaClones
from tests import validateFullTree, validateCloneProps

def main(args):
    np.random.seed(args.s)
    S,C,G = buildAllTrees(args)
    writeAllTrees(args, G, S, C)
    writeAllDOTs(args, G, S, C)
    writeAllClones(args,G)
    validateFullTree(G)


def buildAllTrees(args):
    S = nx.DiGraph()
    S.add_node(0)
    C = nx.DiGraph()
    C.add_node(0)
    G = nx.DiGraph()
    G.add_node((0,0))

    nmuts = args.m + args.d - 2 # number of mutations
    mut_order = np.random.permutation(nmuts) # mutation order
    snv_counter = 0
    cna_counter = 0
    parent_has_child = False

    idx = 0
    while idx < len(mut_order):
        mut = mut_order[idx]

        # choose clone node at random
        parent_node = list(G.nodes)[np.random.randint(len(G.nodes))]
        
        if parent_node == (0,0) and parent_has_child:
            continue
        else:
            parent_has_child = True

        if mut < args.m - 1:
            # SNV mutation
            snv_counter += 1
            G.add_edge(parent_node, (snv_counter, parent_node[1]))
            S.add_edge(parent_node[0], snv_counter)
        else:
            # CNA mutation
            cna_counter += 1
            G.add_edge(parent_node, (parent_node[0], cna_counter))
            C.add_edge(parent_node[1], cna_counter)
        idx += 1
        
    return S, C, G


def writeAllTrees(args, G, S, C):
    writeFullTree(G, args.o)
    writeSnvTree(S, args.o)
    writeCnaTree(C, args.o)

def writeAllDOTs(args, G, S, C):
    writeFullTreeDOT(G, args.o)
    writeSnvTreeDOT(S, args.o)
    writeCnaTreeDOT(C, args.o)

def writeAllClones(args,G):
    df_clones = getClonesDf(args,G)
    writeFullClones(df_clones, args.o)
    writeSnvClones(df_clones, args.o)
    writeCnaClones(df_clones, args.o)
    writeNoisySnvClones(df_clones, args)
    writeNoisyCnaClones(df_clones, args)


def getClonesDf(args, G):
    clone_proportions = getCloneProportions(args)
    validateCloneProps(clone_proportions)
    data_clone = []
    for idx, clone in enumerate(list(G.nodes)):
        data_clone.append([clone, clone[0], clone[1]] + list(clone_proportions[idx, :]))
    df_clone = pd.DataFrame(data_clone, columns=['clone', 'snv', 'cna'] + [f'sample_{idx}' for idx in range(args.n)])
    return df_clone


def getCloneProportions(args):
    rows = args.m+args.d-1
    cols = args.n
    existence = np.random.binomial(1,args.p,(rows, cols))
    for c,s in enumerate(existence.sum(axis=0)):
        if s == 0:
            existence[np.random.randint(rows),c] = 1
    for r,s in enumerate(existence.sum(axis=1)):
        if s == 0:
            existence[r,np.random.randint(cols)] = 1

    clone_props = np.random.dirichlet([1]*rows, cols).transpose()
    clone_props[clone_props < 0.001] = 0.01
    filtered_props = clone_props * existence
    norm_props = filtered_props / filtered_props.sum(axis=0)
    return norm_props



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples [1]', default = 1)
    parser.add_argument('-m', type=int, help='number of SNV genotypes [5]', default = 5)
    parser.add_argument('-d', type=int, help='number of CNA genotypes [4]', default = 4)
    parser.add_argument('-o', type=str, help='output prefix', default='sample')
    parser.add_argument('-s', type=int, help='seed [0]', default = 0)
    parser.add_argument('-t', type=float, help='noise threshold [0.05]', default = 0.05)
    parser.add_argument('-p', type=float, help='probability of keeping each clone', default=0.2)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
