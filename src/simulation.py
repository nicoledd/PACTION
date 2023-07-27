#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 2021

@author: Palash Sashittal
"""

import pandas as pd
import sys
import argparse
import numpy as np
import networkx as nx

from simulate_helpers.simulate_tree_helper import buildTree
from simulate_helpers.simulate_props_helper import getClonesDf
from postprocess.write import writeTree, writeClones
from test_suite.tests import validateTree, validateProps


def main(args):
    np.random.seed(args.s)
    G = simulateTree(args)
    U = simulateProps(args, list(G.nodes))


def simulateTree(args):
    G = buildTree(args)
    validateTree(G)
    writeTree(G, args.o)
    for modeIdx in range(len(args.m)):
        T = extractSubtreeWithMode(modeIdx, G)
        validateTree(T)
        writeTree(T, args.o, '_mode'+str(modeIdx))
    return G

def extractSubtreeWithMode(modeIdx, G):
    S = nx.DiGraph()
    for u,v in G.edges:
        newU, newV = u[modeIdx], v[modeIdx]
        if newU != newV:
            S.add_edge(newU, newV)
    return S


def simulateProps(args, clones):
    dfClones = getClonesDf(args,clones)
    writeClones(dfClones, args.o)
    for modeIdx in range(len(args.m)):
        dfSubclones = extractSubclonesWithMode(modeIdx, dfClones)
        dfSubclones['genotypes']=dfSubclones.index
        cols = dfSubclones.columns.to_list()
        cols = cols[-1:] + cols[:-1]
        dfSubclones = dfSubclones[cols]
        writeClones(dfSubclones, args.o, '_mode'+str(modeIdx))

def extractSubclonesWithMode(modeIdx, dfClones):
    dfSubclones = dfClones.copy()
    dfSubclones['genotypes'] = dfSubclones['clone'].apply(lambda x: x[modeIdx])
    dfSubclones = dfSubclones.groupby('genotypes').sum(numeric_only=True)
    return dfSubclones



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, help='number of samples [1]', default = 1)
    parser.add_argument('-m', type=int, help='number of genotypes per mode [4]', nargs='*')
    parser.add_argument('-o', type=str, help='output prefix', default='sample')
    parser.add_argument('-s', type=int, help='seed [0]', default = 0)
    parser.add_argument('-t', type=float, help='noise threshold [0]', default = 0)
    parser.add_argument('-p', type=float, help='probability of keeping each clone [0.2]', default=0.2)
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])

    main(args)
