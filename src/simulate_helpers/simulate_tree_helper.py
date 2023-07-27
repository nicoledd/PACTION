import itertools
from itertools import combinations
import math
import random
import networkx as nx
import numpy as np


def buildTree(args):
    G = nx.DiGraph()
    numModes = len(args.m)
    G.add_node((0,)*numModes)
    nmuts = sum(args.m) - numModes
    mutDict = getMutDict(args)
    mutOrder = np.random.permutation(nmuts)
    mutCounter = [0 for _ in range(numModes)]
    parentHasChild = False
    idx = 0
    while idx < len(mutOrder):
        mut = mutOrder[idx]
        parentNode = list(G.nodes)[np.random.randint(len(G.nodes))]
        if parentNode == (0,)*numModes and parentHasChild:
            continue
        else:
            parentHasChild = True
        modeIdx = mutDict[mut]
        mutCounter[modeIdx] += 1
        newNode = list(parentNode)[:modeIdx] + [mutCounter[modeIdx]] + list(parentNode)[modeIdx+1:]
        G.add_edge(parentNode, tuple(newNode))
        idx += 1
    return G


def getMutDict(args):
    mutDict = {}
    currMut = 0
    for i in range(len(args.m)):
        for j in range(args.m[i]-1):
            mutDict[currMut] = i
            currMut += 1
    return mutDict

def extractSubtreeWithMode(modeIdx, G):
    S = nx.DiGraph()
    for u,v in G.edges:
        newU, newV = u[modeIdx], v[modeIdx]
        if newU != newV:
            S.add_edge(newU, newV)
    return S