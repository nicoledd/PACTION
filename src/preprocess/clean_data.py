import pandas as pd
import re
from .prufer import enumPrufer
import networkx as nx



def readProportionMatrices(prop):
    return pd.read_csv(prop, sep='\t') #index_col = 'genotypes'


def getAllPossibleTrees(prop):
    trees = enumPrufer(len(prop))
    nxTrees = []
    for tree in trees:
        G = nx.DiGraph()
        G.add_edges_from(tree)
        nxTrees.append(G)
    return nxTrees


def getTreeEdges(edgefile, nodes):
    edgeList = []
    nodes = [str(node) for node in nodes]
    with open(edgefile, 'r') as inp:
        for line in inp:
            data = line.rstrip('\n').split('\t')
            node_out = data[0]
            node_in = data[1]

            if node_out == node_in:
                continue

            idx_out = nodes.index(node_out)
            idx_in = nodes.index(node_in)

            edgeList.append((idx_out, idx_in))

    G = nx.DiGraph()
    G.add_edges_from(edgeList)
    return G


def ifTrueTreesAreFoundIn(minSolutions, args):
    truesnv, truecna = getTrueTrees(args)
    for solution in minSolutions:
        if set(truesnv) == set(solution.snv_edges) and set(truecna) == set(solution.cna_edges):
            return True
    return False

def getTrueTrees(args):
    df_fsnv = readProportionMatrices(args.p[0])
    df_fcna = readProportionMatrices(args.p[1])
    snv_edges = getTreeEdges(args.truesnv, list(df_fsnv['genotypes']))
    cna_edges = getTreeEdges(args.truecna, list(df_fcna['genotypes']))
    return snv_edges, cna_edges

def getIndexOfTrueGraphs(minSolutions, args):
    truesnv, truecna = getTrueTrees(args)
    for i in range(len(minSolutions)):
        if set(minSolutions[i].snv_edges) == set(truesnv) and set(minSolutions[i].cna_edges) == set(truecna):
            return i


def processTreeFile(filename):
    edges = []
    f = open(filename, 'r')
    for line in f:
        nums = [int(s) for s in re.findall(r'\d+', line)]
        edges.append(sorted(((nums[0],nums[1]),(nums[2],nums[3]))))
    f.close()
    return edges

def processSolutionTree(tree):
    edges = []
    for edge in tree:
        edges.append(sorted(edge))
    return edges
