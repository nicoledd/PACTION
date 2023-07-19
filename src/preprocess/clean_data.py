import pandas as pd
import re
from .prufer import enumPrufer



def readProportionMatrices(args):
    df_fsnv = pd.read_csv(args.fsnv, sep=',', index_col = 'genotypes')
    df_fcna = pd.read_csv(args.fcna, sep=',', index_col = 'genotypes')
    return df_fsnv, df_fcna


def getAllPossibleSnvCnaTrees(args):
    df_fsnv, df_fcna = readProportionMatrices(args)
    nsnv = len(df_fsnv)
    ncna = len(df_fcna)
    snvTrees = enumPrufer(nsnv)
    cnaTrees = enumPrufer(ncna)
    return snvTrees, cnaTrees


def getTreeEdges(edgefile, nodes):
    edgeList = []
    nodes = [str(node) for node in nodes]
    with open(edgefile, 'r') as inp:
        for line in inp:
            data = line.rstrip('\n').split(',')
            node_out = data[0]
            node_in = data[1]

            if node_out == node_in:
                continue

            idx_out = nodes.index(node_out)
            idx_in = nodes.index(node_in)

            edgeList.append((idx_out, idx_in))

    return edgeList


def ifTrueTreesAreFoundIn(minSolutions, args):
    truesnv, truecna = getTrueTrees(args)
    for solution in minSolutions:
        if set(truesnv) == set(solution.snv_edges) and set(truecna) == set(solution.cna_edges):
            return True
    return False

def getTrueTrees(args):
    df_fsnv, df_fcna = readProportionMatrices(args)
    snv_edges = getTreeEdges(args.truesnv, list(df_fsnv.index))
    cna_edges = getTreeEdges(args.truecna, list(df_fcna.index))
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
