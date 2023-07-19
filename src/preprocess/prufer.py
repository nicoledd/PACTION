import numpy as np
import itertools
import math

# n: number of nodes in tree
# return: a list of all potential trees
def enumPrufer(n):

    trees = []

    # 1) get a list of all potential prufer sequences
    pruferSequences = []
    for seq in itertools.product([i for i in range(1,n+1)], repeat=n-2):
        pruferSequences.append(list(seq))

    # 2) for each prufer sequence, compute correct tree with seqToTree()
    # 2.1) use preprocessTree() to turn prufer tree into a valid tree
    for seq in pruferSequences:
        pruferTree = seqToTree(seq,len(seq))
        preprocessedTree = preprocessTree(pruferTree)
        trees.append(preprocessedTree)
    return trees




# given a tree, preprocess it to ensure it turns into a valid phylogeny
def preprocessTree(pruferTree):
    stack = [0]
    m = len(pruferTree)
    visited = [False for _ in range(m)]
    while stack:
        currNode = stack.pop() # current parent
        for i in range(m):
            if visited[i] == False:
                u,v = pruferTree[i]
                if u == currNode:
                    stack.append(v)
                    visited[i] = True
                elif v == currNode:
                    stack.append(u)
                    pruferTree[i] = (v,u)
                    visited[i] = True
    return pruferTree




# given a prufer sequence, return the corresponding tree in the format of a list of edges
def seqToTree(prufer, m):
    
    edges = []
    vertices = m + 2
    vertex_set = [0] * vertices
    for i in range(vertices - 2):
        vertex_set[prufer[i] - 1] += 1
    j = 0
    for i in range(vertices - 2):
        for j in range(vertices):
            if (vertex_set[j] == 0):
                vertex_set[j] = -1
                edges.append((j+1-1,prufer[i]-1))
                vertex_set[prufer[i] - 1] -= 1
                break
    j = 0
    tmpa = 0
    tmpb = 0
    for i in range(vertices):
        if (vertex_set[i] == 0 and j == 0):
            tmpa = i+1
            j += 1
        elif (vertex_set[i] == 0 and j == 1):
            tmpb = i+1
    edges.append((tmpa-1,tmpb-1))
    return edges
 