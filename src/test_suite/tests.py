import networkx as nx
import numpy as np
import unittest

def validateTree(G):
    checkIfTree(G)
    checkRootNodeExists(G)
    checkRootNodeHasNoParents(G)
    checkRootNodeHasOneChild(G)
    checkTreeIsConnected(G)

def checkIfTree(G):
    assert G.number_of_nodes()-1==G.number_of_edges(), 'graph is not a tree'

def checkRootNodeExists(G):
    if type(list(G.nodes)[0]) == tuple:
        numModes = len(list(G.nodes)[0])
        assert G.has_node((0,)*numModes), 'root node ' + str((0,)*numModes) + ' does not exist'
    else:
        assert G.has_node(0), 'root node 0 does not exist'


def checkRootNodeHasNoParents(G):
    if type(list(G.nodes)[0]) == tuple:
        numModes = len(list(G.nodes)[0])
        assert len(G.in_edges((0,)*numModes))==0, 'root node ' + str((0,)*numModes) + ' has parent(s)'
    else:
        assert len(G.in_edges(0))==0, 'root node 0 has parent(s)'


def checkRootNodeHasOneChild(G):
    if type(list(G.nodes)[0]) == tuple:
        numModes = len(list(G.nodes)[0])
        assert len(G[(0,)*numModes])==1, 'root node ' + str((0,)*numModes) + ' has more than 1 child'


def checkTreeIsConnected(G):
    visited = set() # Set to keep track of visited nodes of graph.
    def dfs(visited, graph, node):  #function for dfs 
        if node not in visited:
            visited.add(node)
            for neighbour in graph[node]:
                dfs(visited, graph, neighbour)
    if type(list(G.nodes)[0]) == tuple:
        numModes = len(list(G.nodes)[0])
        dfs(visited, G, (0,)*numModes)
    else:
        dfs(visited, G, 0)
    assert len(visited)==len(G), 'graph is not connected'

def validateProps(clone_props):
    checkColumnsSumToOne(clone_props)
    checkNoRowsAreAllZero(clone_props)

def checkColumnsSumToOne(clone_props):
    column_sums = clone_props.sum(axis=0)
    assert abs(sum(column_sums)-column_sums.shape[0])<0.1, 'columns (each representing a different sample) do not sum to 1'

def checkNoRowsAreAllZero(clone_props):
    row_sums = clone_props.sum(axis=1)
    assert abs(np.count_nonzero(row_sums)-row_sums.shape[0])<0.1, 'at least one row (each representing a clone type) is all zeroes'
