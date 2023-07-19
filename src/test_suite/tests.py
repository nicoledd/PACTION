import networkx as nx
import numpy as np
import unittest

def validateFullTree(G):
    checkIfTree(G)
    checkRootNodeExists(G)
    checkRootNodeHasNoParents(G)
    checkRootNodeHasOneChild(G)
    checkTreeIsConnected(G)

def checkIfTree(G):
    assert G.number_of_nodes()-1==G.number_of_edges(), 'graph is not a tree'

def checkRootNodeExists(G):
    assert G.has_node((0,0)), 'root node (0,0) does not exist'

def checkRootNodeHasNoParents(G):
    assert len(G.in_edges((0,0)))==0, 'root node (0,0) has parent(s)'

def checkRootNodeHasOneChild(G):
    assert len(G[(0,0)])==1, 'root node (0,0) has more than 1 child'

def checkTreeIsConnected(G):
    visited = set() # Set to keep track of visited nodes of graph.
    def dfs(visited, graph, node):  #function for dfs 
        if node not in visited:
            visited.add(node)
            for neighbour in graph[node]:
                dfs(visited, graph, neighbour)
    dfs(visited, G, (0,0))
    assert len(visited)==len(G), 'graph is not connected'

def validateCloneProps(clone_props):
    checkColumnsSumToOne(clone_props)
    checkNoRowsAreAllZero(clone_props)

def checkColumnsSumToOne(clone_props):
    column_sums = clone_props.sum(axis=0)
    assert abs(sum(column_sums)-column_sums.shape[0])<0.1, 'columns (each representing a different sample) do not sum to 1'

def checkNoRowsAreAllZero(clone_props):
    row_sums = clone_props.sum(axis=1)
    assert abs(np.count_nonzero(row_sums)-row_sums.shape[0])<0.1, 'at least one row (each representing a clone type) is all zeroes'
