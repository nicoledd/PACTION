import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

def writeFullTree(G, ofile,suffix=''):
    nx.write_edgelist(G, f'{ofile}_clone_tree' + suffix+'.tsv', data=False, delimiter='\t')

def writeSnvTree(G, ofile,suffix=''):
    nx.write_edgelist(G, f'{ofile}_snv_tree'+suffix+'.csv', data=False, delimiter=',')

def writeCnaTree(G, ofile,suffix=''):
    nx.write_edgelist(G, f'{ofile}_cna_tree'+suffix+'.csv', data=False, delimiter=',')

def writeFullTreeDOT(G, ofile, suffix=''):
    filename = f'{ofile}_clone_tree' + suffix + '.dot'
    with open(filename, 'w') as output:
        output.write(f'digraph N {{\n')
        output.write(f"\toverlap=\"false\"\n")
        output.write(f"\trankdir=\"TB\"\n")
        for edge in G.edges:
            u,v = edge
            output.write(f"\t\"{u[0]}, {u[1]}\" -> \"{v[0]}, {v[1]}\" [style=\"bold\"];\n")
        output.write(f'}}')

def writeCnaTreeDOT(G, ofile, suffix=''):
    filename = f'{ofile}_cna_tree' + suffix + '.dot'
    with open(filename, 'w') as output:
        output.write(f'digraph N {{\n')
        output.write(f"\toverlap=\"false\"\n")
        output.write(f"\trankdir=\"TB\"\n")
        for edge in G.edges:
            u,v = edge
            output.write(f"\t\"{u}\" -> \"{v}\" [style=\"bold\"];\n")
        output.write(f'}}')

def writeSnvTreeDOT(G, ofile,suffix=''):
    filename = f'{ofile}_snv_tree'+suffix+'.dot'
    with open(filename, 'w') as output:
        output.write(f'digraph N {{\n')
        output.write(f"\toverlap=\"false\"\n")
        output.write(f"\trankdir=\"TB\"\n")
        for edge in G.edges:
            u,v = edge
            output.write(f"\t\"{u}\" -> \"{v}\" [style=\"bold\"];\n")
        output.write(f'}}')

def writeFullClones(df_clones, ofile,suffix=''):
    df_clones.to_csv(f'{ofile}_clone'+suffix+'.out', sep='\t', index=False)



def writeSnvClones(df_clones, ofile,suffix=''):
    df_snv = df_clones.groupby('snv').sum(numeric_only=True).drop('cna', axis=1)
    df_snv.index.names = ['genotypes']
    df_snv.to_csv(f'{ofile}_snv'+suffix+'.csv')


def writeNoisySnvClones(df_clones, args, suffix=''):
    df_snv = df_clones.groupby('snv').sum(numeric_only=True).drop('cna', axis=1)  
    df_snv.index.names = ['genotypes']
    noisy_values = df_snv
    noisy_values.where(noisy_values != 0, (1 - args.t)*df_snv.values + args.t * np.random.dirichlet([1]*args.m, args.n).transpose())
    df_snv_noisy = pd.DataFrame(noisy_values, index=df_snv.index, columns=df_snv.columns)
    df_snv_noisy.to_csv(f'{args.o}_snv_noisy'+suffix+'.csv')

def writeCnaClones(df_clones, ofile, suffix=''):
    df_cna = df_clones.groupby('cna').sum(numeric_only=True).drop('snv', axis=1)
    df_cna.index.names = ['genotypes']
    df_cna.to_csv(f'{ofile}_cna'+suffix+'.csv')

def writeNoisyCnaClones(df_clones,args, suffix=''):
    df_cna = df_clones.groupby('cna').sum(numeric_only=True).drop('snv', axis=1)
    df_cna.index.names = ['genotypes']
    noisy_values = df_cna
    noisy_values.where(noisy_values != 0, (1 - args.t)*df_cna.values + args.t * np.random.dirichlet([1]*args.d, args.n).transpose())
    df_cna_noisy = pd.DataFrame(noisy_values, index=df_cna.index, columns=df_cna.columns)
    df_cna_noisy.to_csv(f'{args.o}_cna_noisy'+suffix+'.csv')


def writeText(text, ofile):
    filename = f'{ofile}.txt'
    with open(filename, 'w') as output:
        output.write(text)
