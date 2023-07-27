import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

def writeTree(G, ofile, suffix=''):
    nx.write_edgelist(G, f'{ofile}tree' + suffix +'.tsv', data=False, delimiter='\t')
    writeDOT(G, ofile, suffix)

def writeDOT(G, ofile, suffix=''):
    filename = f'{ofile}tree' + suffix + '.dot'
    with open(filename, 'w') as output:
        output.write(f'digraph N {{\n')
        output.write(f"\toverlap=\"false\"\n")
        output.write(f"\trankdir=\"TB\"\n")
        for u,v in G.edges:
            output.write("\t\"" + str(u) + "\" -> \"" + str(v) + "\" [style=\"bold\"];\n")
        output.write(f'}}')


def writeClones(df_clones, ofile,suffix=''):
    df_clones.to_csv(f'{ofile}props'+suffix+'.out', sep='\t', index=False)


def writeNoisySnvClones(df_clones, args, suffix=''):
    df_snv = df_clones.groupby('snv').sum(numeric_only=True).drop('cna', axis=1)  
    df_snv.index.names = ['genotypes']
    noisy_values = df_snv
    noisy_values.where(noisy_values != 0, (1 - args.t)*df_snv.values + args.t * np.random.dirichlet([1]*args.m, args.n).transpose())
    df_snv_noisy = pd.DataFrame(noisy_values, index=df_snv.index, columns=df_snv.columns)
    df_snv_noisy.to_csv(f'{args.o}_snv_noisy'+suffix+'.csv')

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
