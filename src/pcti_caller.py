from linearPro.pcti_solver import PCTIsolver
from preprocess.clean_data import getTreeEdges, readProportionMatrices, getAllPossibleSnvCnaTrees

def solvePCTI(args):
    df_fsnv, df_fcna = readProportionMatrices(args)
    nsnv = len(df_fsnv)
    ncna = len(df_cna)
    snv_edges = getTreeEdges(args.snv_tree, list(df_fsnv.index))
    cna_edges = getTreeEdges(args.cna_tree, list(df_fcna.index))

    solver = PCTIsolver(df_fsnv.values, df_fcna.values, snv_edges, cna_edges)
    solver.solve()

    solver.writeCloneFile(f'{args.o}_clone_prediction.out', snv_clones = list(df_fsnv.index), cna_clones = list(df_fcna.index))
    solver.writeFullTree(f'{args.o}_clone_tree_prediction.tsv', snv_clones = list(df_fsnv.index), cna_clones = list(df_fcna.index))
    solver.writeDOT(f'{args.o}_clone_tree.dot', snv_clones = list(df_fsnv.index), cna_clones = list(df_fcna.index))
