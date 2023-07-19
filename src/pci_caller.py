import csv
from csv import writer

from linearPro.pci_solver import PCIsolver
from linearPro.pcti_solver import PCTIsolver
from preprocess.clean_data import getIndexOfTrueGraphs, processSolutionTree, processTreeFile, ifTrueTreesAreFoundIn, getTreeEdges, readProportionMatrices, getAllPossibleSnvCnaTrees
from postprocess.write import writeText, writeFullTree, writeFullTreeDOT, writeSnvTree, writeSnvTreeDOT, writeCnaTree, writeCnaTreeDOT, writeFullClones, writeSnvClones, writeCnaClones


def solvePCI(args):

    df_fsnv, df_fcna = readProportionMatrices(args)
    snvTrees, cnaTrees = getAllPossibleSnvCnaTrees(args)
    
    minCorrection = 100
    minSolutions = []

    for snv_edges in snvTrees:
        for cna_edges in cnaTrees:

            solver = PCTIsolver(df_fsnv.values, df_fcna.values, snv_edges, cna_edges)
            solver.solve()

            if solver.correction == None:
                continue
            else:
                if solver.correction < minCorrection:
                    minSolutions = [solver]
                    minCorrection = solver.correction
                elif solver.correction == minCorrection:
                    minSolutions.append(solver)
    writeSolutions(minSolutions, args)


def writeSolutions(minSolutions, args):
    writeText(str(len(minSolutions)), args.o+'num_solutions')
    if minSolutions:
        trueTreesFound = ifTrueTreesAreFoundIn(minSolutions, args)
        if trueTreesFound:
            trueTreesIdx = getIndexOfTrueGraphs(minSolutions, args)
            writeTreesInfo(minSolutions[trueTreesIdx],args,'_true')
        maxEdgeRecallTreesIdx = getIndexOfMaxEdgeRecallTrees(minSolutions,args)
        writeTreesInfo(minSolutions[maxEdgeRecallTreesIdx],args,'_maxEdgeRecall')
    writeAllEdgeRecalls(minSolutions,args,trueTreesIdx)



def writeTreesInfo(solution,args,suffix):
    # 1. write correction
    writeText(str(solution.correction), args.o+'correction' + suffix)
    # 2. write edge recall
    trueTree = processTreeFile(args.truetree)
    predictedTree = processSolutionTree(solution.edges)
    edgeRecall = getEdgeRecall(trueTree, predictedTree)
    writeText(str(edgeRecall), args.o+'edge_recall' + suffix)
    # 3. write all trees
    writeFullTreeDOT(solution.G, args.o, suffix)
    writeCnaTreeDOT(solution.C, args.o, suffix)
    writeSnvTreeDOT(solution.S, args.o, suffix)
    writeFullTree(solution.G, args.o, suffix)
    writeSnvTree(solution.S, args.o, suffix)
    writeCnaTree(solution.C, args.o, suffix)
    # 4. write proportions
    writeFullClones(solution.propDf, args.o,suffix)
    writeSnvClones(solution.propDf, args.o,suffix)
    writeCnaClones(solution.propDf, args.o,suffix)



def writeAllEdgeRecalls(minSolutions, args, trueTreesIdx):
    with open(args.o + 'edgeRecalls.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["EdgeRecall", "TrueTreesOrNot"])
    
    for idx in range(len(minSolutions)):
        solution = minSolutions[idx]
        trueTree = processTreeFile(args.truetree)
        predictedTree = processSolutionTree(solution.edges)
        edgeRecall = getEdgeRecall(trueTree, predictedTree)
        with open(args.o + 'edgeRecalls.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            if idx == trueTreesIdx:
                writer.writerow([str(edgeRecall), "T"])
            else:
                writer.writerow([str(edgeRecall), "F"])


def getEdgeRecall(trueTree,predictedTree):
    intersection = 0
    for ele in predictedTree:
        if ele in trueTree:
            intersection += 1
    edgeRecall = intersection/len(trueTree)
    return edgeRecall


def getIndexOfMaxEdgeRecallTrees(minSolutions,args):
    maxEdgeRecall = 0
    maxEdgeRecallIdx = -1
    idx = 0
    for solution in minSolutions:
        trueTree = processTreeFile(args.truetree)
        predictedTree = processSolutionTree(solution.edges)
        edgeRecall = getEdgeRecall(trueTree, predictedTree)
        if edgeRecall > maxEdgeRecall:
            maxEdgeRecallIdx = idx
            maxEdgeRecall = edgeRecall
        idx += 1
    return maxEdgeRecallIdx


def getMaximumRecall(recallDf):
    return max(recallDf['recall'])

def getTrueRecall(trueTreesIdx, recallDf):
    return recallDf['recall'].iloc[trueTreesIdx]
