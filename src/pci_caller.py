import csv
from csv import writer
import pandas as pd
import networkx as nx

from linear_programs.pcti_solver import PCTIsolver
from preprocess.clean_data import getIndexOfTrueGraphs, processSolutionTree, processTreeFile, ifTrueTreesAreFoundIn, getTreeEdges, readProportionMatrices, getAllPossibleTrees
from postprocess.write import writeText, writeTree, writeClones


def getPropBList(args, i):
    return [readProportionMatrices(args.p[i])]

def getTreeBList(args, i):
    propBList = getPropBList(args, i)
    if args.t[i] == "None":
        return getAllPossibleTrees(propBList[0])
    else:
        return args.t[i]

def getNumSamples(args):
    propAList = initPropAList(args)
    return propAList[0].shape[1]-1

def solveProgressivePCI(args):
    propAList = initPropAList(args)
    treeAList = initTreeAList(args)
    clonesAList = initClonesAList(args)
    nsamples = getNumSamples(args)
    d = {} # clone tuple (0,0,0) -> int 0
    e = {} # int 0 -> clone tuple (0,0,0)

    numOptSols=0
    for i in range(1,len(args.p)):
        propBList = getPropBList(args, i)
        treeBList = getTreeBList(args, i)

        propAListNew, treeAListNew, clonesAListNew = [], [], []

        bestSolution = None
        minCorrection = 100
        for propA, treeA, clonesA in zip(propAList, treeAList, clonesAList):
            for treeB in treeBList:
                propB = propBList[0]
                d,e = getProcessedClones(clonesA, d, e)
                propAFake = getProcessedProps(propA, d)
                treeAFake = getProcessedTree(treeA, d)

                solution = solvePCI(propAFake, propB, treeAFake, treeB, args)
                assert solution!=None, 'no solution, progressive paction failed'

                if solution.correction != None and solution.correction < minCorrection:
                    minCorrection = solution.correction

                    bestSolution = solution
                    propSolFake = solution.proportions
                    treeSolFake = solution.edges
                    clonesSolFake = solution.clones

                    propANew = translateProps(clonesSolFake, propSolFake, e, nsamples)
                    treeANew = translateTrees(treeSolFake, e)
                    clonesANew = translateClones(clonesSolFake, e)

                    propAListNew = [propANew]
                    treeAListNew = [treeANew]
                    clonesAListNew = [clonesANew]
                if i==1 and solution.correction != None and solution.correction == minCorrection:
                    numOptSols += 1
        propAList, treeAList, clonesAList = propAListNew, treeAListNew, clonesAListNew

    # L inf norm as obj val
    # constrain copy-number states

    writeTree(treeAList[0], args.o, '')
    writeClones(propAList[0], args.o, '')

def initPropAList(args):
    return [readProportionMatrices(args.p[0])]

def initClonesAList(args):
    propAList = initPropAList(args)
    return [list(propAList[0]['genotypes'])]

def initTreeAList(args):
    treesA = args.t[0]
    propAList = initPropAList(args)
    if args.t[0] == "None":
        return getAllPossibleTrees(propAList[0])
    else:
        return [getTreeEdges(treesA, list(propAList[0]['genotypes']))]



def solvePCI(propA, propB, treeA, treeB, args):
    #minCorrection = 100
    #minSolutions = []

    propA = propA.drop(columns=['genotypes'])
    propB = propB.drop(columns=['genotypes'])

    #print(propA.values)
    #print(propB.values)
    #print(list(treeA.edges), list(treeB.edges))
    solver = PCTIsolver(propA.values, propB.values, list(treeA.edges), list(treeB.edges))
    solver.solve()
    
    return solver

    '''
    if solver.correction == None:
        continue
    else:
        if solver.correction < minCorrection:
            minSolutions = [solver]
            minCorrection = solver.correction
        elif solver.correction == minCorrection:
            minSolutions.append(solver)

    return minSolutions'''



def writeSolutions(minSolutions, args):
    writeText(str(len(minSolutions)), args.o+'num_solutions')
    if minSolutions:
        trueTreesFound=False
        #trueTreesFound = ifTrueTreesAreFoundIn(minSolutions, args)
        trueTreesIdx = -1
        if trueTreesFound:
            trueTreesIdx = getIndexOfTrueGraphs(minSolutions, args)
            writeTreesInfo(minSolutions[trueTreesIdx],args,'_true')
        maxEdgeRecallTreesIdx = getIndexOfMaxEdgeRecallTrees(minSolutions,args)
        writeTreesInfo(minSolutions[maxEdgeRecallTreesIdx],args,'_maxEdgeRecall')
    writeAllEdgeRecalls(minSolutions,args,trueTreesIdx)



def writeTreesInfo(solution,args,suffix, e):

    # 1. write correction
    writeText(str(solution.correction), args.o+'correction' + suffix)
    # 2. write edge recall
    #trueTree = processTreeFile(args.truetree)
    #predictedTree = processSolutionTree(solution.edges)
    #edgeRecall = getEdgeRecall(trueTree, predictedTree)
    #writeText(str(edgeRecall), args.o+'edge_recall' + suffix)

    writeTree(solution.G, args.o, suffix)
    # 4. write proportions
    writeClones(solution.propDf, args.o,suffix)



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


def translateProps(clonesCombined, propsCombined, e, nsamples):
    data = []
    for a,b in clonesCombined:
        # a is to be expanded, b is to stay the same
        if isinstance(e[a], tuple):
            clone=e[a]+(b,)
        else:
            clone = (e[a],b) # how to combine two tuples into one big tuple?
        data.append([clone] + [propsCombined[sample,a,b] for sample in range(nsamples)])
    translatedProps = pd.DataFrame(data, columns=['genotypes'] + [f'sample_{idx}' for idx in range(nsamples)])
    return translatedProps


def translateTrees(treesCombined, e):
    translatedTree = []
    for uv,wx in treesCombined:
        u,v = uv
        w,x = wx
        if isinstance(e[u], tuple):
            edge1 = e[u] + (v,)
        else:
            edge1 = (e[u],v)
        if isinstance(e[w], tuple):
            edge2 = e[w] + (x,)
        else:
            edge2 = (e[w],(x))
        translatedTree.append((edge1, edge2))
    G = nx.DiGraph()
    for u,v in translatedTree:
        G.add_edge(u,v)
    return G

def translateClones(clonesCombined, e):
    translatedClones = []
    for a,b in clonesCombined:
        if isinstance(e[a], tuple):
            clone = e[a] + (b,)
        else:
            clone = (e[a],b)
        translatedClones.append(clone)
    return translatedClones


def getProcessedClones(trueClones, d, e):
    d = {}
    e = {}
    for i in range(len(trueClones)):
        d[trueClones[i]] = i
        e[i] = trueClones[i]
    return d,e

def getProcessedTree(trueTree, d):
    G = nx.DiGraph()
    for u,v in trueTree.edges:
        G.add_edge(d[u], d[v])
    return G

def getProcessedProps(trueProps, d):
    data = []
    nsamples = trueProps.shape[1]
    trueProps['genotypes'].map(d)
    return trueProps
