import sys
import argparse

from pci_caller import solvePCI
from pcti_caller import solvePCTI


def main(args):
    if args.t.count("None") == 0:
        solvePCTI(args)
    else:
        solvePCI(args)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    # debugging arguments - to specify the exact test case to use
    parser.add_argument('--noise', type=str)
    parser.add_argument('--samples', type=str)
    parser.add_argument('--seed', type=str)
    parser.add_argument('--truesnv', type=str)
    parser.add_argument('--truecna', type=str)
    parser.add_argument('--truetree', type=str)

    propString = 'csv file with proportion matrix for each genotype'
    treeString = 'csv file with edges of trees (if N/A use \'None\')'
    outString = 'output prefix'
    parser.add_argument('-p', type=str, help=propString, required=True, nargs='*')
    parser.add_argument('-t', type=str, help=treeString, required=True, nargs='*')
    parser.add_argument('-o', type=str, help=outString, required=True)

    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    main(args)


'''
# may be useful later?

def recordEdgeRecallsInCsv(recallDf, trueTreesIdx, trueTreesFound,args):
    maxRecall = getMaximumRecall(recallDf)
    if trueTreesFound:
        trueRecall = getTrueRecall(trueTreesIdx, recallDf)
    else:
        trueRecall = 0
    rows = [float(args.noise), int(args.samples), int(args.seed), maxRecall, trueRecall]

    filename = args.o + 'edge_recalls.csv'
    with open(filename, 'a') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(rows)

def writeTreeRecall(minSolutions, args):
    trueTree = processTreeFile(args.truetree)
    treeRecalls = []
    for solution in minSolutions:
        predictedTree = processSolutionTree(solution.edges)
        intersection = 0
        for ele in predictedTree:
            if ele in trueTree:
                intersection += 1
        treeRecalls.append(intersection/len(trueTree))
    df = pd.DataFrame(np.array(treeRecalls), columns = ['recall'])
    return df
'''
