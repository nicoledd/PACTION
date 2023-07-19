import pandas as pd
import matplotlib.pyplot as plt


def writeRecallGraph(minSolutions, args, idxTrueGraphs):
    f1 = pd.read_csv('/scratch/data/nsdong2/projectPACTION/simulations/result_'+str(args.seed)+'_'+str(args.samples)+'_'+str(args.noise)+'/_clone.out', sep='\t')
    clones=f1['clone']
    recalls=[]
    for i in range(len(minSolutions)):
        pred_clones = minSolutions[i].clones
        recalls.append((len(set(clones) & set([str(x) for x in pred_clones]))/len(clones)))
    plotGraph(recalls, args, idxTrueGraphs)


def plotGraph(recalls, args, idxTrueGraphs):
    fig, ax = plt.subplots()
    ax.scatter([i for i in range(len(recalls))], recalls, c='b')
    ax.scatter([idxTrueGraphs], [recalls[idxTrueGraphs]], c='r')
    ax.set_xlabel('solution')
    ax.set_ylabel('recall')
    ax.set_title('seed='+str(args.seed)+', samples='+str(args.samples)+', noise='+str(args.noise)+', compare bruteforce pci sols')
    ax.figure.savefig('/scratch/data/nsdong2/projectPACTION/figures/'+str(args.seed)+'_'+str(args.samples)+'_'+str(args.noise)+'_clone_prediction.png')