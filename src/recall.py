def getRecall(filename, samples):
    for noise in ['0','0.05','0.1','0.15']:
        outfile = open("/scratch/data/nsdong2/projectPACTION/results/" + filename + "paction/" + samples + "sample/" + noise + "_clone_prediction.out", "r")
        outclones = set()
        for line in outfile:
            splitted = line.split()
            if splitted[0][0] == '(':
                outclones.add(splitted[0]+splitted[1])
        trueclones = set()
        truefile = open("/scratch/data/nsdong2/paction/data/extensions/samples" + samples + "noise" + noise + "_clone_tree.tsv", "r")
        for line in truefile:
            splitted = line.split()
            for i in range(len(splitted)):
                if splitted[i][0] == '(':
                    trueclones.add(splitted[i]+splitted[i+1])

        recall = len(outclones.intersection(trueclones))/len(trueclones)
        return recall

def main():
    samples = '1'
    oldrecall = getRecall('', samples)
    newrecall = getRecall('old', samples)
    print(oldrecall, newrecall)

main()