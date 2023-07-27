import numpy as np
import random
import math
import pandas as pd
from test_suite.tests import validateTree, validateProps


def getClonesDf(args, clones):
    cloneProportions = getCloneProportions(args)
    validateProps(cloneProportions)
    dataClone = []
    for idx, clone in enumerate(list(clones)):
        dataClone.append([clone] + list(cloneProportions[idx, :]))
    dfClone = pd.DataFrame(dataClone, columns=['clone'] + [f'sample_{idx}' for idx in range(args.n)])
    return dfClone

def getCloneProportions(args):
    rows = sum(args.m) - len(args.m) + 1
    cols = args.n
    existence = np.random.binomial(1,args.p,(rows, cols))
    for c,s in enumerate(existence.sum(axis=0)):
        if s == 0:
            existence[np.random.randint(rows),c] = 1
    for r,s in enumerate(existence.sum(axis=1)):
        if s == 0:
            existence[r,np.random.randint(cols)] = 1

    clone_props = np.random.dirichlet([1]*rows, cols).transpose()
    clone_props[clone_props < 0.001] = 0.01 # at least 0.05, then add this constraint in linear program
    filtered_props = clone_props * existence
    norm_props = filtered_props / filtered_props.sum(axis=0)
    return norm_props