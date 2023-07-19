import gurobipy as gp
import numpy as np
import pandas as pd
import networkx as nx
import itertools


# parsimonious clone reconciliation problem
class PCIsolver:

    def __init__(self, snv_mat, cna_mat, threads = 1, timelimit = None, verbose = True):
        self.snv_mat = snv_mat
        self.cna_mat = cna_mat
        self.threads = threads
        self.timelimit = timelimit
        self.verbose = verbose

        self.nsamples = self.snv_mat.shape[1]
        self.sol_clones = None
        self.sol_props = None

    def solve(self):

        nsamples = self.snv_mat.shape[1]
        assert nsamples == self.cna_mat.shape[1], 'SNV and CNA matrix sizes do not match up.'

        nsnv = self.snv_mat.shape[0]
        ncna = self.cna_mat.shape[0]

        model = gp.Model('PCIsolver')
        x = model.addVars(nsnv, ncna, vtype=gp.GRB.BINARY, name='x')
        w = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'w')
        y = model.addVars(nsamples, nsnv, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'y')
        d_snv = model.addVars(nsamples, nsnv, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_snv')
        d_cna = model.addVars(nsamples, ncna, vtype = gp.GRB.CONTINUOUS, lb = 0, ub = 1, name = 'delta_cna')

        xsum = gp.LinExpr()
        for i in range(nsnv):
            for j in range(ncna):
                xsum += x[i,j]
        model.addConstr(xsum == nsnv+ncna-1)

        # additional constraint -- must have at least one clone from each type
        snvSum = gp.LinExpr()
        for i in range(nsnv):
            for j in range(ncna):
                snvSum += x[i,j]
            model.addConstr(snvSum >= 1)
            snvSum.clear()
        cnaSum = gp.LinExpr()
        for j in range(ncna):
            for i in range(nsnv):
                cnaSum += x[i,j]
            model.addConstr(cnaSum >= 1)
            cnaSum.clear()

        # encode product w[i,j,k] = y[i,j,k] * x[j,k]
        for i in range(nsamples):
            for j in range(nsnv):
                for k in range(ncna):
                    model.addConstr(w[i,j,k] <= y[i,j,k])
                    model.addConstr(w[i,j,k] <= x[j,k])
                    model.addConstr(w[i,j,k] >= x[j,k] + y[i,j,k] - 1)

        # encode abundance constraint for snv with correction
        for i in range(nsamples):
            for j in range(nsnv):
                sum = gp.LinExpr()
                for k in range(ncna):
                    sum += w[i,j,k]
                model.addConstr(self.snv_mat[j,i] - sum <= d_snv[i,j])
                model.addConstr(sum - self.snv_mat[j,i] <= d_snv[i,j])

        # encode abundance constraint for cna with correction
        for i in range(nsamples):
            for k in range(ncna):
                sum = gp.LinExpr()
                for j in range(nsnv):
                    sum += w[i,j,k]
                model.addConstr(self.cna_mat[k,i] - sum <= d_cna[i,k])
                model.addConstr(sum - self.cna_mat[k,i] <= d_cna[i,k])

        # encode total abundance constraint
        for i in range(nsamples):
            sum = gp.LinExpr()
            for j in range(nsnv):
                for k in range(ncna):
                    sum += w[i,j,k]
            model.addConstr(sum == 1)

        # set objective function
        obj_sum = gp.LinExpr()
        for i in range(nsamples):
            for j in range(nsnv):
                obj_sum += d_snv[i,j]
            for k in range(ncna):
                obj_sum += d_cna[i,k]
        model.setObjective(obj_sum, gp.GRB.MINIMIZE)

        model.setParam(gp.GRB.Param.Threads, self.threads)
        model.optimize()

        if model.status == gp.GRB.OPTIMAL:
            solx = model.getAttr('x', x)
            self.sol_clones = [key for key, val in solx.items() if val >= 0.5]
            self.sol_props = model.getAttr('x', w)

    def writeCloneFile(self, clone_file, snv_clones = None, cna_clones = None):
        clone_data = []
        for clone in self.sol_clones:
            #for sample in range(self.nsamples):
            if snv_clones:
                snv_clone = snv_clones[clone[0]]
            else:
                snv_clone = clone[0]
            if cna_clones:
                cna_clone = cna_clones[clone[1]]
            else:
                cna_clone = clone[1]

            clone_data.append([clone, snv_clone, cna_clone] + [self.sol_props[sample, clone[0], clone[1]] for sample in range(self.nsamples)])
        df_clone = pd.DataFrame(clone_data, columns=['clone', 'snv_clone', 'cna_clone'] + [f'sample_{idx}' for idx in range(self.nsamples)])
        df_clone.to_csv(clone_file, sep='\t', index=False)












