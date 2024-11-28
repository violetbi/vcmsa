import os, sys

import numpy as np
import gurobipy as gp
from gurobipy import GRB


def ilp_solver(aa_emb: np.ndarray, intvs: np.ndarray, verbose=0, timeout=600):
    m = gp.Model("min weighted k-clique in k-partite graph")
    # set model parameters
    m.Params.OutputFlag = verbose
    m.Params.TIME_LIMIT = timeout # seconds
    m.Params.LogToConsole = 1
    m.Params.Threads = 8
    m.Params.MIPFocus = 3
    # m.Params.NumericFocus = 0
    # m.Params.PoolSearchMode = 0
    # m.Params.PoolSolutions = 10
    # m.Params.Method = 4

    num_seq = len(intvs)
    intv_sz = list(int(i) for i in intvs[:, 1] - intvs[:, 0])

    aa_vars = []
    for i in range(num_seq):
        aa_var = m.addVars(intv_sz[i], vtype=GRB.BINARY)
        aa_vars.append(aa_var)
    if verbose:
        print(f"Number of aa variable={np.sum(intv_sz)}")
    obj = gp.LinExpr()
    z_vars = []
    for i in range(len(aa_emb)): # number of sequence
        m.addConstr(gp.quicksum(aa_vars[i]) == 1)
        for j in range(len(aa_emb)):
            if i == j:
                continue
            z_var = m.addVars(intv_sz[i], intv_sz[j], vtype=GRB.CONTINUOUS)
            for k1 in range(intvs[i, 1] - intvs[i, 0]):
                for k2 in range(intvs[j, 1] - intvs[j, 0]):
                    m.addConstr(z_var[k1, k2] <= aa_vars[i][k1])
                    m.addConstr(z_var[k1, k2] <= aa_vars[j][k2])
                    m.addConstr(z_var[k1, k2] >= aa_vars[i][k1] + aa_vars[j][k2] - 1)

                    d = np.linalg.norm(aa_emb[i][k1] - aa_emb[j][k2])
                    obj += z_var[k1, k2] * d
            z_vars.append(z_var)
    if verbose:
        print(f"Number of z slack variable={sum(len(z_var) for z_var in z_vars)}")
    
    obj = obj / num_seq # average pairwise distance

    m.setObjective(obj, GRB.MINIMIZE)
    m.update()
    try:
        m.optimize()
    except gp.GurobiError():
        print("Optimize failed due to non-convexity")
        sys.exit(1)
    if verbose:
        print("ILP finished under time constraint: ", timeout)
    if m.SolCount == 0:
        print("No solution found under current time constraint")
        print("Either input is not feasible or increase the time limit")
        sys.exit(1)
    if verbose:
        print("Suboptimal/Optimal solution found")
    midpoints = []
    for i in range(num_seq):
        # print(f"sequence {i}")
        vars = [aa_vars[i][j].x for j in range(intv_sz[i])]
        midpoints.append(vars.index(1) + intvs[i, 0])
    print("Objective: ", obj.getValue())
    return midpoints


def divide_conquer_midpoints(aa_emb: np.ndarray, seq_sz: np.ndarray, verbose: int):
    intvs = np.array([[0, seq_sz[i]] for i in range(len(seq_sz))])
    midpoints_mm = acc_solve(aa_emb, intvs, verbose)
    return midpoints_mm

def acc_solve(aa_emb: np.ndarray, intvs: np.ndarray, layer=0, label="X", verbose=0):
    if verbose:
        print()
        print(f"Layer={layer}, label={label}, interval size={intvs[:, 1] - intvs[:, 0]}")
    midpoints = []
    if any(intvs[:, 1] - intvs[:, 0] == 0):
        return midpoints
    mid = ilp_solver(aa_emb, intvs, verbose=verbose)
    intv_left = np.array([[intvs[i, 0], mid[i]] for i in range(len(intvs))])
    intv_right = np.array([[mid[i] + 1, intvs[i, 1]] for i in range(len(intvs))])
    mid1 = acc_solve(aa_emb, intv_left, layer=layer + 1, label='L', verbose=verbose)
    mid2 = acc_solve(aa_emb, intv_right, layer=layer + 1, label='R', verbose=verbose)
    midpoints.extend(mid1)
    midpoints.append(mid)
    midpoints.extend(mid2)
    return midpoints