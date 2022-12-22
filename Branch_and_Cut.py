"""
Branch and Cut algorithm, implemented by gurobi recall
author: Charles Lee
date: 2022.11.07
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import gurobipy as gp
from gurobipy import GRB
from read_data import read_data, read_data_from_Augerat

def callback_BC(model, where):
    if where == GRB.Callback.MIPSOL:
        vars = model._vars
        vals = model.cbGetSolution(vars)
        graph = model._graph
        Cuts = seperate(vars, vals, graph)
        for cut in Cuts:
            model.cbLazy(cut >= 0)
    if where == GRB.Callback.MIPNODE and model.cbGet(GRB.Callback.MIPNODE_STATUS) == GRB.OPTIMAL:
        vars = model._vars
        vals = model.cbGetNodeRel(vars)
        graph = model._graph
        Cuts = seperate(vars, vals, graph)
        for cut in Cuts:
            model.cbLazy(cut >= 0)
        
def seperate(vars, vals, graph):
    Cuts = []
    nCus = graph.nodeNum
    superVertex = []
    cusLabels = [1]*nCus
    while True:
        # find the open edge with the highest weight
        iMax = -1
        jMax = -1
        wMax = 0.001
        for i in range(1, nCus):
            for j in range(1, nCus):
                if cusLabels[i] == 0 or cusLabels[j] == 0 or i==j:
                    continue
                idx = getIndex(i, j, nCus)
                weight = vals[idx]
                if weight > wMax:
                    wMax = weight
                    iMax = i
                    jMax = j
        if iMax == -1 and jMax == -1:
            break

        superVertex.clear()
        superVertex.append(iMax)
        superVertex.append(jMax)
        cusLabels[iMax] = 0
        cusLabels[jMax] = 0

        # while there is any open customer
        while True:
            cus = identifyMaxWeightCustomer(superVertex, cusLabels, vals, graph)
            if cus != -1:
                cusLabels[cus] = 0
                superVertex.append(cus)
            else:
                break
        
        # check violation and add cut
        lhs1 = lhs2 = rhs = 0
        lhsExpr1 = gp.LinExpr()
        lhsExpr2 = gp.LinExpr()
        for i in range(len(superVertex)):
            for j in range(nCus):
                if j in superVertex:
                    continue
                idx1 = getIndex(j, superVertex[i], nCus)
                lhs1 += vals[idx1]
                lhsExpr1 += vars[idx1]
                idx2 = getIndex(superVertex[i], j, nCus)
                lhs2 += vals[idx2]
                lhsExpr2 += vars[idx2]
            rhs += graph.demand[superVertex[i]]
        rhs = math.ceil(rhs / graph.capacity)
        if lhs1 < rhs:
            Cuts.append(lhsExpr1 - rhs)
        if lhs2 < rhs:
            Cuts.append(lhsExpr2 - rhs)
    return Cuts
        
def getIndex(i, j, nCus):
    assert i>=0 and i<nCus and j>=0 and j<nCus, "illegal node index"
    return i*nCus + j

def identifyMaxWeightCustomer(superVertex, cusLabels, var_vals, graph):
    nCus = graph.nodeNum
    wMax = 0.001
    customer = -1
    for j in range(1, nCus):
        if j in superVertex or cusLabels[j] == 0:
            continue
        weight = 0
        for i in range(len(superVertex)):
            idx1 = getIndex(superVertex[i], j, nCus)
            weight += var_vals[idx1]
            idx2 = getIndex(j, superVertex[i], nCus)
            weight += var_vals[idx2]
        if weight > wMax:
            wMax = weight
            customer = j
    return customer

class BranchAndCut():
    def __init__(self, graph):
        self.graph = graph
        self.model = self.build_seperate_VRP_model(graph)

    def build_seperate_VRP_model(self, graph):
        # building model
        MODEL = gp.Model('VRPTW')

        # data preprocess
        points = list(range(graph.nodeNum))
        A = [(i, j) for i in points for j in points]
        D = graph.disMatrix

        ## add variates
        x = MODEL.addVars(A, vtype=GRB.BINARY, name='x')
        ## set objective
        MODEL.modelSense = GRB.MINIMIZE
        MODEL.setObjective(gp.quicksum(x[i, j] * D[i, j] for i, j in A))
        ## set constraints(flow balance)
        MODEL.addConstrs(gp.quicksum(x[i, j] for j in points if j!=i)==1 for i in points[1:]) # depot not included
        MODEL.addConstrs(gp.quicksum(x[i, j] for i in points if i!=j)==1 for j in points[1:]) # depot not included

        # set model params
        MODEL.setParam('OutputFlag', 1)
        MODEL.setParam('LazyConstraints', 1)

        # update model
        MODEL.update()

        return MODEL

    def run(self):
        # pass data into callback function
        self.model._vars = self.model.getVars()
        self.model._graph = self.graph
        self.model.optimize(callback_BC)

if __name__ == "__main__":
    file_name = "solomon_100\C101.txt"
    graph = read_data(file_name)
    # file_name = "A-n32-k5.txt"
    # graph = read_data_from_Augerat(file_name)
    alg = BranchAndCut(graph)
    alg.run()