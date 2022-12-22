# applying gurobipy to solve VRPTW

import gurobipy as gp
from gurobipy import GRB
from read_data import read_data
import numpy as np
import matplotlib.pyplot as plt
from time import time

def gurobi_VRPTW(problem):
    # read data
    vehicleNum = problem.vehicleNum
    capacity = problem.capacity
    # seperate data
    location = problem.location
    demand = problem.demand
    readyTime = problem.readyTime
    dueTime = problem.dueTime
    serviceTime = problem.serviceTime

    # building model
    MODEL = gp.Model('VRPTW')

    nodeNum = problem.nodeNum
    points = list(range(nodeNum))
    A = [(i, j) for i in points for j in points]
    D = {(i, j): np.linalg.norm(location[i]-location[j]) for i, j in A}

    ## add variates
    x = MODEL.addVars(A, vtype=GRB.BINARY)
    s = MODEL.addVars(points, vtype=GRB.CONTINUOUS)
    c = MODEL.addVars(points, vtype=GRB.CONTINUOUS)
    ## set objective
    MODEL.modelSense = GRB.MINIMIZE
    MODEL.setObjective(gp.quicksum(x[i, j] * D[i, j] for i, j in A))
    ## set constraints
    ### 1. flow balance
    MODEL.addConstrs(gp.quicksum(x[i, j] for j in points if j!=i)==1 for i in points[1:]) # depot not included
    MODEL.addConstrs(gp.quicksum(x[i, j] for i in points if i!=j)==1 for j in points[1:]) # depot not included
    ### 2. avoid subring / self-loop
    M = 1e7
    MODEL.addConstrs(s[i] + D[i, j] + serviceTime[i] - M * (1 - x[i, j]) <= s[j] for i, j in A if j!=0)
    ### 3. time constraints
    MODEL.addConstrs(s[i] >= readyTime[i] for i in points)
    MODEL.addConstrs(s[i] <= dueTime[i] for i in points)
    ### 4. capacity constraints
    MODEL.addConstrs(c[i] - demand[j] + M * (1 - x[i, j])>= c[j] for i, j in A if j!=0)
    MODEL.addConstrs(c[i] <= capacity for i in points)
    MODEL.addConstrs(c[i] >= 0 for i in points)
    ##serviceTime 5. vehicle number constraint
    MODEL.addConstr(gp.quicksum(x[0, j] for j in points) <= vehicleNum)

    # optimize the model
    MODEL.optimize()

    # get the routes
    routes = []
    for j in range(1, nodeNum):
        if round(x[0, j].X) == 1:
            route = [0]
            route.append(j)
            i = j
            while j != 0:
                for j in range(nodeNum):
                    if round(x[i, j].X) == 1:
                        route.append(j)
                        i = j
                        break
            routes.append(route)
 
    return routes, MODEL.ObjVal

def show_routes(location, routes):
    plt.figure()
    plt.scatter(location[1:, 0], location[1:, 1])
    plt.scatter(location[0:1, 0], location[0:1, 1], s = 150, c = 'r', marker='*')
    for route in routes:
        plt.plot(location[route, 0], location[route, 1], c='r')
    plt.show()

if __name__ == "__main__":
    file_name = "solomon_100/C101.txt"
    # 101\102-1s, 101_25-0.02s, 103-200s
    problem = read_data(file_name)
    time1 = time()
    routes, obj = gurobi_VRPTW(problem)
    time2 = time()
    show_routes(problem.location, routes)
    print("optimal obj: {}\ntime consumption: {}".format(obj, time2-time1))

