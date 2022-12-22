"""
Restruct the code: Branch and Bound (VRPTW for example)
Author: Mingzhe Li
Date: 2022.09.10
"""

import numpy as np
import matplotlib.pyplot as plt
import gurobipy as gp
from gurobipy import GRB
import math
import copy
import time
from read_data import read_data

class Node():
    """
    Node class 
    """
    def __init__(self, model, local_LB=np.inf):
        """properties before being solved"""
        self.model = model.copy() # apply gurobi function: copy()
        self.local_LB = local_LB
        self.EPS = 1e-6

        """properties after being solved"""
        self.IP_obj = np.inf
        self.LP_obj = np.inf
        self.x_sol = {} # contain solution. key: variable name; value: variable value
        self.x_int_sol = {} # contain integer version of solution. 
        self.is_integer = False # whether get an integer solution or not 
        self.is_feasible = False # whether LP has solution
        self.branch_var_list = [] # contain non-integer variables for branch

        """display properties"""
        self.way_of_opt = '---' # the way of getting optimal IP, By Simple/Rounding
        self.depth = 0 # the search depth of node
        self.prune_info = '---' # if been pruned, why been pruned, By Inf/Bnd/Opt
        self.has_showed_way_of_opt = 0 # whether showed the way of opt

    def solve_and_update(self):
        """
        solve the model and update properties
        """

        """ Solve and Check feasibility """
        self.model.optimize()
        self.is_feasible = self.model.Status == 2
        if not self.is_feasible:
            return

        """ Update x_sol and Round to x_int_sol """
        p_num2 = 0 # count p_num^2
        for var in self.model.getVars():
            self.x_sol[var.VarName] = var.x
            if var.VarName[0] == 'x': # round for integer variables 
                self.x_int_sol[var.VarName] = round(var.x)
                p_num2 += 1
            else:
                self.x_int_sol[var.VarName] = var.x
        p_num = int(math.sqrt(p_num2))

        """ Check integer solution"""
        self.is_integer = True
        for key in self.x_sol:
            if abs(self.x_sol[key] - self.x_int_sol[key]) > self.EPS:
                # one variable non-integer (EPS as criteria)
                self.is_integer = False
                self.branch_var_list.append(key)
        if self.is_integer:
            self.way_of_opt = 'By Simplex'
        else:
            self.way_of_opt = 'By Rounding'
        
        """ Update LP / IP"""
        self.LP_obj = self.model.ObjVal
        if self.is_integer:
            self.IP_obj = self.model.ObjVal
        # else: # round to get integer solution
        #     UB_value = 0
        #     # 补充针对VRPTW的取整方法(但效果不大) : 某点无入度则补充从depot过去的路程，某点无出度则补充回到depot的路程，最坏情况所有点都单点来回往返
        #     points = list(range(p_num))
        #     in_free_list = {} # 如果有从i出去的路，则不添加从i到0的路
        #     out_free_list = {} # 如果有到达j的路，则不添加从0到j的路
        #     for i in range(1, p_num):
        #         for j in range(1, p_num):
        #             varname = "x[{},{}]".format(i, j)
        #             if self.x_int_sol[varname] == 1:
        #                 varname = "s[{},{}]".format(i, j)
        #                 if self.x_int_sol["s[{}]".format(i)] >= self.x_int_sol["s[{}]".format(j)]: 
        #                     # 若违法子圈约束, 则仍需要添加进出的路
        #                     continue
        #                 in_free_list[i] = 1
        #                 out_free_list[j] = 1

        #     for i in range(1, p_num):
        #         if i not in in_free_list:
        #             varname = "x[{},{}]".format(0, i)
        #             self.x_int_sol[varname] = 1
        #         if i not in out_free_list:
        #             varname = "x[{},{}]".format(i, 0)
        #             self.x_int_sol[varname] = 1

        #     for var in self.model.getVars():
        #         UB_value += var.obj * self.x_int_sol[var.VarName]
        #     self.IP_obj = UB_value

        return
        
class Branch_and_Bound():
    """
    Branch and Bound Algorithm class
    """
    def __init__(self, model):
        model.setParam('OutputFlag', 0) # don't show gurobi information
        self.original_model = model
        self.root_node = Node(model=model.relax()) # node contain the relaxed model only
        self.node_list = [] # the nodes for exploration
        self.global_LB = -np.inf # the best we could get possibly
        self.global_UB = np.inf # the best integer solution we get
        self.incumbent_node = self.root_node # record the best node
        
        """ Strategy part """
        self.branch_strategy = 'first found' # 'min_inf' / 'max_inf' / first found'
        self.search_strategy = 'best LB first'

        """ display part"""
        self.iter_cnt = 0 # record iteration times
        self.Gap = np.inf # Gap between incumbent and LB
        self.fea_sol_cnt = 0 # count the number of feasible solution
        self.BB_tree_size = 0 # size of Branch and Bound tree

    def root_init(self):
        """
        Solve the root node, initialize the global UB/LB and branch
        """
        self.root_node.solve_and_update()
        assert self.root_node.is_feasible == True, "Original problem not feasible"
        self.global_UB = self.root_node.IP_obj
        self.global_LB = self.root_node.LP_obj
        self.incumbent_node = self.root_node
        self.current_node = self.root_node
        if not self.root_node.is_integer:
            self.branch(self.root_node)
        
    def branch(self, node):
        """
        branch according to the branch_var_list
        """
        if self.branch_strategy == 'first found':
            varname = node.branch_var_list[0]
            value = node.x_sol[varname]
        elif self.branch_strategy == 'min_inf':
            min_diff = np.inf
            for vi, varname in enumerate(node.branch_var_list):
                value = node.x_sol[varname]
                diff = min(value-math.floor(value), math.ceil(value)-value)
                if diff < min_diff:
                    best_varname = varname
            varname = best_varname
            value = node.x_sol[varname]
        elif self.branch_strategy == 'max_inf':
            max_inf = -np.inf
            for vi, varname in enumerate(node.branch_var_list):
                value = node.x_sol[varname]
                diff = min(value-math.floor(value), math.ceil(value)-value)
                if diff > max_inf:
                    best_varname = varname
            varname = best_varname
            value = node.x_sol[varname]
        
        self.branch_var_name = varname

        left_node = Node(node.model, node.LP_obj) # copy model and set LP_obj as local_LB
        left_node.depth = node.depth + 1
        self.BB_tree_size += 1
        var_l = left_node.model.getVarByName(varname)
        left_node.model.addConstr(var_l <= math.floor(value))
        self.node_list.append(left_node)

        right_node = Node(node.model, node.LP_obj) # copy model and set LP_obj as local_LB
        right_node.depth = node.depth + 1
        self.BB_tree_size += 1
        var_r = right_node.model.getVarByName(varname)
        right_node.model.addConstr(var_r >= math.ceil(value))
        self.node_list.append(right_node)

    def search(self):
        """
        get a node from the node_list
        """
        if self.search_strategy == 'best LB first':
            best_LB = np.inf
            best_i = 0
            for i in range(len(self.node_list)):
                if self.node_list[i].local_LB < best_LB:
                    best_i = i
                    best_LB = self.node_list[i].local_LB 
            return self.node_list.pop(best_i)
    
    def display_MIP_logging(self):
        """
        Show the MIP logging.

        :param iter_cnt:
        :return:
        """

        if (self.iter_cnt <= 0):
            print('|%6s  |' % 'Iter', end='')
            print(' \t\t %1s \t\t  |' % 'BB tree', end='')
            print('\t %10s \t |' % 'Current Node', end='')
            print('    %11s    |' % 'Best Bounds', end='')
            print(' %8s |' % 'incumbent', end='')
            print(' %5s  |' % 'Gap', end='')
            print(' %5s  |' % 'Time', end='')
            print(' %6s |' % 'Feasible', end='')
            print('     %10s      |' % 'Branch Var', end='')
            print()
            print('| %4s   |' % 'Cnt', end='')
            print(' %5s |' % 'Depth', end='')
            print(' %8s |' % 'ExplNode', end='')
            print(' %10s |' % 'UnexplNode', end='')
            print(' %4s |' % 'InfCnt', end='')
            print('    %3s   |' % 'Obj', end='')
            print('%10s |' % 'PruneInfo', end='')
            print(' %7s |' % 'Best UB', end='')
            print(' %7s |' % 'Best LB', end='')
            print(' %8s |' % 'Objective', end='')
            print(' %5s  |' % '(%)', end='')
            print(' %5s  |' % '(s)', end='')
            print(' %8s |' % ' Sol Cnt', end='')
            print(' %7s  |' % 'Max Inf', end='')
            print(' %7s  |' % 'Max Inf', end='')
            print()
        if(self.incumbent_node == None):
            print('%2s' % ' ', end='')
        elif(self.incumbent_node.way_of_opt == 'By Rouding' and self.incumbent_node.has_showed_way_of_opt == False):
            print('%2s' % 'R ', end='')
            self.incumbent_node.has_showed_way_of_opt = True
        elif (self.incumbent_node.way_of_opt == 'By Simplex' and self.incumbent_node.has_showed_way_of_opt == False):
            print('%2s' % '* ', end='')
            self.incumbent_node.has_showed_way_of_opt = True
        # elif (self.current_node.has_a_int_sol_by_heur == True and incumbent_node.has_showed_heu_int_fea == False):
        #     print('%2s' % 'H ', end='')
        #     self.incumbent_node.has_showed_heu_int_fea = True
        else:
            print('%2s' % ' ', end='')

        print('%3s' % self.iter_cnt, end='')
        print('%10s' % self.current_node.depth, end='')
        print('%9s' % self.iter_cnt, end='')
        print('%12s' % len(self.node_list), end='')
        if (len(self.current_node.branch_var_list) > 0):
            print('%11s' % len(self.current_node.branch_var_list), end='')
        else:
            if (self.current_node.model.status == 2):
                print('%11s' % 'Fea Int',
                        end='')  # indicates that this is a integer feasible solution, no variable is infeasible
            elif (self.current_node.model.status == 3):
                print('%11s' % 'Inf Model',
                        end='')  # indicates that this is a integer feasible solution, no variable is infeasible
            else:
                print('%11s' % '---',
                        end='')  # indicates that this is a integer feasible solution, no variable is infeasible
        if (self.current_node.model.status == 2):
            print('%12s' % round(self.current_node.model.ObjVal, 2), end='')
        else:
            print('%12s' % '---', end='')
        print('%10s' % self.current_node.prune_info, end='')
        print('%12s' % round(self.global_UB, 2), end='')
        print('%10s' % round(self.global_LB, 2), end='')
        if(self.incumbent_node == None):
            print('%11s' % '---', end='')
        else:
            print('%11s' % round(self.incumbent_node.IP_obj, 2), end='')
        if (self.Gap != '---'):
            print('%9s' % round(100 * self.Gap, 2), end='%')
        else:
            print('%8s' % 100 * self.Gap, end='')
        print('%8s' % round(self.end_time - self.start_time, 0), end='s')
        print('%9s' % self.fea_sol_cnt, end=' ')
        print('%14s' % self.branch_strategy, end='')
        print('%9s' % self.branch_var_name, end='')
        print()
    
    def show_result(self):
        self.CPU_time = time.time() - self.start_time
        print('\n')
        if len(self.node_list) == 0:
            print("Unexplored node list empty")
        else:
            print("Global LB and UB meet")
        print("Branch and bound terminates !!!")
        print("\n\n ------------ Summary ------------")
        print("Incumbent Obj: {}".format(self.incumbent_node.IP_obj))
        print("Gap: {}%".format(round(self.Gap * 100) if self.Gap < np.inf else 'inf')  )
        print("BB tree size: {}".format(self.BB_tree_size))
        print("CPU time: {}s".format(self.CPU_time))
        print(" --------- Solution  --------- ")
        for key in self.incumbent_node.x_int_sol.keys():
            if self.incumbent_node.x_int_sol[key] == 1:
                print('{} = {}'.format(key, self.incumbent_node.x_int_sol[key]))

    def run(self):
        self.start_time = time.time()
        """ initalize the node """
        self.root_init() # solve root node and update global_LB/UB and branch
        """ branch and bound """
        while len(self.node_list) > 0 and self.global_LB < self.global_UB:
            """ search part """
            self.current_node = self.search() # get a node from node_list

            """ solve and prune """
            # prune1: By Bnd
            if self.current_node.local_LB >= self.global_UB: 
                self.current_node.prune_info = 'By Bnd'
            else:
                # prune2: By Inf
                self.current_node.solve_and_update()
                if not self.current_node.is_feasible:            
                    self.current_node.prune_info = 'By Inf'
                else:
                    # prune3: By Opt
                    incum_update = False # record whether incumbent updated
                    if self.current_node.is_integer:                 
                        self.fea_sol_cnt += 1
                        if self.current_node.IP_obj < self.global_UB: # update best solution
                            self.global_UB = self.current_node.IP_obj
                            self.incumbent_node = self.current_node
                            incum_update = True
                            if (self.global_UB < np.inf):
                                self.Gap = (self.global_UB - self.global_LB) / self.global_UB
                        self.current_node.prune_info = 'By Opt'

            """ branch part """
            if self.current_node.prune_info == '---': # if not been pruned
                self.branch(self.current_node)
        
            """ display logging """
            if self.iter_cnt % 100 == 0 or incum_update: # display when iter 10 times or when update incumbent
                self.end_time = time.time()
                self.display_MIP_logging()
                
            self.iter_cnt += 1

        if len(self.node_list) == 0:
            self.global_LB = self.global_UB
        self.Gap = (self.global_UB - self.global_LB) / self.global_UB
        self.display_MIP_logging()
        """ show result """
        self.show_result()
        
# define model
def VRPTW(problem, lmt_pnum=1000):
    """
    VRPTW model
    """
    BIG_M = 1e7
    # read data
    vehicleNum = problem.vehicleNum
    capacity = problem.capacity
    customers = problem.customers[:lmt_pnum] # depot index 0, including x, y, demand, ready_time, due_time, serviceTime
    # seperate data
    location = problem.location[:lmt_pnum]
    demand = problem.demand[:lmt_pnum]
    readyTime = problem.readyTime[:lmt_pnum]
    dueTime = problem.dueTime[:lmt_pnum]
    serviceTime = problem.serviceTime[:lmt_pnum]

    # building model
    MODEL = gp.Model('VRPTW')
    # MODEL.setParam('OutputFlag', 0) # don't show gurobi information

    p_num = len(customers)
    points = list(range(p_num))
    A = [(i, j) for i in points for j in points]
    D = {(i, j): np.linalg.norm(location[i]-location[j]) for i, j in A}

    ## add variates
    x = MODEL.addVars(A, vtype=GRB.BINARY, name='x')
    s = MODEL.addVars(points, vtype=GRB.CONTINUOUS, name='s') # service start time
    c = MODEL.addVars(points, vtype=GRB.CONTINUOUS, name='c') # remained capacity
    ## set objective
    MODEL.modelSense = GRB.MINIMIZE
    MODEL.setObjective(gp.quicksum(x[i, j] * D[i, j] for i, j in A))
    ## set constraints
    ### 1. flow balance
    MODEL.addConstrs(gp.quicksum(x[i, j] for j in points if j!=i)==1 for i in points[1:]) # depot not included
    MODEL.addConstrs(gp.quicksum(x[i, j] for i in points if i!=j)==1 for j in points[1:]) # depot not included
    ### 2. avoid subring / self-loop
    MODEL.addConstrs(s[i] + D[i, j] + serviceTime[i] - BIG_M * (1 - x[i, j]) <= s[j] for i, j in A if j!=0)
    ### 3. time constraints
    MODEL.addConstrs(s[i] >= readyTime[i] for i in points)
    MODEL.addConstrs(s[i] <= dueTime[i] for i in points)
    ### 4. capacity constraints
    MODEL.addConstrs(c[i] - demand[j] + BIG_M * (1 - x[i, j])>= c[j] for i, j in A if j!=0)
    MODEL.addConstrs(c[i] <= capacity for i in points)
    MODEL.addConstrs(c[i] >= 0 for i in points)
    ### 5. vehicle number constraint
    MODEL.addConstr(gp.quicksum(x[0, j] for j in points) <= vehicleNum)

    # solve the model
    MODEL.optimize()

    return MODEL, location

def MIP():
    MIP_model = gp.Model('my model')

    # create decision variable
    x = {}
    x[0] = MIP_model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.INTEGER, name='x_0')
    x[1] = MIP_model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.INTEGER, name='x_1')

    # create objective
    # sense >=, <=, ==, MAX, MIN
    MIP_model.setObjective(-100 * x[0] - 150 * x[1], GRB.MINIMIZE) # transform to MINIMIZE 

    # add constraints
    MIP_model.addConstr(2 * x[0] + x[1] <= 10, name='con_0')
    MIP_model.addConstr(3 * x[0] + 6 * x[1] <= 40, name='con_1')

    # solve
    MIP_model.optimize()
    print('obj = {}'.format(MIP_model.ObjVal))
    for var in MIP_model.getVars():
        print('{} = {}'.format(var.VarName, var.x))
    
    return MIP_model

# show result of VRPTW
def draw_VRPTW(x_int_sol, location):
    p_num = len(location)
    
    plt.figure()
    plt.scatter(location[0:1, 0], location[0:1, 1], s = 150, c = 'r', marker='*')
    plt.scatter(location[1:, 0], location[1:, 1], c = 'b')
    connects = []
    for i in range(p_num):
        for j in range(p_num):
            if x_int_sol["x[{},{}]".format(i, j)] == 1:
                connects.append([i, j]) 
    for conn in connects:
        points = location[conn, :]
        plt.plot(points[:, 0], points[:, 1], c='r')
    
    plt.show()

if __name__ == "__main__":
    problem = read_data("Research\LinearProgramming\solomon_100\C101.txt")
    model, location = VRPTW(problem=problem, lmt_pnum=12)
    alg = Branch_and_Bound(model)
    alg.run()
    draw_VRPTW(alg.incumbent_node.x_int_sol, location)

    # model = MIP()
    # alg = Branch_and_Bound(model)
    # alg.run()
