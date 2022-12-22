"""
Branch and Price
author: Charles Lee
date: 2022.11.07
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import time
from copy import deepcopy
import gurobipy as gp
from gurobipy import GRB
from read_data import read_data
from labeling import Labeling
from VRP_heuristics import *

class BPNode():
    def __init__(self, graph, RLMP, SP, global_column_pool):
        self.graph = graph
        self.global_column_pool = global_column_pool
        self.RLMP = RLMP.copy()
        self.SP = SP
        # model part
        self.duals_of_RLMP = {}
        self.SP_must_include = []
        self.SP_cant_include = []
        self.max_not_improve_cnt = 10
        # algorithm part
        self.local_LB = np.inf
        self.IP_obj = np.inf
        self.LP_obj = np.inf
        self.EPS = 1e-6 
        self.is_feasible = False
        self.is_integer = False
        self.inf_var_list = []
        self.x_sol = {}
        self.x_int_sol = {}
        # display part
        self.cg_iter_cnt = 0
        self.depth = 0
        self.has_showed_way_of_opt = False
        self.way_of_opt = "---"
        self.prune_info = "---"
    
    def generate(self):
        subNode = BPNode(graph, RLMP, SP, global_column_pool)
        subNode.local_LB = self.LP_obj
        subNode.depth = self.depth + 1
        subNode.SP_must_include = self.SP_must_include.copy()
        subNode.SP_cant_include = self.SP_cant_include.copy()
        return subNode
    
    def solve_and_update(self):
        """ solve and check feasibility """
        self.is_feasible = self.column_generation()
        if (~self.is_feasible):
            return
        """ update x_sol and round to x_int_sol """
        vars = self.RLMP.getVars()
        self.is_integer = True
        for var in vars:
            var_name = var.VarName
            var_val = var.X
            self.x_sol[var_name] = var_val
            self.x_int_sol[var_name] = round(var_val)
            # check integer
            if (abs(round(var_val) - var_val) > self.EPS):
                self.is_integer = False
                self.inf_var_list.append(var_name)
        """ update LP / IP """
        self.LP_obj = self.RLMP.ObjVal
        if self.is_integer:
            self.way_of_opt = "By Simplex"
            obj = 0
            for var in vars:
                var_name = var.VarName
                obj += self.x_int_sol[var_name] * var.X
            self.IP_obj = obj
        else:
            # solve RMP to get a integer solution as UB
            self.solve_final_RMP_and_update_IPobj()
    
    def column_generation(self):
        self.set_SP()
        best_RLMP_obj = np.inf
        not_improve_cnt = 0
        while True:
            # solve RLMP and get duals
            self.solve_RLMP_and_get_duals()
            if self.RLMP.ObjVal < best_RLMP_obj:
                not_improve_cnt = 0
                best_RLMP_obj = self.RLMP.ObjVal
                self.SP.setParam("PoolSolutions", 1) 
            else:
                not_improve_cnt += 1
                if not_improve_cnt > self.max_not_improve_cnt:
                    self.SP.setParam("PoolSolutions", not_improve_cnt//10) 
            # solve SP
            self.solve_SP()
            if self.SP.Status != 2:
                return 0 # node infeasible
            # break if can't find route to improve RLMP
            if self.SP.ObjVal >= -self.EPS:
                return 1 # node feasible
            # get columns and add into RLMP
            self.get_columns_from_SP_and_add()
            # self.get_columns_from_Labeling_and_add()
            self.cg_iter_cnt += 1
        self.reset_SP()
    
    def solve_RLMP_and_get_duals(self):
        self.RLMP.optimize()
        assert self.RLMP.Status == 2, "RLMP infeasible"
        for cons in self.RLMP.getConstrs():
            cons_name = cons.ConstrName 
            self.duals_of_RLMP[cons_name] = cons.Pi
    
    def solve_SP(self):
        # update objective with duals
        obj = 0.0
        for i in range(self.graph.nodeNum):
            for j in self.graph.feasibleNodeSet[i]:
                var_name = f"x[{i},{j}]"
                cons_name = f"R{i}"
                coef = self.graph.disMatrix[i, j] - self.duals_of_RLMP[cons_name]
                obj += coef * self.SP.getVarByName(var_name)
        self.SP.setObjective(obj, GRB.MINIMIZE) 
        # optimize SP
        self.SP.optimize()
                
    def set_SP(self):
        """ update must_include / cant_include constraints """
        for pair in self.SP_must_include:
            pi, pj = pair
            for j in range(self.graph.nodeNum):
                if j == pj:
                    continue
                var_name = f"x[{i},{j}]"
                self.SP.getVarByName(var_name).setParam("LB", 0) #! check
        for pair in self.SP_cant_include:
            pi, pj = pair
            var_name = f"x[{i},{j}]" 
            self.SP.getVarByName(var_name).setParam("LB", 0) #! check
        self.SP.update()

    def reset_SP(self):
        """ reset must_include / cant_include constraints """
        for pair in self.SP_must_include:
            pi, pj = pair
            for j in range(self.graph.nodeNum):
                if j == pj:
                    continue
                var_name = f"x[{i},{j}]"
                self.SP.getVarByName(var_name).setParam("LB", 1) #! check
        for pair in self.SP_cant_include:
            pi, pj = pair
            var_name = f"x[{i},{j}]" 
            self.SP.getVarByName(var_name).setParam("LB", 1) #! check
        self.SP.update()
    
    def get_columns_from_SP_and_add(self):
        """ get columns from SP and add into RLMP """
        new_route = []
        route_length = 0
        col_num = self.graph.nodeNum
        new_column = np.zeros(self.graph.nodeNum)
        for SolutionNumber in range(self.SP.SolCount): #!check
            # get column from SP
            self.get_column_from_SP(new_route, route_length, new_column, SolutionNumber)
            self.add_column_into_RLMP(new_route, route_length, new_column)
        # display CG iteration info
        print(f"CG_iter {self.cg_iter_cnt}: RMPobj={self.RLMP.ObjVal}, SPobj={self.SP.ObjVal}")
    
    def get_column_from_SP(self, new_route, route_length, new_column, SolutionNumber=0):
        new_route.append(0)
        current_i = 0
        while True:
            for j in self.graph.feasibleNodeSet[current_i]:
                var_name = f"x[{current_i},{j}]"
                self.SP.setParam("SolutionNumber", SolutionNumber)
                var_val = round(self.SP.getVarByName(var_name).X)
                if var_val == 1:
                    route_length += self.graph.disMatrix[current_i, j]
                    new_route.append(j)
                    new_column[j] = 1
                    current_i = j
                    break
            if current_i == 0:
                break

    def add_column_into_RLMP(self, new_route, route_length, new_column):
        # update column pool
        new_column_name = "y_" + str(len(self.global_column_pool))
        self.global_column_pool[new_column_name] = new_route
        # update RLMP
        new_RLMP_column = gp.Column()
        new_RLMP_column.addTerms(new_column, self.RLMP.getConstrs())
        self.RLMP.addVar(obj=route_length, vtype=GRB.CONTINUOUS, column=new_RLMP_column, name=new_column_name)
    
    def solve_final_RMP_and_update_IPobj(self):
        # convert RLMP into RMP
        for var in self.RLMP.getVars():
            var.vtype = 'B'
        # optimize RMP
        self.RLMP.optimize()
        # update IP_obj if feasible
        if self.RLMP.Statuc == 2:
            self.IP_obj = self.RLMP.ObjVal
            self.way_of_opt = "By RMP"
    
    def get_columns_from_Labeling_and_add(self):
        # create tmp_graph, update must_include / cant_include constraints
        tmp_graph = deepcopy(self.graph)
        for pair in self.SP_must_include:
            pi, pj = pair
            for j in range(tmp_graph.nodeNum):
                if j==pj:
                    continue
                tmp_graph.disMatrix[pi][j] = np.inf
        for pair in self.SP_cant_include:
            pi, pj = pair
            tmp_graph.disMatrix[pi][pj] = np.inf
        # create labeling algorithm, set duals and solve
        alg = Labeling(tmp_graph)
        duals = np.zeros(tmp_graph.nodeNum)
        for di in range(tmp_graph.nodeNum):
            cons_name = f"R{di}"
            cons = self.RLMP.getConstrByName(cons_name)
            duals[di] = cons.Pi
        alg.set_dual(duals)
        alg.run()
        routes = alg.best_routes
        objs = alg.best_objs
        min_obj = min(objs)
        # add routes into RLMP
        col_num = self.graph.nodeNum
        for route in routes: #!check
            # calculate route_length
            route_length = 0
            new_column = np.zeros(self.graph.nodeNum)
            for i in range(1, len(route)):
                route_length += tmp_graph.disMatrix[route[i-1], route[i]]
                new_column[i] = 1
            self.add_column_into_RLMP(route, route_length, new_column)
        print(f"CG_iter {self.cg_iter_cnt}: RMPob{self.RLMP.ObjVal}, min_obj={self.SP.ObjVal}")
        return min_obj
        
class BranchAndPrice():
    def __init__(self, graph):
        # graph info
        self.graph = graph
        self.global_column_pool = {}
        # build and set RLMP, SP
        self.RLMP = self.set_RLMP_model(graph)
        self.SP = self.set_SP_model(graph)    
        # create nodes
        self.root_node = BPNode(graph, self.RLMP, self.SP, self.global_column_pool)
        self.incumbent_node = self.root_node
        self.current_node = self.root_node
        # set strategies
        self.branch_strategy = "max_inf"
        self.search_strategy = "best_LB_first"
        # algorithm part
        self.node_list = []
        self.global_LB = np.inf
        self.global_UB = -np.inf
        # display parament
        self.iter_cnt = 0
        self.Gap = np.inf
        self.fea_sol_cnt = 0
        self.BP_tree_size = 0
        self.branch_var_name = ""

    def set_RLMP_model(self, graph):
        RLMP = gp.Model()
        # init solution with Solomon-Insert1
        routes = solomon_insertion(graph)
        routes_length = []
        routes_a = np.zeros((len(routes), graph.nodeNum))
        for ri, route in enumerate(routes):
            length = 0
            for pi in range(1, len(route)):
                length += graph.disMatrix[route[pi-1], route[pi]]
                routes_a[ri, route[pi]] = 1
            routes_length.append(length)
        # add init solution in RLMP
        ## add variables
        y_list = list(range(len(routes)))
        y = RLMP.addVars(y_list, vtype="C", name="y")
        ## set objective
        RLMP.setObjective(gp.quicksum(y[i] * routes_length[i] for i in range(len(routes))), GRB.MINIMIZE)
        ## set constraints
        RLMP.addConstr(gp.quicksum(y[i] for i in range(len(routes))) <= graph.vehicleNum)
        RLMP.addConstrs(gp.quicksum(y[i] * routes_a[i, j] for i in range(len(routes))) >= 1 for j in range(1, graph.nodeNum))

        RLMP.setParam("OutputFlag", 0)
        RLMP.setParam("PoolSearchMode", 2)
        RLMP.setParam("PoolSolutions", 1)
        RLMP.update()
        return RLMP
    
    def set_SP_model(self, graph):
        SP = gp.Model()
        ## add variables
        points = list(range(graph.nodeNum))
        A_list = [(i, j) for i in points for j in graph.feasibleNodeSet[i]]
        x = SP.addVars(A_list, vtype="B", name="x")
        t = SP.addVars(points, vtype="C", name="t")
        ## set objective 
        SP.setObjective(gp.quicksum(x[i, j] * graph.disMatrix[i][j] \
            for i in points for j in graph.feasibleNodeSet[i]))
        ## set constraints
        ### 1. flow balance
        SP.addConstrs(gp.quicksum(x[i, j] for j in graph.feasibleNodeSet[i] if j!=i)==1 for i in points[1:]) # depot not included
        SP.addConstrs(gp.quicksum(x[i, j] for i in graph.availableNodeSet[j] if i!=j)==1 for j in points[1:]) # depot not included
        ### 2. time window & sub-ring
        M = 1e7
        SP.addConstrs(t[i] + graph.disMatrix[i, j] + graph.serviceTime[i] - M * (1 - x[i, j]) <= t[j] for i, j in A_list if j!=0)
        SP.addConstrs(t[i] >= graph.readyTime[i] for i in points)
        SP.addConstrs(t[i] <= graph.dueTime[i] for i in points)

        # set model params
        SP.setParam("OutputFlag", 0)
        SP.update()
        return SP
    
    def root_init(self):
        self.root_node.solve_and_update()
        self.global_LB = self.root_node.LP_obj
        self.global_UB = self.root_node.IP_obj
        self.incumbent_node = self.root_node
        self.current_node = self.root_node
        if self.root_node.is_integer == False:
            self.branch(self.root_node)
    
    def search(self):
        """ best_LB_first: choose the node with best LB to search """
        best_node_i = 0
        if search_strategy == "best_LB_first":
            min_LB = np.inf
            for node_i in range(len(self.node_list)):
                LB = self.node_list[node_i].local_LB
                if LB < min_LB:
                    min_LB = LB
                    best_node_i = node_i
        best_node = self.node_list.pop(best_node_i)
        return best_node
    
    def branch(self, node):
        # get flow of each arc
        flow_matrix = np.zeros((self.graph.nodeNum, self.graph.nodeNum))
        vars = self.RLMP.getVars()
        for var in vars:
            var_val = var.X
            var_name = var.VarName
            if var_val > 0:
                route = self.global_column_pool[var_name]
                for j in range(1, len(route)):
                    i = j - 1
                    flow_matrix[i][j] += var_val
        # max_inf: choose the arc farthest to integer to branch 
        best_i = best_j = 0
        if self.branch_strategy == "max_inf":
            max_flow = -np.inf
            for i in range(self.graph.nodeNum):
                for j in range(self.graph.nodeNum):
                    if flow_matrix[i, j] > max_flow:
                        max_flow = flow_matrix[i, j]
                        best_i = i
                        best_j = j
        # branch on the chosen variable
        self.branch_var_name = f"x[{i},{j}]"
        ## 1. branch left, must include xij
        leftNode = node.generate()
        ### delete routes constains i or j, but xij != 1
        for var in vars:
            var_name = var.VarName
            route = self.global_column_pool[var_name]
            if best_i not in route and best_j not in route:
                continue
            elif best_i in route and best_j in route and \
                route.index(best_i)+1 == route.index(best_j):
                continue
            else:
                var.set("UB", 0) # make var_i invalid
        ### add into must include list
        leftNode.SP_must_include.append([best_i, best_j])
        self.node_list.append(leftNode)
        self.BP_tree_size+=1
        ## 2. branch right, must include xij
        rightNode = node.generate()
        ### delete routes with xij == 1
        for var in vars:
            var_name = var.VarName
            route = self.global_column_pool[var_name]
            if best_i in route and best_j in route and \
                route.index(best_i)+1 == route.index(best_j):
                var.set("UB", 0) # make var_i invalid
        ### add into must include list
        rightNode.SP_cant_include.append([best_i, best_j])
        self.node_list.append(rightNode)
        self.BP_tree_size+=1

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

if __name__ == "__main__":
    file_name = "solomon_100\C101.txt"
    graph = read_data(file_name)
    alg = BranchAndPrice(graph)
    alg.run()






