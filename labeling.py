# labeling algorithm for solving VRPTW subModel
# author: Charles Lee
# date: 2022.09.28

import numpy as np
import matplotlib.pyplot as plt
import math
from read_data import read_data

class Label():
    def __init__(self, path, q, t, obj):
        self.path = path # current route
        self.q = q # current load of route
        self.t = t # start time of node
        self.obj = obj # the check number
    
    @staticmethod
    def if_dominate(l1, l2, feasibleNodeSet):
        """check if l1 dominates l2 or on contrary

        Args:
            l1 (Label): label one 
            l2 (Label): label two
        Return:
            res (int): 0 stands for non-dominate, 1 for l1 dominate l2, 2 for l2 dominate l1
        """
        dominate_num = 0 
        if l1.obj <= l2.obj:
            dominate_num += 1
        if l1.q <= l2.q:
            dominate_num += 1
        if l1.t <= l2.t:
            dominate_num += 1
        
        if dominate_num == 3:
        # if dominate_num == 3 and set(l1.path).issubset(l2.path):
        # if dominate_num == 3 and len([pi for pi in l1.path if pi in feasibleNodeSet[l1.path[-1]]])==0:
            return 1
        elif dominate_num == 0:
        # elif dominate_num == 0 and set(l2.path).issubset(l1.path):
        # elif dominate_num == 0 and len([pi for pi in l2.path if pi in feasibleNodeSet[l2.path[-1]]])==0:
            return 2
        else:
            return 0

class Labeling():
    def __init__(self, graph, select_num=1):
        self.graph = graph
        self.select_num = select_num # routes number generated at one time
        self.Q = [[] for i in range(self.graph.nodeNum)] # queue for each points, containing part-routes that ends in the point

        self.nodeNum = self.graph.nodeNum
        self.dualValue = np.zeros(self.nodeNum)
        self.vehicleNum = self.graph.vehicleNum
        self.capacity = self.graph.capacity
        self.location = np.zeros((self.nodeNum, 2))
        for i in range(self.nodeNum):
            self.location[i] = np.array(self.graph.location[i])
        self.demand = self.graph.demand
        self.disMatrix = np.zeros((self.nodeNum, self.nodeNum))
        for i in range(self.nodeNum):
            for j in range(self.nodeNum):
                self.disMatrix[i, j] = self.graph.disMatrix[i][j]
        self.readyTime = self.graph.readyTime
        self.dueTime = self.graph.dueTime
        self.serviceTime = self.graph.serviceTime
        
    def set_dual(self, Dual):
        self.dualValue = Dual
    
    def dominant_add(self, label, node):
        """
        add label to node, while checking dominance
        input:
            label (Label): label to add
            node (int): idx of the node
        update:
            self.Q (dict[int:List]): queue for each points
        """
        li = 0
        while li < len(self.Q[node]):
            labeli = self.Q[node][li]
            flag = Label.if_dominate(label, labeli, self.graph.feasibleNodeSet)
            # if l1 dominates l2, pop(l2)
            if flag == 1:
                self.Q[node].pop(li)
            # if l2 dominates l1, not add l1
            elif flag == 2:
                return 
            li += 1
        self.Q[node].append(label)
    
    def node_expand(self, node):
        """
        expand each labels in the node
        input:
            node (int): idx of node to expand
        update:
            self.Q (dict[int:List]): queue of node 
        """
        availabel_list = self.graph.feasibleNodeSet[node]
        while self.Q[node]: # while not empty
            label = self.Q[node].pop()
            for next_node in availabel_list: # next_node: the next node
                if node == 0 and next_node == 0: #!
                    continue
                if next_node in label.path[1:]: # terminus 0 not included #!
                    continue
                q_ = label.q + self.demand[next_node]
                t_arrive = label.t + self.serviceTime[node] + self.disMatrix[node, next_node]
                if q_ > self.capacity or t_arrive > self.dueTime[next_node]: # check feasibility
                    continue
                t_ = max(self.readyTime[next_node], t_arrive)
                # the correlation formula
                obj_ = label.obj + self.disMatrix[node, next_node] - self.dualValue[next_node]
                path_ = label.path.copy()
                path_.append(next_node)
                new_label = Label(path_, q_, t_, obj_)
                self.dominant_add(new_label, next_node) # add node and check dominance
    
    def select_best(self):
        pareto_labels = self.Q[0] #!
        pareto_labels.sort(key=lambda label:label.obj)
        routes = [label.path for label in pareto_labels]
        objs = [label.obj for label in pareto_labels]
        return routes[:self.select_num], objs[:self.select_num]

    def run(self):
        label0 = Label([0], 0, 0, 0)
        self.Q[0].append(label0)
        self.node_expand(0) 
        while 1:
            # break if all queues empty
            break_flag = 1
            for node in range(1, self.nodeNum): #!
                queue = self.Q[node]
                if len(queue) > 0: # if not empty
                    break_flag = 0 # 0 as not break
            if break_flag:
                break
            # expand each node
            for node in range(1, self.nodeNum):
                self.node_expand(node)
        routes, objs = self.select_best()
        return routes, objs
        
if __name__ == "__main__":
    # graph = Data()
    # graph.read_solomon(path=r"labeling20220921/r101.txt", customerNum=100)
    # graph.get_euclidean_distance_matrix()
    graph = read_data("solomon_100/c101.txt")
    alg = Labeling(graph=graph, select_num=100)
    Dual = np.arange(0, alg.nodeNum)
    # Dual = [0] * alg.nodeNum
    Dual = np.arange(alg.nodeNum)
    alg.set_dual(Dual)
    routes, objs = alg.run()
    for ri, route in enumerate(routes):
        print("{} obj: {}, route: {}".format(ri+1, objs[ri], route))
    


