# some heuristrics for VRP
# author: Charles Lee
# date: 2022.09.23

import numpy as np
import matplotlib.pyplot as plt
import math
import copy
from read_data import read_data
import time

# constructive heuristics
def solomon_insertion(problem):
    """solomon insertion algorithm to get an initial solution for VRP

    Args:
        problem (Problem): all information needed in VRPTW
            positions (ndarray[N, 2]): positions of all points, depot as index 0
            demands (ndarray[N]): demands of all points, depot as 0
            capacity (int): capacity of each car

    Returns:
        routes (List): routes consist of idx of points
    """

    """ set paraments """
    miu = 1
    lamda = 1 # ps: lambda is key word
    alpha1 = 1
    alpha2 = 0

    """ read data and preprocess """
    vehicle_num = problem.vehicleNum
    capacity = problem.capacity
    positions = problem.location
    demands = problem.demand
    ready_times = problem.readyTime
    due_times = problem.dueTime
    service_time = problem.serviceTime
    p_num = problem.nodeNum
    dist_m = problem.disMatrix
    unassigned_points = list(range(1, p_num)) 

    """ construct a route each circulation """
    routes = []
    while len(unassigned_points) > 0: 
        # initiate load, point_list
        load = 0
        point_list = unassigned_points.copy() # the candidate point set
        route_start_time_list = [0] # contains time when service started each point
        # choose the farthest point as s
        max_d = 0
        for p in point_list:
            dist = dist_m[0, p]
            # start_time = max(dist, ready_times[p])
            # if start_time > due_times[p]: # exclude point break time constraint
            #     continue
            if dist > max_d:
                max_d = dist
                max_p = p # farthest point as max_pi
                # best_start_time = start_time
        route = [0, max_p] # route contains depot andcustomer points 
        # route_start_time_list.append(best_start_time) 
        point_list.remove(max_p) # ? change data structure to save time cost
        unassigned_points.remove(max_p)
        load += demands[max_p]

        """ add a point each circulation """
        while len(point_list) > 0:
            c2_list = [] # contains the best c1 value
            best_insert_list = [] # contains the best insert position
            # ? add data structure to save time cost
            resi_load = capacity - load
            # find the insert position with lowest additional distance
            pi = 0
            while pi < len(point_list):
                u = point_list[pi]
                # remove if over load
                if demands[u] > resi_load: 
                    point_list.pop(pi)
                    continue
                
                best_c1 = np.inf 
                for ri in range(len(route)):
                    i = route[ri]
                    if ri == len(route)-1:
                        rj = 0
                    else:
                        rj = ri+1
                    j = route[rj]
                    # c11 = diu + dui - miu*dij
                    c11 = dist_m[i, u] + dist_m[u, j] - miu * dist_m[i, j]
                    # c12 = bju - bj # ! seems sth wrong?
                    # bj = route_start_time_list[rj]
                    # bu = max(route_start_time_list[ri] + service_time[i] + dist_m[i, u], ready_times[u])
                    # bju = bu + service_time[u] + dist_m[u, j]
                    # c12 = bju - bj
                    c12 = 0

                    # remove if over time window
                    # if bu > due_times[u] or bju > due_times[j]:
                    #     continue
                    # PF = c12
                    # pf_rj = rj
                    # overtime_flag = 0
                    # while PF > 0 and pf_rj < len(route)-1:
                    #     pf_rj += 1
                    #     bju = max(bju + service_time[route[pf_rj-1]] + dist_m[route[pf_rj-1], route[pf_rj]], \
                    #         ready_times[route[pf_rj]]) # start time of pf_rj
                    #     if bju > due_times[route[pf_rj]]:
                    #         overtime_flag = 1
                    #         break
                    #     PF = bju - route_start_time_list[pf_rj] # time delay
                    # if overtime_flag == 1:
                    #     continue

                    # c1 = alpha1*c11(i,u,j) + alpha2*c12(i,u,j)
                    c1 = alpha1*c11 + alpha2*c12
                    # find the insert pos with best c1
                    if c1 < best_c1:
                        best_c1 = c1
                        best_insert = ri+1
                # remove if over time (in all insert pos)
                # if best_c1 == np.inf:
                #     point_list.pop(pi)
                #     continue
                c2 = lamda * dist_m[0, u] - best_c1
                c2_list.append(c2)
                best_insert_list.append(best_insert)
                pi += 1
            if len(point_list) == 0:
                break
            # choose the best point
            best_pi = np.argmax(c2_list)
            best_u = point_list[best_pi]
            best_u_insert = best_insert_list[best_pi] 
            # update route
            route.insert(best_u_insert, best_u)
            point_list.remove(best_u)
            unassigned_points.remove(best_u) # when point is assigned, remove from unassigned_points
            load += demands[best_u]
            # update start_time
            # start_time = max(route_start_time_list[best_u_insert-1] + service_time[route[best_u_insert-1]] + dist_m[i, u], ready_times[best_u])
            # route_start_time_list.insert(best_u_insert, start_time)
            # for ri in range(best_u_insert+1, len(route)):
            #     start_time = max(route_start_time_list[best_u_insert-1] + service_time[route[best_u_insert-1]] + dist_m[i, u], ready_times[best_u])
            #     route_start_time_list[ri] = start_time
        route.append(0)
        routes.append(route) 
    return routes
                    
def nearest_neighbour(problem):
    """nearest neighbour algorithm to get an initial solution for VRP

    Args:
        problem (Problem): all information needed in VRPTW
            positions (ndarray[N, 2]): positions of all points, depot as index 0
            demands (ndarray[N]): demands of all points, depot as 0
            capacity (int): capacity of each car

    Returns:
        routes (List): routes consist of idx of points
    """

    """ read data and preprocess """
    vehicle_num = problem.vehicleNum
    capacity = problem.capacity
    positions = problem.location
    demands = problem.demand
    ready_times = problem.readyTime
    due_times = problem.dueTime
    service_time = problem.serviceTime
    p_num = problem.nodeNum
    dist_m = problem.disMatrix
    unassigned_points = list(range(1, p_num)) 
    
    routes = []
    """ construct a route each circulation """
    while len(unassigned_points) > 0:
        points_list = unassigned_points.copy()
        route = [0]
        cur_p = 0
        load = 0
        """ add a point each circulation """
        while len(points_list) > 0:
            min_d = np.inf
            pi = 0
            while pi < len(points_list):
                p = points_list[pi]
                if load + demands[p] > capacity:
                    points_list.remove(p)
                    continue
                dist = dist_m[cur_p, p]
                if dist < min_d:
                    min_d = dist
                    best_p = p
                pi += 1
            if len(points_list) == 0:
                break
            route.append(best_p)
            points_list.remove(best_p)
            unassigned_points.remove(best_p)
            load += demands[best_p]
            cur_p = best_p
        route.append(0)
        routes.append(route)
    return routes

def nearest_addition(problem):
    """nearest addition algorithm to get an initial solution for VRP

    Args:
        problem (Problem): all information needed in VRPTW
            positions (ndarray[N, 2]): positions of all points, depot as index 0
            demands (ndarray[N]): demands of all points, depot as 0
            capacity (int): capacity of each car

    Returns:
        routes (List): routes consist of idx of points
    """

    """ read data and preprocess """
    vehicle_num = problem.vehicleNum
    capacity = problem.capacity
    positions = problem.location
    demands = problem.demand
    ready_times = problem.readyTime
    due_times = problem.dueTime
    service_time = problem.serviceTime
    p_num = problem.nodeNum
    dist_m = problem.disMatrix
    unassigned_points = list(range(1, p_num)) 
    
    routes = []
    """ construct a route each circulation """
    while len(unassigned_points) > 0:
        points_list = unassigned_points.copy()
        route = [0]
        cur_p = 0
        load = 0
        """ add a point each circulation """
        while len(points_list) > 0:
            min_addition = np.inf
            pi = 0
            while pi < len(points_list):
                p = points_list[pi]
                if load + demands[p] > capacity:
                    points_list.remove(p)
                    continue
                # calculate addition
                min_d = np.inf
                for ri in range(len(route)):
                    if ri == len(route)-1:
                        rj = 0
                    else:
                        rj = ri + 1
                    i, j = route[ri], route[rj]
                    dist_add = dist_m[i, p] + dist_m[p, j] - dist_m[i, j]
                    if dist_add < min_d:
                        min_d = dist_add
                        best_insert = ri+1
                if min_d < min_addition:
                    min_addition = min_d
                    best_p = p
                    best_p_insert = best_insert
                pi += 1
            if len(points_list) == 0:
                break
            route.insert(best_p_insert, best_p)
            points_list.remove(best_p)
            unassigned_points.remove(best_p)
            load += demands[best_p]
            cur_p = best_p
        route.append(0)
        routes.append(route)
    return routes

def farthest_addition(problem):
    """farthest addition algorithm to get an initial solution for VRP

    Args:
        problem (Problem): all information needed in VRPTW
            positions (ndarray[N, 2]): positions of all points, depot as index 0
            demands (ndarray[N]): demands of all points, depot as 0
            capacity (int): capacity of each car

    Returns:
        routes (List): routes consist of idx of points
    """

    """ read data and preprocess """
    vehicle_num = problem.vehicleNum
    capacity = problem.capacity
    positions = problem.location
    demands = problem.demand
    ready_times = problem.readyTime
    due_times = problem.dueTime
    service_time = problem.serviceTime
    p_num = problem.nodeNum
    dist_m = problem.disMatrix
    unassigned_points = list(range(1, p_num)) 
    
    routes = []
    """ construct a route each circulation """
    while len(unassigned_points) > 0:
        points_list = unassigned_points.copy()
        route = [0]
        cur_p = 0
        load = 0
        """ add a point each circulation """
        while len(points_list) > 0:
            max_addition = 0
            pi = 0
            if len(route) == 83:
                print("")
            while pi < len(points_list):
                p = points_list[pi]
                if load + demands[p] > capacity:
                    points_list.remove(p)
                    continue
                # calculate addition
                min_d = np.inf
                for ri in range(len(route)):
                    if ri == len(route)-1:
                        rj = 0
                    else:
                        rj = ri + 1
                    i, j = route[ri], route[rj]
                    dist_add = dist_m[i, p] + dist_m[p, j] - dist_m[i, j]
                    if dist_add < min_d:
                        min_d = dist_add
                        best_insert = ri+1
                if min_d >= max_addition:
                    max_addition = min_d
                    best_p = p
                    best_p_insert = best_insert
                pi += 1
            if len(points_list) == 0:
                break
            route.insert(best_p_insert, best_p)
            points_list.remove(best_p)
            unassigned_points.remove(best_p)
            load += demands[best_p]
            cur_p = best_p
        route.append(0)
        routes.append(route)
    return routes

def CW_saving(problem):
    """Clarke-Wright Saving Algorithm to get an initial solution for VRP

    Args:
        problem (Problem): all information needed in VRPTW
            positions (ndarray[N, 2]): positions of all points, depot as index 0
            demands (ndarray[N]): demands of all points, depot as 0
            capacity (int): capacity of each car

    Returns:
        routes (List): routes consist of idx of points
    """

    """ read data and preprocess """
    vehicle_num = problem.vehicleNum
    capacity = problem.capacity
    positions = problem.location
    demands = problem.demand
    ready_times = problem.readyTime
    due_times = problem.dueTime
    service_time = problem.serviceTime
    p_num = problem.nodeNum
    dist_m = problem.disMatrix
    unassigned_points = list(range(1, p_num)) 
    
    """ initial allocation of one vehicle to each customer """
    X = np.zeros((p_num, p_num)) # the connection matrix, X[i, j]=1 shows i to j
    for p in range(1, p_num):
        X[p, 0] = 1
        X[0, p] = 1
    
    """ calculate saving sij and order """
    S = []
    for i in range(1, p_num):
        for j in range(i+1, p_num):
            sij = dist_m[0, i] + dist_m[j, 0] - dist_m[i, j]
            S.append([i, j, sij])
    S.sort(key=lambda s:s[2]) # sort by sij in increasing order

    """ each step find the largest sij and link them """ 
    out_map = {} # save points already out to other point
    in_map = {} # save points already in by other point
    while len(S) > 0:
        ss = S.pop()
        i, j = ss[:2]
        # exclude if already been connected
        if i in out_map or j in in_map:
            continue
        # exclude if overload
        load_l = demands[i]
        load_r = demands[j]
        i_ = i
        j_ = j
        while 1: # find the previous point until 0
            for i_pre in range(p_num): 
                if X[i_pre, i_] == 1:
                    break
            if i_pre == 0:
                break
            load_l += demands[i_pre]
            i_ = i_pre
        while 1: # find the next point until 0
            for j_next in range(p_num): 
                if X[j_, j_next] == 1:
                    break
            if j_next == 0:
                break
            load_r += demands[j_next]
            j_ = j_next
        total_load = load_l + load_r
        if total_load > capacity: # exclude
            continue
        # link i and j
        X[i, 0] = 0
        X[i, j] = 1
        X[0, j] = 0
        out_map[i] = 1
        in_map[j] = 1
    
    """ translate X to route """
    routes = []
    for j in range(1, p_num):
        if X[0, j] == 1:
            route = [0]
            route.append(j)
            i = j
            while j != 0:
                for j in range(p_num):
                    if X[i, j] == 1:
                        route.append(j)
                        i = j
                        break
            routes.append(route)
    
    return routes
                
def sweep_algorithm(problem):
    """ sweep algorithm to get an initial solution for VRP

    Args:
        problem (Problem): all information needed in VRPTW
            positions (ndarray[N, 2]): positions of all points, depot as index 0
            demands (ndarray[N]): demands of all points, depot as 0
            capacity (int): capacity of each car

    Returns:
        routes (List): routes consist of idx of points
    """

    """ read data and preprocess """
    vehicle_num = problem.vehicleNum
    capacity = problem.capacity
    positions = problem.location
    demands = problem.demand
    ready_times = problem.readyTime
    due_times = problem.dueTime
    service_time = problem.serviceTime
    p_num = problem.nodeNum
    dist_m = problem.disMatrix
    unassigned_points = list(range(1, p_num)) 
    
    """ sort unassigned points by angle """
    points_angles = np.zeros(p_num)
    for i in range(1, p_num):
        y_axis = positions[i, 1] - positions[0, 1]
        x_axis = positions[i, 0] - positions[0, 0]
        r = np.sqrt(x_axis**2 + y_axis**2)
        cospi = x_axis / r
        angle = math.acos(cospi)
        if y_axis < 0:
            angle = 2*np.pi - angle 
        points_angles[i] = angle
    sort_idxs = np.argsort(-points_angles) # sort by angle in decrease order
    unassigned_points = sort_idxs.tolist()

    """ construct a route each circulation """
    routes = [[0]]
    routes_load = [0]
    while len(unassigned_points) > 0:
        p = unassigned_points.pop()
        if routes_load[-1] + demands[p] < capacity:
            routes[-1].append(p)
            routes_load[-1] += demands[p]
        else:
            routes[-1].append(0)
            routes.append([0, p])
            routes_load.append(demands[p])
    routes[-1].append(0)
    
    return routes

def cluster_routing(problem):
    """ two-phase (cluster first, routing second) algorithm to get an initial solution for VRP

    Args:
        problem (Problem): all information needed in VRPTW
            positions (ndarray[N, 2]): positions of all points, depot as index 0
            demands (ndarray[N]): demands of all points, depot as 0
            capacity (int): capacity of each car

    Returns:
        routes (List): routes consist of idx of points
    """

    """ set paraments """
    cluster_num = 4
    diff_eps = 1e-2

    """ read data and preprocess """
    vehicle_num = problem.vehicleNum
    capacity = problem.capacity
    positions = problem.location
    demands = problem.demand
    ready_times = problem.readyTime
    due_times = problem.dueTime
    service_time = problem.serviceTime
    p_num = problem.nodeNum
    dist_m = problem.disMatrix
    unassigned_points = list(range(1, p_num)) 

    """ cluster first """
    cluster_centers = positions[1:1+cluster_num].copy() # initiate cluster_centers with first points
    while 1:
        # find best cluster for each point
        clusters = [[] for _ in range(cluster_num)] # contains point_idxs of each cluster
        np.random.shuffle(unassigned_points) # shuffle to make randomness
        for ui in range(len(unassigned_points)):
            i = unassigned_points[ui]
            min_c = np.inf
            for k in range(cluster_num):
                jk = cluster_centers[k]
                d0i = dist_m[0, i]
                dijk = np.linalg.norm(positions[i] - cluster_centers[k])
                djk0 = np.linalg.norm(cluster_centers[k] - positions[0])
                cki = (d0i + dijk +djk0) - 2*djk0 # ? is the second part of formula needed?
                if cki < min_c:
                    min_c = cki
                    best_k = k
            clusters[best_k].append(i)
        # update cluster_centers, until nearly no change 
        diff = 0
        for k in range(cluster_num):
            assert len(clusters[k]) > 0, "cluster empty, maybe cluster number too high"
            center = np.mean(positions[clusters[k]], 0) 
            diff += sum(abs(cluster_centers[k] - center))
            cluster_centers[k] = center
        if diff < diff_eps:
            break
    
    """ show cluster result (optional) """
    show = True
    if show:
        plt.scatter(positions[:1, 0], positions[:1, 1], s=200, marker='*')
        plt.scatter(cluster_centers[:, 0], cluster_centers[:, 1], c = 'r', s = 100, marker='+')
        for cluster in clusters:
            plt.scatter(positions[cluster, 0], positions[cluster, 1])
    plt.show()

    """ routing second """
    routes = []
    for cluster in clusters:
        cluster.insert(0, 0) # add depot
        sub_problem = copy.deepcopy(problem)
        sub_problem.customers = sub_problem.customers[cluster]
        # apply other algorithm to do subrouting
        sub_routes = solomon_insertion(sub_problem)
        # translate sub_points to points
        for route in sub_routes:
            for ri in range(len(route)):
                route[ri] = cluster[route[ri]]
        routes += sub_routes
    
    return routes

# neighbour stuctures (operators)
class Relocate():
    def __init__(self, k=2):
        self.k = k # how many points relocate together, k=1:relocate, k>1:Or-Opt

    def run(self, solution):
        """relocate point and the point next to it randomly inter/inner route (capacity not considered)

        Args:
            solution (List[int]): idxs of points of each route (route seperate with idx 0)

        Returns:
            neighbours (List[List[int]]): idxs of points of each route (seperate with idx 0) of each neighbour 
        """
        neighbours = []
        # 1. choose a point to relocate
        for pi in range(len(solution)-1):
            if math.prod(solution[pi:pi+self.k]) == 0: # don't relocate 0
                continue
            # 2. choose a position to put
            for li in range(1, len(solution)-self.k): # can't relocate to start/end
                neighbour = solution.copy()
                points = []
                for _ in range(self.k):
                    points.append(neighbour.pop(pi))
                for p in points[::-1]:
                    neighbour.insert(li, p)
                neighbours.append(neighbour)
        return neighbours     

class Exchange():
    def __init__(self, k=1):
        self.k = k # how many points exchange together

    def run(self, solution):
        """exchange two points randomly inter/inner route (capacity not considered)
        ps: Exchange operator won't change the points number of each vehicle

        Args:
            solution (List[int]): idxs of points of each route (route seperate with idx 0)

        Returns:
            neighbours (List[List[int]]): idxs of points of each route (seperate with idx 0) of each neighbour 
        """
        neighbours = []
        # 1. choose point i
        for pi in range(len(solution)-self.k):
            # 2. choose point j
            for pj in range(pi+self.k, len(solution)-self.k): 
                if math.prod(solution[pi:pi+self.k]) == 0 or math.prod(solution[pj:pj+self.k]) == 0: # don't exchange 0
                    continue
                neighbour = solution.copy()
                tmp = neighbour[pi:pi+self.k].copy()
                neighbour[pi:pi+self.k] = neighbour[pj:pj+self.k]
                neighbour[pj:pj+self.k] = tmp
                neighbours.append(neighbour)
        return neighbours    

# tools
def evaluate(problem, routes):
    """evaluate the objective value and feasibility of route

    Args:
        problem (Problem): informations of VRPTW
        routes (List): solution of problem, to evaluate
    Return:
        obj (double): objective value of the route (total distance)
    """
    # read data
    vehicle_num = problem.vehicleNum
    capacity = problem.capacity
    positions = problem.location
    demands = problem.demand
    ready_times = problem.readyTime
    due_times = problem.dueTime
    service_time = problem.serviceTime
    p_num = problem.nodeNum
    dist_m = problem.disMatrix

    obj = 0
    # calculate total routes length
    total_dist = 0
    for route in routes:
        route_dist = 0
        for ri in range(len(route)):
            p1 = route[ri]
            if ri == len(route)-1:
                p2 = route[0]
            else:
                p2 = route[ri+1]
            dist = dist_m[p1, p2]
            route_dist += dist
        total_dist += route_dist

    # check capacity constraint
    overload_cnt = 0
    for route in routes:
        route_load = 0
        for ri in range(len(route)):
            route_load += demands[route[ri]]
        if route_load > capacity:
            overload_cnt += 1
            # obj += np.inf
    print('overload: {}routes'.format(overload_cnt))
    
    # check time window constraint
    overtime_cnt = 0
    for route in routes:
        cur_time = 0
        for ri in range(len(route)):
            p1 = route[ri]
            if ri == len(route)-1:
                p2 = route[0]
            else:
                p2 = route[ri+1]
            cur_time += np.linalg.norm(positions[p1] - positions[p2])
            if cur_time < ready_times[p2]:
                cur_time = ready_times[p2]
            if cur_time > due_times[p2]: # compare start_time with due_time
                overtime_cnt += 1
                # obj += np.inf
                break
            cur_time += service_time[p2]
    print('overtime: {}routes'.format(overtime_cnt))

    obj += total_dist
    return obj

def show_routes(positions, routes):
    for ri, route in enumerate(routes):
        print("route {}: {}".format(ri, route))
    plt.figure()
    plt.scatter(positions[1:, 0], positions[1:, 1])
    plt.scatter(positions[0:1, 0], positions[0:1, 1], s = 150, c = 'r', marker='*')
    for route in routes:
        plt.plot(positions[route, 0], positions[route, 1], c='r')
    plt.show()

if __name__ == "__main__":
    file_name = "solomon_100/C101.txt"
    problem = read_data(file_name)
    time1 = time.time()
    # routes = nearest_neighbour(problem)
    # routes = nearest_addition(problem)
    # routes = farthest_addition(problem)
    routes = solomon_insertion(problem)
    # routes = CW_saving(problem)
    # routes = sweep_algorithm(problem)
    # routes = cluster_routing(problem)
    time2 = time.time()
    obj = evaluate(problem, routes)
    show_routes(problem.customers[:, :2], routes)
    print("vehicel_num: {}, obj: {}, time consumption: {}".format(len(routes), obj, time2-time1))


