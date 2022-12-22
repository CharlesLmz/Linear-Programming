# Simplex Method
# by Charles Lee
# on 2022.07.30

import math
import numpy as np
import matplotlib.pyplot as plt

class Result():
    """contrainer for optimization result"""
    def __init__(self, obj, x_values, x_b=None):
        self.obj = obj
        self.X = x_values
        self.x_b = x_b

def simplex_method(c, A, b, show=True):
    """solving LP in canonical form
    input: 
        c : value coefficients
        A : restraint coefficients
        b : resource coefficients
    output:
        x : optimal solution
    """
    c = np.array(c)
    A = np.array(A)
    b = np.array(b)

    assert len(c) <= len(A[0]), "error: number of variate varies"
    assert len(b) == len(A), "error: number of restraint varies"
    if len(c) < len(A[0]): # add 0 for relaxation variates if not enough
        c = np.hstack((c, [0]*(len(A[0])-len(c))))

    var_num = len(c)
    base_num = len(b)

    # select base variants
    x = list(range(var_num))
    x_b = x[-base_num:] # select the last base_num variates
    x_n = x[:-base_num] # others as nonbase variates

    while 1:
        B = A[:,x_b] # the coefficient of base variates
        N = A[:,x_n] # the coefficient of nonbase variates
        c_b = c[x_b] # the value coefficient of base variates
        c_n = c[x_n] # the value coefficient of nonbase variates
        
        # calculate inverse matrix of B (Gauss transform)
        if np.linalg.det(B) == 0: # if singular matrix, problem has no solution under current x_b
            ts_res = two_stage_method(A, b)
            if ts_res is None:
                if show:
                    print("problem has no solution")
                return Result(obj=None, x_values=None, x_b=x_b)
            x_b, x_n = ts_res
            continue

        B_inv = np.linalg.inv(B) # calculate the inverse matrix of B
        N_inv = np.dot(B_inv, N) # calculate the transformed N
        x_b_value = np.dot(B_inv,b) # calculate the value of base variates / omega_0j
        r_n = c_n - sum((c_b * N_inv.T).T) # calculate the check number r_n

        # if x_b_value not all positive, not feasible solution, apply two-stage method
        if min(x_b_value) < 0:
            ts_res = two_stage_method(A, b)
            if ts_res is None:
                if show:
                    print("problem has no solution")
                return Result(obj=None, x_values=None, x_b=x_b)
            x_b, x_n = ts_res
            continue

        # select variates in 
        if (r_n <= 0).all(): # end alg and return answer when r_n all negative
            obj = sum(c_b * x_b_value)
            x_values = np.zeros(var_num)
            x_values[x_b] = x_b_value
            return Result(obj=obj, x_values=x_values, x_b=x_b)
        in_N_index = np.argmax(r_n)

        # select variates out 
        min_epsilon = float('inf')
        for i in range(len(b)):
            if N_inv[i,in_N_index] <= 0:
                continue
            epsilon = x_b_value[i] / N_inv[i,in_N_index]
            if epsilon < min_epsilon:
                min_epsilon = epsilon
                out_B_index = i
        if min_epsilon == float('inf'):# if epsilon all negative, obj unbounded
            x_values = np.zeros(var_num)
            x_values[x_n[in_N_index]] = float('inf')
            obj = sum(c * x_values)
            if show:
                print("obj unbounded")
            return Result(obj=obj, x_values=x_values, x_b=x_b)
        
        # variates in and out
        x_n.append(x_b.pop(out_B_index))
        x_b.append(x_n.pop(in_N_index))

def dual_simplex_method(c, A, b, show=True):
    """solving LP in canonical form using dual simplex method (min default)
    input: 
        c : value coefficients
        A : restraint coefficients
        b : resource coefficients
    output:
        x : optimal solution
    """
    assert len(c) == len(A[0]), "error: number of variate varies"
    assert len(b) == len(A), "error: number of restraint varies"
    if len(c) < len(A[0]): # add 0 for relaxation variates if not enough
        c = np.hstack((c, [0]*len(A[0]-len(c))))

    c = np.array(c)
    A = np.array(A)
    b = np.array(b)
    
    var_num = len(c)
    base_num = len(b)

    # select base variants
    x = list(range(var_num))
    x_b = x[-base_num:] # select the last base_num variates
    x_n = x[:-base_num] # others as nonbase variates

    while 1:
        B = A[:,x_b] # the coefficient of base variates
        N = A[:,x_n] # the coefficient of nonbase variates
        c_b = c[x_b] # the value coefficient of base variates
        c_n = c[x_n] # the value coefficient of nonbase variates
        
        # calculate inverse matrix of B (Gauss transform)
        if np.linalg.det(B) == 0: # if singular matrix, problem has no solution
            if show:
                print("problem has no solution")
            return Result(obj=None, x_values=None, x_b=x_b)
        B_inv = np.linalg.inv(B)  # calculate the inverse matrix of B
        N_inv = np.dot(B_inv, N)  # calculate the transformed N

        x_b_value = np.dot(B_inv,b) # calculate the value of base variates / omega_0j
        r_n = c_n - sum((c_b * N_inv.T).T) # calculate the check number r_n

        # if r_n not all negative, not dual feasible solution
        assert max(r_n) <= 0, "initial dual feasible solution needed"
        # select variates out
        if min(x_b_value) >= 0:
            obj = sum(c_b * x_b_value)
            x_values = np.zeros(var_num)
            x_values[x_b] = x_b_value
            return Result(obj=obj, x_values=x_values, x_b=x_b)
        out_B_index = np.argmin(x_b_value) # select the variate out with lowest value

        # select variates in 
        min_epsilon = float('inf')
        for j in range(len(x_n)):
            if N_inv[out_B_index, j] >= 0:
                continue
            epsilon = -r_n[j] / N_inv[out_B_index, j]
            if epsilon < min_epsilon:
                min_epsilon = epsilon
                in_N_index = j
        if min_epsilon == float('inf'):# if epsilon all negative, obj unbounded
            x_values = np.zeros(var_num)
            x_values[x_b[out_B_index]] = float('inf')
            obj = sum(c * x_values)
            if show:
                print("obj unbounded")
            return Result(obj=obj, x_values=x_values, x_b=x_b)

        # variates in and out
        x_n.append(x_b.pop(out_B_index))
        x_b.append(x_n.pop(in_N_index))

def two_stage_method(A, b):
    '''two-stage method for initial feasible solution
    input:
        A : restraint matrix, np-array
        b : resource coefficients, np-array
    output:
        x_b, x_n : the indexs of base variates
    '''
    init_var_num = len(A[0]) # initial variates number
    A_ = np.hstack((A, np.eye(len(A)))) # add artifitial variates
    c = np.zeros(len(A_[0]))
    c[-len(A_):] = -1 # minimize artifitial variates
    res =  simplex_method(c, A_, b, show=False)
    if res.obj is None or res.obj != 0: # if artifitial variates not all 0, no solution
        return None
    x_b = res.x_b
    x_n = [n for n in list(range(init_var_num)) if n not in x_b]
    return x_b, x_n

if __name__ == "__main__":
    # problem1 (obj 36)
    # c = [3,4,5]
    # A = [[1,1,1,1,0,0],[1,0,3,0,1,0],[0,4,2,0,0,1]]
    # b = [10,12,12]
    # problem2 (obj 1.69)
    # c = [-0.8, -0.5, -0.9, -1.5]
    # A = [[1000,1500,1750,3250,-1,0,0],[0.6,0.27,0.68,0.3,0,-1,0],[17.5,7.7,0,30,0,0,-1]]
    # b = [4000,1,30]
    # problem3 
    # c = [1,2]
    # A = [[3,1,1,0],[3,4,0,1]]
    # b = [6,12]
    # problem4
    # c = [1/2,3/2]
    # A = [[1,3,1,0],[1,1,0,-1]]
    # b = [6,4]
    # problem5
    c = [5,4]
    A = [[1,1,1,0,0],[10,6,0,1,0],[1,0,0,0,-1]]
    b = [5,45,4]

    res = simplex_method(c, A, b)
    print('X:{}'.format(res.X))
    print('obj:{}'.format(res.obj))