# read data from dateset

import numpy as np

class Problem():
    def __init__(self, vehicleNum, capacity, location, demand,\
        readyTime=None, dueTime=None, serviceTime=None):
        self.vehicleNum = vehicleNum
        self.capacity = capacity
        self.location = location
        self.demand = demand
        self.readyTime = readyTime
        self.dueTime = dueTime
        self.serviceTime = serviceTime
        self.nodeNum = len(self.location)

        if readyTime is not None:
            self.cal_disMatrix()
            self.cal_feasibleNodeSet()
        else:
            self.cal_disMatrix(roundFlag=True) 
    
    def cal_disMatrix(self, roundFlag=False):
        self.disMatrix = np.zeros((self.nodeNum, self.nodeNum))
        for i in range(self.nodeNum):
            for j in range(self.nodeNum):
                if roundFlag:
                    self.disMatrix[i, j] = round(np.linalg.norm(self.location[i] - self.location[j]))
                else:
                    self.disMatrix[i, j] = np.linalg.norm(self.location[i] - self.location[j])        

    def cal_feasibleNodeSet(self):
        self.feasibleNodeSet = [[] for _ in range(self.nodeNum)]
        self.availableNodeSet = [[] for _ in range(self.nodeNum)]
        for i in range(self.nodeNum):
            for j in range(self.nodeNum):
                if self.readyTime[i] + self.serviceTime[i] + self.disMatrix[i, j] <= self.dueTime[j] and i!=j:
                    self.feasibleNodeSet[i].append(j)
                    self.availableNodeSet[j].append(i)

def read_data(file_name):
    """
    read VRPTW data from dataset
    input: file_name
    output: problem object (including (int)vehicleNum, (int capacity, (numpy-array[25, 6])customers)
            ps:customers include x, y, demand, ready_time, due_time, service_time
    """
    with open(file_name) as file_object:
        lines = file_object.readlines()
    
    # load vehicle setting
    vehicle = list(map(int, lines[4].split()))
    vehicleNum, capacity = vehicle

    # load customers setting
    location = []
    demand = []
    readyTime = []
    dueTime = []
    serviceTime = []
    for line in lines[9:]:
        cust = list(map(int, line.split()))
        if cust == []:
            continue
        location.append(cust[1:3])
        demand.append(cust[3])
        readyTime.append(cust[4])
        dueTime.append(cust[5])
        serviceTime.append(cust[6])
    location = np.array(location)
    demand = np.array(demand)
    readyTime = np.array(readyTime)
    dueTime = np.array(dueTime)
    serviceTime = np.array(serviceTime)
    prob = Problem(vehicleNum, capacity, location, demand, readyTime, dueTime, serviceTime)
    return prob
        
def read_data_from_Augerat(file_name):
    """
    read VRPTW data from Augerat dataset
    input: file_name
    output: problem object (including (int)vehicleNum, (int capacity, (numpy-array[25, 6])customers)
            ps:customers include x, y, demand, ready_time, due_time, service_time
    """
    with open(file_name) as file_object:
        lines = file_object.readlines()
    
    # load vehicle setting
    vehicleNum = nodeNum = int(lines[3].split()[2])
    capacity = int(lines[5].split()[2])

    # load customers setting
    location = []
    demand = []
    for line in lines[7:7+nodeNum]:
        cust = list(map(int, line.split()))
        location.append(cust[1:3])
    for line in lines[8+nodeNum:8+2*nodeNum]:
        cust = list(map(int, line.split()))
        demand.append(cust[1])
    location = np.array(location)
    demand = np.array(demand)
    prob = Problem(vehicleNum, capacity, location, demand)
    return prob

if __name__ == "__main__":
    file_name = "solomon_100/C106.txt"
    prob = read_data(file_name)
    file_name = "A-n32-k5.txt"
    prob = read_data_from_Augerat(file_name)
