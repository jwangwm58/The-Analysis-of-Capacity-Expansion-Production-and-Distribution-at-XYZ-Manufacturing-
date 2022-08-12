"""
Created on Fri Dec 10 15:24:19 2021

"""

import gurobipy as gp
from gurobipy import GRB
from gurobipy import *
#i
U =[16000, 12000, 14000, 10000, 13000]
CO =[2000, 1600, 1800, 900, 1500]
OP =[420, 380, 460, 280, 340]
R =[190, 150, 160, 100, 130]
SD =[170, 120, 130, 80, 110] #shutdown

C = [[0.12, 0.13, 0.08, 0.05],
[0.1, 0.03, 0.1, 0.09],
[0.05, 0.07, 0.06, 0.03],
[0.06, 0.03, 0.07, 0.07],
[0.06, 0.02, 0.04, 0.08]] #Cij


D = [1000,1200,1800,1200,1000,1400,1600,1000] #Dk
DR = [0.20,0.20,0.25,0.20,0.20,0.25,0.25,0.2] #Dr_k


S = [[0.09,0.1,0.06,0.05,0.08,0.09,0.02,0.12],
[0.05,0.07,0.12,0.04,0.03,0.09,0.03,0.08],
[0.06,0.09,0.07,0.09,0.09,0.04,0.11,0.07],
[0.07,0.08,0.09,0.06,0.1,0.07,0.06,0.09]]

M = 1000000

T = 10 #10 years  t
I = len(U) #5  i
J = len(C[0]) #4  j
K = len(D) #8  k



m = gp.Model("Model")
m.modelSense = GRB.MINIMIZE

X = m.addVars(I, J, T, vtype=GRB.CONTINUOUS, name = "x")
V = m.addVars(J, K, T, vtype=GRB.CONTINUOUS, name = "V") #Vjkt
IV = m.addVars(J, T, vtype=GRB.CONTINUOUS, name = "IV") #IVjt
P = m.addVars(I, T, vtype=GRB.BINARY, name = "IV") #Pit
Z = m.addVars(I, T, vtype=GRB.BINARY, name = "IV") #Zit constructed
SZ = m.addVars(I, T, vtype=GRB.INTEGER,ub=1, name = "IV") #SZit

m.addConstrs((SZ[i,t] == quicksum(Z[i,t1] for t1 in range(t+1))
              for i in range(I) for t in range(T)), "cons1")
m.addConstrs((P[i,t] <= SZ[i,t]
              for i in range(I) for t in range(T)), "cons2")

PR = m.addVars(I, T, vtype=GRB.BINARY, name = "PR") #PRit
PS = m.addVars(I, T, vtype=GRB.BINARY, name = "PS") #PSit
W = m.addVars(I, T, vtype=GRB.CONTINUOUS, name = "W") #Wit
over = m.addVars(I, T, vtype=GRB.BINARY, name = "over") #Over_it

m.addConstrs((P[i,t] - P[i,t-1] <= PR[i,t] - PS[i,t]
              for i in range(I) for t in range(1,T)), "cons3")

m.addConstrs((PR[i,t] <= P[i,t]
              for i in range(I) for t in range(T)), "cons19")

m.addConstrs((PS[i,t] <= P[i,t-1]
              for i in range(I) for t in range(1,T)), "cons20")

m.addConstrs((PR[i,t] <= (1-P[i,t-1])
              for i in range(I) for t in range(1,T)), "cons21")

#very important
m.addConstrs( ((P[i,t] - P[i,t-1]) + M*PS[i,t] >=0
              for i in range(I) for t in range(1,T)), "cons22")

m.addConstrs((Z[i,t] <= PR[i,t] 
              for i in range(I) for t in range(T)), "cons4")

m.addConstrs((W[i,t] <= X.sum(i,'*',t) * 3
              for i in range(I) for t in range(T)), "cons5")
m.addConstrs(((W[i,t] - 9000) * (over[i,t] - 0.5) >= 0
              for i in range(I) for t in range(T)), "cons6")

cost_W = quicksum( ((1 - over[i,t])*W[i,t] * 0.15 + over[i,t]* (9000 * 0.15 + (W[i,t] - 9000)*0.12))* (1+0.3)**(t-1)
                  for i in range(I) for t in range(T))
cost_A = quicksum(W[i,t] * 4.7 * 0.02 * (1+0.3)**(t-1) for i in range(I) for t in range(T))




m.setObjective(quicksum( (CO[i] * Z[i,t] + OP[i] * P[i,t] + R[i] * PR[i,t] + SD[i]*PS[i,t]) * (1+0.3)**(t-1)
                        for i in range(I) for t in range(T))
                + quicksum( C[i][j]*X[i,j,t]*(1+0.3)**(t-1) for i in range(I) for j in range(J) for t in range(T))
                + quicksum( S[j][k] * V[j,k,t] * (1+0.3)**(t-1) for j in range(J) for t in range(T) for k in range(K))
                + cost_W + cost_A)


m.addConstrs(( X.sum(i,'*',t) <= U[i] * P[i,t]
              for i in range(I) for t in range(T)), "cons7")

m.addConstrs(( IV[j,0] == X.sum('*',j,0) - V.sum(j,'*',0)
              for j in range(J)), "cons8")

m.addConstrs(( IV[j,t-1] + quicksum(X[i,j,t] for i in range(I)) - quicksum(V[j,k,t] for k in range(K))
              == IV[j,t] for t in range(1,T) for j in range(J)),"cons9")
                                    
m.addConstrs(( V.sum("*",k,t) == D[k] * (1+ DR[k] * (t-1)) for t in range(T) for k in range(K) ),"cons10")

m.addConstrs(( IV[j,t] + IV[j,t-1]/2 <= 4000 for t in range(1,T) for j in range(J) ),"cons11")

m.addConstrs(( IV[j,0] / 2 <= 4000 for j in range(J) ),"cons12")

m.addConstrs(( X.sum("*", j,t) <= 12000 for j in range(J) for t in range(T)),"cons13")

m.addConstrs(( V.sum(j, "*",t) <= 12000 for j in range(J) for t in range(T)),"cons14")

m.addConstrs(( W[i,t] * 4.7 <= 60000 for t in range(T) for i in range(I)),"cons18")

# Prints model in a form to look at it line by line (to ensure constraints and objective function written correctly)
m.write('Model.lp')

# Optimize Model
m.optimize()

#printSolution()
# Objective Function Value
print('\nTotal Costs: %g' % m.objVal)


# Solution (Variable Values):
print('SOLUTION:')
for i in range(I):
    for t in range(T):
        if Z[i,t].x > 0.99:
            # print('******************************')
            print('Plant %s is constructed in year %s' % ((i+1), (t+1)))
            for t2 in range(t,T):
                if PR[i,t2].x > 0.99:
                    print('Plant %s is reopened in year %s' % ((i+1), (t2+1)))
                if P[i,t2].x > 0.99:
                    print('Plant %s is operating in year %s' % ((i+1), (t2+1)))
                if PS[i,t2].x > 0.99:
                    print('Plant %s is closed in year %s' % ((i+1), (t2+1)))
                
            # for t1 in range(t,T):
            #     for j in range(J):
            #         if X[i,j,t1].x > 0:
            #             print('Transport %g units to warehouse %s in year %s' % ((X[i,j,t1].x), (j+1),(t1+1)))
                        
for i in range(I):
    for j in range(J):
        for t in range(T):
            if X[i,j,t].x > 0:
                print('Transport %g units to warehouse %s in year %s' % ((X[i,j,t].x), (j+1),(t+1)))
            
            
for j in range(J):
    print('Warehouse %s:' % (j+1))
    for k in range(K):
        for t in range(T):
            if V[j,k,t].x > 0:
                print('Transport %g units to Retail_Center %s in year %s' % ((V[j,k,t].x), (k+1),(t+1)))
