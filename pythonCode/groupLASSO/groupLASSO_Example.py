#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 15:53:05 2021

@author: cassiebuhler
% Group LASSO - Example
% Code is adapted from Group LASSO code by Boyd https://web.stanford.edu/~boyd/papers/admm/ and converted to Python
% This code generates random problems that are solved with Conjugate
% Gradient Method (CGM) with and without Cubic Regularization.
"""
import random
import math
from scipy.sparse import spdiags
import numpy as np
from groupLASSO_cg_hybridCubic import groupLASSO_cg_hybridCubic
from groupLASSO_cg_powellRestarts import groupLASSO_cg_powellRestarts

random.seed(1)
np.random.seed(1)

model = 2 # 0 = CG without Cubic Reg, 1 = CG with cubic Reg, 2 = BOTH MODELS
names = ["CG WITH POWELL RESTARTS","CG WITH HYBRID CUBIC REGULARIZATION"]


# Group lasso example with random data
# Generate problem data
m = 1500      #amount of data
K = 4       # number of blocks
numProbs = 100 #number of problems
upper = 1000 #upperbound for partition 
alpha = 1.0 #over-relaxation parameter
lambdas = np.zeros(K)

time = np.zeros((numProbs,2))
iters = np.zeros((numProbs,2))
status = np.zeros((numProbs,2))

if model != 2: #print out solver name before solving problem(s)
    print("-"*40)
    print(names[model])
    print("-"*40)
    
for rr in range(numProbs):
    print("Problem %d"% rr)
    partition = np.random.randint(low = 0, high = upper,size = K)
    n = sum(partition)# number of features
    p = 100/n#  sparsity density
    
    # generate block sparse solution vector
    x = np.zeros(n)
    start_ind = 0
    cum_part = np.cumsum(partition)
    for i in range(K):
        sel = np.arange(start_ind,cum_part[i])
        x[sel] = 0
        if( random.random() < p):
            # fill nonzeros
            sel = np.arange(start_ind,cum_part[i])
            x[sel] = np.random.randn(partition[i],)
        start_ind = cum_part[i]
        
    # generate random data matrix
    A = np.random.randn(m,n)
    # normalize columns of A
    A = A*spdiags(1./np.sqrt(sum(A**2)),0,n,n)
    # generate measurement b with noise
    b = np.dot(A,x) + np.dot(math.sqrt(0.001),np.random.randn(m,))
    # lambda max
    start_ind = 0
    for i in range (K):
        sel = np.arange(start_ind,cum_part[i])
        lambdas[i] = np.linalg.norm(A[:,sel].T*b)
        start_ind = cum_part[i] 
    lambda_max = max(lambdas);
    
    # regularization parameter
    lamb = 0.01*lambda_max
    xtrue = x   # save solution
    #Solve problem
    if model == 0:
        x1, history1 = groupLASSO_cg_powellRestarts(A, b, lamb, partition, alpha)
    elif model == 1:
        x2, history2 = groupLASSO_cg_hybridCubic(A, b, lamb, partition, alpha)
    else:
        print("-"*40)
        print('--'+names[0]+'--')
        x1, history1 = groupLASSO_cg_powellRestarts(A, b, lamb, partition, alpha)
        print('--'+names[1]+'--')
        x2, history2 = groupLASSO_cg_hybridCubic(A, b, lamb, partition, alpha)
        time[rr,:] = [history1['time'],history2['time']] #saving time and iterations for performance profiles
        iters[rr,:] = [history1['iters'],history2['iters']]
        status[rr,:] = [history1['status'],history2['status']]
    print("-"*40)

