#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 15:53:05 2021
Last modified on March 11, 2022

@author: cassiebuhler
"""
import random
import math
from scipy.sparse import spdiags
import numpy as np
from groupLASSO_cg_withCubic import groupLASSO_cg_withCubic

random.seed(1)
np.random.seed(1)

# Group lasso example with random data
# Generate problem data
m = 1000      #amount of data
K = 3       # number of blocks
alpha = 1.0 #over-relaxation parameter
lambdas = np.zeros(K)
for rr in range(100):
    print("Problem %d"% rr)
    partition = np.random.randint(low = 0, high = 5000,size = K)
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
    x, history = groupLASSO_cg_withCubic(A, b, lamb, partition, alpha)
