#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 15:20:20 2021
Last modified on March 11, 2022
@author: cassiebuhler
"""
from scipy.sparse import spdiags
import numpy as np
from huber_cubic import huber_cubic
import scipy.sparse as sparse

# Generate problem data
np.random.seed(22)


m = 1000 #number of examples
n = 500 # number of features

for rr in range(100):
    print("Problem %d "% rr)
    x0 = np.random.randn(n,)
    A = np.random.randn(m,n)
    A = A*spdiags(1./np.sqrt(sum(A**2)).T,0,n,n)#normalize columns
    b = np.dot(A,x0) + np.sqrt(0.01)*np.random.randn(m,)
    noise = 10*sparse.random(m,1,density = 200/m) # add sparse, large noise
    b = b + np.squeeze(noise.toarray())
    lambda_max = np.linalg.norm( A.T@b, np.inf)
    lamb = 0.001*lambda_max;
    # Solve problem
    x, history = huber_cubic(A, b, lamb, 1.0, 1.0)


