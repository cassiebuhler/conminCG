#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 15:20:20 2021
Last modified on March 11, 2022
@author: cassiebuhler
"""
from scipy.sparse import spdiags
import numpy as np
from huber_cg_withCubic import huber_cg_withCubic
import scipy.sparse as sparse

# Generate problem data
np.random.seed(0)


m = 1000 #number of examples
n = 500 # number of features
alpha = 1.0 #over-relaxation parameter 
for rr in range(100):
    print("Problem %d "% rr)
    x0 = np.random.randn(n,)
    A = np.random.randn(m,n)
    A = A*spdiags(1./np.sqrt(sum(A**2)).T,0,n,n)#normalize columns
    b = np.dot(A,x0) + np.sqrt(0.01)*np.random.randn(m,)
    noise = 10*sparse.random(m,1,density = 200/m) # add sparse, large noise
    b = b + np.squeeze(noise.toarray())
    # Solve problem
    x, history = huber_cg_withCubic(A, b, alpha)


