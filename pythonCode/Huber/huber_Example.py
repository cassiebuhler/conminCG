#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 15:20:20 2021
@author: cassiebuhler

% Huber Function - Example
% Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/ and converted to Python
% This code generates random problems that are solved with Conjugate
% Gradient Method (CGM) with and without Cubic Regularization.
"""
from scipy.sparse import spdiags
import numpy as np
import scipy.sparse as sparse

from huber_cg_withCubic import huber_cg_withCubic
from huber_cg_noCubic import huber_cg_noCubic
from getPerformanceProfiles import getPerformanceProfiles


model = 2 # 0 = CG without Cubic Reg, 1 = CG with cubic Reg, 2 = BOTH MODELS
perfProf = True # outputs a performance profile comparing two methods. 

names = ["CG WITHOUT CUBIC REGULARIZATION","CG WITH CUBIC REGULARIZATION"]


# Generate problem data
np.random.seed(0)


m = 5000 #number of examples
n = 2000 # number of features
numProbs = 100 #number of problems
alpha = 1.0 #over-relaxation parameter 


time = np.zeros((numProbs,2))
iters = np.zeros((numProbs,2))
status = np.zeros((numProbs,2))

if model != 2:
    print("-"*40)
    print(names[model])
    print("-"*40)


for rr in range(numProbs):
    print("Problem %d "% rr)
    x0 = np.random.randn(n,)
    A = np.random.randn(m,n)
    A = A*spdiags(1./np.sqrt(sum(A**2)).T,0,n,n)#normalize columns
    b = np.dot(A,x0) + np.sqrt(0.01)*np.random.randn(m,)
    noise = 10*sparse.random(m,1,density = 200/m) # add sparse, large noise
    b = b + np.squeeze(noise.toarray())
    # Solve problem
    if model == 0:
        x, history = huber_cg_noCubic(A, b, alpha)
    elif model == 1:
        x, history = huber_cg_withCubic(A, b, alpha)
    else:
        print("-"*40)
        print('--'+names[0]+'--')
        x1, history1 = huber_cg_noCubic(A, b, alpha)
        print('--'+names[1]+'--')
        x2, history2 = huber_cg_withCubic(A, b, alpha)
        time[rr,:] = [history1['time'],history2['time']] #saving time and iterations for performance profiles
        iters[rr,:] = [history1['iters'],history2['iters']]
        status[rr,:] = [history1['status'],history2['status']]
    print("-"*40)

if model == 2 and perfProf == True:
    getPerformanceProfiles(iters,time,status) #output performance profiles

