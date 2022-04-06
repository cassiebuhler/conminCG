#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 28, 2022
@author: cassiebuhler


Solves the following problem via Conjugate Gradient WITHOUT Cubic regularization:

minimize h( Az-b ) where h is the Huber loss funciton

The solution is returned in the vector x.

history is a dictionary that contains the objective values, l2 norm of gradients, time elapsed,
number of iterations, solution status (0 = solved, 1 = Search direction is not descent direction, 
2 = Iterations limit reached, 3 = search direction is undefined), and if a powell restart was needed (TRUE/FALSE)
"""
import numpy as np
import time as time
from scipy import optimize
def huber_cg_noCubic(A, b, alpha):

    
    #initialize
    t_start = time.time()

    QUIET    = 0
    MAX_ITER = 1000
    ABSTOL   = 1e-4
    RELTOL   = 1e-2
    
    # Data preprocessing
    m,n = A.shape
    
    # CG solver
    x = 10*np.ones(n)
    c = grad(A, b, x, m, 1.0)
    
    nrst = n
    restart = 0
    inPowell = False
    
    # initialize to None
    c0 = None
    x0 = None
    
    status = None
    objs = []
    gradNorms = []
    history = {}
    history['objective'] = []
    history['gradNorm'] = []
    

    for k in range(MAX_ITER):
        xTx = np.dot(x,x)
        cTc = np.dot(c,c)

		# Check for convergence
        if ( np.sqrt(cTc) <= np.sqrt(n)*ABSTOL + RELTOL*np.sqrt(xTx) ):
            status = 0 #solved
            break

		# Compute step direction
        if ( restart == 0 ):
            dx = -c
        else:
            if ( abs(np.dot(c,c0)/cTc) > 0.2) and (nrst != n ):
			# Test to see if the Powell restart criterion holds
                nrst = n 
                inPowell = True;

			# If performing a restart, update the Beale restart vectors
            if ( nrst == n ):
                pt = alpha*dx
                yt = c - c0
                ytTyt = np.dot(yt,yt)
                cTyt = np.dot(pt,yt)
            
            p = alpha*dx
            y = c - c0
            pTc = np.dot(pt,c)
            yTc = np.dot(yt,c)

            u1 = -pTc/ytTyt
            u2 = 2*pTc/cTyt - yTc/ytTyt
            u3 = cTyt/ytTyt
            dx = -u3*c - u1*yt - u2*pt
            
            if ( nrst != n ):
                u1 = -np.dot(y,pt)/ytTyt
                u2 = -np.dot(y,yt)/ytTyt + 2*np.dot(y,pt)/cTyt
                u3 = np.dot(p,y)
                temp = cTyt/ytTyt*y + u1*yt + u2*pt
                u4 = np.dot(temp, y)
            
                u1 = -np.dot(p,c)/u3
                u2 = (u4/u3 + 1)*np.dot(p,c)/u3 - np.dot(c,temp)/u3
                dx = dx - u1*temp - u2*p

		# Check that the search direction is a descent direction
        dxTc = np.dot(dx, c)
        if ( dxTc > 0 ):
            status = 1 
            print('Search direction is not a descent direction.')
            break

		# Save the current point
        x0 = x
        c0 = c


        if ( restart == 0 ):
            restart = 1
        else:

            if ( nrst == n ):
                nrst = 0
            nrst +=  1
            restart = 2


        
        afind = lambda a: objective(A, b, x + a*dx)
        alpha = optimize.fminbound(afind, 0, 10)
        
        # Take the step and update function value and gradient
        x = x0 + alpha*dx
        c = grad(A, b, x, m, 1.0)
        objs.append(objective(A, b, x))
        gradNorms.append(np.linalg.norm(c))

    if not QUIET:
        elapsedTime = time.time() - t_start
        print('time elapsed: ' + str(elapsedTime))
        print('n = %d, Iters = %d, powellRestart = %s\n'% (n, k, inPowell == 1))

        
    if k == MAX_ITER:
        status = 2 
        print('Iterations limit reached.')
        
    z = x
    history['objective'] = objs
    history['gradNorm'] = gradNorms
    history['status'] = status
    history['time'] = elapsedTime
    history['iters'] = k
    history['powellRestart'] = inPowell == 1
    return z, history

def objective(A, b,x):
    p =  1/2*huber(A@x - b, 1.0)
    return p 

def huber(x,gamma):
    x1 = -x - 0.5
    x2 = 0.5*x*x
    x3 = x - 0.5
    d = sum(x1*(x <= -gamma) + x2*(-gamma < x)*(x < gamma ) + x3*(x >= gamma))
    return d

def huber_grad(x,n,gamma):
    c = -np.ones(n)*(x <= -gamma) + x*( -gamma < x)*(x < gamma ) + np.ones(n)*(x >= gamma)
    return c


def grad(A, b, x, n, gamma):
    c = A.T@huber_grad((A@x-b),n,gamma)
    return c

