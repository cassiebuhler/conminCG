#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on March 28, 2022
@author: cassiebuhler
  
Solves the following problem via Conjugate Gradient WITHOUT Cubic regularization:

  minimize 1/2*|| Ax - b ||_2^2 + \lambda sum(norm(x_i))

The input p is a K-element vector giving the block sizes n_i, so that x_i
is in R^{n_i}.

The solution is returned in the vector x.

history is a dictionary that contains the objective values, l2 norm of gradients, time elapsed,
number of iterations, solution status (0 = solved, 1 = Search direction is not descent direction, 
2 = Iterations limit reached, 3 = search direction is undefined), and if a powell restart was needed (TRUE/FALSE)
"""
import numpy as np
import time as time
from scipy import optimize

def groupLASSO_cg_noCubic(A, b, lamb, p, alpha):
    #initialize
    t_start = time.time()
    QUIET    = 0
    MAX_ITER = 1000
    ABSTOL   = 1e-4
    RELTOL   = 1e-2

    
    #data preprocessing
    [m,n] = A.shape

    # check that sum(p) = total number of elements in x
    if sum(p)!= n:
        raise Exception('invalid partition')
    
    #cumulative partition
    cum_part = np.cumsum(p)
    # CG solver
    x = 0.1*np.ones(n)
    c = grad(A,b,lamb,x,cum_part)
    nrst = n
    restart = False
    inPowell = False
    
    #set to None because we will define later
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
        if (np.sqrt(cTc) <= np.sqrt(n)*ABSTOL + RELTOL*np.sqrt(xTx) ):
            status = 0 #solved
            break
        
        # Compute step direction
        if restart == False:
            dx = -c
        else:
            # test to see if powell restart criterion holds
            if (nrst != n) and (restart >1) and (abs(np.dot(c,c0)/cTc) > 0.2):
                nrst = n
                inPowell = True
                
            # If performing a restart, update the Beale restart vectors
            if ( nrst == n ):
                pt = alpha*dx
                yt = c - c0
                ytTyt = np.dot(yt,yt)
                cTyt = np.dot(pt,yt)
                    
            p = alpha*dx
            y = c - c0
            u1 = -np.dot(pt,c)/ytTyt;
            u2 = 2*np.dot(pt,c)/cTyt - np.dot(yt,c)/ytTyt;
            u3 = cTyt/ytTyt;
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
            if ( nrst == n ) :
                nrst = 0
            nrst +=  1
            restart = 2
  


        afind = lambda a: objective(A, b, lamb, cum_part,  x + a*dx, x)
        alpha = optimize.fminbound(afind, 0, 10)
  

        
        # Take the step and update function value and gradient
        x = x0 + alpha*dx
        c = grad(A, b, lamb, x, cum_part)
        objs.append(objective(A, b, lamb, cum_part,  x, x))
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

    return z,history



	
def objective(A, b, lamb, cum_part, x, z):
    obj = 0
    start_ind = 0
    for i in range(len(cum_part)):
        sel = range(start_ind,cum_part[i])
        obj = obj + np.linalg.norm(z[sel])
        start_ind = cum_part[i]
    p = ( 1/2*sum((A@x - b)**2) + lamb*obj )
    return p


def grad(A, b, lamb, x, cum_part):
    start_ind = 0
    c =np.dot(A.T,(np.dot(A,x)-b))
    for i in range(len(cum_part)):
        sel = range(start_ind,cum_part[i])
        c[sel] = c[sel] +lamb*x[sel]/(np.linalg.norm(x[sel]))
        start_ind = cum_part[i] 
    return c

