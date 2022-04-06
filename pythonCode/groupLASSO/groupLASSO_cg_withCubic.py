#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 14:12:02 2021
@author: cassiebuhler

Solves the following problem via Conjugate Gradient WITH Cubic regularization:

  minimize 1/2*|| Ax - b ||_2^2 + \lambda sum(norm(x_i))

The input p is a K-element vector giving the block sizes n_i, so that x_i
is in R^{n_i}.

The solution is returned in the vector x.

history is a dictionary that contains the objective values, l2 norm of gradients, time elapsed,
number of iterations, solution status (0 = solved, 1 = Search direction is not descent direction, 
2 = Iterations limit reached, 3 = search direction is undefined), and if cubic regularization was invoked (TRUE/FALSE)
"""
import numpy as np
import time as time
from scipy import optimize

def groupLASSO_cg_withCubic(A, b, lamb, p, alpha):
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
    pertcnt = 0 
    nrst = n
    restart = False
    
    inPowell = False
    dprat = 0.0
    ddprat = 0.0
    lam = 0
    k = 0
    
    #set to None because we will define later
    c0 = None
    dprat0 = None
    prat = None
    prat0 = None
    olamb = None   
    oolamb = None   
    alpha0 = None
    restart0 = None
    x0 = None
    dx0 = None
    c00 = None
    
    status = None
    objs = []
    gradNorms = []
    history = {}
    history['objective'] = []
    history['gradNorm'] = []

    while (k < MAX_ITER):
        k += 1
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
                if nrst == 0:
                    nrst = n
                inPowell = True
                if lam == 0.0:
                    if dprat0 != 0.0:
                        ddprat = ( abs(np.dot(c,c0)/cTc) - 2*prat + prat0 ) / ( (lam - olamb)*(oolamb-olamb) )
                    dprat0 = dprat
                    if lam >0:
                        dprat = ( abs(np.dot(c,c0)/cTc) - prat ) / ( lam - olamb )
                    if prat is not None:
                        prat0 = prat
                        
                    prat = abs(np.dot(c,c0)/cTc)
                    alpha = alpha0
                    restart = restart0
                    x = x0
                    dx = dx0
                    c = c0
                    c0 = c00
                    xTx = np.dot(x,x)
                    cTc = np.dot(c,c)
                    
                    if olamb is not None:
                        oolamb = olamb
                    olamb = lam
                    if lam == 0.0:
                        lam = prat/0.2
                    else:
                        if pertcnt <= 0.2:
                            lam = olamb - prat/dprat
                        else:
                            if ( dprat*dprat - 2*ddprat*prat > 0) and (abs(ddprat) > 1e-8 ):
                                lam = olamb + max( (-dprat + np.sqrt(dprat*dprat - 2*ddprat*prat))/ddprat, (-dprat - np.sqrt(dprat*dprat - 2*ddprat*prat))/ddprat )
                            else:
                                lam = 2*olamb
                
                        if lam < 1e-12:
                            lam = 2*olamb

                else:
                    lam = 2*lam
                pertcnt += 1
                k -= 1
                
            else:
                dx0 = dx
                if ( pertcnt > 0 ):
                    lam = lam/2.0
                else:
                    lam = 0.0
                dprat = 0.0
                dprat0 = 0.0
                prat0 = 0.0
                pertcnt = 0
            # If performing a restart, update the Beale restart vectors
            if ( nrst == n ):
                pt = alpha*dx
                yt = c - c0
                ytTyt = np.dot(yt,yt)
                cTyt = np.dot(pt,yt)
                cTct = np.dot(pt,pt)
                    
            p = alpha*dx
            y = c - c0
            pTc = np.dot(pt,c)
            yTc = np.dot(yt,c)
    
            u1 = -pTc/ytTyt
            u2 = 2*pTc/cTyt - yTc/ytTyt
            u3 = cTyt/ytTyt;

            if ( pertcnt == 0 ) :
                dx = -u3*c - u1*yt - u2*pt
            else:
                bracket = lam*lam + 2*lam*ytTyt/cTyt + ytTyt/cTct
                a = -cTyt/(lam*cTyt + ytTyt)
                b1 = (-lam*cTyt - 2*ytTyt)*ytTyt/(cTyt*cTct*bracket*(lam*cTyt + ytTyt))
                d = lam/(bracket*(lam*cTyt + ytTyt))
                e = ytTyt/(cTct*bracket*(lam*cTyt + ytTyt))

                bracket = lam*lam+lam/u3+lam*ytTyt/cTyt+cTyt/(u3*cTct)
                denom = u3*u3*cTct*lam*bracket + u3*cTct*bracket
                d = u3*u3*lam*(cTct/cTyt) / denom

                dx = a*c + b1*pTc*pt + d*yTc*yt + e*yTc*pt + e*pTc*yt
                
            if ( nrst != n ):
                if ( pertcnt == 0 ):
                    u1 = -np.dot(y,pt)/ytTyt
                    u2 = -np.dot(y,yt)/ytTyt + 2*np.dot(y,pt)/cTyt
                    u3 = np.dot(p,y)
                    temp = cTyt/ytTyt*y + u1*yt + u2*pt
                    u4 = np.dot(temp, y)
                    
                    u1 = -np.dot(p,c)/u3
                    u2 = (u4/u3 + 1)*np.dot(p,c)/u3 - np.dot(c,temp)/u3
                    dx = dx - u1*temp - u2*p
                else:
                    ptTy = np.dot(pt,y)
                    ytTy = np.dot(yt,y)
                    temp1 = -(a*y + b1*ptTy*pt + d*ytTy*yt + e*ytTy*pt + e*ptTy*yt)

                    a2 = 1.0/(lam*u3+1.0)
                    b2 = lam*b1
                    d2 = lam*d
                    e2 = lam*e
                    ptTp = np.dot(pt,p)
                    ytTp = np.dot(yt,p)
                    temp2 = a2*p + b2*ptTp*pt + d2*ytTp*yt + e2*ytTp*pt + e2*ptTp*yt

                    u10 = np.dot(temp1, c)
                    u11 = np.dot(temp2, c)
                    u12 = np.dot(temp1, y)
                    u13 = np.dot(temp2, y)
                    u14 = np.dot(temp2, p)
                    u15 = np.dot(p, y)
                    denom = lam*u14*(u15 + u12) + u13*u13
                    if denom < np.finfo(float).eps:
                        status = 3 
                        print('Search direction is undefined.')
                        break
                    dx = dx + lam*u14*u10*temp1/denom - (u15 + u12)*u11*temp2/denom + u13*u10*temp2/denom + u13*u11*temp1/denom

		# Check that the search direction is a descent direction
        dxTc = np.dot(dx, c)
        if ( dxTc > 0 ):
            status = 1 
            print('Search direction is not a descent direction.')
            break
		# Save the current point
        x0 = x
        if c0 is not None:
            c00 = c0
        c0 = c
        alpha0 = alpha
        restart0 = restart

        if ( restart == 0 ):
            restart = 1
        else:
            if ( pertcnt == 0 ):
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
        print('n = %d, Iters = %d, invokedCubic = %s\n'% (n, k, inPowell == 1))
    
    if k == MAX_ITER:
        status = 2 
        print('Iterations limit reached.')
        
    z = x
    history['objective'] = objs
    history['gradNorm'] = gradNorms
    history['status'] = status
    history['time'] = elapsedTime
    history['iters'] = k
    history['invokedCubic'] = inPowell == 1

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

