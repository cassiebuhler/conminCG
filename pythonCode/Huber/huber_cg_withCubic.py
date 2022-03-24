#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 14:53:11 2021
Last modified on March 11, 2022
@author: cassiebuhler
"""
import numpy as np
import time as time
from scipy import optimize
def huber_cg_withCubic(A, b, alpha):

    
    #initialize
    t_start = time.time()

    QUIET    = 0
    MAX_ITER = 200
    ABSTOL   = 1e-4
    RELTOL   = 1e-2
    
    # Data preprocessing
    m,n = A.shape
    
    # CG solver
    x = 10*np.ones(n)
    c = grad(A, b, x, m, 1.0)
    
    pertcnt = 0
    nrst = n
    restart = 0
    
    inPowell = False
    dprat = 0.0
    ddprat = 0.0
    lam = 0
    k = 0
    
    # initialize to None
    c0 = None
    dprat0 = None
    prat = None
    prat0 = None
    olambda = None
    oolambda = None
    alpha0 = None
    restart0 = None
    x0 = None
    dx0 = None
    c00 = None
    
    

    while (k < MAX_ITER):
        k += 1
        xTx = np.dot(x,x)
        cTc = np.dot(c,c)

		# Check for convergence
        if ( np.sqrt(cTc) <= np.sqrt(n)*ABSTOL + RELTOL*np.sqrt(xTx) ):
            break

		# Compute step direction
        if ( restart == 0 ):
            dx = -c
            # print("%4d :\t%14.6e\t %14.6e\t | %d \n", k, f, sqrt(cTc/max(1.0, xTx)), pertcnt);
        else:
            if ( abs(np.dot(c,c0)/cTc) > 0.2) and (restart > 1) and (nrst != n ):
			# Test to see if the Powell restart criterion holds
                if (nrst == 0):
                    nrst = n 
                inPowell = True;
                if ( lam == 0.0) or (lam > 1e6):
                    if ( dprat0 != 0.0 ):
                        ddprat = ( abs(np.dot(c,c0)/cTc) - 2*prat + prat0 ) / ( (lam - olambda)*(oolambda-olambda) )
                    dprat0 = dprat
                    if ( lam > 0 ):
                        dprat = ( abs(np.dot(c,c0)/cTc) - prat ) / ( lam - olambda )
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
                    if olambda is not None:
                        oolambda = olambda
                    olambda = lam
                    if ( lam == 0.0 ):
                        lam = prat/0.2
                    else:
                        if ( pertcnt <= 2 ):
                            lam = olambda - prat/dprat
                        else:
                            if ( dprat*dprat - 2*ddprat*prat > 0) and (abs(ddprat) > 1e-8 ):
                                lam = olambda + max( (-dprat + np.sqrt(dprat*dprat - 2*ddprat*prat))/ddprat, (-dprat - np.sqrt(dprat*dprat - 2*ddprat*prat))/ddprat )
                            else:
                                lam = 2*olambda
                        if ( lam < 1e-12 ) :
                            lam = 2*olambda

                else:
                    lam = 2*lam
                pertcnt +=   1
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
				# print("%4d :\t%14.6e\t %14.6e\t | %d \n", k, f, sqrt(cTc/max(1.0, xTx)), pertcnt);
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
            u3 = cTyt/ytTyt
            
            if ( pertcnt == 0 ):
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
                    u15 = np.dot(p, y);	
                    
                    denom = lam*u14*(u15 + u12) + u13*u13
                    
                    dx = dx + lam*u14*u10*temp1/denom - (u15 + u12)*u11*temp2/denom + u13*u10*temp2/denom + u13*u11*temp1/denom;


		# Check that the search direction is a descent direction
        dxTc = np.dot(dx, c)
        if ( dxTc > 0 ):
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
            restart = 2
            if ( pertcnt == 0 ):
                if ( nrst == n ):
                    nrst = 0
                nrst +=  1
                restart = 2


        
        afind = lambda a: objective(A, b, x + a*dx)
        alpha = optimize.fminbound(afind, 0, 10)
        
        # Take the step and update function value and gradient
        x = x0 + alpha*dx
        c = grad(A, b, x, m, 1.0)
        



    
    if not QUIET:
        print('time elapsed: ' + str(time.time() - t_start))
    
    if k == MAX_ITER:
        print('MAX ITERATION REACHED')
        
    z = x
    history = x
    print('Iters = %d, invokedCubic = %s\n'% (k, inPowell == 1))
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

