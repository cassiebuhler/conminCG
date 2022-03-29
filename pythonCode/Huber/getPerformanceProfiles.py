#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 12:18:58 2022

@author: cassiebuhler
"""
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt


def getRatio(mat, isIter):
    best = mat.min(axis=1)
    ratios = mat/best[:,None]
    maxIter = 1000
    if isIter:
        ratios[mat == maxIter] = np.nan; #if it doesn't solve, set ratio to nan 
    return ratios

def performanceProf(ratios):
  ratio = sorted(ratios)
  counts = Counter(ratio)
  counts = sorted(counts.items(), key=lambda x: x[0])
  probs = []
  taus = []
  prob = 0
  for count in counts:
    prob += count[1]/len(ratio)
    probs.append(prob)
    taus.append(count[0])
  return np.vstack((taus,probs))


def getPerformanceProfiles(iters,time):
    ratios_time = getRatio(time,0) 
    ratios_iters = getRatio(iters,1)

    ppIters_cgNoCubic = performanceProf(ratios_iters[:,0])
    ppIters_cgCubic = performanceProf(ratios_iters[:,1])

    ppTime_cgNoCubic = performanceProf(ratios_time[:,0])
    ppTime_cgCubic = performanceProf(ratios_time[:,1])
    plt.figure(0)
    plt.title("Performance Profile: Iterations") 
    plt.xlabel(r"Tau")
    plt.ylabel("Probability")
    plt.plot(ppIters_cgNoCubic[0,:],ppIters_cgNoCubic[1,:], label = 'No Cubic', marker='o', linestyle='dashed', markersize = 8,linewidth = 2)
    plt.plot(ppIters_cgCubic[0,:],ppIters_cgCubic[1,:], label = 'With Cubic', marker='o', linestyle='dashed', markersize = 8,linewidth = 2)
    plt.legend()
    plt.ylim([0,1.1])

    
    plt.figure(1)
    plt.title("Performance Profile: Time") 
    plt.xlabel(r"Tau")
    plt.ylabel("Probability")
    plt.plot(ppTime_cgNoCubic[0,:],ppTime_cgNoCubic[1,:], label = 'No Cubic', marker='o', linestyle='dashed', markersize = 8,linewidth = 2)
    plt.plot(ppTime_cgCubic[0,:],ppTime_cgCubic[1,:], label = 'With Cubic', marker='o', linestyle='dashed', markersize = 8,linewidth = 2)
    plt.legend()
    plt.ylim([0,1.1])
    return

