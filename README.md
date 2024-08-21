# Conmin-CG
Conjugate gradient minimization methods (CGM) and their accelerated variants are widely used in unconstrained optimization. We propose the hybrid cubic regularization of conjugate gradient minimization methods (CGM). This repository contains implementations in C, Python, and MATLAB. 

- Authors: Cassidy K. Buhler, Hande Y. Benson, and David F. Shanno
- Code is adapted from [Conmin](https://dl.acm.org/doi/pdf/10.1145/355921.355933) by David F. Shanno.

## Overview

### The Issue with CGM 

CGM requires "restarting" the step direction in the following scenarios: 1) every $n$ iterations (Beale Restart) and 2) when conjugacy is lost (Powell Restart). On problems that frequently lose conjugacy, the latter becomes an issue as continuously resetting the step direction reduces local convergence from superlinear to linear. 

### Improving CGM Step Quality

In our paper, ["Regularized Step Directions in Nonlinear Conjugate Gradient Methods"](https://arxiv.org/abs/2110.06308), we proposed **CGM with hybrid cubic regularization** to improve the step quality of CGM. 
- Using Shanno's reformulation of CGM as a memoryless BFGS method, we derive new formulas for the regularized step direction.
- CGM with hybrid cubic regularization uses a regularized step in lieu of a Powell Restart. 
  
#### Results: 
- The regularized step direction uses the same order of computational burden per iteration as its non-regularized version.
- CGM with hybrid cubic regularization exhibits global convergence with fewer assumptions.
- In numerical experiments (refer to C implementation), the new step directions are shown to require fewer iteration counts, improve runtime, and reduce the need to reset the step direction. 

#### Takeaway:
- CGM with hybrid cubic regularization exhibits the same memoryless and matrix-free properties as Shanno's CGM reformulation while outperforming it in both iterations and runtime.


### Extensions to Machine Learning 
- We are extending this work to include machine learning problems and have implemented Conmin-CG in Python and MATLAB.
- These implementations solve large instances of machine learning problems (Huber and Group LASSO) and compare it to other state-of-the-art solvers (ADMM). 


## Contents

#### C Code:
- This code solves unconstrained problems from the [CUTEst test set](https://github.com/ralna/CUTEst) where the models were converted from SIF to an [AMPL formulation](https://vanderbei.princeton.edu/ampl/nlmodels/cute/index.html).
- We solve each instance with the following solvers for comparison:
  1. CGM with hybrid cubic regularization.
  2. CGM without hybrid cubic regularization (in our paper, we refer to this method as “CGM with Powell Restarts”)


#### Matlab & Python Code:
- We solve large instances of machine learning problems (Huber and Group LASSO) with the following solvers:
  1. CGM with hybrid cubic regularization
  2. CGM without hybrid cubic regularization (aka CGM with Powell Restarts)
  3. Stephen Boyd's [ADMM](https://web.stanford.edu/~boyd/papers/admm/).

## Contact
Email: cassidy.buhler@gmail.com
