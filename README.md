# conminCG
We propose the hybrid cubic regularization of conjugate gradient minimization methods (CGM). This repository contains implementations in C, Python, and MATLAB. 

- Authors: Cassidy K. Buhler, Hande Y. Benson, and David F. Shanno
- Code is adapted from [Conmin]((https://dl.acm.org/doi/pdf/10.1145/355921.355933)) by David F. Shanno
- The C implementation is described in our paper, *Regularized Step Directions in Nonlinear Conjugate Gradient Methods*. [Preprint](https://arxiv.org/abs/2110.06308). 

## Contents

C Code:
- This code solves unconstrained problems from the [CUTEst test set](https://github.com/ralna/CUTEst) where the models were converted from SIF to an [AMPL formulation](https://vanderbei.princeton.edu/ampl/nlmodels/cute/index.html).
- We solve each instance with the following solvers for comparison:
  1. CGM with hybrid cubic regularization.
  2. CGM without hybrid cubic regularization (in our paper, we refer to this method as “CGM with Powell Restarts”)



Matlab & Python Code:
- We solve large instances of machine learning problems (Huber and Group LASSO) with the following solvers:
  1. CGM with hybrid cubic regularization
  2. CGM without hybrid cubic regularization (aka CGM with Powell Restarts)
  3. And for comparison, Stephen Boyd's [ADMM](https://web.stanford.edu/~boyd/papers/admm/).

## Abstract 

Conjugate gradient minimization methods (CGM) and their accelerated variants are widely used. We focus on the use of cubic regularization to improve the CGM direction independent of the steplength computation. In this paper, we propose the *Hybrid Cubic Regularization of CGM*, where regularized steps are used selectively. Using Shanno's reformulation of CGM as a memoryless BFGS method, we derive new formulas for the regularized step direction. We show that the regularized step direction uses the same order of computational burden per iteration as its non-regularized version. Moreover, the *Hybrid Cubic Regularization of CGM* exhibits global convergence with fewer assumptions. In numerical experiments, the new step directions are shown to require fewer iteration counts, improve runtime, and reduce the need to reset the step direction. Overall, the *Hybrid Cubic Regularization of CGM* exhibits the same memoryless and matrix-free properties, while outperforming CGM as a memoryless BFGS method in iterations and runtime.

Paper: https://arxiv.org/abs/2110.06308

## Contact
Email: cassidy.buhler@gmail.com
