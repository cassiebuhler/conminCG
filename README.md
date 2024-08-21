# conminCG - C Implementation

- This repository contains the code from our paper, *Regularized Step Directions in Nonlinear Conjugate Gradient* Methods. [Preprint](https://arxiv.org/abs/2110.06308)
- Authors: Cassidy K. Buhler, Hande Y. Benson, and David F. Shanno
- Code is adapted from Conmin by David F. Shanno (https://dl.acm.org/doi/pdf/10.1145/355921.355933)


## Abstract 

Conjugate gradient minimization methods (CGM) and their accelerated variants are widely used. We focus on the use of cubic regularization to improve the CGM direction independent of the steplength computation. In this paper, we propose the Hybrid Cubic Regularization of CGM, where regularized steps are used selectively. Using Shanno's reformulation of CGM as a memoryless BFGS method, we derive new formulas for the regularized step direction. We show that the regularized step direction uses the same order of computational burden per iteration as its non-regularized version. Moreover, the Hybrid Cubic Regularization of CGM exhibits global convergence with fewer assumptions. In numerical experiments, the new step directions are shown to require fewer iteration counts, improve runtime, and reduce the need to reset the step direction. Overall, the Hybrid Cubic Regularization of CGM exhibits the same memoryless and matrix-free properties, while outperforming CGM as a memoryless BFGS method in iterations and runtime.


## Description
We provide following solvers: 

- CGM with hybrid cubic regularization. 
- CGM without hybrid cubic regularization. In the paper, this non-regularized method is referred to as “CGM with Powell Restarts”. 

This code solves unconstrained nonlinear optimization problems with and without hybrid cubic regularization. The problems must be formulated in AMPL, and the user needs to have AMPL installed.

In our paper, we used the AMPL models available [here](https://vanderbei.princeton.edu/ampl/nlmodels/cute/index.html). These have been converted to AMPL from SIF and originated from the [CUTEst repository](https://github.com/ralna/CUTEst).


## Installation

The software must be compiled prior to use.  Makefiles are provided for MacOS/Linux/Unix (makefile) and Windows (makefile.vc), and the user should edit the makefile for their platform with information on the locations of their AMPL and C libraries.  Once the edits are completed, simply run
```
make
```
at the command prompt for MacOS/Linux/Unix, or
```
make -f makefile.vc
```
at a developer command prompt (Visual Studio) for Windows.  This will create the executable conmin (or conmin.exe for Windows).


## How to Use

In order to use the solver, it needs to be specified as the solver of choice for AMPL.  This can be done inside the AMPL model or at the AMPL prompt with the line
```
option solver conmin;
```

There are several solver options which can be specified using conmin_options inside the AMPL model or at the AMPL prompt:
```
option conmin_options "verbose=2";
```
will set the output level to 2, which means one line per iteration.  The verbose option can be set to any integer from 0 to 5.  

Other options that can be set by the user are as follows:
- `timing`: Toggle between 0 and 1 to provide runtime information
- `iterlim`: Maximum number of iterations
- `inftol`: Convergence tolerance

Once the solver and its options are chosen, the user can solve the AMPL model at the command prompt as 
```
ampl modelname.mod
```
or at the AMPL prompt as
```
model modelname.mod; solve;
```

## Contact
Email: cassidy.buhler@gmail.com





