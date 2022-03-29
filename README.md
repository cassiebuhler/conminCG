# ConminCG
## Constrained Minimization: Cubic Regularization for Conjugate Gradient Methods
Conjugate gradient minimization methods (CGM) and their accelerated variants are
widely used in machine learning applications. We focus on the use of cubic regularization to
improve the CGM direction independent of the steplength (learning rate) computation. Using
Shanno’s reformulation of CGM as a memoryless BFGS method, we derive new formulas for the
regularized step direction, which can be evaluated without additional computational effort. The
new step directions are shown to improve iteration counts and runtimes and reduce the need to
restart the CGM.


Paper: https://arxiv.org/abs/2110.06308

## Contents
We have implemented cubic regularization to CGM in C, Python, and Matlab. 

C code:
- The C implementation is connected to AMPL and includes CGM with and without cubic regularization. 

Matlab code:
- CGM with and without cubic regularization
- Stephen Boyd's ADMM code (https://web.stanford.edu/~boyd/papers/admm/) for comparison. 
- Two machine learning examples: Huber and Group LASSO.

Python code:
- CGM with and without cubic regularization
- Two machine learning examples: Huber and Group LASSO.

## Contact
Email: cb3452@drexel.edu 
