# ConminCG
Conjugate gradient minimization methods (CGM) and their accelerated variants are widely used. We focus on the use of cubic regularization to improve the CGM direction independent of the steplength computation. In this paper, we propose the Hybrid Cubic Regularization of CGM, where regularized steps are used selectively. Using Shanno's reformulation of CGM as a memoryless BFGS method, we derive new formulas for the regularized step direction. We show that the regularized step direction uses the same order of computational burden per iteration as its non-regularized version. Moreover, the Hybrid Cubic Regularization of CGM exhibits global convergence with fewer assumptions. In numerical experiments, the new step directions are shown to require fewer iteration counts, improve runtime, and reduce the need to reset the step direction. Overall, the Hybrid Cubic Regularization of CGM exhibits the same memoryless and matrix-free properties, while outperforming CGM as a memoryless BFGS method in iterations and runtime.

Paper: https://arxiv.org/abs/2110.06308

## Contents
We have implemented hybrid cubic regularization to CGM in C, Python, and Matlab. 

C code:
- The C implementation is connected to AMPL and includes CGM with and without hybrid cubic regularization. 

Matlab code:
- CGM with and without cubic regularization
- Stephen Boyd's ADMM code (https://web.stanford.edu/~boyd/papers/admm/) for comparison. 
- Two machine learning examples: Huber and Group LASSO.

Python code:
- CGM with and without cubic regularization
- Two machine learning examples: Huber and Group LASSO.

## Contact
Email: cassidy.buhler@gmail.com
