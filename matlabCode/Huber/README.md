# ABOUT
## Matlab Code - Huber Examples

- Authors: Cassidy Buhler and Hande Benson
- Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/

## WHAT THE CODE DOES 

This code minimizes the Huber loss function with Conjugate Gradient Method (CGM) with and without Cubic Regularization.
We also include Boyd's ADMM for comparison. The data is randomly generated. 

## HOW TO USE IT

Variables can be modified in *huber_Example.m* file. This is the driver code.

Depending on the solver you pick, it will call *huber_cg_hybridCubic.m*, *huber_cg_powellRestarts.m*, and/or *huber_admm.m*. The latter is the code downloaded from Boyd's website linked above. 

To reproduce results of our paper, *huber_Example_moreCases.m* runs 4 cases of varying (m,n) sizes.   

1. **MODIFY SOLVER**

The model option be set to integers from 1 to 4 to control the optimizer used:

- `model = 1` CG without cubic regularization 
- `model = 2` CG with cubic regularization 
- `model = 3` ADMM (Boyd's code)
- `model = 4` solve using all the above. In other words, solve each problem 3 times, using a different solver each time. 

2. **MODIFY PROBLEM SIZE**

 - `m`: Number of rows in feature matrix
 - `n`: Number of columns in feature matrix

3. **MODIFY NUMBER OF PROBLEMS**

 - `numProbs`: Number of randomly generated problems 

## CONTACT 

Email: cb3452@drexel.edu 

