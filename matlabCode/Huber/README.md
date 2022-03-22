# ABOUT
## Matlab Code - Huber Examples

- Authors: Cassidy Buhler and Hande Benson
- Date code was last modified: March 11, 2022
- Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/

## WHAT THE CODE DOES 

This code generates random problems that are solved with Conjugate Gradient Method (CGM) with and without Cubic Regularization.
We also include Boyd's ADMM for comparison.

## HOW TO USE IT

1. **MODIFY SOLVER**

There are 3 solvers that can be called by the toggling the variable called *model*. 

*model = 1;* -> solve using CG without cubic regularization 

*model = 2;* -> solve using CG with cubic regularization 

*model = 3;* -> solve using ADMM (Boyd's code)

*model = 4;* -> solve using all the above. In other words, solve each problem 3 times, using a different solver each time. 

2. **MODIFY PROBLEM SIZE**

To change the problem size, you can modify *m* and *n*. 

3. **MODIFY NUMBER OF PROBLEMS**

We run the problem 100 times by using a for loop. To change the number of problems, modify the *rr* variable. 


## CONTACT 

Email: cb3452@drexel.edu 

