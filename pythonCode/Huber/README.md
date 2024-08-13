
# ABOUT
## Python Code - Huber Example

- Authors: Cassidy Buhler and Hande Benson
- Code is adapted from Huber fitting code by Boyd https://web.stanford.edu/~boyd/papers/admm/

## WHAT THE CODE DOES 

This code generates random problems that are solved with Conjugate Gradient Method (CGM) with Cubic Regularization.

## HOW TO USE IT

The driver code is the file *huber_Example.py*. Run this file. 

1. **MODIFY SOLVER**

There are 3 solvers that can be called by the toggling the variable called *model*. 

*model = 0* -> solve using CG without cubic regularization 

*model = 1* -> solve using CG with cubic regularization 

*model = 2* -> solve using all the above. In other words, solve each problem twice, using a different solver each time. 

2. **MODIFY PROBLEM SIZE**

To change the problem size, you can modify *m* and *n*. 

3. **MODIFY NUMBER OF PROBLEMS**

We run the problem 100 times by using a for loop. To change the number of problems, modify the numProbs variable. 


## CONTACT 

Email: cassidy.buhler@gmail.com
