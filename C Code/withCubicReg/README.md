# ABOUT

## C Code - With Cubic Regularization

- Authors: Cassidy Buhler, Hande Y. Benson, David F. Shanno
- Code is adapted from Conmin by David F. Shanno (https://dl.acm.org/doi/pdf/10.1145/355921.355933)
- Incorporates cubic regularization as described in http://www.optimization-online.org/DB_HTML/2021/09/8610.html

## WHAT THE CODE DOES

It solves unconstrained nonlinear optimization problems using the conjugate gradient method *with cubic regularization*.  The problems must be formulated in AMPL, and the user needs to have AMPL installed.

In our paper, we used the AMPL models available [here](https://vanderbei.princeton.edu/ampl/nlmodels/cute/index.html). These have been converted to AMPL from SIF and originated from the [CUTEst repository](https://github.com/ralna/CUTEst).

## INSTALLATION

The software must be compiled prior to use.  Makefiles are provided for MacOS/Linux/Unix (makefile) and Windows (makefile.vc), and the user should edit the makefile for their platform with information on the locations of their AMPL and C libraries.  Once the edits are completed, simply run
```
make
```
at the command prompt for MacOS/Linux/Unix, or
```
make -f makefile.vc

```
at a developer command prompt (Visual Studio) for Windows.  This will create the executable conmin (or conmin.exe for Windows).


## HOW TO USE IT

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

## CONTACT
Email: benson@drexel.edu

