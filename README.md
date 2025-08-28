# MG_FDQGMRES_method
This repository contain the codes (written in C++) for the multigrid preconditioned Flexible Direct Quasi Generalized Minimal Residual method (FDQGMRES) for variable-coefficient Poisson equations. Finite Volume Method is used to discretize the Poisson equation.

The code is designed to solve the second test problem (Problem 2) of Tatebe (1993).

# Reference
Tatebe, O., 1993, November. The multigrid preconditioned conjugate gradient method. In NASA. Langley Research Center, The Sixth Copper Mountain Conference on Multigrid Methods, Part 2.

# Software requirements
This solver needs:

- gcc

# How to install the required packages (on a Linux system)

To install gcc (compiler for C++ code)

```bash
sudo apt install build-essential
```

# How to compile and run the code

To compile the code

```bash
g++ ht.cpp -o output
```
To run this code

```bash
./output
```
