# MG_FDQGMRES_method
This repository contain the codes (written in C++) for the multigrid preconditioned Flexible Direct Quasi Generalized Minimal Residual method (FDQGMRES) for variable-coefficient Poisson equations. Finite Difference Method is used to discretize the Poisson equation.

The codes are designed following the works of Saad (2003) and Wessling (1992).

My hand-written notes on various Krylov subspace-based iterative methods and the geometric multigrid method are provided in the following Google Drive links.

Krylov subspace-based methods: https://drive.google.com/file/d/1uav4r3Ulcez8wGPhFXVXYzOZBUtvQJDR/view?usp=drive_link

Geometric multigrid method: https://drive.google.com/file/d/1ZEe2J0fHEkPguUGeVJhy_0zrihB5x2V8/view?usp=drive_link

The code is designed to solve the second test problem (Problem 2) of Tatebe (1993).

# Reference
Saad, Y., 2003. Iterative methods for sparse linear systems. Society for Industrial and Applied Mathematics.
Wesseling, P., 1992. An introduction to multigrid methods. John Wiley & Sons.
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
