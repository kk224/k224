# A Practical Attack on kk224

The codes for attacking kk224 via an allocating approach are provided. 

The constructed polynomial system for Block i is contained in 'blocki/data_processing/keccakEquations/' where i = 1 or 2. The strategy for solving these quadratic polynomial systems is given below.

Let Q be a quadratic polynomial system. We enumerate the values of K selected linear polynomials, such that a solvable linear system L is obtained from equations in Q. Let sol_L be the solution(s) to the linear system solved via Gaussian Elimination. Then we verify this solution by the other equations in Q afterwards. If sol_L is a solution to Q, then it is returned; otherwise, other values of the K selected linear polynomials are tried. Particularly, not all equations should be satisfied for the system in Block 1, and we only  expect the equations hold as many as possible.

To speed up the computations, polynomials are converted to matrices in our implementation. We also provide an efficient algorithm for solving the systems using GPU. The packages CUDD and M4RI are used to process data.


Structures
-----------------------------------------------------------
    .
    ├── block1                      # codes for solving systems in Block 1
        ├── data_processing         # processing data
            ├── keccakEquations     # polynomial equations generated from the 1st block
            ├── b1P2M               # codes for converting polynomials to matrices            
            └── Mat                 # matrices converted from polynomial equations
        └── run                     # source and executable files for solving nonlinear systems by GPU
    ├── block2                      # codes for solving systems in Block 2
        ├── data_processing         # processing data
            ├── keccakEquations     # polynomial equations generated from the 2nd block
            ├── b1P2M               # codes for converting polynomials to matrices            
            └── Mat                 # matrices converted from polynomial equations
        └── run                     # source and executable files for solving nonlinear systems by GPU
    └── README.md

There are some 'Readme.md' files in some folders describing the detailed procedures.
