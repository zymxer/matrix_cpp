# Linear System Solvers

## Project Overview
The goal of this project is the implementation and analysis of the following methods for solving linear systems of equations:

- **Jacobi method**  
- **Gauss-Seidel method**  
- **LU factorization**  

## Assumptions
1. Tests for each method are conducted on linear systems that may arise from the discretization of differential equations.  
2. Implemented matrix operations and solution methods are **not parallelized**.  

## Features
- Fully implemented `Matrix` class with support for basic arithmetic operations, norms, and residual calculations.  
- Iterative methods (Jacobi and Gauss-Seidel) with configurable convergence threshold and maximum iterations.  
- Direct method (LU factorization) for exact solutions.  
- Solution tracking with iteration counts, residual norms, and execution time.  

## Usage
Include `Matrix.h` in your project and call the desired method:

```cpp
matrix A(size, size);
matrix b(size, 1);

// Fill A and b with data
solution sol = matrix::jacobi(A, b, 1e-6);
solution sol2 = matrix::gaussSeidel(A, b, 1e-6);
solution sol3 = matrix::lu(A, b);
```
## Performance

Example execution time (in seconds) for different matrix sizes from tests conducted in `main.cpp`:

| Matrix Size | Jacobi  | Gauss-Seidel | LU Factorization |
|------------:|--------:|-------------:|----------------:|
| 1000        | 0.0427  | 0.0413       | 0.1810          |
| 2000        | 0.1700  | 0.1703       | 1.4969          |
| 3000        | 0.4191  | 0.4266       | 6.8084          |
| 4000        | 0.8067  | 0.7657       | 17.1684         |
| 5000        | 1.1678  | 1.1649       | 32.5792         |
