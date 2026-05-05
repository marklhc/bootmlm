# Get the square root of a matrix using eigenvalue decomposition, and solve for a linear system

Solve for \\x\\ in the linear system \\Ax = b\\, where \\A\\ is a
symmetric matrix square root of \\M\\ with eigenvalue decomposition such
that \\AA = M\\.

## Usage

``` r
solve_eigen_sqrt(M, b)
```

## Arguments

- M:

  a symmetric positive definite/semi-definite matrix

- b:

  a numeric vector or matrix
