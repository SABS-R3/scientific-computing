---
title: "LDL decomposition"
date: "2014-04-01"
markup: "mmark"
weight: 5 
---

## $LDL$ decomposition

It is often very benificial when solving linear systems to consider and take advantage 
of any special structure that the matrix $A$ might possesses. The $LDL$ decomposition is 
a varient on LU decomposition which is only applicable to a symmetric matrix $A$ (i.e. 
$A = A^T$). The advantage of using this decomposition is that it takes advantage of the 
redundent entries in the matrix to reduce the amount of computation to $n^3/3$, which is 
about a half that required for the $LU$ decomposition.

### Other reading


- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 4.1
- https://en.wikipedia.org/wiki/Cholesky_decomposition#LDL_decomposition_2

### Software

[`scipy.linalg.ldl`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.ldl.html#scipy.linalg.ldl)
