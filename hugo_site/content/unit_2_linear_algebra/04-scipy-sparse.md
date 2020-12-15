---
title: "Scipy.sparse"
weight: 4 
markup: "mmark"
---

There are seven available sparse matrix types in `scipy.sparse`:

- `csc_matrix`: Compressed Sparse Column format
- `csr_matrix`: Compressed Sparse Row format
- `bsr_matrix`: Block Sparse Row format
- `lil_matrix`: List of Lists format
- `dok_matrix`: Dictionary of Keys format
- `coo_matrix`: COOrdinate format (aka IJV, triplet format)
- `dia_matrix`: DIAgonal format

As indicated by the excellent 
[documentation](https://docs.scipy.org/doc/scipy/reference/sparse.html), the 
`dok_matrix` or `lil_matrix` formats are preferable to construct matrices as they 
support basic slicing and indexing similar to a standard NumPy array.

You will notice that the FD matrix we have constructed for the Poisson problem is 
composed entirely of diagonal elements, as is often the case. If you were constructing a 
similar matrix in MATLAB, you would use the 
[`spdiags`](https://uk.mathworks.com/help/matlab/ref/spdiags.html) function, and 
`scipy.sparse` has its own 
[equivalent](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.spdiags.html). 
However, all the `scipy.sparse` formats also have special methods 
[`setdiag`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.lil_matrix.setdiag.html) 
which provide a more object-orientated method of doing the same thing.

Scipy has a few different direct solvers for sparse matrics, given below:
 
[`spsolve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve.html#scipy.sparse.linalg.spsolve): 
This solves $Ax=b$ where $A$ is converted into CSC or CSR form
 
[`spsolve_triangular`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.spsolve_triangular.html#scipy.sparse.linalg.spsolve_triangular): 
Solves $Ax=b$, where $A$ is assumed to be triangular.

 
[`factorized`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.factorized.html#scipy.sparse.linalg.factorized): 
This computes the $LU$ decomposition of the input matrix $A$, returning a Python 
function that can be called to solve $Ax = b$

[`splu`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.splu.html#scipy.sparse.linalg.splu): 
This computes the $LU$ decomposition of the input matrix $A$ using the popular SuperLU 
library. It returns a Python object of class
[`SuperLU`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.SuperLU.html#scipy.sparse.linalg.SuperLU), 
that has a `solve` method you can use to solve $Ax = b$

Note, `scipy.sparse.linalg` also has many iterative solvers, which we will investigate 
further in the next chapter.

### Problem: construction in `scipy.sparse`

Your goal for this problem is to construct the FD matrix $A$ given in the Finite 
difference section, using `scipy.sparse`, and:

1. Visualise the matrix $A$ using the Matplotlib 
   [`spy`](https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.spy.html) plot
2. Solve the Poisson problem using $f(x) = 2 \cos(x) / e^x$ and $g(x) = \sin(x) / e^x$. 
   Check your answer using the analytical solution $u_{a}(x) = \sin(x) / e^x$.
3. Vary the number of discretisation points $N$ and calculate $AA$ using both sparse and 
   dense matrices. For each $N$ calculate the time to calculate the matix 
   multiplicatiion using Python's 
   [`time.perf_counter`](https://docs.python.org/3/library/time.html#time.perf_counter), 
   and plot time verus $N$ for dense and sparse matrix multiplicatiion. Comment on how 
   the time varies with $N$.
5. Vary the number of discretisation points $N$ and solve the Poisson problem with 
   varying $N$, and with using both the sparse and direct $LU$ solvers. For each $N$ 
   record the time taken for both the dense and sparse solvers, and record the numerical 
   error $||\mathbf{v} - \mathbf{v}_a||_2$. Generate plots of both error and time versus 
   $N$, and comment on how they vary with $N$

