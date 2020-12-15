---
title: "Sparse matrices"
weight: 1 
markup: "mmark"
---

## Why sparse matrices

Taking advantage of any special structure in the matrix of interest is always of great 
importance when designing a linear algebra algorithm/solver. Thus far we have discussed 
special structures such as symmetric or positive definite matrices, but one of the most 
common matrix structures in scientific computing is that of a *sparse matrix*, or a 
matrix containing many zero elements. Since zeros can be ignored in many computations, 
for example multiplication or addition, a sparse matrix will specify a special data 
structure so that only the *non-zero* elements of the matrix are actually stored and 
used in computation. Note also that a sparse matrix can itself be symmetric or positive 
definite, and it is often neccessary to take all these properties into account when 
designing your algorithm.

The most obvious benifit of a sparse matrix is low memory usage. A dense matrix needs to 
store all the elements of a matrix of size $n$, and thus the memory requirements scale 
as $\mathcal{O}(n^2)$. For example, the following show the memory requirements of a 
matrix of double precision numbers (taken from the excellent 
[scipy-lectures](http://scipy-lectures.org/advanced/scipy_sparse/introduction.html#why-sparse-matrices)) 


![least squares problem](/scientific-computing/images/unit_02/sparse_versus_dense.svg)

A sparse matrix only stores non-zero elements, and in many different applications this 
represents a huge memory saving as matrices are often very sparse, holding only a few 
non-zero elements. Some typical applications are:

- solution of partial differential equations (PDEs), such as the finite difference 
  method illustrated below
- applications of graphs or networks (e.g. electrical networks, website links), where 
  non-zero elements of the matrix represent edges between nodes


Note that while a sparse matrix has obvious benifits in terms of matrix multiplication, 
where the zero elements can simply be ignored, direct solver algorithms such as $LU$ 
decodecomposition for the problem $Ax = b$, where $A$ is sparse, need considerably more 
thought as the zeros in $A$ can have propogating effects, and there is no guarentee that 
the decomposition of a $A$ or its inverse will be itself sparse, there can be a 
significant amount of what is known as *fill-in*. This fact motivates a separate class 
of *iterative* (as opposed to *direct*) solvers that only rely on the matrix 
multiplication of $A$ with a vector, ignoring the internal sparsity structure of $A$ and 
only taking advantage of the increased speed of the matrix multiplication itself. These 
iterative solvers will be covered in the following chapter, but in this chapter we will 
focus on the practical requirements of constructing and using sparse matrices using the 
`scipy.sparse` library.

## An example sparse format

The COOrdinate format (COO) is also known as the "ijv" or "triplet" format, and stores 
the non-zero elements in three arrays, `row`, `col`, and `data`. The `data[i]` value is 
the non-zero entry in row `row[i]` and column `col[i]` of the matrix. The advantages of 
this format are:

- fast format for contructing sparse matrices
- fast conversions to/from the CSR and CSC formats
- fast matrix-vector multiplication
- fast elementwise operations (e.g. multiply each element by 2 is just `data * 2`)

However, slicing using this format is difficult. 

Here are some examples of the COO matrix format using 
[`scipy.sparse.coo_matrix`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.coo_matrix.html). 
Again, these have been taken from 
[scipy-lectures](http://scipy-lectures.org/advanced/scipy_sparse/introduction.html#why-sparse-matrices), 
which is an excellent resource and contains examples of the other sparse matrix formats 
implemented in Scipy.

### create empty COO matrix:

```python
mtx = sparse.coo_matrix((3, 4), dtype=np.int8)
mtx.todense()
```

Output:
```
matrix([[0, 0, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 0]], dtype=int8)
```

### create using (data, ij) tuple:

```python
row = np.array([0, 3, 1, 0])
col = np.array([0, 3, 1, 2])
data = np.array([4, 5, 7, 9])
mtx = sparse.coo_matrix((data, (row, col)), shape=(4, 4))
mtx
mtx.todense()
```

Output:

```
>>> mtx
<4x4 sparse matrix of type '<class 'numpy.int64'>'
        with 4 stored elements in COOrdinate format>
>>> mtx.todense()
matrix([[4, 0, 9, 0],
        [0, 7, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 5]])
```


### duplicates entries are summed together:

```python
row = np.array([0, 0, 1, 3, 1, 0, 0])
col = np.array([0, 2, 1, 3, 1, 0, 0])
data = np.array([1, 1, 1, 1, 1, 1, 1])
mtx = sparse.coo_matrix((data, (row, col)), shape=(4, 4))
mtx.todense()
```

Output:

```
>>> mtx.todense()
matrix([[3, 0, 1, 0],
        [0, 2, 0, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 1]])
```

### no slicing…:

```python
mtx[2, 3]
```

Output:

```
>>> mtx[2, 3]
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
TypeError: 'coo_matrix' object is not subscriptable
```


## An example sparse matrix arising from a finite difference discretistaion

Many matrices in scientific computing contain mostly zeros, particularly those arising 
from the discretistaion of partial differential equations (PDEs). Here we will construct 
a sparse matrix using `scipy.sparse` that is derived from the finite difference 
discretistaion of the Poisson equation. In 1D, Poisson equation is

$$u_{xx} = f(x)\text{ for }0 \le x \le 1$$

The central FD approximation of $u_{xx}$ is:

$$u_{xx} \approx \frac{u(x + h) - 2u(x) + u(x-h)}{h^2}$$

We will discretise $u_{xx} = 0$ at $N$ regular points along $x$ from 0 to 1, given by 
$x_1$, $x_2$:

              +----+----+----------+----+> x
              0   x_1  x_2    ... x_N   1

Using this set of point and the discretised eqution, this gives a set of $N$ equations 
at each interior point on the domain:

$$\frac{v_{i+1} - 2v_i + v_{i-1}}{h^2} = 0 \text{ for } i = 1...N$$

where $v_i \approx u(x_i)$.

To solve these equations we will need additional equations at $x=0$ and $x=1$, known as 
the *boundary conditions*. For this example we will use $u(x) = g(x)$ at $x=0$ and $x=1$ 
(also known as a non-homogenous Dirichlet bc), so that $v_0 = g(0)$, and $v\_{N+1} = 
g(1)$, and the equation at $x_1$ becomes:

$$\frac{v_{i+1} - 2v_i + g(0)}{h^2} = 0$$

and the equation at $x_N$ becomes:

$$\frac{g(1) - 2v_i + v_{i-1}}{h^2} = 0$$

We can therefore represent the final $N$ equations in matrix form like so:

$$
\frac{1}{h^2}
\begin{bmatrix} -2      & 1      &         &   &     \\
 1      & -2     & 1       &       & \\
&\ddots & \ddots  &  \ddots &\\
&        & 1      &  -2     &  1     \\
&        &        &   1     & -2     \end{bmatrix}
\begin{bmatrix} v_1    \\
v_2    \\
\vdots \\
v_{N-1}\\
v_{N}  
\end{bmatrix}
= \begin{bmatrix} -g(0)    \\
0    \\
\vdots \\
0    \\
-g(1)
\end{bmatrix}
$$

The relevent sparse matrix here is $A$, given by


$$
A = \begin{bmatrix} -2      & 1      &         &   &     \\
 1      & -2     & 1       &       & \\
&\ddots & \ddots  &  \ddots &\\
&        & 1      &  -2     &  1     \\
&        &        &   1     & -2     \end{bmatrix}
$$

As you can see, the number of non-zero elements grows linearly with the size $N$, so a 
sparse matrix format is much prefered over a dense matrix holding all $N^2$ elements!

## Additional Reading

For more on the Finite Difference Method for solving PDEs:

K. W. Morton and D. F. Mayers. Numerical Solution of Partial Differential Equations: An
Introduction. Cambridge University Press, 2005.

## Software

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

Your goal for this problem is to construct the FD matrix $A$ given above, using 
`scipy.sparse`, and:

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

