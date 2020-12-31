---
title: "Jacobi and Relaxation Methods"
weight: 6 
markup: "mmark"
---


## Jacobi Method

The Jacobi method is the simplest of the iterative methods, and relies on the fact that 
the matrix is *diagonally dominant*. Starting from the problem definition:

$$
A\mathbf{x} = \mathbf{b}
$$

we decompose $A$ in to $A = L + D + U$, where $L$ is lower triangular, $D$ is diagonal, 
$U$ is upper triangular. 

$$
A\mathbf{x} = L\mathbf{x} + D\mathbf{x} + U\mathbf{x} =  \mathbf{b}
$$

We then assume that we have an initial guess at the solution $\mathbf{x}^0$, and try to 
find a new estimate $\mathbf{x}^1$. Assuming that the diagonal $D$ dominates over $L$ 
and $U$, a sensible choice would be to insert $x^0$ and the unknown $x^1$ into the 
equation like so:

$$
L\mathbf{x}^0 + D\mathbf{x}^1 + U\mathbf{x}^0 =  \mathbf{b}
$$

we can rearrange to get an equation for $x^1$. This is easily solved as we can take the 
inverse of the diagonal matrix by simply inverting each diagonal element individually:

$$
D\mathbf{x}_1 =  \mathbf{b} - (L+U)\mathbf{x}_0
$$

Thus we end up with the general Jacobi iteration:

$$
\mathbf{x}_{k+1} =  D^{-1}(\mathbf{b} - (L+U)\mathbf{x}_k)
$$

## Relaxation methods

The Jacobi method is an example of a relaxation method, where the matrix $A$ is split 
into a dominant part $M$ (which is easy to solve), and the remainder $N$. That is, $A = 
M - N$

$$M\mathbf{x}_{k+1} = N\mathbf{x}_k + \mathbf{b}$$
$$\mathbf{x}_{k+1} = M^{-1}N\mathbf{x}_k + M^{-1}\mathbf{b}$$

For the Jacobi method $M = D$ and $N = -(L + U)$. Other relaxation methods include 
Gauss-Seidel, where $M = (D + L)$ and $N = -U$, and successive over-relaxation (SOR), 
where $M = \frac{1}{\omega} D + L$ and $N = -\frac{1 - \omega}{\omega} D - U$, where 
$\omega$ is the *relaxation* parameter.

For any relaxation method to converge we need $\rho(M^{-1}N) < 1$, where $\rho()$ is the 
*spectral radius* of $M^{-1} N$, which is defined as the largest eigenvalue $\lambda$ of 
a a given matrix $G$:

$$
\rho(G) = \max{|\lambda|: \lambda \in \lambda(G)}
$$

For the SOR method, the relaxation parameter $\omega$ is generally chosen to minimise 
$\rho(M^{-1}N)$, so that the speed of convergence is maximised. In some cases this 
optimal $\omega$ is known, for example for finite difference discretisation of the 
[Poisson equation](https://www.sciencedirect.com/science/article/pii/S0893965908001523).
However, in many cases sophisticated eigenvalue analysis is required to determine the 
optimal $\omega$. 

### Other Reading

- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 10 
- Barrett, R., Berry, M., Chan, T. F., Demmel, J., Donato, J., Dongarra, J., ... & Van 
  der Vorst, H. (1994). Templates for the solution of linear systems: building blocks 
  for iterative methods. Society for Industrial and Applied Mathematics.

### Problems


This exercise involves the manipulation and solution of the linear system resulting from 
the finite difference solution to Poisson's equation in *two* dimensions. Let $A$ be a 
sparse symmetric positive definite matrix of dimension $(N-1)^2 \times (N-1)^2$ created 
using `scipy.sparse` (for a given $N$) by the function
`buildA` as follows:
```python
import numpy as np
import scipy.sparse as sp

def buildA(N):
    dx = 1 / N
    nvar = (N - 1)**2;
    e1 = np.ones((nvar), dtype=float);
    e2 = np.copy(e1)
    e2[:N-1:] = 0
    e3 = np.copy(e1)
    e3[N-1:N-1:] = 0
    A = sp.spdiags(
        (-e1, -e3, 4*e1, -e2, -e1),
        (-N-1, -1, 0, 1, N-1), nvar, nvar
    )
    A = A / dx**2;
    return A
```

and let $\mathbf{f}_1$ and $\mathbf{f}_2$ be the vectors defined in
`buildf1` and `buildf2`

```python
def buildf1(N):
    x = np.arange(0, 1, 1/N).reshape(N, 1)
    y = x.T
    f = np.dot(np.sin(np.pi*x), np.sin(np.pi*y))
    return f[1:,1:].reshape(-1,1)
```

```python
def buildf2(N):
    x = np.arange(0, 1, 1/N).reshape(N, 1)
    y = x.T
    f = np.dot(np.maximum(x,1-x), np.maximum(y,1-y))
    return f[1:,1:].reshape(-1, 1)
```

We will consider manipulation of the matrix $A$ and solution of the linear
systems $A\mathbf{U}_i=\mathbf{f}_i$. The solution to this linear system
corresponds to a finite difference solution to Poisson's equation $-\nabla^2 u
= f$ on the unit square with zero Dirichlet boundary conditions where $f$ is
either $\sin(\pi x) \sin (\pi y)$ or $\max(x,1-x) \max(y,1-y)$. PDEs of this type occur 
(usually with some additional reaction and or convection terms) very frequently
in mathematical modelling of physiological processes, and even in image
analysis. 

1. Write a function to solve a linear system using the Jacobi method. In
  terms of $N$, how many iterations does it take to converge? (Try
  $N=4,8,16,32,64$.)
2. Write a function to solve a linear system using the SOR method. For
  $N=64$ and right-hand-side $\mathbf{f}_2$ determine numerically the best
  choice of the relaxation parameter to 2 decimal places and compare this
  with theory.
