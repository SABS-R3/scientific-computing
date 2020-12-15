---
title: "Finite Difference Matrix"
weight: 3 
markup: "mmark"
---

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
