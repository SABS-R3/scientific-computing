---
title: "Matrix form of equations"
weight: 1 
markup: "mmark"
---

## Systems of linear equations

Linear algebra is largely concerned with representing, manipulating and solving large 
systems of linear equations. Consider the following 2 linear equations:

$$
\begin{aligned}
a_1x+b_1y &= c_1, \quad (1)\\
a_2x+b_2y &= c_2, \quad (2)
\end{aligned}
$$

where the values $\;x\;$ and $\;y\;$ are to be found, and $\;a_1, \;b_1, \;a_2, \;b_2, 
\;c_1\;$ and $\;c_2\;$ are given constants. We know that we can use linear combinations 
of these two equations to solve this sytem for $\;x\;$ and $\;y\;$, like so:

$$
\begin{aligned}
 (1) \times b_2:~~~~~~~~~~~~~~~ b_2a_1x+b_2b_1y &=& b_2c_1, \quad (3)\\
  (2) \times b_1:~~~~~~~~~~~~~~~ b_1a_2x+b_1b_2y &=& b_1c_2, \quad (4)\\
  (3) - (4):~~~~~~~~~~~~~~~ b_2a_1x-b_1a_2x &=& b_2c_1-b_1c_2.
\end{aligned}
$$

## The Matrix

We can also write these linear equations using a matrix. Matrices are structures that 
allow us to more easily manipulate linear systems. While not particularly useful for 
just 2 equations, using a matrix representation allows us to generalise to, say $N$ 
equations and $M$ unknowns, or to solve large systems of equations on a computer.

Consider the original system:

$$
\begin{aligned}
a_1x+b_1y &= c_1, \\
a_2x+b_2y &= c_2. 
\end{aligned}
$$

We rewrite this, in the form of a matrix as:

$$
\left(\begin{matrix}a_1&b_1\\ a_2&b_2\end{matrix}\right)
\left(\begin{matrix}x\\y\end{matrix}\right)
=\left(\begin{matrix}c_1\\ c_2 \end{matrix}\right).
$$
 
Think about how this form relates to the original linear system. 

## Geometry of linear equations

Consider the following system of equations

$$
\begin{aligned}
 x + -2y &= -1, \\
-x +  3y &=  3. 
\end{aligned}
$$

That is,

$$A = \left(\begin{matrix} 1 & -2 \\ -1 & 3\end{matrix}\right)$$

Plotting these two linear equations on a graph shows graphically the solution to this 
equation given by

$$\left(\begin{matrix} x \\ y \end{matrix}\right) = A^{-1} \left(\begin{matrix} -1 \\ 3 
\end{matrix}\right) = \left(\begin{matrix} -1 \\ 3 \end{matrix}\right)$$

![single solution](/scientific-computing/images/unit_01/01-sim1.svg)

Now lets consider two different system of equations represented by the matrix:

$$
A = \left(\begin{matrix} 1 & -2 \\ -1 & 2\end{matrix}\right)
$$

$$
\begin{aligned}
 x + -2y &= -1, \\
-x +  2y &=  3. \end{aligned}
$$

and 

$$
\begin{aligned}
 x + -2y &= -1, \\
-x +  2y &=  1. \end{aligned}
$$

![(left) infinite solutions, (right) no 
solutions](/scientific-computing/images/unit_01/01-sim2.svg)

The first gives the plot on the left, while the second, which has a different vector of 
constants on the RHS, gives the plot on the right. You can see that depending on the 
constants, the system of equations can have an infinite number of solutions, or no 
solutions at all. 

The matrix $A$ in this case is singular, and therefore does not have an inverse. Looking 
at the equations again, you can see that the two rows of the matrix $A$ are multiples of 
the other, and thus there is only *one* independent row. That is, the *rank* of the 
matrix is one.

## Singular matrices

The *rank* of an $\;n\,\times\,n\;$ matrix $\;A\;$ is the number of linearly independent rows in $\;A\;$ (rows not combinations of other rows).

When $\;\text{rank}(A) < n\;$ then

- The matrix is said to be 'rank deficient'
- The system $\;A\textbf{x} = \textbf{b}\;$ has *fewer* equations than unknowns
- The matrix is said to be singular
- The matrix is said to be underdetermined
- $A\;$ has no inverse
- The determinant of $\;A\;$ is 0
- The equation $\;A\textbf{u} = \textbf{0}\;$ has non-trivial solutions ($\textbf{u} \neq \textbf{0}$)

"## Null space\n",
    "When a matrix is singular we can find non-trivial solutions to $\\;A\\textbf{u} = \\textbf{0}$.\n",
    "\n",
    "These are vectors which form a *null space* for $\\;A$.\n",
    "\n",
    "These vectors  make no difference to the effect that $A$ is having:\n",
    "\n",
    "$$\n",
    " A(\\textbf{x} + \\textbf{u}) =  A\\textbf{x} + A\\textbf{u} =\n",
    " A\\textbf{x} + \\textbf{0} =   A\\textbf{x}.\n",
    "$$ "

### The determinant

One way of solving a system of equations represented by $A x = b$ is to calculate the 
*inverse* of A, giving the solution as $x = A^{-1} b$. This can be done by calculating 
what is known as the *determinant* of $A$.

If

$$A = \left(\begin{matrix} p & q \\ r & s\end{matrix}\right)$$

then the **determinant** of A is:

$$|A| = ps-qr$$

The inverse of $A$ can be found using the determinant:

$$
A^{-1} = \frac{1}{ps-qr} \left(\begin{matrix} s & -q \\ -r & p\end{matrix}\right)
$$

Calculating the inverse of a matrix using its determinant can be very costly for larger 
matrices, therefore other algorithms are used (e.g. Gaussian Elimination, which is 
introduced in the next section)

If $|A| = 0$, A is said to be **singular** (have no inverse). Graphically, this is 
represented by the parallel or non-intersecting lines in the figure above. 

### Using Python to calculate the inverse

To find $A^{-1}$  for

$$
A = \left(\begin{matrix}3&0&2\\ 3&0&-3\\ 0&1&-1\end{matrix}\right)
$$
 
you can using numpy like so:
 
```python
A = np.array([[3, 0, 2], [3, 0, -3], [0, 1, 1]])
np.linalg.inv(A)
```

Output:

```
array([[ 0.2       ,  0.13333333,  0.        ],
       [-0.2       ,  0.2       ,  1.        ],
       [ 0.2       , -0.2       , -0.        ]])
```

It doesn't always work. Consider the following system

$$
A = \left(\begin{matrix}1&1&1\\ 2&4&2\\ 7&10&7\end{matrix}\right)
$$

```python
A = np.array([[1, 1, 1], [2, 5, 2], [7, 10, 7]])
np.linalg.inv(A)
```

Output:

```
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
  File "<__array_function__ internals>", line 6, in inv
  File "/home/mrobins/git/scientific-computing/env/lib/python3.6/site-packages/numpy/linalg/linalg.py", line 546, in inv
    ainv = _umath_linalg.inv(a, signature=signature, extobj=extobj)
  File "/home/mrobins/git/scientific-computing/env/lib/python3.6/site-packages/numpy/linalg/linalg.py", line 88, in _raise_linalgerror_singular
    raise LinAlgError("Singular matrix")
numpy.linalg.LinAlgError: Singular matrix
```

### Null space

When a matrix is singular we can find non-trivial solutions to  $Ax=0$.

These are vectors which form a null space for $A$.

These vectors make no difference to the effect that $A$ is having:

$$
A(x+u) = Ax + Au = Ax+0 =Ax.
$$
 
Note that any combination or scaling of vectors in the null space is also in the null 
space. That is, if $Au=0$ and $Av=0$ then

$$
A(\alpha u+ \beta v)=0
$$

The number of linearly independent vectors in the null space is denoted $\text(null)(A)$ 
and

$$
\text{null}(A)+\text{rank}(𝐴)=\text{order}(A).
$$

### Null space example

Previous example of a singular system:

$$
A = \left(\begin{matrix}1&1&1\\ 2&4&2\\ 7&10&7\end{matrix}\right)
$$


```python
A = np.array([[1, 1, 1], [2, 5, 2], [7, 10, 7]])
np.linalg.matrix_rank(A)
```

Output:

```
2
```


```python
import scipy.linalg
scipy.linalg.null_space(A)
```

Output:

```
array([[ 0.70710678],
       [ 0.        ],
       [-0.70710678]])
```


remember, scaled vectors in the null space are also in the null space, for example, 
$x=1,y=0,z=−1$ is in the null space.

Try it:

$$
\left(\begin{matrix} 1& 1& 1\\ 2& 5& 2\\ 7& 11&7 \end{matrix}\right)
  \left(\begin{matrix} -1000\\ 0 \\ 1000 \end{matrix}\right) = \quad ?
$$

```python
np.matmul(A,np.array([-1000,0,1000]))
```

Output:

```
array([0, 0, 0])
```

## Column space 

Consider the following matrix $A$

$$
A = \left(\begin{matrix} 1 & -2 \\ -1 & 2\end{matrix}\right)
$$

Multiplying this matrix by a 2d vector $x$ give another 2d vector $b$

$$
\left(\begin{matrix} 1 & -2 \\ -1 & 2\end{matrix}\right)
\left(\begin{matrix} x_1 \\ x_2\end{matrix}\right)
=
\left(\begin{matrix} b_1 \\ b_2\end{matrix}\right)
$$

If the equation above is *consistent* (i.e. true), then we can say that the vector $b$ 
lies in the *column space* of $A$. The column space of $A$ is all the vectors that can 
be made by linear combinations of the columns of $A$. Here the columns of $A$ are $c_1 = 
(1, -1)^T$ and $c_2=  (-2, 2)^T$, and you can see the $b$ is equal to the linear 
combination of $x_1 c_1 + x_2 c_2$.

The matrix $A$ given above is singular with rank 1. Above we gave two examples with 
different $b$ vectors

$$
\left(\begin{matrix} 1 & -2 \\ -1 & 2\end{matrix}\right)
\left(\begin{matrix} x_1 \\ x_2\end{matrix}\right)
=
\left(\begin{matrix} -1 \\ 3\end{matrix}\right)
$$

and 

$$
\left(\begin{matrix} 1 & -2 \\ -1 & 2\end{matrix}\right)
\left(\begin{matrix} x_1 \\ x_2\end{matrix}\right)
=
\left(\begin{matrix} -1 \\ 1\end{matrix}\right)
$$

We showed graphically that $b = (-1, 3)^T$ had an infinite number of solutions $x$, 
whereas $b=(-1, 1)$ had no solutions. Therefore, $b = (-1, 3)^T$ is *within* the column 
space of $A$, whereas $b=(-1, 1)$ is not.

Naturally, if the matrix $A$ is *not* singular, then all possible vectors $b$ lie in the 
column space of $A$.

## Eigenvectors: motivation

The *eigenvalues* and *eigenvectors* give an indication of how much effect the matrix has and in what direction. 

- $A=\left(\begin{matrix} cos(45)&-sin(45)\\ sin(45)&cos(45)\\\end{matrix}\right)$ has no scaling effect.


- $B=\left(\begin{matrix} 2& 0 \\ 0&\frac{1}{2}\\\end{matrix}\right)$ doubles in the $x$-direction, but halves in the $y$-direction.


Repeated applications of $\;A\;$ stay the same distance from the origin, but repeated applications of $\;B\;$ move towards $\;(\infty, 0).$

Eigenvalues and eigenvectors are useful in the following applications:

- Transitions with probability
- Markov chains
- Google Page ranks
- Solution of systems of linear ODEs
- Stability of systems of nonlinear ODEs


Mathematically, let $A\;$ be a matrix, $\;\textbf{v}\;$ be a *non-zero* vector, and
$\;\lambda\;$ be a scalar, 

If,

$$A \textbf{v} = \lambda \textbf{v}$$

then $\;\textbf{v}\;$ is called an *eigenvector* and $\;\lambda\;$ is the corresponding *eigenvalue*.

Note that if $\;\textbf{v}\;$ is a solution, then so is a scaling $\;a\textbf{v}$:

$$A (a \textbf{v}) = \lambda (a \textbf{v}).$$

## Calculating Eigenvalues/Eigenvectors

You are probably already familiar with calculating eigenvalues and eigenvectors 
analytically. If we write:

$$
\begin{eqnarray*}
A \textbf{v} &=& \lambda \textbf{v},\\
A \textbf{v} -  \lambda I \textbf{v}&=& \textbf{0},\\
(A  -  \lambda I) \textbf{v}&=& \textbf{0}.
\end{eqnarray*}
$$

Since $v$ is non-zero, $(A - \lambda I)$ must be singular, and so we find the 
eigenvalues be solving 

$$|A-\lambda I|=0.$$

For example, 

$$A=\left(\begin{matrix}-2&-2\\ 1&-5\\\end{matrix}\right)$$

$$|A-\lambda I|=\left\vert\begin{matrix}-2-\lambda&-2\\ 1&-5-\lambda\end{matrix}\right\vert=(-2-\lambda)(-5-\lambda)-(-2)$$

$$=10+5\lambda+\lambda^2+2\lambda+2=\lambda^2+7\lambda+12=(\lambda+3)(\lambda+4)=0$$

Finding the same eigenvalues, and associated eigenvectors using Python can be done using 
Numpy:

```python
A = np.array([[-2, -2], [1, -5]])
np.linalg.eig(A)
```

Output:

```
(array([-3., -4.]), array([[0.89442719, 0.70710678],
       [0.4472136 , 0.70710678]]))
```


## Other Reading

Linear algebra by Ward Cheney
Linear algebra and its applications by David C. Lay.

Strang, G. (2016). Introduction to linear algebra (Fifth ed.). Wellesley.

Linear algebra and its applications by Gilbert Strang

lots of supplimentary material for this last two via MIT course page here:
https://github.com/mitmath/1806/blob/master/summaries.md


- LA from an ODE perspective
Kapitula, T. (2015). Ordinary Differential Equations and Linear Algebra. Society for 
Industrial and Applied Mathematics.

## Problems

1. Describe the intersection of the three planes u+v+w+z = 6 and u+w+z = 4 and u+w = 2 
   (all in four-dimensional space). Is it a line or a point or a fourth equation that 
   leaves us with no solution. an empty set? What is the intersection if the fourth 
   plane u = −1 is included? Find a fourth equation that leaves us with no solution.

2. Sketch these three lines and decide if the equations are solvable: 3 by 2 system
  x + 2y = 2 x − y = 2 y = 1.
What happens if all right-hand sides are zero? Is there any nonzero choice of right- 
hand sides that allows the three lines to intersect at the same point?

3. Write a Python function that takes in a triangular matrix $A$ represented as an 
   `ndarray`, and a rhs vector $b$, and solves the equation $A x = b$. i.e. the function 
   will solve the following triangular system for $x = (x_1, x_2, x_3)$:

$$
\begin{aligned}
A_{11} x_1 + A_{12} x_2 + A_{13} x_3 &= b_1, \\
A_{22} x_2 + A_{23} x_3 &= b_2, \\
A_{33} x_3 &= b_3
\end{aligned}
$$
