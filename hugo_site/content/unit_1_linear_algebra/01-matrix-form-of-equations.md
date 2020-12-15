---
title: "Matrix form of equations"
date: "2014-04-02"
questions:
- "What is the relationship between matrices and systems of linear equations?"
- "What is a singular matrix and when does it occur?"
- "What is Gaussian Elimination and why is it useful?"
objectives:
- "Understand the main useful concepts for the solution of systems of linear equations"
- "Understand singular matrices and the rank of a matrix"
- "Understand and be able to implement Gaussian Elimination"
weight: 1 
markup: "mmark"
---


## Simultaneous equations

Consider 2 simultaneous equations:

$$
\begin{aligned}
a_1x+b_1y &= c_1, \quad (1)\\
a_2x+b_2y &= c_2, \quad (2)
\end{aligned}
$$

where the values $\;x\;$ and $\;y\;$ are to be found, and $\;a_1, \;b_1, \;a_2, \;b_2, \;c_1\;$ and $\;c_2\;$ are given constants.

$$
\begin{aligned}
 (1) \times b_2:~~~~~~~~~~~~~~~ b_2a_1x+b_2b_1y &=& b_2c_1, \quad (3)\\
  (2) \times b_1:~~~~~~~~~~~~~~~ b_1a_2x+b_1b_2y &=& b_1c_2, \quad (4)\\
  (3) - (4):~~~~~~~~~~~~~~~ b_2a_1x-b_1a_2x &=& b_2c_1-b_1c_2.
\end{aligned}
$$

## The Matrix

Matrices are a structure that allow us to more easily manipulate linear systems. 

Consider the original system

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

### The determinant

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
introduced below)

If $|A| = 0$, A is said to be **singular** (have no inverse).

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
What happens if all right-hand sides are zero? Is there any nonzero choice of right- hand sides that allows the three lines to intersect at the same point?

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
