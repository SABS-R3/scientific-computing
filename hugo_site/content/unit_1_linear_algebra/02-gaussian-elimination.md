---
title: "Gaussian Elimination"
date: "2014-04-02"
questions:
- "What is the relationship between matrices and systems of linear equations?"
- "What is a singular matrix and when does it occur?"
- "What is Gaussian Elimination and why is it useful?"
objectives:
- "Understand the main useful concepts for the solution of systems of linear equations"
- "Understand singular matrices and the rank of a matrix"
- "Understand and be able to implement Gaussian Elimination"
weight: 2 
markup: "mmark"
---

## Gaussian Elimination

Consider the problem: Find $x, y,z$ such that 

$$
\begin{aligned}
eq1: &  2x & + & y &  -&  z & = & 3 \\
eq2: & x & &  &+ &5z & = & 6 \\
eq3: & -x &+& 3y& -& 2z & = &  3
\end{aligned}
$$

**Gaussian elimination -- step 1** reduce the above system of equations so that the 
unknown $x$ is removed from the last two equations:

$$
\begin{aligned}
eq1: & 2x & +& y& -& z & =  &3 \\
eq2 ~\rightarrow eq1  - 2 \times eq2: & & &  y & - &11 z & =& -9 \\
eq3~ \rightarrow ~ eq1 + 2 \times eq3: & & &  7y &- &5z& = & 9 
\end{aligned}
$$

In this case the 2 coefficient for x in the first row is the *pivot*, and we are using 
this pivot value to convert the other two x coefficients in the rows below to zeros.

**Gaussian elimination -- step 2** remove the unknown $y$ from the last equation:

$$
\begin{aligned}
eq1: & 2x + y - z & =  3 \\
eq2: & ~~~ y - 11 z & = -9 \\
eq3 ~ \rightarrow ~ 7 \times eq2 - eq3: & ~~~ -72  z = -72 
\end{aligned}
$$

Now the pivot is the 1 coefficient for $y$ in eq2.

$$
\begin{aligned}
eq1: & 2x + y - z & =  3 \\
eq2: & ~~~ y - 11 z & = -9 \\
eq3:  &  ~~~ -72  z & = -72 \end{aligned} 
$$

This system is said to be *upper triangular*. It is also known as *row echelon* form, 
and the leading coefficients ([2, 1, -72] in this case) are known as the *pivots*.

**Gaussian elimination -- step 3** We can now use back substitution to obtain $x,y,z$.  
In this case 

$$ z = 1,$$
$$ eq2: ~~ y - 11 = -9 ~~ \implies ~~ y = 2,$$
$$ eq1: ~~ 2x +2 -1 = 3 , ~~ \implies ~~ x = 1.$$ 

### Pivoting

Consider the following system 

$$
\begin{aligned}
eq1: &  x & + & y & + & z & =  & a \\
eq2: & 2x & + & 2y & + &  5z &  = & b \\
eq3: & 4x & +&  6y  & + & 8z &  = &  c 
\end{aligned} $$

This firstly reduces to 

$$
\begin{aligned}
eq1: & x  & + & y & + & z & =  & a \\
eq2:& & & & &  3z &   = & b' \\
eq3: & & & 2y &  + & 4z &  = & c' 
\end{aligned}
$$

The problem here is that we have zero for the pivot in eq2. This can easily be switched 
into upper triangular form by switching rows two and three.

**Partial pivoting**: In general, we should be worried about both zero and very small 
pivot values, as in the latter case they will lead to division by a small value, which 
can cause large roundoff errors. So common practice is to select a row/pivot value such 
that the pivot value is as large as possible 

### Singular matrices in Gaussian Elimination

$$
\begin{aligned}
eq1: &  x & + & y & + & z & =  & a \\
eq2: & 2x & + & 2y & + &  5z &  = & b \\
eq3: & 4x & +&  4y  & + & 8z &  = &  c 
\end{aligned} $$

This firstly reduces to 

$$
\begin{aligned}
eq1: & x  & + & y & + & z & =  & a \\
eq2:& & & & &  3z &   = & b' \\
eq3:& & & & &  4z &  = & c' \end{aligned}
$$

We cannot solve this by switching rows in this case, which means that the matrix is 
singular and has no inverse

### Gaussian Elimination Rules

1. Operate on LHS and RHS (or RHSs) at the same time
2. Replace row with a sum/combination of rows
3. Work on one column at a time, choosing a pivot (leading non-zero entry in a chosen 
   row), and eliminating all other non-zero values below that
3. Switch rows to avoid zeros on the diagonal (*pivoting*)
4. If (3) does not work, zeros on the diagonal (*pivots*) indicate a singular matrix

**Computational cost**: If the number of equations $n$ is large, then a number of 
operations for gaussian elimination is $\mathcal{O}(n^3)$.

### Pseudocode

[Wikipedia](https://en.wikipedia.org/wiki/Gaussian_elimination#Pseudocode) has a 
pseudocode implementation of the gaussian elimination algorithm which is helpful to 
understand how it works:

    h := 1 /* Initialization of the pivot row */
    k := 1 /* Initialization of the pivot column */

    while h ≤ m and k ≤ n
        /* Find the k-th pivot: */
        i_max := argmax (i = h ... m, abs(A[i, k]))
        if A[i_max, k] = 0
            /* No pivot in this column, pass to next column */
            k := k+1
        else
             swap rows(h, i_max)
             /* Do for all rows below pivot: */
             for i = h + 1 ... m:
                    f := A[i, k] / A[h, k]
                    /* Fill with zeros the lower part of pivot column: */
                    A[i, k] := 0
                    /* Do for all remaining elements in current row: */
                    for j = k + 1 ... n:
                         A[i, j] := A[i, j] - A[h, j] * f
             /* Increase pivot row and column */
             h := h + 1
             k := k + 1

## Problems

1. Code a Python function that takes an 2D numpy array representing a matrix $A$, and a 
   1D array representing a RHS $b$, and returns the solution of the linear equation $Ax 
   = b$. If you wish you can assume that the matrix has an inverse. Try it out on a few 
   test matrices and check your answer using 
   [`scipy.linalg.solve`](https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.linalg.solve.html).

## Condition Number

Gaussian Elimination might still fail if $A$ is close to being singular, if a slight 
change to its values causes it to be singular. In this case simple round-off error in 
the floating point calculations can lead to zeros in the pivot positions. 

Even if the pivot value is not exactly zero, a pivot value close to zero can lead to 
large differences in the final result. In this case the matrix would be *nearly 
singular*, or *ill-conditioned*. Most linear algebra packages will include a method of 
calculating the *condition number* of a matrix, which evaluates how sensitive the 
solution is to the input values of the matrix or rhs vector. An identity matrix has a 
condition number of 1, while an exactly singular matrix has a condition number of 
infinity.

## Problems

Suppose an experiment leads to the following system of equations:

$$
\begin{aligned}
4.5 x_1 + 3.1 x_2 &= 19.249, (1)\\
1.6 x_1 + 1.1 x_2 &= 6.843.
\end{aligned}
$$

  1. Solve system (1), and then solve system (2), below, in which the data on the right 
     have been rounded to two decimal places. In each case, find the exact solution. 

$$
\begin{aligned}
4.5 x_1 + 3.1 x_2 &= 19.25, (2)\\
1.6 x_1 + 1.1 x_2 &= 6.84.
\end{aligned}
$$

  2. The entries in (2) differ from those in (1) by less than .05%. Find the percentage 
     error when using the solution of (2) as an approximation for the solution of (2).

  3. Use `numpy.linalg.cond` to produce the condition number of the coefficient matrix 
     in (1).

