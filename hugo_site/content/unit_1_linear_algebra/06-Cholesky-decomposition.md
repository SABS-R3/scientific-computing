---
title: "Cholesky decomposition"
date: "2014-04-01"
markup: "mmark"
weight: 6 
---

## Cholesky decomposition

*Symmetric positive definite* matrices are a very special type of matrix that often 
arise in practice. From a computational point of view, this class of matrix is very 
attractive because it is possible to decompose a symmetic positive definite matrix $A$ 
very efficiently into a single lower triangular matrix $G$ so that $A = GG^T$. 

A matrix $A$ is positive definite if $x^T A x > 0$  for any nonzero $x \in \mathbb{R}$. 
This statement by itself is not terribly intuitive, so lets look at also look at an 
example of a $2 \times 2$ matrix

$$
A = \left(\begin{matrix}
a_{11} & a_{12} \\
a_{21} & a_{22}
\end{matrix}\right)
$$

If $A$ is symmetic positive definite (SPD) then

$$
\begin{aligned}
x &= (1, 0)^T \Rightarrow x^T A x = a_{11} > 0 \\
x &= (0, 1)^T \Rightarrow x^T A x = a_{22} > 0 \\
x &= (1, 1)^T \Rightarrow x^T A x = a_{11} + 2a_{12} + a_{22} > 0 \\
x &= (1,-1)^T \Rightarrow x^T A x = a_{11} - 2a_{12} + a_{22} > 0 \\
\end{aligned}
$$

The first two equations show that the diagonal entries of $A$ must be positive, and 
combining the last two equations imply $|a_{12}| \le (a_{11} + a_{22}) / 2$, that is 
that the matrix has much of its "mass" on the diagonal (note: this is *not* the same as 
the matrix being diagonally dominant, where $|a_{ii}| > \sum_{i=1...n,j \ne i} 
|a_{ij}|$). These two observations for our $2 \times 2$ matrix also apply for a general 
$n \times n$ SPD matrix. One of the very nice consequences of this "weighty" diagonal 
for SPD matrices is that it precludes the need for pivoting.

It can be shown that if $A$ is a SPD matrix, then the $LDL^T$ decomposition exists and 
that $D = \text{diag}(d_1, ..., d_n)$ has positive diagonal entries. Therefore, it is 
straightforward to see that $LDL^T$ = $GG^T$, where $G = L \text{diag}(\sqrt{d_1}, ..., 
\sqrt{d_n})$. The decomposition $A = GG^T$ is known as the cholesky decomposition and 
can be efficiently constructed in $n^3 / 3$ flops. There are a number of algorithms to 
construct this decomposition, and both the [wikipedia 
entry](https://en.wikipedia.org/wiki/Cholesky_decomposition) and Chapter 4.2 of the 
Matrix Computations textbook by Golub and Van Loan gives a number of different varients.

Note that a $LDL$ decomposition can also be used to calculate a cholesky decomposition, 
and this could be more efficient approach since (a) the SPD structure means that we can 
neglect pivoting in the $LDL$ decomposition, and (b) the $LDL$ decomposition does not 
requiring taking the square root of the diagonal elements. 

### Other Reading

- Golub, G. H. & Van Loan, C. F. Matrix Computations, 3rd Ed. (Johns Hopkins University 
  Press, 1996). Chapter 4.2
- https://en.wikipedia.org/wiki/Cholesky_decomposition

### Software

- 
  [`scipy.linalg.cholesky`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cholesky.html)
- 
[`scipy.linalg.cho_factor`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_factor.html)
- 
[`scipy.linalg.cho_solve`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve.html)
- 
[`scipy.linalg.cholesky_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cholesky_banded.html)
- 
[`scipy.linalg.cho_solve_banded`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.cho_solve_banded.html)

## Problems

Imagine that we wanted to sample an array of values $x_i$, for $i = 1...n$, where each 
value is sampled from an independent Normal distribution with standard deviation 
$\sigma$

 $$x_i \sim \mathcal{N}(0, \sigma)$$

 This could be achieved, for example, by sampling from a Normal distribution with unit 
 standard deviation, a function that typically exists in any computer language, then 
 multiplying by $\sigma$

 $$x_i = \sigma \eta$$

 where $\eta \sim \mathcal{N}(0, 1)$

 Now imagine that instead of an independent Normal distribution you wish to sample 
 $\mathbf{x} = [x_1, x_2, ..., x_n]$ from a multivariate Normal distribution with some 
 covariance matrix $\Sigma$

 $$\mathbf{x} \sim \mathcal{N}(\mathbf{0}, \Sigma)$$

 We can achive this in practice by using the Cholesky decomposition. A covariance 
 matrix is a symmetic positive semidefinite matrix (i.e. $x^T \Sigma x \ge 0$}, and 
 therefore can be decomposed into  $\Sigma = LL^T$. We can then draw a sample from 
 $\mathcal{N}(\mathbf{0}, \Sigma)$ by scaling an independently generated random vector 
 by $L$

 $$\mathbf{x} = L \mathbf{\eta}$$

 where each element of the vector $\eta$ is $\eta_i \sim \mathcal{N}(0, 1)$.

 Write Python code to randomly sample an n-dimensional vector $x$ from 
 
 1. an independent Normal distribution with variance $\sigma^2$

 2. a multivariate normal distribution using a covariance matrix $\Sigma_{ij} = \exp[(i- 
    j)^2/ \sigma^2]$

 3. a multivariate normal distribution with $\Sigma = \sigma^2 I$. Show that this 
    algorithm reduces to that used for (1).
