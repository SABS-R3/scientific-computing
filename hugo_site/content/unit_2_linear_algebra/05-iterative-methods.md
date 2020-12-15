---
title: "Iterative methods"
weight: 5 
markup: "mmark"
---

Previously we have discussed *direct* linear algebra solvers based on decompositions of 
the original matrix $A$. The amount of computational effort required to achieve these 
decomposisions is $\mathcal{O}(n^3)$, where $n$ is the number of rows of a square 
matrix. They are therefore unsuitable for the large, sparse systems of equations that 
are typically encountered in scientific applications. An alternate class of linear 
algebra solvers are the *iterative* methods, which produce a series of *approximate* 
solutions $x_k$ to the $A x = b$ problem. The performance of each algorithm is then 
based on how quickly, or how many iterations $k$ are required, for the solution $x_k$ to 
converge to within a set tolerance of the true solution $x$.

