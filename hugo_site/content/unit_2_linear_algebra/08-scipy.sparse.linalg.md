---
title: "Iterative solvers in Scipy"
weight: 8 
markup: "mmark"
---

Once again the best resource for Python is the [`scipi.sparse.linalg` 
documentation](https://docs.scipy.org/doc/scipy/reference/sparse.linalg.html). The 
available iterative solvers in Scipy are:

- [BIConjugate Gradient iteration 
  (BiCG)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.bicg.html#scipy.sparse.linalg.bicg)
  - Applicable to non-symmetric problems. Requires the matrix-vector product of $A$ 
    and its transpose $A^T$.
- [Quasi-Minimal Residual iteration 
  (QMR)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.qmr.html#scipy.sparse.linalg.qmr)
  - Applicable to non-symmetric $A$
  - Designed as an improvement of BiCG, avoids one of the two failure situations of 
    BiCG
  - Computational costs slightly higher than BiCG, still requires the transpose 
    $A^T$.
- [Conjugate Gradient Squared iteration 
  (CGS)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cgs.html#scipy.sparse.linalg.cgs)
  - Applicable to non-symmetric $A$
  - Often converges twice as fast as BiCG, but is often irregular and can diverge if 
    starting guess is close to solution.
  - Unlike BiCG, the two matrix-vector products cannot be parallelized.
- [BIConjugate Gradient STABilized iteration 
  (BiCGSTAB)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.bicgstab.html#scipy.sparse.linalg.bicgstab)
  - Applicable to non-symmetric $A$
  - Computational cost similar to BiCG, but does not require the transpose of $A$.
  - Simliar convergence speed as CGS, but avoids the irregular convergence properties 
    of this method
- [Conjugate Gradient iteration 
  (CG)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cg.html#scipy.sparse.linalg.cg)
  - Applicable only to symmetric positive definite $A$.
  - Speed of convergences depends on condition number

- [Generalized Minimal RESidual iteration 
  (GMRES)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html#scipy.sparse.linalg.gmres)
  - Applicable non-symmetric $A$
  - Best convergence properties, but each additional iteration becomes increasingly 
    expensive, with large storage costs.
  - To limit the increasing cost with additional iterations, it is necessary to 
    periodically *restart* the method. When to do this is highly dependence on the 
    properties of $A$
  - Requires only matrix-vector products with $A$
- 
  [LGMRES](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lgmres.html#scipy.sparse.linalg.lgmres)
  - Modification to GMRES that uses alternating residual vectors to improve 
    convergence.
  - It is possible to supply the algorithm with "guess" vectors used to augment the 
    Krylov subspace, which is useful if you are solving several very similar 
    matrices one after another.
- [MINimum RESidual iteration 
  (MINRES)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.minres.html#scipy.sparse.linalg.minres)
  - Applicable to symmetric $A$
- 
[GCROT(m,k)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gcrotmk.html#scipy.sparse.linalg.gcrotmk)

`scipy.sparse.linalg` also contains two iterative solvers for least-squares problems, 
[`lsqr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsqr.html#scipy.sparse.linalg.lsqr) 
and 
[`lsmr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsmr.html#scipy.sparse.linalg.lsmr)

### Problems

{{% notice question %}}
For this problem we are going to compare many of the different types of solvers, both 
direct and iterative, that we've looked at thus far.

Note: We will be using the Finite Difference matrix $A$ based on the two-dimensional 
finite difference approximation of the Poisson problem that we developed in the previous 
exercise.

For $N=4,8,16,32,64,128$ try the following:
1. Solve the linear systems using $\mathbf{U}_i=A^{-1} \mathbf{f}_i$ (see 
   [`scipy.linalg.inv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.inv.html) 
   and record the time this takes on a $\log$-$\log$ graph. (Omit the case $N=128$
  and note this may take a while for $N=64$.)
2. Solve the linear systems using the $\text{LU}$ and Cholesky decompositions. Plot the 
   time this takes on the same graph.
3. Now solve the systems iteratively using a conjugate gradients solver (you can use the 
   one in `scipy.linalg.sparse`, or you can code up your own). How many iterations 
   are needed for each problem? Explain the results for the right-hand-side 
   $\mathbf{f}_1$. For the right-hand-side $\mathbf{f}_2$ what is the relationship 
   between the number of iterations and $N$. How long do the computations take?
4. Repeat using the `scipy.sparse.linalg` BICGSTAB and GMRES solvers.
{{% /notice %}}
