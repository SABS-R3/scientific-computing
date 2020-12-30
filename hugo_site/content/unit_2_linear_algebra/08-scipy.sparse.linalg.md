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
- [BIConjugate Gradient STABilized iteration 
  (BiCGSTAB)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.bicgstab.html#scipy.sparse.linalg.bicgstab)
- [Conjugate Gradient iteration 
  (CG)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cg.html#scipy.sparse.linalg.cg)
- [Conjugate Gradient Squared iteration 
  (CGS)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cgs.html#scipy.sparse.linalg.cgs)
- [Generalized Minimal RESidual iteration 
  (GMRES)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gmres.html#scipy.sparse.linalg.gmres)
- 
  [LGMRES](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lgmres.html#scipy.sparse.linalg.lgmres)
- [MINimum RESidual iteration 
  (MINRES)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.minres.html#scipy.sparse.linalg.minres)
- [Quasi-Minimal Residual iteration 
  (QMR)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.qmr.html#scipy.sparse.linalg.qmr)
- 
  [GCROT(m,k)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.gcrotmk.html#scipy.sparse.linalg.gcrotmk)

`scipy.sparse.linalg` also contains two iterative solvers for least-squares problems, 
[`lsqr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsqr.html#scipy.sparse.linalg.lsqr) 
and 
[`lsmr`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.lsmr.html#scipy.sparse.linalg.lsmr)

### Problems

Note: based on the Finite Difference matrix $A$ in a previous exercise:

For $N=4,8,16,32,64,128$ try the following:
1. Solve the linear systems using $\mathbf{U}_i=A^{-1} \mathbf{f}_i$ (see 
   [`scipy.linalg.inv`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.inv.html) 
   and record the time this takes on a $\log$-$\log$ graph. (Omit the case $N=128$
  and note this may take a while for $N=64$.)
2. Solve the linear systems using the $\text{LU}$ and Cholesky decompositions. Plot the 
   time this takes on the same graph.
3. Now solve the systems iteratively using a conjugate gradients solver (you can use the 
   one in `scipy.linalg.sparse`, or you can code up your own). How many iterations are 
   needed for each problem? Explain the results for the right-hand-side $\mathbf{f}_1$. 
   For the right-hand-side $\mathbf{f}_2$ what is the relationship between the number of 
   iterations and $N$. How long do the computations take?
4. Repeat using the `scipy.sparse.linalg` BICGSTAB and GMRES solvers.
