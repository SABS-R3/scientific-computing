---
title: "Course Overview"
date: 2020-11-26T16:52:22Z
draft: false
pre: "1. "
---

Scientific computing is a broad term, and encompasses everything from a small script to 
plot a dataset for publication, to a high performance numerical PDE solver that is 
running on a supercomputer. Underlying and supporting all these implementations are two 
important elements, without which the field of scientific computing could not exist. The 
first is the development and publishing of **general purpose algorithms** for solving 
common problems, such as algorithms for the solution of ordinary differential equations 
(e.g. Euler or Runge-Kutta methods), or to minimise a non-linear function (e.g. Newtons 
method or Nelder–Mead). The second is the development of software **frameworks** for 
ease of implementation of these algorithms, on both desktop and high performance 
computers. The most important category of software frameworks is the ubiquitous Linear 
Algebra libraries for representing and operating on matrices and vectors. The BLAS APIs 
that were developed in the 1970s were hugely important in developing a common set of 
functions for linear algebra that could be backed by processor-specific implementations 
for high execution speed. After these came higher level libraries such as ATLAS and 
LAPACK, leading eventually to mass market linear algebra applications like MATLAB.

This course attempts to give an introduction to a useful set of algorithms for 
scientific computing, along with practical skills in implementation of these algorithms 
using linear algebra and other scientific libraries. On the algorithm side, we will 
focus on:

- Linear algebra direct solvers and matrix decompositions
- Sparse matrices and iterative solvers
- Solution of Ordinary Differential Equations
- Non-linear Optimisation
- Bayesian Inference and MCMC sampling

For a two week course it is impractical to give a proper treatment of all of these 
topics. Therefor we will attempt to touch on some of the more important algorithms 
within each topic and provide links to further reading to broaden your knowledge.

On the software side, we will focus exclusively on implementation in 
[Python](https://www.python.org/), using the [Numpy](https://numpy.org/) library for 
multidimensional arrays, and the [Scipy](https://www.scipy.org/) library for both sparse 
matrices and many, many other general purpose scientific computing algorithms. 


