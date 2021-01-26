---
title: "Project description"
weight: 1 
markup: "mmark"
---

The development of efficient and accurate general purpose algorithms for scientific 
computing is one of the primary motivators for researchers specialising in this field. 
However, an often overlooked aspect of scientific computing research is the proper 
*implementation* and *publication* of these algorithms in order for other researchers to 
make use of an algorithm. This is even more important for those working in more applied 
fields, whether in biomedicine or engineering, who might not have sufficient specialised 
technical knowledge to even read, for example, the latest literature in advanced MCMC 
sampling methods.

Taking an even broader view to research in general. While it is important to pursue a 
rigorous research agenda in any specialist subject and generate important research 
outcomes, of equal importance is the proper *publication* of this work. How is it 
properly communicated to both those working inside and outside your field, and/or to the 
general public.

In the last three days of this course, you will pick one (or more than one) scientific 
computing algorithm of your choice and publish both an implementation of the algorithm 
as well as an explanation of its theoretical and/or practical properties (i.e. how it 
works). Which algorithm you choose, your implementation and how you publish it is 
entirely up to you.

For example, you might choose to:
- implement the algorithm as a Jupyter notebook, interleaving markdown text and code to 
  explain the algorithm.
- implement the algorithm as a Python library on Github, with a README.md, Sphinx 
  documentation and/or separate pdf document describing both how to use the library and 
  a mathematical description of the algorithm.
- implement the algorithm in a fork of an existing package like Scipy, or 
  [Pints](https://github.com/pints-team/pints)
- Produce an explanatory video or animation describing the algorithm, or which is used 
  within a written explanation. Note that Matplotlib is able to generate animations, or 
  you can [use](https://observablehq.com/@d3/learn-d3-animation) a more purpose built 
  like [D3](https://d3js.org/).

For inspiration, here is a link to the online journal 
[Distill](https://distill.pub/about/), a machine learning journal that specialises in 
articles publishing clear explanations of existing techniques using web technologies 
like D3.

You can pick any scientific computing algorithm you wish, you do not have to be confined 
to the topics discussed during this course (e.g. perhaps you want to focus on a
technique you might use for your DPhil). Below are a few suggestions to get you started. 
Note that you can also choose one of the algorithms already presented during this course 
like QR decomposition, and explain it better than I did!

### Linear algebra

- Eigenvalue solvers like Power iteration method. See list of iterative algorithms on 
  [wikipedia](https://en.wikipedia.org/wiki/Eigenvalue_algorithm#Iterative_algorithms)
- Singular Value Decomposition
- Matrix-matrix multiplication (and numerous methods to speed it up)
- Rank-revealing QR decomposition

### ODE solvers 

- Explicit
  - Single step methods such as Runge-Kutta 
  - Multistep methods like Adams-Bashforth

- Implicit
  - Single step like implicit runge-kutta
  - Multistep like Backwards Difference Formula or Adams–Moulton

### Non-linear Optimisation
- Lots! https://en.wikipedia.org/wiki/List_of_algorithms#Optimization_algorithms

Some examples:

- Nonlinear optimisation methods like Gauss–Newton, or BFGS
- Stochastic gradient descent popular in ML (extensions AdaGrad, RMSProp or Adams)
- Monte Carlo methods like simulated annealing

### MCMC samplers

- [Here](https://pints.readthedocs.io/en/stable/mcmc_samplers/index.html) is the list of 
  samplers available in the Pints library, for example:
  - Adaptive Covariance MCMC
  - Hamiltonian Monte Carlo (HMC)
  - Population MCMC


{{% notice question %}}

Pick a scientific computing algorithm and provide both an re-usable implementation of 
the algorithm, as well as some form of publication that provides a clear explanation of 
how it works, and/or any other theoretical or practical properties it has (e.g. its 
execution time/memory complexity, error bounds etc.). 

Please hand in your project, in whatever format is most suitable, by 17:30 on Friday 
12th Feb 2021.

{{% /notice %}}
