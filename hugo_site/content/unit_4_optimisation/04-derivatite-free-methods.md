---
title: "Derivative-free methods"
weight: 4 
markup: "mmark"
---

The line search and trust region methods introduced in the previous lesson all required 
that the user be able to calculate the gradient of the function $\nabla f$. However, in 
many cases the gradient is either not available or too error-prone to be of use. For 
example, the function $f$ might be only available as a compiled executable or the result 
of a physical experiment. The model might be stochastic, or the gradient evaluation 
might be noisy due to numerical innacuracies, or of sufficiently complexity that the 
gradient is either unknown or too expensive to compute. 

Here we describe two of the most common methods for derivative-free optimisation, using 
a finite difference approximation to approximate the derivative, and the [Nelder-Mead 
algorithm](https://doi.org/10.1093/comjnl/7.4.308), which is a Simplex search method. 
However, there are a large number of derivative-free methods, ranging from the classical  
[*Direct Search 
methods*](https://www.sciencedirect.com/science/article/pii/S0377042700004234) like 
*Pattern search*, *Simplex search*, *Rosenbrock'* or *Powell's* methods. Then there are 
emulator or model -based methods that build up a model of the function $f$ and minimise 
that using a gradient-based method, a powerful example of this class of methods is 
[Bayesian 
Optimisation](http://papers.nips.cc/paper/4522-practical-bayesian-optimization). Many 
global optimsiation algorithms are derivative-free, including *randomised algorithms* 
such as [Simulated Annealing](https://science.sciencemag.org/content/220/4598/671), and 
*evolution-based* strategies such as the popular [Covariance matrix adaptation evolution 
strategy (CMA-ES)](https://arxiv.org/abs/1604.00772), or *swarm algorithms* inspired 
from bees/ants like [Particle Swarm 
Optimisation](https://doi.org/10.1109/ICNN.1995.488968).
