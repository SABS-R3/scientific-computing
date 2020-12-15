---
title: "Derivative-free Methods"
teaching: 10
exercises: 0
questions:
- "What ?"
- "What is the Nelder-Mead method for non-linear optimisation?"

objectives:
- "Implement the Nelder-Mead algorithm in Python"

katex: true
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

## Finite Difference approximation

The simplest way of converting a gradient-based optimisation algorithm to a derivative 
free one is to approximate the gradient of the function using finite differences.

The Finite Difference (FD) method is based on taking a Taylor series expansion of either 
$f(x+h)$ and $f(x-h)$ (and others) for a small parameter $f$ about $x$. Consider a 
smooth function $f(x)$ then its Taylor expansion is

$$f(x+h) = f(x) + h f'(x) + \frac{h^2}{2} f''(x) + \frac{h^3}{6} f'''(x) + \frac{h^4}{24} f'''''(x) + \ldots $$

$$f(x-h) = f(x) - h f'(x) + \frac{h^2}{2} f''(x) - \frac{h^3}{6} f'''(x) + \frac{h^4}{24} f'''''(x) - \ldots $$

From this, we can compute three different \emph{schemes} (approximations) to $u'(x)$:

**Forward difference**:
$$u'(x) = \frac{u(x+h)-u(x)}{h} + O(h)$$

**Backward difference**:
$$u'(x) = \frac{u(x)-u(x-h)}{h} + O(h)$$

**Centered difference**:
$$u'(x) = \frac{u(x+h)-u(x-h)}{2h} + O(h^2)$$

Finite difference approximations are easily computed, but suffer from innacuracies which 
can cause optimisation algorithms to fail or perform poorely. As well as the error in 
the FD approximation itself (e.g. $O(h^2)$ for centered difference), the function 
evaluation itself might have some numerical or stochastic "noise". If this noise 
dominates over the (small) step size $h$, then it is entirely probable that the 
calculated steepest descent $-\nabla f(x)$ will **not** be a direction of descent for 
$f$.

### Software

It is very common that optimisation libraries provide a finite difference approximation 
to the Jacobian $\nabla f$ if it is not supplied, as is done for the gradient-based 
methods in [`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html).

More dedicated libraries can give superior approximations to the gradient, like the 
[`numdifftools`](https://numdifftools.readthedocs.io/en/latest/index.html) package. This 
library provides higher order FD approximations and *Richardson extrapolation* to 
evaluate the limit of $h \rightarrow 0$, and can calculate Jacobians and Hessians of 
user-supplied functions. 

### Problems

1. Implement a simple 2D quadratic function $f(x, y) = x^2 + y^2$, the 2D [Rosenbrock
function](https://en.wikipedia.org/wiki/Rosenbrock_function), and the 2D [Rastrigin 
function](https://en.wikipedia.org/wiki/Rastrigin_function), and visualise these 
functions over the domain $-5 < x < 5$ and $-5 < y < 5$

2. Use [`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html) to
   minimise each of these functions, starting at $x = y = 2.5$, using the following 
   methods:
   - Nelder-Mead Simplex
   - Conjugate Gradient
   - BFGS Quasi-Newton
   - Newton-CG
   - SHG Global Optimisation

   In each case perform the optimisation with and without a user-supplied jacobian and 
   evaluate the effect on the number of evaluations of the function $f$ required to 
   converge to the opimium.

3. Try using [`numdifftools`](https://numdifftools.readthedocs.io/en/latest/index.html) 
   to calculate the jacobian. See if you can improve on the number of function 
   evaluations required to converge to an minimum.


## Nelder-Mead method

The [Nelder-Mead method](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) is 
popular and implementations exist in many optimisation software libraries. It is based 
on the idea of a simplex in parameter space of dimension $n$, which is formed from the  
convex hull of $n + 1$ points in $\mathcal{R}^n$. These points $x_i$ are ordered 
according to their function value so that

$$
f(x_1) \le f(x_2) \le \cdots \le f(x_{n+1})
$$

For each iteration of the algorithm, there are five different points of interest, the 
first of which is the centroid of the $n$ points with the lowest $f$ values

$$
\bar{x} = \sum_{i=1}^n x_i
$$

The other four points are defined by considering the line joining $\bar{x}$ and the 
point with the highest $f$ value $x_{n+1}$

$$
\bar{x}(t) = \bar{x} + t(x_{n+1} - \bar{x}
$$

The four points are the *reflection*, *expanding*, the *outside contraction* and *inside 
constraction* points, given by $\bar{x}(-1)$, $\bar{x}(-2)$, $\bar{x}(1/2)$, and 
$\bar{x}(-1/2)$ respectivly.

The Nelder-Mead algorithm tries to replace $x_{n+1}$ by reflecting, expanding, or 
contracting the simplex to one of these points. If it cannot find a better point, then 
all the vertices on the simplex are shrunk towards the best vertex $x_1$. 

[Scholarpedia](http://www.scholarpedia.org/article/Nelder-Mead_algorithm) and 
[Wikipedia](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) provide diagrams 
and pseudocode of the Nelder-Mead algorithm, as does Chapter 9.5 of the Nocedal and 
Write textbook given below

### Problems

- Code up the Nelder-Mead algorithm and compare its performance against the steepest 
  descent, Newton and dogleg algorithms you did in the last lesson. You can evaluate 
  them on the 2D quadratic function $f(x, y) = x^2 + y^2$, the 2D [Rosenbrock
function](https://en.wikipedia.org/wiki/Rosenbrock_function) or on one of many different 
[optimisation test 
functions](https://en.wikipedia.org/wiki/Test_functions_for_optimization)

## Other Reading

- Nelder, John A.; R. Mead (1965). "A simplex method for function minimization". 
  Computer Journal. 7 (4): 308–313. doi:10.1093/comjnl/7.4.308.
- Numerical optimization by Nocedal, Jorge; Wright, Stephen J., 1960-, Chapter 9
