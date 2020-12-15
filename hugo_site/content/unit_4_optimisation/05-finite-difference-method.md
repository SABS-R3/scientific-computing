---
title: "Finite difference method"
weight: 4 
markup: "mmark"
---

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


