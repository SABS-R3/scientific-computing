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

{{% notice question %}}
Use the following methods from 
[`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html) to minimize 
the 2D [Rosenbrock 
function](https://en.wikipedia.org/wiki/Rosenbrock_function):
  - Nelder-Mead Simplex
  - Conjugate Gradient
  - BFGS Quasi-Newton
  - Newton-CG
  - SHG Global Optimisation

Either use $x_0 = (−1.2, 1)^T$ as the starting point, or experiment with your own.

In each case perform the optimisation with and without a user-supplied jacobian and 
evaluate the effect on the number of evaluations of the function $f$ required to 
converge to the optimum. Optional: You can take the derivatives by hand, or use 
automatic differentiation via the [`autograd`](https://github.com/HIPS/autograd) or 
[`JAX`](https://github.com/google/jax) packages
{{% /notice %}}


{{% expand "Expand for solution" %}}
{{% notice solution %}}
```python
import autograd.numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, shgo
from autograd import grad

def convex_function(x):
    return np.sum(np.array(x)**2, axis=0)

def rosenbrock(x):
    a = 1.0
    b = 100.0
    return (a-x[0])**2 + b*(x[1]-x[0]**2)**2

def rastrigin(x):
    A = 10.0
    n = len(x)
    ret = A*n
    for i in range(n):
        ret += x[i]**2 - A * np.cos(2*np.pi*x[i])
    return ret

def visualise_functions():
    nx, ny = (100, 100)
    x = np.linspace(-5, 5, nx)
    y = np.linspace(-5, 5, ny)
    xv, yv = np.meshgrid(x, y)
    eval_convex = convex([xv, yv])
    eval_rastrigin = rastrigin([xv, yv])
    eval_rosenbrock = rosenbrock([xv, yv])

    plt.contourf(x, y, eval_convex)
    plt.colorbar()
    plt.show()
    plt.clf()
    plt.contourf(x, y, eval_rosenbrock)
    plt.colorbar()
    plt.show()
    plt.clf()
    plt.contourf(x, y, eval_rastrigin)
    plt.colorbar()
    plt.show()

def optimize(function, method, autodiff):
    eval_points_x = []
    eval_points_y = []
    def fill_eval_points(xk):
        eval_points_x.append(xk[0])
        eval_points_y.append(xk[1])

    x0 = np.array([2.5, 2.5])
    if autodiff:
        jac = grad(function)
    else:
        jac = None

    if method == 'shgo':
        bounds = [(-10, 10), (-10.0, 10.0)]
        res = shgo(function, bounds, callback=fill_eval_points,
                    options={'disp': True})
    else:
        res = minimize(function, x0, method=method, callback=fill_eval_points,
                        jac = jac,
                    options={'disp': True})

    nx, ny = (100, 100)
    x = np.linspace(-5, 5, nx)
    y = np.linspace(-5, 5, ny)
    xv, yv = np.meshgrid(x, y)
    eval_function = np.empty_like(xv)
    eval_function = function([xv, yv])

    print(function.__name__)
    plt.clf()
    plt.contourf(x, y, eval_function)
    plt.plot(eval_points_x, eval_points_y, 'x-k')
    plt.colorbar()
    neval = res.nfev
    try:
        neval += res.njev
    except:
        print('no njev')
    try:
        neval += res.nhev
    except:
        print('no nhev')
    plt.title('iterations: {} evaluations: {}'.format(res.nit, neval))

    if autodiff:
        ext = '-auto.pdf'
    else:
        ext = '.pdf'

    plt.savefig(function.__name__ + '-' + method + ext)


if __name__ == '__main__':
    for f in [convex_function, rosenbrock, rastrigin]:
        for m in ['shgo','nelder-mead', 'cg', 'bfgs', 'newton-cg']:
            for a in [False, True]:
                if m == 'newton-cg' and a == False:
                    continue
                if m == 'shgo' and a == True:
                    continue
                optimize(f, m, a)
```
{{% /notice %}}
{{% /expand %}}

