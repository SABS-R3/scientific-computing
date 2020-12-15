---
title: "Nonlinear Optimisation"
teaching: 10
exercises: 0
questions:
- "What is non-linear optimisation?"
- "What are the challenges in non-linear optimisation"
- "What algorithms exist for solving these problems?"
- "What are the line search and trust region methods?"

objectives:
- "Explain the main concepts of non-linear optimisation"
- "Describe two main classes of optimisation methods, line searcha and trust region"
katex: true
markup: "mmark"
---

## Mathematical formulation

Optimisation aims to find the minimum (or equivilently the maximum) of some *objective*, 
or *loss* function $f$, given a set of $n$ parameters $\theta$

$$
\min_{\theta \in \mathcal{R}^n} f(\theta)
$$

We might also have a set of *constraints*, for example a parameter might be required to 
be non-negative (e.g. a concentration or population number). These are often written as 
a set of equality $\mathcal{E}$ and inequality $\mathcal{I}$ constraints

$$
\min_{\theta \in \mathcal{R}^n} f(\theta) \text{ subject to } \begin{cases}
c_i(x) = 0, & i \in \mathcal{E} \\
c_i(x) \ge 0, & i \in \mathcal{I} \end{cases}
$$

Or these might simply be defined as *bounds* in parameter space that restrict the 
minimisation to a given domain $\Omega \in \mathcal{R}^n$

$$
\min_{\theta \in \Omega} f(\theta)
$$


## Useful terms 

*Modelling* is the process of defining the objective function $f$, the parameters of 
interest $\theta$, and the constraints. The algorithms for performing the minimisation 
fall under the field of optimisation. Sub-fields of this are concerned with the 
minimisation of discrete function, often called *integer programming*. Confusingly, it 
is common to see the terms "optimisation" and "programming" used interchangably, as thie 
latter term was coined before the 1940s, and does not refer to computer software 
programming at all.

If the function $f$ is linear, then there are specific algorithms for this class of 
problem that fall under the topic of *linear programming*, or *linear optimisation*. The 
more general problem of a non-linear $f$ is termed *non-linear programming*, or 
*non-linear optimisation*. If a set of equality and/or inequality constraints are needed 
then algorithms that deal with these fall under the topic of *constrained* optimisation.

An important distinction when looking at optimisation problems is the notion of *global* 
versus *local* optimisation. The latter finds a point in parameter space $\theta_m$ that 
has a function value $f(\theta_m)$ greater than the surrounding points, but might not 
neccessarily by the global minimum. These algorithms are often initialised to a point 
that is near to the minima of interest. The more general problem of global optimisation 
is significantly more difficult as it requires that the optimisation be robust to 
finding and rejecting such local minima. For a function that is *convex*, then local and 
global minimisation are the same, which is very advantagous since local minimisation 
algorithms are often both faster and often have more guarentees of convergence. The 
function $f$ is a convex function if its domain $\Omega$ is a convex set, and for any 
two points $\theta_x$ and $\theta_y$:

$$
f(\alpha \theta_x + (1 - \alpha) \theta_y ) \le \alpha f(\theta_x) + (1 - \alpha) 
f(\theta_y), \text{ for all } \alpha \in [0, 1].
$$

The term *convex programming* is used to describe the case of contrained optimisation 
where $f$ is convex, the equality constraints are linear and the inequality contraints 
are concave.

## Non-linear optimisation and Root-Finding

Non-linear optimisation is closely related to finding the roots, or zeros, of a 
function. This can be seen easily by considering the fact that at each local minima or 
maxima of the function the value of the gradient of $f$ is zero, i.e. $\nabla f = 0$. 
Therefore finding a local minima or maxima of $f$ corresponds to finding the zeros of 
the function $g = \nabla f$.

### Gradient descent and line search methods

One of the simplest local optimisation algoriths is *gradient descent*. It is 
initialised at some point in parameter space $a_0$, and at each iteration the function 
$f(x)$ is reduced by following the direction of *steepest descent* $-\nabla f(a)$ 

$$
a_{n+1} = a_n - \gamma \nabla f(a_n)
$$

This is an example of an imporant class of algorithms called the *line search* methods. 
These algorithms choose a *search direction*  $p_k$ at each iteration $k$, and search 
along the 1D line from the initial point $a_k$ to a new point 

$$
a_{k+1} = a_k + \alpha p_k
$$ 

with a lower function value. The problem at each iteration becomes a one-dimensional 
optimisation problem along $p_k$ to find the optimal value of $\alpha$. Each line search 
algorithm is thus defined on how it chooses both the search direction $p_k$ and the 
optimal $\alpha$.

### Plataus with low gradient 

An obvious downside to simple gradient descent can be seen for functions which have 
regions of zero or small gradients, or plataus. Here a gradient descent algorithm with a 
constant $\gamma$ will proceed very slowly, if at all. This motivates another important 
line search algorithm, *Newtons method*.

The Newtons direction $p^N_k$ can be derived by considering the second-order Taylor 
expansion of the function $f(x)$

$$
f(a_k + p) \approx f(a_k) + p^T \nabla f(a_k) + \frac{1}{2} p^T \nabla^2 f(a_k) = 
m_k(p).
$$

We find the value of $p$ that minimises $m_k(p)$ by setting the derivative of $m_k$ to 
zero, leading to 

$$
p_k^N = - (\nabla^2 f(a_k))^{-1} \nabla f(a_k)
$$

Unlike the steepest descent, Newtons method has a natural step length $\alpha \approx 
1$, which is suitable for a wide variety of problems and can quickly cross areas of low 
gradient. Natually, since the algorithm is based on a *second-order* approximation of 
the function $f$, it works better if this approximation is reasonably accurate.

Newtons method can be used as long as the inverse of the second derivative of the 
function $(\nabla^2 f(a_k))^{-1}$, exists (e.g. it will always exist for a positive 
definite $\nabla^2 f$). However, even when this inverse does exist it is possible that 
the direction $p^N_k$ does not satisfy the descent condition $f(a_k + \alpha p^N_k) < 
f(a_k)$ (or equivilently $\nabla f(a_k)^T p^N < 0$), so many modifications to Newtons 
methods, falling under a class of methods called *Quasi-Newton* methods, have been 
proposed to satisfy this descent condition. 

Quasi-Newton methods do not require the (often onerous) calculation of the hession 
$\nabla^2 f(x)$ like Newtons, instead they form an approximation to the hessian $B_k 
\approx \nabla^2 f(a_k)$ that is updated at each step using the information given by the 
gradient evaluations $\nabla f(a_k)$. Two popular methods of performing this update are 
the *symmetric-rank-one* (SR1), and the *Broyden, Fletcher, Goldfarb, and Shanno, 
(BFGS)* formula. Once the approximation $B_k$ is formed then the search direction is 
calculated via

$$
p_k = -B_k^{-1} \nabla f(a_k)
$$

For more details of other line search methods, and other concepts such as the *Wolfe 
conditions* for calculating the step length $\alpha$, please see Chapter 3 of the 
Nocedal and Wright textbook, or in the other textbooks listed at the end of this 
lession. Finally, it should be noted that the *conjugate gradient* method can also be 
used for non-linear optimisation, where the search direction is given by

$$
p_k = -\nabla f(a_k) + \beta_k p\_{k-1}
$$

### Step length

In line search methods, choosing the step length $\alpha_k$ is a non-trivial task. 
Ideally we would want to chose $\alpha_k$ to minimise the function along the 
one-dimensional search direction $p_k$. That is, we wish to minimise

$$
\phi(\alpha_k) = f(a_k + \alpha_k p_k),\text{ }\alpha_k > 0.
$$

In general this is too expensive to do this minimisation exactly, so approximate methods 
are used so that multiple trial $\alpha_k$ values are trialed, stopping when a candidate 
is found that satisfies a set of *conditions*. There are two main conditions used, the 
*Wolfe conditions* and the *Goldstein* conditions.

The two Wolfe conditions are the *sufficient decrease* condition, which ensures that the 
reduction in the function value is proportional to the step length $\alpha_k$ and the 
gradient in the direction of the step 

$$
f(a_k + \alpha_k p_k) \le f(a_k) + c_1 \alpha_k \nabla f(a_k)^T p_k.
$$

The second Wolfe condition is the *curvature* condition, which prevents unacceptibly 
short steps by ensuring that the slope of $\phi$ is greater than some constant $c_2$ 
times the initial slope $\phi'(0)$

$$
\nabla f(a_k + \alpha_k p_k)^T p_k \ge c_2 \nabla f(a_k)^T p_k,
$$

where $c_2 \in (c_1, 1)$. Typical values are $c_1 = 10^{-4}$ and $c_2 = 0.9$. The 
*strong Wolfe* conditions restrict the gradient $\phi'$ to be small, so as to exclude 
points that are too far from stationary points of $\phi$

$$
f(a_k + \alpha_k p_k) \le f(a_k) + c_1 \alpha_k \nabla f(a_k)^T p_k.
$$
$$
|\nabla f(a_k + \alpha_k p_k)^T p_k| \ge c_2 |\nabla f(a_k)^T p_k|,
$$

The Goldstein conditions are similar in spirit to the Wolfe condtitions, and are formed 
from the two inequalities

$$
f(a_k) + (1 - c) \alpha_k \nabla f(a_k)^T p_k \le f(a_k + \alpha_k p_k) \le f(a_k) + c 
\alpha_k \nabla f(a_k)^T p_k.
$$

with $0 < c < 1/2$. The first inequality prevents small step sizes while the second is 
the same sufficient decrease condition as in the Wolfe conditions. The Goldstein 
conditions are often used in Newton-type methods but for quasi-Newton methods the Wolfe 
conditions are prefered. The diagrams from the text by Nocedal and Wright
illustrate the two conditions

![](/conditions.jpg)

Algorithms for choosing candidate step size values $\alpha_k$ can be complicated, so we 
will only mention here one of the simplest, which is the *backtracking* method. This 
approach implicitly satisfies the condition on too small $\alpha_k$, and only repeatedly 
test for the common sufficient decrease condition that appears in both the Wolfe and 
Goldstein condtitions.

Choose $\bar{\alpha} > 0$, $\rho \in (0, 1)$, $c \in (0, 1)$ \\
$\alpha_k := \bar{\alpha}$ \\
**repeat** until $f(a_k + \alpha_k p_k) \le f(a_k) + c \alpha_k \nabla f(a_k)^T p_k$ \\
&nbsp;&nbsp; $\alpha_k := \rho \alpha_k$ \\
**end repeat**


### Problems 

1. Program the steepest descent and Newton algorithms using the backtracking line 
   search, Algorithm 3.1. Use them to minimize the [Rosenbrock 
   function](https://en.wikipedia.org/wiki/Rosenbrock_function). Set the initial step 
   length $\alpha_0 = 1$ and print the step length used by each method at each 
   iteration. First try the initial point $x_0 = (1.2, 1.2)^T$ and then the more 
   difficult starting point $x_0 = (−1.2, 1)^T$. 
2. Plot the function surface using `matplotlib` and overlay the line search segments so 
   you can visualise the progress of your algorithm. 
3. Repeat (1) and (2) above using the line seach implemented in Scipy 
   [`scipy.optimize.line_search`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.line_search.html), 
   which uses the strong Wolfe conditions.

### Saddle points

Saddle point pose a particular challenge in non-linear optimisation, particularly in 
higher dimensions. The plots below show two examples of saddle points in two dimensions. 
Like local minima and maxima, these are stationary points where the gradient of the 
function is zero $\nabla f = 0$, but where the value of the function rises along certain 
directions and reduces along others (left plot). An alternative type of saddle point 
arises when the hessian is singular, and are characterised by a plateau around the 
stationary point, like the [monkey saddle](https://en.wikipedia.org/wiki/Monkey_saddle) 
depicted in the plot to the right. 

![Two examples of saddle points](/saddle.svg)

Gradient descent perform poorly, with very slow convergence near saddle points, and 
Newtons methods tend to be 
[attracted](https://www.offconvex.org/2016/03/22/saddlepoints/) and 
[trapped](https://arxiv.org/abs/1405.4604) in saddle points, due to the presence of 
negative eigenvalues in the hessian matrix. This can be avoided by approximating the 
hessian with a positive definite matrix $B_k \approx \nabla^2 f$, like some quasi-Newton 
methods (i.e. BFGS). Alternativly, the class of *trust region methods* restate the 
optimisation problem as an sequence of optimisations of a second order approximation to 
$f$ in a local *trust-region* surrounding the current point $a_k$. The exact solution to 
each of these subproblems can be shown to be $(\nabla^2 f(a_k) + \lambda I)^{-1} \nabla 
f(a_k)$, where $\lambda$ is chosen to be large enought to make $(\nabla^2 f(a_k) + 
\lambda I)$ positive definite. Therefore, by design the trust-region methods aim avoid 
this problem of negative eigenvalues. 

### Trust region methods

Like many line search methods, trust region methods also use the second order Taylor 
expansion of $f$ around $a_k$

$$
f(a_k + p) \approx f(a_k) + g_k^T p + \frac{1}{2} p^T B_k p = m_k(p)
$$

where $g_k = \nabla f(a_k)$, $B_k$ is an approximation to the hessian matrix $B_k 
\approx \nabla^2 f(a_k)$ or the hessian itself $B_k = \nabla^2 f(a_k)$. Trust region 
methods aim to find the $p$ that minimises $m_k$ in a local trust region  $||p|| < 
\Delta_k$ around the current point $a_k$, where $\Delta_k$ is the trust region radius. 

Solving the minimisation given above is normally done approximately, with different 
trust region methods varying how the approximation is achieved. Choosing the 
trust-region radius is fundamental to this class of methods, and is done by comparing 
the actual to the predicted reduction in the function value

$$
\rho_k = \frac{f(a_k) - f(a_k + p_k)}{m_k(0) - m_k(p_k)}.
$$

Since $m_k(0) - m_k(p_k)$ is always positive, if $\rho_k$ is negative then the the 
actual function value is increasing, the step is rejected and the trust region radius 
$\Delta_k$ is decreased in order to improve the approximate model $m_k$. If $\rho_k$ is 
positive but much smaller than one then we do not alter $\Delta_k$. If $\rho_k$ is close 
to or greater than 1 we can be confident in our model and thus increase $\Delta_k$. The 
general algorithm for a trust region method (reproduced from the text by Nocedal and 
Wright cited below) is:

Given $a_0$, $\hat{\Delta} > 0$, $\Delta_0 \in (0, \hat{\Delta})$, and $\nu \in [0, 
\frac{1}{4})$: \\
**for** $k = 0, 1, 2, ...$ \\
&nbsp;&nbsp; Obtain $p_k$ by (approximatly) minimising $m_k(p)$ where $||p|| < \Delta_k$ 
\\
&nbsp;&nbsp; $\rho_k := \frac{f(a_k) - f(a_k + p_k)}{m_k(0) - m_k(p_k)}$ \\
&nbsp;&nbsp; **if** $\rho_k < \frac{1}{4}$ \\
&nbsp;&nbsp;&nbsp;&nbsp; $\Delta\_{k+1} := \frac{1}{4} \Delta_k$ \\
&nbsp;&nbsp; **else** \\
&nbsp;&nbsp; &nbsp;&nbsp;**if** $\rho_k > \frac{3}{4}$ and $||p_k|| = \Delta_k$ \\
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp; $\Delta\_{k+1} := \min(2 \Delta_k, \hat{\Delta})$ 
\\
&nbsp;&nbsp; &nbsp;&nbsp;**else** \\
&nbsp;&nbsp;&nbsp;&nbsp; &nbsp;&nbsp; $\Delta\_{k+1} := \Delta_k$ \\
&nbsp;&nbsp; **if** $\rho\_k > \nu$ \\
&nbsp;&nbsp; &nbsp;&nbsp; $a\_{k+1} := a_k + p_k$\\
&nbsp;&nbsp;**else** \\
&nbsp;&nbsp; &nbsp;&nbsp; $a\_{k+1} := a_k$\\
**end for** 

We will describe two algorithms for minimising $m_k(p)$, the *Cauchy point* and the 
*dogleg* methods. The Cauchy point first solves a linear version of $m_k$ defined as

$$
p^s_k = \min_{p \in \mathcal{R}^n} f(a_k) + g_k^T p \text{ for }||p|| \le \Delta_k
$$

Subsequently, $p^s_k$ is used to find the scalar $\tau_k > 0$ such that

$$
\tau_k = \min_{\tau \ge 0} m_k(\tau p_k^s) \text{ for }||\tau p_k^s|| \le \Delta_k
$$

Finally, the Cauchy point is given as $p_k^C = \tau_k p_k^s$.

The solution to this problem can be shown to be

$$
p_k^C = -\tau_k \frac{\Delta_k}{|| g_k ||} g_k,
$$

where 

$$
\tau_k = \begin{cases}
1 & \text{if }g_k^T B_k g_k \le 0 \\
\min (||g_k||^3 / (\Delta_k g_k^T B_k g_k), 1) & \text{otherwise}.
\end{cases}
$$

The second method we describe is the *dogleg* method, which is applicable when $B_k$ is 
a positive definite matrix. If the original hessian is positive definite then this 
method is directly applicable, or one of the quasi-Newton positive definite 
approximation to the hessian could also be used. The dogleg method is derived by 
considering the path of the $p$ that minimises $m_k(p)$ with increasing $\Delta_k$, 
which forms a curved path in parameter space. The method approximates this path with two 
straight line segments. The first segment follows the steepest descent direction and is 
given by

$$
p_k^U = -\frac{g_k^T g_k}{g_k^T B_k g_k} g_k
$$

The second step is along the path between $p_k^U$ and $p^B_k = -B_k^{-1} g_k$. In the 
case where $p_k^B$ is *inside* the trust region $||p_k^B|| \le \Delta_k$ then $p_k^B$ 
can be used without modification. Otherwise the point of intersection with the 
trust-region radius must be calculated, which can be done by solving the following 
quadratic equation


$$
||p_k^U + (\tau - 1)(p_k^B - p_k^U)||^2 = \Delta_k^2
$$

with the second segment being defined by

$$
\tilde{p}_k(\tau) = \begin{cases}
\tau p_k^U, & 0 \le \tau 1, \\
p_k^U + (\tau - 1)(p_k^B - p_k^U), & 1 \le \tau \le 2.
\end{cases}
$$


### Problems

1. Let $f(x) = 10 (x_2 − x^2_1)^2 + (1 − x_1)^2$. At $x = (0,−1)$ draw the contour lines 
   of the quadratic model 
   
   $$
   m_k(p) = f(a_k) + g_k^T p + \frac{1}{2} p^T B_k p
   $$
   
   assuming that $B_k$ is the Hessian of $f$. Draw the family of solutions of $\min_{p 
   \in \mathcal{R}^n}m_k(p)$ so that $||p|| \le \Delta_k$  as the trust region radius 
   varies from $\Delta_k = 0$ to $\Delta_k = 2$. Repeat this at $x = (0, 0.5)$.
2. Write a program that implements the dogleg method. Choose $B_k$ to be the exact 
   Hessian. Apply it to solve [Rosenbrock’s 
   function](https://en.wikipedia.org/wiki/Rosenbrock_function). Experiment with the 
   update rule for the trust region by changing the constants in the trust region 
   algorithm given above, or by designing your own rules.

### Other reading

- Numerical optimization by Nocedal, Jorge; Wright, Stephen J., 1960-
- Bazaraa, Sherali, and Shetty, Nonlinear Optimization, 2/e, Wiley
- Griva, Nash, and Sofer, Linear and Nonlinear Optimization, 2/e, SIAM Press
- Luenberger, Linear and Nonlinear Programming, 3/e, Springer
- Bertsekas, Nonlinear Programming, Athena
- Ruszczynski, Nonlinear Optimization, Princeton University Press
- Broyden, C. G. (1972). "Quasi-Newton Methods". In Murray, W. (ed.). Numerical Methods 
  for Unconstrained Optimization. London: Academic Press. pp. 87–106. ISBN 
  0-12-512250-0.

### Software

- Scipy has a wide variety of (mostly) line search and trust region algorithms in 
  [`scipy.optimize`](https://docs.scipy.org/doc/scipy/reference/optimize.html)
- [Pints](https://pints.readthedocs.io/en/latest/optimisers/index.html) (Parameter 
  Inference on Time Series problems) has a number of global optimisation algorithms
- [NLOpt](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/) provides a common 
  interface for a large number of routines for non-linear optimisation
- [pygmo](https://github.com/esa/pygmo2) is a scientific Python library for massively 
  parallel optimisation


### Global Optimization

simplest is multi-start
- evolutionary strategies: Covariance matrix adaptation evolution strategy (CMA-ES) is a 
  popular example: 

[1]	The CMA Evolution Strategy: A Tutorial Nikolaus Hanse, arxiv 
    https://arxiv.org/abs/1604.00772
[2]	Hansen, Mueller, Koumoutsakos (2006) “Reducing the time complexity of the derandomized evolution strategy with covariance matrix adaptation (CMA-ES)”. Evolutionary Computation https://doi.org/10.1162/106365603321828970

- Others: , swarm algorithms inspired from bees/ants etc. e.g. Particle Swarm 
  Optimisation: 
  
  [1]	Kennedy, Eberhart (1995) Particle Swarm Optimization. IEEE International 
      Conference on Neural Networks https://doi.org/10.1109/ICNN.1995.488968



a wide variety of metaheuristics algorithms (wikipedia has a large list: 
https://en.wikipedia.org/wiki/Metaheuristic)

  - See a (Non-exhaustive) list of methods here: 

https://github.com/pints-team/pints/issues/684


## Textbooks



Gould, Nicholas I. M.; Toint, Philippe L. (2000). "A Quadratic Programming Bibliography" 
(PDF). RAL Numerical Analysis Group Internal Report 2000-1. 
ftp://ftp.numerical.rl.ac.uk/pub/qpbook/qp.pdf


## Software

{% include links.md %}
