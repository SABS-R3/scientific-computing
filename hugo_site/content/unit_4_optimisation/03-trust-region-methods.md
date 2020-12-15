---
title: "Trust Region Methods"
weight: 3 
markup: "mmark"
---

### Saddle points

Saddle point pose a particular challenge in non-linear optimisation, particularly in 
higher dimensions. The plots below show two examples of saddle points in two dimensions. 
Like local minima and maxima, these are stationary points where the gradient of the 
function is zero $\nabla f = 0$, but where the value of the function rises along certain 
directions and reduces along others (left plot). An alternative type of saddle point 
arises when the hessian is singular, and are characterised by a plateau around the 
stationary point, like the [monkey saddle](https://en.wikipedia.org/wiki/Monkey_saddle) 
depicted in the plot to the right. 

![Two examples of saddle points](/scientific-computing/images/unit_04/saddle.svg)

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

1. Let $f(x) = 10 \left( x_2 − x^2_1 \right)^2 + (1 − x_1)^2$. At $x = (0,−1)$ draw the 
   contour lines of the quadratic model 
   
   $$
   m_k(p) = f(a_k) + g_k^T p + \frac{1}{2} p^T B_k p,
   $$
   
   assuming that $B\_k$ is the Hessian of $f$. Draw the family of solutions of $\min_{p 
   \in \mathcal{R}^n}m_k(p)$ so that $||p|| \le \Delta_k$  as the trust region radius 
   varies from $\Delta_k = 0$ to $\Delta_k = 2$. Repeat this at $x = (0, 0.5)$.

2. Write a program that implements the dogleg method. Choose $B_k$ to be the exact 
   Hessian. Apply it to solve [Rosenbrock’s 
   function](https://en.wikipedia.org/wiki/Rosenbrock_function). Experiment with the 
   update rule for the trust region by changing the constants in the trust region 
   algorithm given above, or by designing your own rules.
