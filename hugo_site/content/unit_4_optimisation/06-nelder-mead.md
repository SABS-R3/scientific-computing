---
title: "Nelder-Mead method"
weight: 6 
markup: "mmark"
---

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

![Nelder–Mead simplex search over the Rosenbrock banana 
function](/scientific-computing/images/unit_04/Nelder-Mead_Rosenbrock.gif) 

[Scholarpedia](http://www.scholarpedia.org/article/Nelder-Mead_algorithm) and 
[Wikipedia](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) provide diagrams 
and pseudocode of the Nelder-Mead algorithm, as does Chapter 9.5 of the Nocedal and 
Write textbook given below

## Other Reading

- Nelder, John A.; R. Mead (1965). "A simplex method for function minimization". 
  Computer Journal. 7 (4): 308–313. doi:10.1093/comjnl/7.4.308.
- Numerical optimization by Nocedal, Jorge; Wright, Stephen J., 1960-, Chapter 9

### Problems

- Code up the Nelder-Mead algorithm and compare its performance against the steepest 
  descent, Newton and dogleg algorithms you did in the last lesson. You can evaluate 
  them on the 2D quadratic function $f(x, y) = x^2 + y^2$, the 2D [Rosenbrock
function](https://en.wikipedia.org/wiki/Rosenbrock_function) or on one of many different 
[optimisation test 
functions](https://en.wikipedia.org/wiki/Test_functions_for_optimization)


