import matplotlib.pylab as plt
import numpy as np
from matplotlib import cm
import scipy.optimize

def steepest_descent(x, f, grad_f, hessian_f):
    return -grad_f(x)

def newton(x, f, grad_f, hessian_f):
    A = hessian_f(x)
    return -np.linalg.solve(A, grad_f(x))

def backtracking(f, grad_f, x, p):
    alpha = 1.0
    c = 1e-4
    rho = 0.5
    f_x = f(x)
    c_grad_f_x_p = c*grad_f(x).dot(p)
    while f(x + alpha*p) > f_x + alpha * c_grad_f_x_p:
        alpha = rho * alpha
    return alpha

def scipy_line_search(f, grad_f, x, p):
    return scipy.optimize.line_search(f, grad_f, x, p)[0]

def minimise(f, grad_f, hessian_f, x0, tol=1e-5, max_iter=1000, search_direction=steepest_descent, line_search=backtracking):
    x = np.copy(x0)
    alpha = 1.0
    p = np.zeros_like(x)
    p[0] = 1.
    i = 0
    array_x = np.empty((max_iter+1, len(x)), dtype=float)
    array_alpha = np.empty(max_iter, dtype=float)
    array_x[0, :] = x
    while np.linalg.norm(alpha * p) > tol and i < max_iter:
        p = search_direction(x, f, grad_f, hessian_f)
        alpha = line_search(f, grad_f, x, p)
        x += alpha * p
        array_x[i+1, :] = x
        array_alpha[i] = alpha
        i += 1
    return array_x[:i,:], array_alpha[:i]

def rosenbrock(x, y):
    return (1 - x)**2 + (y - x**2)**2

def rosenbrock_2d(x):
    return (1 - x[0])**2 + (x[1] - x[0]**2)**2

def grad_rosenbrock_2d(x):
    return np.array([
        4*x[0] * (x[0]**2 - x[1]) + 2*x[0] - 2,
        -2*x[0]**2 + 2*x[1]
    ])

def hessian_rosenbrock_2d(x):
    return np.array([
        [12 * x[0]**2 - 4*x[1] + 2, -4*x[0]],
        [-4*x[0], 2],
    ])

# verify grad and hessian using finite differences
h = 1e-5
x0 = np.array([0., 0.])
x0_x = np.array([h, 0.])
x0_y = np.array([0., h])
fd_grad_rosenbrock = np.array([
    (rosenbrock_2d(x0_x) - rosenbrock_2d(x0)) / h,
    (rosenbrock_2d(x0_y) - rosenbrock_2d(x0)) / h,
])
np.testing.assert_almost_equal(fd_grad_rosenbrock, grad_rosenbrock_2d(x0), decimal=4)
fd_hessian_rosenbrock = np.vstack((
    (grad_rosenbrock_2d(x0_x) - grad_rosenbrock_2d(x0)) / h,
    (grad_rosenbrock_2d(x0_y) - grad_rosenbrock_2d(x0)) / h,
))
np.testing.assert_almost_equal(fd_hessian_rosenbrock, hessian_rosenbrock_2d(x0),
                               decimal=4)

# plot rosenbrock surface
x = np.linspace(-1.5, 2.0, 100)
y = np.linspace(0.2, 3.0, 100)
X, Y = np.meshgrid(x, y)
Z = np.vectorize(rosenbrock)(X, Y)

fig = plt.figure()
ax = fig.gca(projection='3d')
fig2 = plt.figure()
ax2 = fig2.add_subplot(1,1,1)

surf = ax.plot_surface(X, Y, Z, alpha=0.5, cmap=cm.coolwarm)

# pick two initial starting points
x01 = np.array([1.2, 1.2])
x02 = np.array([-1.2, 1])

# run all methods and plot results
for x0 in [x01, x02]:
    for p_method in [steepest_descent, newton]:
        for ls_method in [backtracking, scipy_line_search]:
            label = '{}-{}'.format(p_method.__name__, ls_method.__name__)
            pos, alpha = \
                minimise(rosenbrock_2d, grad_rosenbrock_2d, hessian_rosenbrock_2d, x0,
                         search_direction=p_method, line_search=ls_method)
            print('method', label, 'finished in', len(alpha), 'iterations')
            xs = pos[:, 0]
            ys = pos[:, 1]
            ax.plot(xs, ys, rosenbrock(xs, ys), label=label, linewidth=2)
            ax2.plot(alpha, label=label)
ax2.legend()
plt.show()
