import numpy as np
import matplotlib.pylab as plt

def f(x):
    return 10*(x[1] - x[0]**2)**2 + (1-x[0])**2

def grad_f(x):
    return np.array([
        40*x[0] * (x[0]**2 - x[1]) + 2*x[0] - 2,
        -20*x[0]**2 + 20*x[1]
    ])

def hessian_f(x):
    return np.array([
        [120 * x[0]**2 - 40*x[1] + 2, -40*x[0]],
        [-40*x[0], 20],
    ])

# plot contours
x = np.linspace(-0.5, 2.0, 100)
y = np.linspace(-1.5, 2.0, 100)
X, Y = np.meshgrid(x, y)
Z = f([X, Y])
plt.contour(X, Y, Z, levels=30)
min_x = np.array([1., 1.])

for x0 in [np.array([0., -1.]), np.array([0., 0.5])]:
    N = 20
    x = np.empty((N, 2), dtype=float)
    plt.plot(x0[0], x0[1], '.')
    plt.plot(min_x[0], min_x[1], '+')
    for i, delta in enumerate(np.linspace(1e-3, 2, N)):
        # if minimum of f is within radius, this is the min point, otherwise it is
        # somewhere on the border of the search radius
        if np.sum((min_x - x0)**2) < delta**2:
            x[i, :] = min_x
        else:
            min_f = 1e10
            for omega in np.linspace(0, 2*np.pi, 200):
                test_x = x0 + delta * np.array([np.cos(omega), np.sin(omega)])
                test_f = f(test_x)
                if test_f < min_f:
                    x[i, :] = test_x
                    min_f = test_f
    plt.plot(x[:, 0], x[:, 1], '-.', linewidth=2)

plt.show()

def line_sphere_intersect(o, u, r2):
    """
    algorithm from https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection

    here we assume that the line intersects the sphere, and we are only interested in
    the most positive solution
    """
    u_dot_o = u.dot(o)
    det = u_dot_o**2 - (np.dot(o, o) - r2)
    return -u_dot_o + np.sqrt(det)

def dogleg(delta, g, B, x0):
    """
    implements the dogleg method of minimising m(p)

    returns min p and a boolean that is True if min p lies on the trust region border
    """
    delta2 = delta**2
    pB = - np.linalg.solve(B, g)

    # if pB is within trust region then use that
    if np.dot(pB, pB) < delta2:
        return pB, False

    pU = - np.dot(g, g) / np.dot(g, B @ g) * g

    # find intersection of segment x0 -> pU -> pB with trust region
    if np.dot(pU, pU) > delta2:
        u = (pU - x0) / np.linalg.norm(pU - x0)
        intersect = line_sphere_intersect(x0, u, delta2)
        return x0 + intersect * u, True
    else:
        u = (pB - pU) / np.linalg.norm(pB - pU)
        intersect = line_sphere_intersect(pU, u, delta2)
        return pU + intersect * u , True


def minimise(f, grad_f, hessian_f, x0, tol=1e-5, max_iter=1000):
    nu = 0.2
    delta_hat = 1.0
    delta = delta_hat
    a = np.copy(x0)
    a_store = np.empty((max_iter+1, len(x0)), dtype=float)
    a_store[0, :] = a
    j = 1
    for i in range(max_iter):
        g = grad_f(a)
        B = hessian_f(a)
        f_a = f(a)
        p, p_equal_delta = dogleg(delta, g, B, a)
        rho = (f_a - f(a + p)) / (-g.dot(p) - 0.5*p.dot(B @ p))
        if rho < 0.25:
            delta *= 0.25
        else:
            if rho > 0.75 and p_equal_delta:
                delta = min(2*delta, delta_hat)
        if rho > nu:
            a += p
            a_store[j, :] = a
            j += 1
            if np.linalg.norm(p) < tol:
                break
    return a_store[:j, :]

# plot contours
x = np.linspace(-0.5, 2.0, 100)
y = np.linspace(-1.5, 2.0, 100)
X, Y = np.meshgrid(x, y)
Z = f([X, Y])
plt.contour(X, Y, Z, levels=30)
min_x = np.array([1., 1.])

# do the minimisation for both points and plot all iterations
for x0 in [np.array([0., -1.]), np.array([0., 0.5])]:
    plt.plot(x0[0], x0[1], '.')
    plt.plot(min_x[0], min_x[1], '+')
    a_store = minimise(f, grad_f, hessian_f, x0)
    print('done in', a_store.shape[0], 'iterations')
    plt.plot(a_store[:, 0], a_store[:, 1])

plt.show()

