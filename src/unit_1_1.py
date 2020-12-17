import numpy as np
import matplotlib.pylab as plt
import scipy
import scipy.linalg
import sys


x = np.linspace(-1, 4, 100)

def plot_lines(b1, b2, b3):
    y1 = (b1 - x) / 2
    y2 = x - b2
    y3 = b3 * np.ones_like(x)
    plt.plot(x, y1)
    plt.plot(x, y2)
    plt.plot(x, y3)
    plt.show()


plot_lines(2, 2, 1)
plot_lines(0, 0, 0)
plot_lines(2, -1, 1)


def solve_triangular(A, b):
    n = len(b)
    x = np.empty_like(b)
    for i in range(n-1, -1, -1):
        x[i] = b[i]
        for j in range(n-1, i, -1):
            x[i] -= A[i, j] * x[j]
        x[i] /= A[i, i]
    return x

def random_upper_triangular(n):
    R = np.random.rand(n, n)
    A = np.zeros((n, n))
    triu = np.triu_indices(n)
    A[triu] = R[triu]
    return A

As = [
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    random_upper_triangular(3),
    random_upper_triangular(4),
    random_upper_triangular(5),
    random_upper_triangular(6),
]

bs = [
    np.array([1, 2, 3]),
    np.random.rand(3),
    np.random.rand(4),
    np.random.rand(5),
    np.random.rand(6),
]

for A, b in zip(As, bs):
    x_scipy = scipy.linalg.solve_triangular(A, b)
    x_mine = solve_triangular(A, b)
    np.testing.assert_almost_equal(x_scipy, x_mine)


