import numpy as np
import matplotlib.pylab as plt
import scipy
import scipy.linalg
import sys

def solve_triangular(A, b):
    n = len(b)
    x = np.empty_like(b)
    for i in range(n-1, -1, -1):
        x[i] = b[i]
        for j in range(n-1, i, -1):
            x[i] -= A[i, j] * x[j]
        x[i] /= A[i, i]
    return x


def gaussian_elimination(A):
    m, n = A.shape

    # initialise the pivot row and column
    h = 0
    k = 0
    while h < m and k < n:
        # Find the k-th pivot:
        i_max = np.argmax(A[h:, k]) + h
        if A[i_max, k] == 0:
            # No pivot in this column, pass to next column
            k = k+1
        else:
            # swap rows
            A[[h, i_max], :] = A[[i_max, h], :]
            # Do for all rows below pivot:
            for i in range(h+1, m):
                f = A[i, k] / A[h, k]
                # Fill with zeros the lower part of pivot column:
                A[i, k] = 0
                # Do for all remaining elements in current row:
                for j in range(k + 1, n):
                    A[i, j] = A[i, j] - A[h, j] * f
            # Increase pivot row and column
            h = h + 1
            k = k + 1
    return A

def solve_gaussian_elimination(A, b):
    augmented_system = np.concatenate((A, b.reshape(-1, 1)), axis=1)
    gaussian_elimination(augmented_system)
    return solve_triangular(augmented_system[:, :-1], augmented_system[:, -1])

def random_matrix(n):
    R = np.random.rand(n, n)
    A = np.zeros((n, n))
    triu = np.triu_indices(n)
    A[triu] = R[triu]
    return A

def random_non_singular_matrix(n):
    A = np.random.rand(n, n)
    while np.linalg.cond(A) > 1/sys.float_info.epsilon:
        A = np.random.rand(n, n)
    return A

As = [
    np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
    random_non_singular_matrix(3),
    random_non_singular_matrix(4),
    random_non_singular_matrix(5),
    random_non_singular_matrix(6),
]

bs = [
    np.array([1, 2, 3]),
    np.random.rand(3),
    np.random.rand(4),
    np.random.rand(5),
    np.random.rand(6),
]

for A, b in zip(As, bs):
    x_scipy = scipy.linalg.solve(A, b)
    x_mine = solve_gaussian_elimination(A, b)
    np.testing.assert_almost_equal(x_scipy, x_mine)

A = np.array([
    [4.5, 3.1],
    [1.6, 1.1],
])
x1 = scipy.linalg.solve(A, np.array([19.249, 6.843]))
x2 = scipy.linalg.solve(A, np.array([19.25, 6.84]))
print('percentage error is ', np.linalg.norm(x2 - x1) / np.linalg.norm(x1) * 100)
print('condition number is ', np.linalg.cond(A))


