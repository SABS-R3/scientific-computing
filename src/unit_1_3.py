import numpy as np
import matplotlib.pylab as plt
import scipy
import scipy.linalg
import sys

def lu_decomposition(A):
    m, n = A.shape

    LU = np.copy(A)
    pivots = np.empty(n, dtype=int)
    # initialise the pivot row and column
    h = 0
    k = 0
    while h < m and k < n:
        # Find the k-th pivot:
        pivots[k] = np.argmax(LU[h:, k]) + h
        if LU[pivots[k], k] == 0:
            # No pivot in this column, pass to next column
            k = k+1
        else:
            # swap rows
            LU[[h, pivots[k]], :] = LU[[pivots[k], h], :]
            # Do for all rows below pivot:
            for i in range(h+1, m):
                f = LU[i, k] / LU[h, k]
                # Store f as the new L column values
                LU[i, k] = f
                # Do for all remaining elements in current row:
                for j in range(k + 1, n):
                    LU[i, j] = LU[i, j] - LU[h, j] * f
            # Increase pivot row and column
            h = h + 1
            k = k + 1
    return LU, pivots

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

def pivots_to_row_indices(pivots):
    n = len(pivots)
    indices = np.array(range(0, n))
    for i, p in enumerate(pivots):
        indices[i], indices[p] = indices[p], indices[i]
    return indices

def calculate_L_mult_U(LU):
    L = np.tril(LU)
    np.fill_diagonal(L, 1)
    U = np.triu(LU)
    return L @ U

for A in As:
    LU_scipy, pivots_scipy = scipy.linalg.lu_factor(A)
    row_indices_scipy = pivots_to_row_indices(pivots_scipy)
    LU_mine, pivots_mine = lu_decomposition(A)
    row_indices_mine = pivots_to_row_indices(pivots_mine)

    np.testing.assert_almost_equal(calculate_L_mult_U(LU_scipy),  A[row_indices_scipy])
    np.testing.assert_almost_equal(calculate_L_mult_U(LU_mine),  A[row_indices_mine])

