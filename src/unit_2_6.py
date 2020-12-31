import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import scipy.optimize
import matplotlib.pylab as plt

def buildA(N):
    dx = 1 / N
    nvar = (N - 1)**2;
    e1 = np.ones((nvar), dtype=float);
    e2 = np.copy(e1)
    e2[:N-1:] = 0
    e3 = np.copy(e1)
    e3[N-2:N-1:] = 0
    A = sp.spdiags(
        (-e1, -e3, 4*e1, -e2, -e1),
        (-(N-1), -1, 0, 1, N-1), nvar, nvar
    )
    A = A / dx**2;
    return A


def buildf1(N):
    x = np.arange(0, 1, 1/N).reshape(N, 1)
    y = x.T
    f = np.dot(np.sin(np.pi*x), np.sin(np.pi*y))
    return f[1:,1:].reshape(-1,1)

def buildf2(N):
    x = np.arange(0, 1, 1/N).reshape(N, 1)
    y = x.T
    f = np.dot(np.maximum(x,1-x), np.maximum(y,1-y))
    return f[1:,1:].reshape(-1, 1)

def jacobi(A, b, x0=None, tol=1e-5, max_iter=1000):
    if x0 is None:
        x0 = np.zeros_like(b)
    x = np.copy(x0)
    b_norm = np.linalg.norm(b)

    # jacobi method: M = D
    M = A.diagonal().reshape(-1, 1)
    invM = 1/M

    # main relaxation iteration
    for i in range(max_iter):
        r = b - A @ x
        error = np.linalg.norm(r) / b_norm
        if error < tol:
            break
        x += invM * r
    return x, i

def SOR(A, b, omega, x0=None, tol=1e-5, max_iter=300):
    if x0 is None:
        x0 = np.zeros_like(b)
    x = np.copy(x0)
    b_norm = np.linalg.norm(b)

    # SOR method
    D = sp.spdiags((A.diagonal()), (0), *A.shape)
    L = sp.tril(A, k=-1)
    M = (1/omega) * D + L

    # main relaxation iteration
    for i in range(max_iter):
        r = b - A @ x
        error = np.linalg.norm(r) / b_norm
        if error < tol:
            break
        x += sp.linalg.spsolve_triangular(M, r)
    return x, i


num = 20
iterations = np.empty((num, 2), dtype=int)
iterations[:] = np.nan
Ns = np.logspace(0.5, 1.5, num=num, dtype=int)
for j, buildf in enumerate((buildf1, buildf2)):
    for i, N in enumerate(Ns):
        A = buildA(N)
        f = buildf(N)
        max_iter = 10*N
        x, iters = jacobi(A, f, max_iter=max_iter)
        if i < max_iter:
            iterations[i, j] = iters

plt.plot(Ns, iterations)
plt.xlabel('N')
plt.ylabel('iterations')
plt.show()

N = 64
A = buildA(N)
f = buildf2(N)

def SOR_iterations(omega):
    x, i = SOR(A, f, omega, max_iter=10, tol=1e-32)
    return np.linalg.norm(A @ x - f)

res = scipy.optimize.minimize_scalar(SOR_iterations, bracket=[0.1, 1.0, 1.99], tol=1e-2)
print('ideal omega is', res.x, 'versus analytic value of', 2 / (1 + np.sin(np.pi/N)))


