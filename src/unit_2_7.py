import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg
import scipy.optimize
import matplotlib.pylab as plt
import time

def buildA(N):
    dx = 1 / N
    nvar = (N - 1)**2;
    e1 = np.ones((nvar), dtype=float);
    e2 = np.copy(e1)
    e2[::N-1] = 0
    e3 = np.copy(e1)
    e3[N-2::N-1] = 0
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




num = 20
times = np.empty((num, 2, 6), dtype=float)
iterations = np.empty((num, 2, 3), dtype=int)
times[:] = np.nan
iterations[:] = np.nan
Ns = np.logspace(0.5, 2.0, num=num, dtype=int)
for j, buildf in enumerate((buildf1, buildf2)):
    for i, N in enumerate(Ns):
        print('doing j=',j,' and N=',N)
        A = buildA(N)
        Adense = A.toarray()
        f = buildf(N)

        t0 = time.perf_counter()
        lu_A = scipy.linalg.lu_factor(Adense)
        x_using_lu = scipy.linalg.lu_solve(lu_A, f)
        t1 = time.perf_counter()
        times[i, j, 0] = t1 - t0
        np.testing.assert_almost_equal(A @ x_using_lu, f)

        t0 = time.perf_counter()
        cho_A = scipy.linalg.cho_factor(A.toarray())
        x_using_cho = scipy.linalg.cho_solve(cho_A, f)
        t1 = time.perf_counter()
        np.testing.assert_almost_equal(A @ x_using_cho, f)
        np.testing.assert_almost_equal(x_using_cho, x_using_lu)
        times[i, j, 1] = t1 - t0

        if N <= 64:
            t0 = time.perf_counter()
            invA = scipy.linalg.inv(Adense)
            x_using_inv = invA @ f
            t1 = time.perf_counter()
            np.testing.assert_almost_equal(x_using_inv, x_using_lu)
            times[i, j, 2] = t1 - t0

        global iters
        def count_iters(x):
            global iters
            iters += 1

        iters = 0
        t0 = time.perf_counter()
        x_using_cg, info = sp.linalg.cg(A, f, callback=count_iters)
        t1 = time.perf_counter()
        np.testing.assert_almost_equal(x_using_cg.reshape(-1, 1), x_using_lu,
                                       decimal=5)
        iterations[i, j, 0] = iters
        times[i, j, 3] = t1 - t0

        iters = 0
        t0 = time.perf_counter()
        x_using_bicgstab, info = sp.linalg.bicgstab(A, f, callback=count_iters)
        t1 = time.perf_counter()
        np.testing.assert_almost_equal(x_using_bicgstab.reshape(-1, 1), x_using_lu,
                                       decimal=5)
        iterations[i, j, 1] = iters
        times[i, j, 4] = t1 - t0

        iters = 0
        t0 = time.perf_counter()
        x_using_gmres, info = sp.linalg.gmres(A, f, callback=count_iters)
        t1 = time.perf_counter()
        np.testing.assert_almost_equal(x_using_gmres.reshape(-1, 1), x_using_lu,
                                       decimal=5)
        iterations[i, j, 2] = iters
        times[i, j, 5] = t1 - t0

plt.loglog(Ns, times[:,0,:], ls='-')
plt.gca().set_prop_cycle(None)
plt.loglog(Ns, times[:,1,:], ls='--')
plt.xlabel('N')
plt.ylabel('time')
plt.legend(['lu', 'cholesky', 'inv', 'cg', 'bicgstab', 'gmres'])
plt.show()

plt.plot(Ns, iterations[:,0,:], ls='-')
plt.gca().set_prop_cycle(None)
plt.plot(Ns, iterations[:,1,:], ls='--')
plt.xlabel('N')
plt.ylabel('iterations')
plt.legend(['cg', 'bicgstab', 'gmres'])
plt.show()

print(iterations)

# Krylov subspace solvers only take 1 iteration to solve with b = f1 because x is a
# scalar multiple of f (i.e. f1 is an eigenvector for A, and
A = buildA(10)
f = buildf1(10)
np.testing.assert_almost_equal(A@f, 19.57739348*f)

