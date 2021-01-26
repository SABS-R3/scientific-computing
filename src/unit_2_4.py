import scipy.sparse.linalg
import scipy.sparse as sp
import matplotlib.pylab as plt
import numpy as np

N = 100
e = np.ones(N, dtype=float)
A = sp.spdiags([e, -2*e, e], [-1, 0, 1], N, N, format='csc')

plt.spy(A)
plt.show()

L = 2*np.pi
x = np.linspace(0, L, N+2)
h = x[1] - x[0]
fcos = 2 * np.cos(x) / np.exp(x)
analytical_cos = -np.sin(x) / np.exp(x)
fsin = 2 * np.sin(x) / np.exp(x)
analytical_sin = np.cos(x) / np.exp(x)

# Use sparse lu decomposition, matrix is tridiagonal so no fill-in
splu = sp.linalg.splu(A / h**2)
fig, axs = plt.subplots(1, 2)
axs[0].spy(splu.L)
axs[0].set_title('L')
axs[1].spy(splu.U)
axs[1].set_title('U')
plt.show()

# Use decomposition to solve both problems
for f, a in zip([fcos, fsin], [analytical_cos, analytical_sin]):
    b = f[1:-1]
    b[0] -= a[0] / h**2
    b[-1] -= a[-1] / h**2
    v = splu.solve(b)
    plt.plot(x[1:-1], v, label='solution')
    plt.plot(x, a, label='analytical')
    plt.legend()
    plt.show()

import time
num = 100
times = np.empty(num, dtype=float)
times_dense = np.empty(num, dtype=float)
times_dense[:] = np.nan
Ns = np.logspace(0.5, 6, num=num, dtype=int)
for i, N in enumerate(Ns):
    e = np.ones(N, dtype=float)
    A = sp.spdiags([e, -2*e, e], [-1, 0, 1], N, N, format='csc')

    t0 = time.perf_counter()
    AA = A @ A
    t1 = time.perf_counter()
    times[i] = t1 - t0

    if N < 2000:
        Adense = A.toarray()
        t0 = time.perf_counter()
        AA = Adense @ Adense
        t1 = time.perf_counter()
        times_dense[i] = t1 - t0

plt.clf()
plt.loglog(Ns, times, label='sparse @')
plt.loglog(Ns, times_dense, label='dense @')
plt.xlabel('N')
plt.ylabel('time taken')
plt.legend()
plt.show()

times = np.empty(num, dtype=float)
times_dense = np.empty(num, dtype=float)
times_dense[:] = np.nan

for i, N in enumerate(Ns):
    e = np.ones(N, dtype=float)
    A = sp.spdiags([e, -2*e, e], [-1, 0, 1], N, N, format='csc')

    x = np.linspace(0, L, N+2)
    h = x[1] - x[0]
    fcos = 2 * np.cos(x) / np.exp(x)
    analytical_cos = -np.sin(x) / np.exp(x)

    A /= h**2

    b = fcos[1:-1]
    b[0] -= analytical_cos[0] / h**2
    b[-1] -= analytical_cos[-1] / h**2

    t0 = time.perf_counter()
    splu = sp.linalg.splu(A)
    v = splu.solve(b)
    t1 = time.perf_counter()
    times[i] = t1 - t0
    if N < 2000:
        Adense = A.toarray()
        t0 = time.perf_counter()
        lu = scipy.linalg.lu_factor(Adense)
        v = scipy.linalg.lu_solve(lu, b)
        t1 = time.perf_counter()
        times_dense[i] = t1 - t0

plt.clf()
plt.loglog(Ns, times, label='sparse LU')
plt.loglog(Ns, times_dense, label='dense LU')
plt.xlabel('N')
plt.ylabel('time taken')
plt.legend()
plt.show()


