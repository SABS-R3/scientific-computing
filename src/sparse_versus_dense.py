import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(0, 1e6, 10)
plt.plot(x, 8.0 * (x**2) / 1e6, lw=5)

plt.xlabel(r'size $n$')

plt.ylabel('memory [MB]')
plt.savefig('sparse_versus_dense.svg')
