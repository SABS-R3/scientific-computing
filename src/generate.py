import matplotlib.pylab as plt
import numpy as np
import seaborn
seaborn.set()

fig, ax = plt.subplots(figsize=(6,4))
x = np.linspace(-2, 4, 100)
y = (1 + x) / 2
ax.plot(x, y)
y = (3 + x) / 3
ax.plot(x, y)
ax.set_aspect('equal')
ax.grid(True, which='both')
fig.savefig('01-sim1.svg')
plt.close(fig)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12,4))
y = (1 + x) / 2
ax1.plot(x, y)
y = (3 + x) / 2
ax1.plot(x, y)
ax1.set_aspect('equal')
ax1.grid(True, which='both')

y = (1 + x) / 2
ax2.plot(x, y)
y = (1 + x) / 2
ax2.plot(x, y)
ax2.set_aspect('equal')
ax2.grid(True, which='both')
fig.savefig('01-sim2.svg')
plt.close(fig)

