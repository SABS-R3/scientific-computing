import matplotlib.pyplot as plt
import numpy as np

x = np.outer(np.linspace(-2, 2, 30), np.ones(30))
y = x.copy().T
z = x**2 - y**2

fig = plt.figure(figsize=plt.figaspect(0.5))
ax1 = fig.add_subplot(1, 2, 1, projection='3d')

ax1.plot_surface(x, y, z, cmap='viridis', edgecolor='none')
ax1.scatter3D([0], [0], [0], color='red', marker='o', alpha=1.0)
ax1.set_xlabel(r'$x$')
ax1.set_ylabel(r'$y$')
ax1.set_zlabel(r'$x^2 - y^2$')

z = x**3 - 3*x*y**2

ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.plot_surface(x, y, z, cmap='viridis', edgecolor='none')
ax2.scatter3D([0], [0], [0], color='red', marker='o', alpha=1.0)
ax2.set_xlabel(r'$x$')
ax2.set_ylabel(r'$y$')
ax2.set_zlabel(r'$x^3 - 3xy^2$')
plt.savefig('saddle.svg')
plt.show()




