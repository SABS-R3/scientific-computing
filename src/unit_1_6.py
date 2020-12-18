import numpy as np
import matplotlib.pylab as plt
import scipy
import scipy.linalg

def sample_zero_mean_random_field(covariance_matrix):
    N = covariance_matrix.shape[0]
    L, _ = scipy.linalg.cho_factor(covariance_matrix)
    indep_random_sample = np.random.normal(size=(N, 1))
    return np.triu(L) @ indep_random_sample

n = 100
sigma1 = 1.0
sigma2 = 10.0
xs = np.arange(0, n)
K1 = sigma1**2 * np.eye(n)
K2 = sigma1**2 * np.exp(
    -(xs[:, np.newaxis] - xs[:, np.newaxis].T)**2 / sigma2**2
)
K2 += 1e-5 * np.eye(n)

plt.plot(sample_zero_mean_random_field(K1))
plt.plot(sample_zero_mean_random_field(K2))
plt.xlabel('x')
plt.ylabel('y')
plt.show()
