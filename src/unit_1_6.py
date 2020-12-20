import numpy as np
import matplotlib.pylab as plt
import scipy
import scipy.linalg

n = 100
xs = np.arange(0, n)

def construct_covariance(sigma1, sigma2):
    K = sigma1**2 * np.exp(
        -(xs[:, np.newaxis] - xs[:, np.newaxis].T)**2 / sigma2**2
    )
    K += 1e-5 * np.eye(n)
    return K

def sample_zero_mean_random_field(covariance_matrix):
    L, _ = scipy.linalg.cho_factor(covariance_matrix, lower=True)
    indep_random_sample = np.random.normal(size=n)
    return np.dot(np.tril(L), indep_random_sample)

sigma1 = 1.0
sigma2 = 10.0
K1 = sigma1**2 * np.eye(n)
K2 = construct_covariance(sigma1, sigma2)

plt.plot(sample_zero_mean_random_field(K1))
plt.plot(sample_zero_mean_random_field(K2))
plt.show()

def log_likelihood(sigma1, sigma2, x):
    K = construct_covariance(sigma1, sigma2)
    L, lower = scipy.linalg.cho_factor(K, lower=True)

    logdet = 2 * np.sum(np.log(np.diag(L)))
    xinvSx = x.dot(scipy.linalg.cho_solve((L, lower), x))

    return -0.5 * (logdet + xinvSx)

sigma1 = np.linspace(0.5, 1.5, 100)
sigma2 = np.linspace(5.0, 15.0, 100)
Sigma1, Sigma2 = np.meshgrid(sigma1, sigma2)
x = sample_zero_mean_random_field(K2)

L = np.vectorize(log_likelihood, excluded=['x'])(Sigma1, Sigma2, x=x)
max_log_likelihood = log_likelihood(1, 10, x)
levels = np.linspace(0.9 * max_log_likelihood, max_log_likelihood, 5)

plt.clf()
contours = plt.contour(Sigma1, Sigma2, L, levels=levels, colors='black')
plt.clabel(contours, inline=True, fontsize=8)

plt.imshow(L, extent=[0.5, 1.5, 5.0, 15.0], origin='lower',
           cmap='RdGy', alpha=0.5, vmin=levels[0], aspect='auto')
c = plt.colorbar();
c.set_label(r'$\mathcal{L}$')
plt.xlabel(r'$\sigma_1$')
plt.ylabel(r'$\sigma_2$')
plt.show()

