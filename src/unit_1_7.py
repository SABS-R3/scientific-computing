import pandas as pd
import matplotlib.pylab as plt
import numpy as np
import scipy.linalg

names = ['year', 'month', 'maxTemp', 'minTemp', 'hoursFrost', 'rain', 'hoursSun']
df    = pd.read_csv('OxfordWeather.txt',
                     delim_whitespace=True, header=None, names=names)

x = df.month.values.reshape(-1,1)
y = df.hoursSun.values.reshape(-1,1)

# we have to sort the data according to x if we want a nice plot
i = np.argsort(x,axis=0).reshape(-1)
x = x[i]
y = y[i]

# polynomial model $y = a x^2 + b x + c$,
# or $M beta = y$ where $beta = (a, b, c)$
M = np.concatenate([np.ones_like(x), x, x**2], axis=1)
Q, R = scipy.linalg.qr(M, mode='economic')
beta = np.linalg.solve(R, Q.T @ y)

# check against lstsq
np.testing.assert_almost_equal(
    beta, np.linalg.lstsq(M, y, rcond=None)[0]
)

plt.plot(x, y, 'o')
plt.plot(x, M @ beta, 'r')
plt.xlabel('month')
plt.ylabel('hoursSun')
plt.show()
