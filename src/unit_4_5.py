import autograd.numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, shgo
from autograd import grad


def convex_function(x):
    return np.sum(np.array(x)**2, axis=0)

def rosenbrock(x):
    a = 1.0
    b = 100.0
    return (a-x[0])**2 + b*(x[1]-x[0]**2)**2

def rastrigin(x):
    A = 10.0
    n = len(x)
    ret = A*n
    for i in range(n):
        ret += x[i]**2 - A * np.cos(2*np.pi*x[i])
    return ret

def visualise_functions():
    nx, ny = (100, 100)
    x = np.linspace(-5, 5, nx)
    y = np.linspace(-5, 5, ny)
    xv, yv = np.meshgrid(x, y)
    eval_convex = convex([xv, yv])
    eval_rastrigin = rastrigin([xv, yv])
    eval_rosenbrock = rosenbrock([xv, yv])

    plt.contourf(x, y, eval_convex)
    plt.colorbar()
    plt.show()
    plt.clf()
    plt.contourf(x, y, eval_rosenbrock)
    plt.colorbar()
    plt.show()
    plt.clf()
    plt.contourf(x, y, eval_rastrigin)
    plt.colorbar()
    plt.show()

def optimize(function, method, autodiff):
    eval_points_x = []
    eval_points_y = []
    def fill_eval_points(xk):
        eval_points_x.append(xk[0])
        eval_points_y.append(xk[1])

    x0 = np.array([2.5, 2.5])
    if autodiff:
        jac = grad(function)
    else:
        jac = None

    if method == 'shgo':
        bounds = [(-10, 10), (-10.0, 10.0)]
        res = shgo(function, bounds, callback=fill_eval_points,
                    options={'disp': True})
    else:
        res = minimize(function, x0, method=method, callback=fill_eval_points,
                        jac = jac,
                    options={'disp': True})

    nx, ny = (100, 100)
    x = np.linspace(-5, 5, nx)
    y = np.linspace(-5, 5, ny)
    xv, yv = np.meshgrid(x, y)
    eval_function = function([xv, yv])

    print(function.__name__)
    plt.clf()
    plt.contourf(x, y, eval_function)
    plt.plot(eval_points_x, eval_points_y, 'x-k')
    plt.colorbar()
    neval = res.nfev
    try:
        neval += res.njev
    except:
        print('no njev')
    try:
        neval += res.nhev
    except:
        print('no nhev')
    plt.title('iterations: {} evaluations: {}'.format(res.nit, neval))

    if autodiff:
        ext = '-auto.pdf'
    else:
        ext = '.pdf'

    plt.savefig(function.__name__ + '-' + method + ext)


if __name__ == '__main__':
    for f in [convex_function, rosenbrock, rastrigin]:
        for m in ['shgo','nelder-mead', 'cg', 'bfgs', 'newton-cg']:
            for a in [False, True]:
                if m == 'newton-cg' and a == False:
                    continue
                if m == 'shgo' and a == True:
                    continue
                optimize(f, m, a)
