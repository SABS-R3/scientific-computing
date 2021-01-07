from sympy import symbols, simplify, pprint
import sympy as sp

pU = symbols('pU')
pB = symbols('pB')
pUx = symbols('pUx')
pUy = symbols('pUy')
pBx = symbols('pBx')
pBy = symbols('pBy')
tau = symbols('tau')
delta = symbols('delta')

quadratic = \
    (pUx + (tau - 1)*(pBx - pUx))**2 + \
    (pUy + (tau - 1)*(pBy - pUy))**2 + \
    delta**2

quadratic = sp.collect(sp.expand(quadratic), tau)
pprint(quadratic)
tau1, tau2 = sp.solve(quadratic, tau)
pprint(simplify(tau1))
pprint(simplify(tau2))

